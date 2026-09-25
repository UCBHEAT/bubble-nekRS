// PID controller that keeps the bubble (gas-phase centroid) at the domain
// center, by applying a spatially uniform acceleration to the whole domain.
//
// The controller output F_pid is a force per unit mass (an acceleration, in the
// same nondimensional units as the -1/Fr^2 gravity term in buoyancySource) that
// is added to every velocity point. In a fully periodic domain a uniform
// acceleration is exactly an acceleration of the reference frame: it shifts the
// whole velocity field uniformly and does not change the bubble/liquid relative
// motion (the net gravity + pressure gradient force on the domain is zero, as
// are the net viscous and surface tension forces, so in the continuum F_pid is
// the only term that changes the domain's mean momentum; the discretization
// does not conserve it exactly, and in the t = 0-10 run the mean velocity
// drifts from integral(F_pid dt) by ~5% of the frame velocity). Once the
// bubble reaches terminal velocity, the frame rises with the bubble and the
// liquid flows down past it.
//
//   e      = x_c - x_0          gas centroid offset from the domain center
//   u_gas                       gas-phase mean velocity, used as de/dt
//                               (derivative on measurement)
//   F_pid  = -(kp*e + ki*integral(e dt) + kd*u_gas)
//
// u_gas is the material velocity of the gas, which is exactly de/dt between
// level set reinitializations. Each CLSR reinit also shifts the centroid
// slightly (upward for a rising bubble, by ~3e-4 per reinit with dt = 0.0025
// and reinit every 10 steps, a drift of ~0.012-0.017 per unit time), which
// u_gas does not see. At steady state the integral term therefore settles at
// ki*integral = kd*(reinit drift rate) rather than 0, and slip_v (from u_gas)
// is a few percent lower than the rise rate of the level set centroid.
//
// Gains are read from the [PID] section of the .par file. The plant (F_pid ->
// x_c) is a double integrator, so the closed-loop characteristic polynomial is
// s^3 + kd*s^2 + kp*s + ki, which is stable (Routh-Hurwitz) for kd > 0, ki > 0
// and kd*kp > ki (which implies kp > 0); with ki = 0 it needs kd > 0, kp > 0.

#include <algorithm>
#include <fstream>
#include <sstream>

typedef struct pidState {
    // Gains, from the [PID] par section.
    dfloat kp;
    dfloat ki;
    dfloat kd;
    // Domain center (volume centroid of the mesh).
    dfloat center[3];
    // Gas phase volume and centroid offset from the domain center (the error e).
    dfloat gas_volume;
    dfloat error[3];
    // Time integral of the error. Written to data.csv so it can be restored on
    // restart.
    dfloat integral[3];
    // Gas and liquid phase mean velocities. u_gas is used as de/dt;
    // u_gas - u_liquid is the bubble slip (rise) velocity.
    dfloat u_gas[3];
    dfloat u_liquid[3];
    // PID force per unit mass, applied uniformly in the next time step.
    dfloat force[3];
    // Simulation time of the last update, for integrating the error.
    double time;
} pidState_t;

static pidState_t pid;

/**
 * Restore the integral term on restart from the data.csv row written at the
 * restart time (the rest of the controller state is re-measured from the
 * restart fields). If there is no such row, the integral starts from zero.
 *
 * @param restart_time simulation time of the restart file
 */
void pidRestoreIntegral(double restart_time)
{
    // data.csv rows are written with 4 decimals of time, so accept the row
    // closest to the restart time within this tolerance.
    const double tolerance = 1e-3;
    double state[4] = {0, 0, 0, 0}; // found, integral x, y, z

    if (platform->comm.mpiRank() == 0) {
        std::ifstream f("data.csv");
        std::string line;
        std::vector<std::string> header;
        double best = tolerance;
        while (std::getline(f, line)) {
            std::vector<std::string> cols;
            std::stringstream ss(line);
            std::string col;
            while (std::getline(ss, col, ',')) cols.push_back(col);
            if (!cols.empty() && cols[0] == "Time") {
                // Header. bubble.udf writes it only when it creates data.csv, so
                // a data.csv started with an older column layout (no pid_int_*
                // columns) is never restored from: remove data.csv when the
                // column layout changes.
                header = cols;
                continue;
            }
            const auto column = [&](const std::string& name) {
                const auto it = std::find(header.begin(), header.end(), name);
                return (it == header.end()) ? -1 : int(it - header.begin());
            };
            const int itime = column("Time");
            const int ix = column("pid_int_x"), iy = column("pid_int_y"), iz = column("pid_int_z");
            if (itime < 0 || ix < 0 || iy < 0 || iz < 0 || int(cols.size()) != int(header.size())) continue;
            try {
                const double t = std::stod(cols[itime]);
                // Keep the last row within the tolerance, so a row from a later
                // rerun supersedes an earlier one at the same time.
                if (std::abs(t - restart_time) <= best) {
                    const double integral[3] = {std::stod(cols[ix]), std::stod(cols[iy]), std::stod(cols[iz])};
                    best = std::max(std::abs(t - restart_time), 1e-12);
                    state[0] = 1;
                    for (int d = 0; d < 3; d++) state[1 + d] = integral[d];
                }
            } catch (const std::exception&) {
                // Skip malformed rows (e.g. a partially written last line).
            }
        }
    }
    MPI_Bcast(state, 4, MPI_DOUBLE, 0, platform->comm.mpiComm());

    for (int d = 0; d < 3; d++) pid.integral[d] = state[1 + d];
    if (platform->comm.mpiRank() == 0) {
        if (state[0] > 0) {
            printf("PID: restored integral=(%g, %g, %g) from data.csv at t=%g\n",
                    pid.integral[0], pid.integral[1], pid.integral[2], restart_time);
        } else {
            printf("PID: WARNING no data.csv row at restart time t=%g, integral starts from 0\n",
                    restart_time);
        }
    }
}

/**
 * Read gains from the [PID] par section and compute the domain center, and on
 * restart restore the integral term. Call once from UDF_Setup, after the mesh
 * is available.
 */
void pidSetup()
{
    mesh_t* mesh = nrs->meshV;

    pid = {};
    platform->par->extract("pid", "kp", pid.kp);
    platform->par->extract("pid", "ki", pid.ki);
    platform->par->extract("pid", "kd", pid.kd);

    // Domain center = integral(x dV)/V (exact for the box mesh).
    const occa::memory o_xyz[3] = {mesh->o_x, mesh->o_y, mesh->o_z};
    for (int d = 0; d < 3; d++) {
        pid.center[d] = platform->linAlg->innerProd(mesh->Nlocal, o_xyz[d], mesh->o_Jw,
                platform->comm.mpiComm()) / mesh->volume;
    }

    if (platform->comm.mpiRank() == 0) {
        printf("PID: kp=%g ki=%g kd=%g center=(%g, %g, %g)\n", pid.kp, pid.ki, pid.kd,
                pid.center[0], pid.center[1], pid.center[2]);
        if ((pid.kp != 0 || pid.ki != 0 || pid.kd < 0) &&
                !(pid.kd > 0 && pid.ki >= 0 && pid.kd*pid.kp > pid.ki)) {
            printf("PID: WARNING gains do not satisfy kd > 0, ki >= 0 and kd*kp > ki, the loop is unstable\n");
        }
    }

    if (!platform->options.getArgs("RESTART FILE NAME").empty()) {
        double restart_time = 0;
        platform->options.getArgs("START TIME", restart_time);
        pidRestoreIntegral(restart_time);
    }
}

/**
 * Measure the gas phase volume, centroid offset, and mean velocity, and the
 * liquid phase mean velocity.
 */
void pidMeasure()
{
    mesh_t* mesh = nrs->meshV;
    fluidSolver_t* fluid = nrs->fluid.get();
    const occa::memory o_psi = nrs->scalar->o_solution("cls");
    const occa::memory o_xyz[3] = {mesh->o_x, mesh->o_y, mesh->o_z};
    const occa::memory o_uvw[3] = {fluid->o_solution("x"), fluid->o_solution("y"), fluid->o_solution("z")};
    MPI_Comm comm = platform->comm.mpiComm();

    // Gas phase quadrature weights o_gasJw = (1-psi)*Jw, so that
    // innerProd(f, o_gasJw) = integral(f*(1-psi) dV).
    auto o_gasJw = platform->deviceMemoryPool.reserve<dfloat>(mesh->Nlocal);
    platform->linAlg->axmyz(mesh->Nlocal, 1.0, o_psi, mesh->o_Jw, o_gasJw);
    platform->linAlg->axpby(mesh->Nlocal, 1.0, mesh->o_Jw, -1.0, o_gasJw);

    pid.gas_volume = platform->linAlg->sum(mesh->Nlocal, o_gasJw, comm);
    const dfloat liquid_volume = mesh->volume - pid.gas_volume;
    for (int d = 0; d < 3; d++) {
        pid.error[d] = platform->linAlg->innerProd(mesh->Nlocal, o_xyz[d], o_gasJw, comm)
                / pid.gas_volume - pid.center[d];
        const dfloat gas_momentum = platform->linAlg->innerProd(mesh->Nlocal, o_uvw[d], o_gasJw, comm);
        const dfloat total_momentum = platform->linAlg->innerProd(mesh->Nlocal, o_uvw[d], mesh->o_Jw, comm);
        pid.u_gas[d] = gas_momentum / pid.gas_volume;
        pid.u_liquid[d] = (total_momentum - gas_momentum) / liquid_volume;
    }
}

/**
 * Update the controller from the current solution at the end of a time step,
 * and compute the PID force for the next time step. Call once per step from
 * UDF_ExecuteStep, after lvlSet::solve (including at tstep=0 during setup).
 *
 * @param time simulation time of the current solution
 * @param tstep time step index (0 during setup)
 */
void pidUpdate(double time, int tstep)
{
    pidMeasure();

    if (tstep > 0) {
        const double dt = time - pid.time;
        for (int d = 0; d < 3; d++) {
            pid.integral[d] += pid.error[d] * dt;
        }
    }
    pid.time = time;

    for (int d = 0; d < 3; d++) {
        pid.force[d] = -(pid.kp*pid.error[d] + pid.ki*pid.integral[d] + pid.kd*pid.u_gas[d]);
    }

    if (platform->comm.mpiRank() == 0) {
        printf("PID: step=%d t=%.8e e=(%+.4e, %+.4e, %+.4e) u_gas=(%+.4e, %+.4e, %+.4e) "
                "F=(%+.4e, %+.4e, %+.4e) Vgas=%.8e\n", tstep, time,
                pid.error[0], pid.error[1], pid.error[2],
                pid.u_gas[0], pid.u_gas[1], pid.u_gas[2],
                pid.force[0], pid.force[1], pid.force[2], pid.gas_volume);
    }
}

/**
 * Add the PID force (per unit mass) to the fluid explicit terms.
 *
 * @param o_uSource fluid explicit terms (all 3 components), which must already
 *     contain the surface tension acceleration as that overwrites the buffer
 */
void pidApply(occa::memory& o_uSource)
{
    for (int d = 0; d < 3; d++) {
        platform->linAlg->add(nrs->meshV->Nlocal, pid.force[d], o_uSource, d*nrs->fieldOffset);
    }
}
