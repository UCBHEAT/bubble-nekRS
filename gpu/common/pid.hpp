// PID controller that keeps the bubble (gas-phase centroid) at a fixed point in
// the mesh, by accelerating the reference frame (moving reference frame).
//
// Shared by the gpu/ cases (gpu/pid-centering, gpu/pid-centering-3D): include
// it from the case .udf as "../common/pid.hpp". A case using it needs
//   - the conservative level set scalar "cls" (psi = 1 in the liquid),
//   - a [PID] par section (kp, ki, kd, and optionally the setpoint x0, y0, z0),
//   - pidSetup() in UDF_Setup, pidUpdate(time, tstep) in UDF_ExecuteStep after
//     the level set solve, and pidApply(o_uSource, dt) in the user source,
//   - its inflow Dirichlet condition to take the velocity from bc->usrwrk[0..2]
//     (the frame velocity, see pidApply), and
//   - for restarts, the pid_int_* and frame_* columns in data.csv.
//
// The controller output F_pid is a force per unit mass (an acceleration, in the
// same nondimensional units as the -1/Fr^2 gravity term in buoyancySource) that
// is added to every velocity point: the fictitious force of a frame
// accelerating at -F_pid. The liquid far from the bubble is at rest in the lab
// frame, so in this frame it moves with the frame velocity
//
//   U_frame(t) = integral(F_pid dt),
//
// which is imposed as the inflow velocity. The uniform acceleration then shifts
// the whole velocity field uniformly and does not change the bubble/liquid
// relative motion, up to O(dt) differences between the forward-Euler inflow
// update and the BDF/EXT integration of the forcing, which a uniform pressure
// gradient absorbs (~0.1% of gravity on average here). Once the bubble reaches
// terminal velocity V, U_frame = -V: the bubble stays put and the liquid flows
// down past it.
//
//   e      = x_c - x_0          gas centroid offset from the setpoint x_0
//   u_gas                       gas-phase mean velocity, used as de/dt
//                               (derivative on measurement)
//   F_pid  = -(kp*e + ki*integral(e dt) + kd*u_gas)
//
// u_gas is the material velocity of the gas, which is exactly de/dt between
// level set reinitializations. Each CLSR reinit also shifts the centroid
// slightly (upward for a rising bubble, by ~3e-4 per reinit with dt = 0.0025
// and reinit every 10 steps, a drift of ~0.012-0.017 per unit time), which
// u_gas does not see. At steady state the integral term therefore settles at
// ki*integral = kd*(reinit drift rate) rather than 0.
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
    // Setpoint x_0 ([PID] x0, y0, z0; default: volume centroid of the mesh).
    dfloat center[3];
    // Gas phase volume and centroid offset from the setpoint (the error e).
    dfloat gas_volume;
    dfloat error[3];
    // Time integral of the error. Written to data.csv so it can be restored on
    // restart.
    dfloat integral[3];
    // Gas and liquid phase mean velocities; u_gas is used as de/dt.
    dfloat u_gas[3];
    dfloat u_liquid[3];
    // PID force per unit mass, applied uniformly in the next time step.
    dfloat force[3];
    // Frame velocity U_frame = integral(F_pid dt) at the last update: the
    // velocity of the quiescent far-field liquid in this frame. Written to
    // data.csv so it can be restored on restart.
    dfloat frame_velocity[3];
    // Simulation time of the last update, for integrating the error.
    double time;
} pidState_t;

static pidState_t pid;

/**
 * Restore the integral term and frame velocity on restart from the data.csv
 * row written at the restart time (the rest of the controller state is
 * re-measured from the restart fields). If there is no such row, the integral
 * starts from zero and the frame velocity is taken from the restart field's
 * (uniform) inflow velocity.
 *
 * @param restart_time simulation time of the restart file
 */
void pidRestoreState(double restart_time)
{
    // data.csv rows are written with 4 decimals of time, so accept the row
    // closest to the restart time within this tolerance.
    const double tolerance = 1e-3;
    const std::string names[6] = {"pid_int_x", "pid_int_y", "pid_int_z",
                                  "frame_u", "frame_v", "frame_w"};
    double state[7] = {0, 0, 0, 0, 0, 0, 0}; // found, integral x/y/z, frame u/v/w

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
                // a data.csv started with an older column layout is never
                // restored from: remove data.csv when the column layout changes.
                header = cols;
                continue;
            }
            const auto column = [&](const std::string& name) {
                const auto it = std::find(header.begin(), header.end(), name);
                return (it == header.end()) ? -1 : int(it - header.begin());
            };
            const int itime = column("Time");
            int index[6];
            bool complete = (itime >= 0);
            for (int k = 0; k < 6; k++) {
                index[k] = column(names[k]);
                complete = complete && (index[k] >= 0);
            }
            if (!complete || int(cols.size()) != int(header.size())) continue;
            try {
                const double t = std::stod(cols[itime]);
                // Keep the last row within the tolerance, so a row from a later
                // rerun supersedes an earlier one at the same time.
                if (std::abs(t - restart_time) <= best) {
                    double values[6];
                    for (int k = 0; k < 6; k++) values[k] = std::stod(cols[index[k]]);
                    best = std::max(std::abs(t - restart_time), 1e-12);
                    state[0] = 1;
                    for (int k = 0; k < 6; k++) state[1 + k] = values[k];
                }
            } catch (const std::exception&) {
                // Skip malformed rows (e.g. a partially written last line).
            }
        }
    }
    MPI_Bcast(state, 7, MPI_DOUBLE, 0, platform->comm.mpiComm());

    if (state[0] == 0) {
        // No data.csv row: take U_frame from the restart field instead, since it
        // also sets the inflow velocity (resetting it to 0 would stop the liquid
        // in one step). The inflow face (boundary ID 1, y = ymax) is uniform and
        // holds U_frame, and no Dirichlet condition has been applied yet.
        mesh_t* mesh = nrs->meshV;
        fluidSolver_t* fluid = nrs->fluid.get();
        MPI_Comm comm = platform->comm.mpiComm();
        auto [x, y, z] = mesh->xyzHost();
        double ymax = -1e300;
        for (dlong n = 0; n < mesh->Nlocal; n++) ymax = std::max(ymax, double(y[n]));
        MPI_Allreduce(MPI_IN_PLACE, &ymax, 1, MPI_DOUBLE, MPI_MAX, comm);
        double sum[4] = {0, 0, 0, 0}; // u, v, w, node count
        const std::string comp[3] = {"x", "y", "z"};
        std::vector<dfloat> u(mesh->Nlocal);
        for (int d = 0; d < 3; d++) {
            fluid->o_solution(comp[d]).copyTo(u, mesh->Nlocal);
            for (dlong n = 0; n < mesh->Nlocal; n++) {
                if (std::abs(y[n] - ymax) < 1e-6) {
                    sum[d] += u[n];
                    if (d == 0) sum[3] += 1;
                }
            }
        }
        MPI_Allreduce(MPI_IN_PLACE, sum, 4, MPI_DOUBLE, MPI_SUM, comm);
        for (int d = 0; d < 3; d++) state[4 + d] = sum[d] / sum[3];
    }

    for (int d = 0; d < 3; d++) {
        pid.integral[d] = state[1 + d];
        pid.frame_velocity[d] = state[4 + d];
    }
    if (platform->comm.mpiRank() == 0) {
        if (state[0] > 0) {
            printf("PID: restored integral=(%g, %g, %g) frame velocity=(%g, %g, %g) from data.csv at t=%g\n",
                    pid.integral[0], pid.integral[1], pid.integral[2],
                    pid.frame_velocity[0], pid.frame_velocity[1], pid.frame_velocity[2], restart_time);
        } else {
            printf("PID: WARNING no data.csv row at restart time t=%g, integral starts from 0, "
                    "frame velocity (%g, %g, %g) taken from the restart field inflow\n", restart_time,
                    pid.frame_velocity[0], pid.frame_velocity[1], pid.frame_velocity[2]);
        }
    }
}

/**
 * Read gains and setpoint from the [PID] par section, and on restart restore
 * the integral term and frame velocity. Call once from UDF_Setup, after the
 * mesh is available.
 */
void pidSetup()
{
    mesh_t* mesh = nrs->meshV;

    pid = {};
    platform->par->extract("pid", "kp", pid.kp);
    platform->par->extract("pid", "ki", pid.ki);
    platform->par->extract("pid", "kd", pid.kd);

    // Setpoint: the volume centroid of the mesh, unless given in [PID].
    const occa::memory o_xyz[3] = {mesh->o_x, mesh->o_y, mesh->o_z};
    const std::string keys[3] = {"x0", "y0", "z0"};
    for (int d = 0; d < 3; d++) {
        pid.center[d] = platform->linAlg->innerProd(mesh->Nlocal, o_xyz[d], mesh->o_Jw,
                platform->comm.mpiComm()) / mesh->volume;
        platform->par->extract("pid", keys[d], pid.center[d]);
    }

    if (platform->comm.mpiRank() == 0) {
        printf("PID: kp=%g ki=%g kd=%g setpoint=(%g, %g, %g)\n", pid.kp, pid.ki, pid.kd,
                pid.center[0], pid.center[1], pid.center[2]);
        if ((pid.kp != 0 || pid.ki != 0 || pid.kd < 0) &&
                !(pid.kd > 0 && pid.ki >= 0 && pid.kd*pid.kp > pid.ki)) {
            printf("PID: WARNING gains do not satisfy kd > 0, ki >= 0 and kd*kp > ki, the loop is unstable\n");
        }
    }

    if (!platform->options.getArgs("RESTART FILE NAME").empty()) {
        double restart_time = 0;
        platform->options.getArgs("START TIME", restart_time);
        pidRestoreState(restart_time);
    }

    // Device copy of the inflow velocity, read by udfDirichlet as bc->usrwrk[0..2].
    platform->app->bc->o_usrwrk.resize(3);
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
            // pid.force still holds the force applied during this step.
            pid.frame_velocity[d] += pid.force[d] * dt;
        }
    }
    pid.time = time;

    for (int d = 0; d < 3; d++) {
        pid.force[d] = -(pid.kp*pid.error[d] + pid.ki*pid.integral[d] + pid.kd*pid.u_gas[d]);
    }

    if (platform->comm.mpiRank() == 0) {
        printf("PID: step=%d t=%.8e e=(%+.4e, %+.4e, %+.4e) u_gas=(%+.4e, %+.4e, %+.4e) "
                "F=(%+.4e, %+.4e, %+.4e) Vgas=%.8e U_frame=(%+.4e, %+.4e, %+.4e)\n", tstep, time,
                pid.error[0], pid.error[1], pid.error[2],
                pid.u_gas[0], pid.u_gas[1], pid.u_gas[2],
                pid.force[0], pid.force[1], pid.force[2], pid.gas_volume,
                pid.frame_velocity[0], pid.frame_velocity[1], pid.frame_velocity[2]);
    }
}

/**
 * Add the PID force (per unit mass) to the fluid explicit terms, and set the
 * inflow velocity for the step, U_frame at the end of the step. Call from the
 * userSource hook, which runs before the step's Dirichlet conditions.
 *
 * @param o_uSource fluid explicit terms (all 3 components), which must already
 *     contain the surface tension acceleration as that overwrites the buffer
 * @param dt time step size of the step being set up
 */
void pidApply(occa::memory& o_uSource, dfloat dt)
{
    for (int d = 0; d < 3; d++) {
        platform->linAlg->add(nrs->meshV->Nlocal, pid.force[d], o_uSource, d*nrs->fieldOffset);
    }

    dfloat inflow[3];
    for (int d = 0; d < 3; d++) inflow[d] = pid.frame_velocity[d] + pid.force[d]*dt;
    platform->app->bc->o_usrwrk.copyFrom(inflow, 3);
}
