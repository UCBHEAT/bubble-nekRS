// PID controller that keeps the bubble (gas-phase centroid) at the domain
// center, by applying a spatially uniform acceleration to the whole domain.
//
// The controller output F_pid is a force per unit mass (an acceleration, in the
// same nondimensional units as the -1/Fr^2 gravity term in buoyancySource) that
// is added to every velocity point. In a fully periodic domain a uniform
// acceleration is exactly an acceleration of the reference frame: it shifts the
// whole velocity field uniformly and does not change the bubble/liquid relative
// motion (the net gravity + pressure gradient force on the domain is zero, so
// F_pid is the only term that changes the domain's mean momentum). Once the
// bubble reaches terminal velocity, the frame rises with the bubble and the
// liquid flows down past it.
//
//   e      = x_c - x_0          gas centroid offset from the domain center
//   de/dt  = u_gas              gas-phase mean velocity (derivative on
//                               measurement: avoids differencing x_c, which
//                               jumps slightly at each level set reinit)
//   F_pid  = -(kp*e + ki*integral(e dt) + kd*de/dt)
//
// Gains are read from the [PID] section of the .par file. The plant (F_pid ->
// x_c) is a double integrator, so kd > 0 is required for stability; with all
// three gains, the closed loop is stable for kd*kp > ki.

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
    // Time integral of the error.
    dfloat integral[3];
    // Gas and liquid phase mean velocities. u_gas is de/dt; u_gas - u_liquid
    // is the bubble slip (rise) velocity.
    dfloat u_gas[3];
    dfloat u_liquid[3];
    // PID force per unit mass, applied uniformly in the next time step.
    dfloat force[3];
    // Simulation time of the last update, for integrating the error.
    double time;
} pidState_t;

static pidState_t pid;

/**
 * Read gains from the [PID] par section and compute the domain center.
 * Call once from UDF_Setup, after the mesh is available.
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
        printf("PID: step=%d t=%.8e e=(%+.4e, %+.4e, %+.4e) de/dt=(%+.4e, %+.4e, %+.4e) "
                "F=(%+.4e, %+.4e, %+.4e)\n", tstep, time,
                pid.error[0], pid.error[1], pid.error[2],
                pid.u_gas[0], pid.u_gas[1], pid.u_gas[2],
                pid.force[0], pid.force[1], pid.force[2]);
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
