// Species: c, at the Schmidt number of case.hpp, and any further scalars
// named c<Sc> in the .par (e.g. scalars = tls, cls, c, c20, c1150): passive
// scalars at the Schmidt number in their name, carried by the same flow (c
// is passive, so the flow does not depend on Sc). Each has its own
// diffusivity 1/Pe, hard sink, species budget and output: data.csv for c,
// data-c<Sc>.csv with the same columns for the others (see bubble.udf).
struct species_t {
    std::string name;            // scalar name
    double Pe;                   // Re*Sc
    std::string csv;             // output file
    // Species budget since the last sample (see UDF_ExecuteStep).
    double budget_inflow = 0.0;  // time integral of the net boundary inflow
    double budget_C0 = 0.0;      // integral of c dV at the window start
    double budget_flux = 0.0;    // net boundary inflow rate at the last step
    double sink_removed = 0.0;   // c zeroed by the hard sink in the window
};
static std::vector<species_t> species;

void customProperties(double t)
{
    mesh_t* mesh = nrs->meshV;
    fluidSolver_t* fluid = nrs->fluid.get();
    scalar_t* scalar = nrs->scalar.get();

    // Properties for momentum transport equation
    const occa::memory o_psi = scalar->o_solution("cls");
    // mu = ((1.0-psi)/muratio + psi)/Re
    weightedMixing(info, o_psi, 1.0/muratio, Re, fluid->o_diffusionCoeff());
    // rho = (1.0-psi)/rhoratio + psi
    weightedMixing(info, o_psi, 1.0/rhoratio, 1.0, fluid->o_transportCoeff());

    // Properties for the species transport equations: CST is off, so each
    // species is a passive scalar with the liquid diffusivity 1/Pe everywhere
    // (the hard sink in UDF_ExecuteStep removes it in the gas).
    for (const auto& s : species) {
        platform->linAlg->fill(mesh->Nlocal, 1.0/s.Pe, scalar->o_diffusionCoeff(s.name));
        // No transport coefficient in the species equations.
        platform->linAlg->fill(mesh->Nlocal, 1.0, scalar->o_transportCoeff(s.name));
    }

    // nekRS-LS only computes the scalar SVV viscosity when userProperties is
    // unset (nrs_t::evaluateProperties), so compute it here or the [SCALAR *]
    // svv regularization has no effect (gpu/mesh-sensitivity does not, so its
    // scalars ran without SVV).
    scalar->mueSVV();
}

void customSource(double t)
{
    scalar_t* scalar = nrs->scalar.get();
    fluidSolver_t* fluid = nrs->fluid.get();
    const occa::memory o_psi = scalar->o_solution("cls");
    const occa::memory o_rho = fluid->o_transportCoeff();
    occa::memory o_uSource = fluid->o_explicitTerms();
    occa::memory o_uSourceY = o_uSource.slice(1*nrs->fieldOffset, nrs->fieldOffset);

    // Surface tension source term for the U equation.
    lvlSet::applySurfaceTensionAcc(We, o_uSource);

    // Buoyancy source terms for the U equation.
    buoyancySource(info, o_psi, o_rho, Fr, o_uSourceY);

    // PID bubble centering force (uniform acceleration) for the U equation,
    // and the matching inflow velocity for this step.
    pidApply(o_uSource, nrs->dt[0]);
}
