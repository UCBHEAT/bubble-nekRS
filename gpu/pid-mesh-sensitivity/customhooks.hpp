// [CASEDATA] scalarSVV = false skips scalar->mueSVV() below (for tests).
static bool scalarSVV = true;

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

    // Properties for species transport equation: CST is off, so c is a passive
    // scalar with the liquid diffusivity 1/Pe everywhere (the hard sink in
    // UDF_ExecuteStep removes it in the gas).
    platform->linAlg->fill(mesh->Nlocal, 1.0/Pe, scalar->o_diffusionCoeff("c"));
    // No transport coefficient in c equation.
    platform->linAlg->fill(mesh->Nlocal, 1.0, scalar->o_transportCoeff("c"));

    // nekRS-LS only computes the scalar SVV viscosity when userProperties is
    // unset (nrs_t::evaluateProperties), so compute it here or the [SCALAR *]
    // svv regularization has no effect (gpu/mesh-sensitivity does not, so its
    // scalars ran without SVV).
    if (scalarSVV) {
        scalar->mueSVV();
    }
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
