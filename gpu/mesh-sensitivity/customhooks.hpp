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

    // Use D1 everywhere.
    platform->linAlg->fill(mesh->Nlocal, 1.0/Pe, scalar->o_diffusionCoeff("c"));
    // No transport coefficient in c equation.
    platform->linAlg->fill(mesh->Nlocal, 1.0, scalar->o_transportCoeff("c"));
}

void customSource(double t)
{
    mesh_t* mesh = nrs->meshV;
    scalar_t* scalar = nrs->scalar.get();
    fluidSolver_t* fluid = nrs->fluid.get();
    const occa::memory o_psi = scalar->o_solution("cls");
    const occa::memory o_phi = scalar->o_solution("tls");
    const occa::memory o_c = scalar->o_solution("c");
    occa::memory o_cSource = scalar->o_explicitTerms("c");
    const occa::memory o_rho = fluid->o_transportCoeff();
    occa::memory o_uSource = fluid->o_explicitTerms();
    occa::memory o_uSourceY = o_uSource.slice(1*nrs->fieldOffset, nrs->fieldOffset);

    // Soft source term for the C equation.
    softSource(info, o_psi, o_c, o_cSource);

    // Surface tension source term for the U equation.
    lvlSet::applySurfaceTensionAcc(We, o_uSource);

    // Buoyancy source terms for the U equation.
    buoyancySource(info, o_psi, o_rho, Fr, o_uSourceY);
}
