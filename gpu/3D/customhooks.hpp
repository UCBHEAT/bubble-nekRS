void customProperties(double t)
{
    mesh_t* mesh = nrs->meshV;
    fluidSolver_t* fluid = nrs->fluid.get();
    scalar_t* scalar = nrs->scalar.get();

    // Properties for momentum transport equation
    const occa::memory o_psi = scalar->o_solution("cls");
    // mu = ((1.0-psi)*muratio + psi)/Re
    weightedMixing(info, o_psi, muratio, Re, fluid->o_diffusionCoeff());
    // rho = (1.0-psi)*rhoratio + psi
    weightedMixing(info, o_psi, rhoratio, 1.0, fluid->o_transportCoeff());

    // Properties for species transport equation
    // D = (psi/(psi+(1.0-psi)*solubilityratio) +
    //     diffratio*(1.0-psi)/((1.0-psi)+psi/solubilityratio))/Pe
    speciesDiff(info, o_psi, solubilityratio, diffratio, Pe,
            scalar->o_diffusionCoeff("c"));
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
    occa::memory o_cstVector = platform->device.malloc<dfloat>(3*(nrs->fieldOffset));
    occa::memory o_cstVectorX = o_cstVector.slice(0*nrs->fieldOffset, nrs->fieldOffset);
    occa::memory o_cstVectorY = o_cstVector.slice(1*nrs->fieldOffset, nrs->fieldOffset);
    occa::memory o_cstVectorZ = o_cstVector.slice(2*nrs->fieldOffset, nrs->fieldOffset);
    const occa::memory o_rho = fluid->o_transportCoeff();
    occa::memory o_uSource = fluid->o_explicitTerms();
    occa::memory o_uSourceY = o_uSource.slice(1*nrs->fieldOffset, nrs->fieldOffset);

    // Calculate interface unit normals.
    opSEM::strongGrad(mesh, nrs->fieldOffset, o_phi, o_cstVector);
    interfaceNormals(info, o_phi, lvlSet::getDeltaFunction(), o_cstVectorX, o_cstVectorY, o_cstVectorZ);

    // Calculate CST vector field.
    speciesSource(info, o_c, o_psi, solubilityratio, diffratio, Pe,
        o_cstVectorX, o_cstVectorY, o_cstVectorZ);

    // Source term is the divergence of the above vector field.
    opSEM::strongDivergence(mesh, nrs->fieldOffset, o_cstVector, o_cSource);

    // Surface tension source term for the U equation.
    lvlSet::applySurfaceTensionAcc(We, o_uSource);

    // Buoyancy source terms for the U equation.
    buoyancySource(info, o_psi, o_rho, Fr, o_uSourceY);
}
