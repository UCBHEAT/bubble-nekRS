void customProperties(double t)
{
    mesh_t* mesh = nrs->meshV;
    fluidSolver_t* fluid = nrs->fluid.get();
    scalar_t* scalar = nrs->scalar.get();

    // Properties for momentum transport equation
    const occa::memory o_psi = scalar->o_solution("psi");
    // mu = ((1.0-psi)*muratio + psi)/Re
    weightedMixing(mesh->Nelements, o_psi, muratio, Re, fluid->o_diffusionCoeff());
    // rho = (1.0-psi)*rhoratio + psi
    weightedMixing(mesh->Nelements, o_psi, rhoratio, 1.0, fluid->o_transportCoeff());

    // Properties for species transport equation
    // D = (psi/(psi+(1.0-psi)*solubilityratio) +
    //     diffratio*(1.0-psi)/((1.0-psi)+psi/solubilityratio))/Pe
    speciesDiff(mesh->Nelements, o_psi, solubilityratio, diffratio, Pe,
            scalar->o_diffusionCoeff("c"));
    // No transport coefficient in c equation.
    constantFill(mesh->Nelements, 1.0, scalar->o_transportCoeff("c"));
}

void customSource(double t)
{
    mesh_t* mesh = nrs->meshV;
    scalar_t* scalar = nrs->scalar.get();
    const occa::memory o_psi = scalar->o_solution("psi");
    const occa::memory o_phi = scalar->o_solution("phi");
    const occa::memory o_c = scalar->o_solution("c");
    occa::memory o_cstVector = platform->device.malloc<dfloat>(3*(nrs->fieldOffset));
    occa::memory o_cstVectorX = o_cstVector.slice(0*nrs->fieldOffset, nrs->fieldOffset);
    occa::memory o_cstVectorY = o_cstVector.slice(1*nrs->fieldOffset, nrs->fieldOffset);
    occa::memory o_cstVectorZ = o_cstVector.slice(2*nrs->fieldOffset, nrs->fieldOffset);

    // Calculate interface unit normals.
    opSEM::strongGrad(mesh, nrs->fieldOffset, o_phi, o_cstVector);
    interfaceUnitNormals(mesh->Nelements, o_phi, o_cstVectorX, o_cstVectorY, o_cstVectorZ);
    // // Workaround for test purposes -- set interface normal vector field to constant (0, 1, 0) vector field.
    // constantFill(mesh->Nelements, 0.0, o_cstVectorX);
    // constantFill(mesh->Nelements, 1.0, o_cstVectorY);
    // constantFill(mesh->Nelements, 0.0, o_cstVectorZ);

    // Calculate CST vector field.
    speciesSource(mesh->Nelements, o_c, o_psi, solubilityratio, diffratio, Pe,
        o_cstVectorX, o_cstVectorY, o_cstVectorZ);

    // Source term is the divergence of the above vector field.
    opSEM::strongDivergence(mesh, nrs->fieldOffset, o_cstVector, scalar->o_explicitTerms("c"));
}
