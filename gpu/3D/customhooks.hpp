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

void myGetCurvature(const occa::memory& o_normals, occa::memory& o_curvature)
{
  auto mesh = nrs->scalar->mesh(nrs->scalar->nameToIndex.find("tls")->second);

  opSEM::strongDivergence(mesh, nrs->scalar->fieldOffset(), o_normals, o_curvature);
  // seems to help, but not by reducing boxiness... it reduces magnitude of the result??
  //opSEM::divergence(mesh, nrs->scalar->fieldOffset(), o_normals, o_curvature);
}

void myApplySurfaceTensionAcc(const dfloat& We, occa::memory &o_sforce)
{
    scalar_t* scalar = nrs->scalar.get();
    auto meshV = nrs->scalar->meshV;

    auto o_delta = lvlSet::getDeltaFunction();
    // this looks good
    scalar->o_solution("debug1").copyFrom(o_delta, o_delta.size());

    auto o_phi = nrs->scalar->o_solution("tls");
    bool avg = false;
    if(platform->options.compareArgs("LVLSET NORMAL AVERAGING", "TRUE")) {
        avg = true;
    }
    lvlSet::normalVector(o_phi, o_sforce, avg);
    // this is boxy (nearly locally constant on each element), leading the gradient
    // to be concentrated on element edges
    //scalar->o_solution("debug1").copyFrom(o_sforce, nrs->fieldOffset);
    //scalar->o_solution("debug2").copyFrom(o_sforce.slice(1*nrs->fieldOffset, nrs->fieldOffset), nrs->fieldOffset);
    //scalar->o_solution("debug3").copyFrom(o_sforce.slice(1*nrs->fieldOffset, nrs->fieldOffset), nrs->fieldOffset);

    auto o_curvDeltabyRho = platform->device.malloc<dfloat>(meshV->Nlocal);
    myGetCurvature(o_sforce, o_curvDeltabyRho);
    // Laplacian is a bit cleaner than double derivative. We ignore the normalization here because
    // mag(grad(phi)) is supposed to be 1 anyways. Turn off averaging.
    //opSEM::strongLaplacian(meshV, nrs->scalar->fieldOffset(), o_phi, o_curvDeltabyRho, false);
    // didn't seem to do anything, with or without averaging
    //scalar->o_solution("debug1").copyFrom(o_curvDeltabyRho, o_curvDeltabyRho.size()); // curvature
    platform->linAlg->axmy(meshV->Nlocal, 1.0, o_delta, o_curvDeltabyRho);
    scalar->o_solution("debug2").copyFrom(o_curvDeltabyRho, o_curvDeltabyRho.size()); // curvature*area

    // Divide by density
    auto o_rho = nrs->fluid->o_prop + 1 * nrs->fluid->fieldOffset;
    platform->linAlg->aydx(meshV->Nlocal, 1.0, o_rho, o_curvDeltabyRho);
    //scalar->o_solution("debug3").copyFrom(o_curvDeltabyRho, o_curvDeltabyRho.size()); // curvature*area/rho

    // There should be no curvature inside the bubble. Currently rho_g being << rho_l
    // amplifies the high noise in o_curvature inside the bubble. Multiplying by delta
    // (surface area density, but also a Dirac delta for the interface location) does
    // not do a good enough job at killing this massive far-from-interface curvature,
    // so we do our own custom cleanup.
    auto o_psi = nrs->scalar->o_solution("cls");
    cleanupCurvature(info, o_psi, o_curvDeltabyRho);
    scalar->o_solution("debug3").copyFrom(o_curvDeltabyRho, o_curvDeltabyRho.size()); // cleanup(curvature*area/rho)

    platform->linAlg->axmyVector(meshV->Nlocal,
                                nrs->scalar->vFieldOffset,
                                0,
                                -1.0/We, //reverse sign (see Nek5000)
                                o_curvDeltabyRho,
                                o_sforce);
    //scalar->o_solution("debug1").copyFrom(o_sforce, nrs->fieldOffset);
    //scalar->o_solution("debug2").copyFrom(o_sforce.slice(1*nrs->fieldOffset, nrs->fieldOffset), nrs->fieldOffset);
    //scalar->o_solution("debug3").copyFrom(o_sforce.slice(1*nrs->fieldOffset, nrs->fieldOffset), nrs->fieldOffset);
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
    //scalar->o_solution("debug1").copyFrom(o_cstVectorY, o_cstVectorY.size());

    // Calculate CST vector field.
    speciesSource(info, o_c, o_psi, solubilityratio, diffratio, Pe,
        o_cstVectorX, o_cstVectorY, o_cstVectorZ);
    //scalar->o_solution("debug2").copyFrom(o_cstVectorY, o_cstVectorY.size());

    // Source term is the divergence of the above vector field.
    opSEM::strongDivergence(mesh, nrs->fieldOffset, o_cstVector, o_cSource);
    //scalar->o_solution("debug3").copyFrom(o_cSource, o_cSource.size());

    // Surface tension source term for the U equation.
    //lvlSet::applySurfaceTensionAcc(We, o_uSource);
    myApplySurfaceTensionAcc(We, o_uSource);

    // Buoyancy source terms for the U equation.
    buoyancySource(info, o_psi, o_rho, Fr, o_uSourceY);
}
