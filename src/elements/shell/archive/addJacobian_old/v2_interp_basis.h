template <class quadrature, class basis, class director, class model>
void TACSShellElement<quadrature, basis, director, model>::addJacobian(
    int elemIndex, double time, TacsScalar alpha, TacsScalar beta,
    TacsScalar gamma, const TacsScalar Xpts[], const TacsScalar vars[],
    const TacsScalar dvars[], const TacsScalar ddvars[], TacsScalar res[],
    TacsScalar mat[]) {
  // Compute the number of quadrature points
  const int nquad = quadrature::getNumQuadraturePoints();

  printf("Begin addJacobian.. \n");

  // midplane displacements & time derivatives
  A2D::Mat<TacsScalar,3,num_nodes> u_temp, udot_temp, uddot_temp;
  for (int inode = 0; inode < num_nodes; inode++) {
    for (int ivar = 0; ivar < 3; ivar++) {
      u_temp(ivar,inode) = vars[vars_per_node * inode + ivar];
      udot_temp(ivar,inode) = dvars[vars_per_node * inode + ivar];
      uddot_temp(ivar,inode) = ddvars[vars_per_node * inode + ivar];
    }
  }
  A2D::A2DObj<A2D::Mat<TacsScalar,3,num_nodes>> u(u_temp), udot(udot_temp), uddot(uddot_temp);
  delete u_temp; delete udot_temp; delete uddot_temp;

  // Compute the node normal directions
  A2D::Mat<TacsScalar, 3, num_nodes> fn, Xcoords(Xpts);
  A2D::Mat<TacsScalar, 9, num_nodes> Xdn;
  // TODO : add A2D later for this guy
  TacsShellComputeNodeNormals<basis>(Xcoords.get_data(), fn.get_data(), Xdn.get_data());

  // Store information about the transformation and derivatives at each node for
  // the drilling degrees of freedom
  A2D::A2DObj<Mat<TacsScalar, 1, num_nodes>> etn;
  A2D::Mat<TacsScalar, 9, num_nodes> XdinvTn, Tn, u0xn, Ctn;
  TacsShellComputeDrillStrain<vars_per_node, offset, basis, director, model>(
      transform, Xdn.get_data(), fn.get_data(), vars, XdinvTn.get_data(), Tn.get_data(), u0xn.get_data(), Ctn.get_data(), 
      etn.get_data());

  A2D::A2DObj<A2D::Mat<TacsScalar,3,num_nodes>> d, ddot, dddot;
  director::template computeDirectorRates<vars_per_node, offset, num_nodes>(
      vars, dvars, ddvars, fn, d.get_data(), ddot.get_data(), dddot.get_data());

  // Set the total number of tying points needed for this element
  int num_tying = basis::NUM_TYING_POINTS;
  A2D::A2DObj<Mat<TacsScalar, num_tying>> ety;
  A2D::Mat<TacsScalar, num_tying, num_tying> d2ety;
  A2D::Mat<TacsScalar, num_tying, usize> d2etyu;
  A2D::Mat<TacsScalar, num_tying, dsize> d2etyd;

  model::template computeTyingStrain<vars_per_node, basis>(Xpts, fn.get_data(), vars, d.get_data(), ety.get_data());

  // Loop over each quadrature point and add the residual contribution
  for (int quad_index = 0; quad_index < nquad; quad_index++) {
    
    // Get the quadrature weight
    double pt[3];
    double weight = quadrature::getQuadraturePoint(quad_index, pt);

    // printf("iquad = %d\n", quad_index);

    // interpolation section
    // ----------------------------------------

    // passive A2D Objs used in interpolation
    A2D::Vec<TacsScalar,3> X, n0;
    A2D::Mat<TacsScalar,3,2> Xxi, nxi;
    A2D::Mat<TacsScalar,3,3> T;
    // passive interpolation vectors
    A2D::Vec<TacsScalar, num_nodes> interp;
    A2D::Mat<TacsScalar, num_nodes, 2> interpGrad;
    
    // active A2D objs used in interpolation
    A2D::A2DObj<A2D::Vec<TacsScalar,1>> et;
    A2D::A2DObj<A2D::Vec<TacsScalar,3>> d0;
    A2D::A2DObj<A2D::Mat<TacsScalar,3,2>> d0xi, u0xi;
    A2D::A2DObj<A2D::SymMat<TacsScalar,3>> e0ty, e0ty_tmp, gty;    
    
    // get the constant interpolation vectors
    basis::template getInterpFields(pt, interp.get_data()); // TODO : make routines for these
    basis::template getInterpFieldsGrad(pt, interp3Grad.get_data());
    basis::template getInterpTyingStrain


    auto interp_stack = A2D::MakeStack(
      MatVecMult(etn, interp1, et),
      MatVecMult(Xpts, interp3, X),
      MatVecMult(fn, interp3, n0),
      MatVecMult(d, intepr3, d0),
      MatMatMult(Xpts, interp3Grad, Xxi),
      MatMatMult(fn, interp3Grad, nxi),
      MatMatMult(d, interp3Grad, d0xi),
      MatMatMult(vars, interp3Grad, u0xi)
      // TODO : Tying strain one ?
    )
    

    // interpolate coordinates, director, midplane displacements with the basis
    basis::template interpFields<3, 3>(pt, Xpts, X.get_data());
    basis::template interpFields<3, 3>(pt, fn, n0.get_data());
    basis::template interpFields<1, 1>(pt, etn, et.value().get_data());
    basis::template interpFields<3, 3>(pt, d, d0.value().get_data());

    basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi.get_data());
    basis::template interpFieldsGrad<3, 3>(pt, fn, nxi.get_data());
    basis::template interpFieldsGrad<3, 3>(pt, d, d0xi.value().get_data());
    basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars, u0xi.value().get_data());

    basis::interpTyingStrain(pt, ety, gty.value().get_data());

    // setup before A2D strain energy stack
    // ------------------------------------

    // Compute the transformation at the quadrature point
    // want to implement transform somehow inside the stack now..
    transform->computeTransform(Xxi.get_data(), n0.get_data(), T.get_data()); 

    // compute ABD matrix from shell theory (prospective)
    A2D::SymMat<TacsScalar,9> ABD; // normally ABD is 6x6, but this one includes transverse shear and drill strains
    con->getABDmatrix(0, pt, X.get_data(), ABD.get_data()); // TODO make this routine

    // passive variables for strain energy stack
    A2D::Vec<TacsScalar, 3> zero;

    // active variables for strain energy stack
    A2D::A2DObj<TacsScalar> detXd, ES_dot, Uelem;
    A2D::A2DObj<A2D::Mat<TacsScalar, 3, 3>> Xd, Xdz, Xdinv, XdinvT;
    A2D::A2DObj<A2D::Vec<TacsScalar,9>> E, S;
    A2D::A2DObj<A2D::Mat<TacsScalar, 3, 3>> u0x_tmp, u1x_tmp1, u1x_tmp2, u1x_tmp3, u1x_term1, u1x_term2, u1x_sum; // temp variables
    A2D::A2DObj<A2D::Mat<TacsScalar, 3, 3>> u0xi_frame, u1xi_frame, u0x, u1x;

    const A2D::MatOp NORMAL = A2D::MatOp::NORMAL, TRANSPOSE = A2D::MatOp::TRANSPOSE;
    const A2D::ShellStrainType STRAIN_TYPE = A2D::ShellStrainType::LINEAR; // if condition on type of model here..

    // printf("Pre strain energy stack\n");

    // compute the strain energy from d0, d0xi, u0xi
    auto strain_energy_stack = A2D::MakeStack(
      // part 1 - compute shell basis and transform matrices (passive portion)
      A2D::ShellAssembleFrame(Xxi, n0, Xd), 
      A2D::ShellAssembleFrame(nxi, zero, Xdz), 
      A2D::ShellAssembleFrame(u0xi, d0, u0xi_frame),
      A2D::ShellAssembleFrame(d0xi, zero, u1xi_frame),
      A2D::MatInv(Xd, Xdinv),
      A2D::MatDet(Xd, detXd),
      A2D::MatMatMult(Xdinv, T, XdinvT),
      // part 2 - compute u0x midplane disp gradient
      A2D::MatMatMult(u0xi_frame, Xdinv, u0x_tmp),
      A2D::MatRotateFrame(u0x_tmp, T, u0x),
      // part 3 - compute u1x director disp gradient
      A2D::MatMatMult(u1xi_frame, Xdinv, u1x_term1),
      // computes u0xi_frame * Xdinv * Xdz * Xdinv => u1x_term2 
      A2D::MatMatMult(u0xi_frame, Xdinv, u1x_tmp1), 
      A2D::MatMatMult(u1x_tmp1, Xdz, u1x_tmp2), 
      A2D::MatMatMult(u1x_tmp2, Xdinv, u1x_term2),
      // compute final u1x = T^T * (u0x * Xdinv - u1x * Xdinv * Xdz * Xdinv) * T
      A2D::MatSum(1.0, u1x_term1, -1.0, u1x_term2, u1x_sum), // for some reason this entry has no hzero?
      A2D::MatRotateFrame(u1x_sum, T, u1x),
      // part 4 - compute transformed tying strain e0ty
      A2D::MatRotateFrame(gty, XdinvT, e0ty),
      // part 5 - compute strains, stresses and then strain energy
      A2D::ShellStrain<STRAIN_TYPE>(u0x, u1x, e0ty, et, E),
      A2D::MatMatMult(ABD, E, S),
      // part 6 - compute strain energy
      A2D::VecDot(E, S, ES_dot),
      A2D::Eval(0.5 * weight * detXd * ES_dot, Uelem)
    );

    // printf("Post strain energy stack defn\n");

    // reverse mode 1st order AD for the strain energy stack
    // -------------------------------------------------
    Uelem.bvalue() = 1.0;
    Uelem.hvalue() = 0.0;
    strain_energy_stack.reverse();

    // printf("Post strain energy stack.reverse\n");

    // reverse through the basis back to the director class, drill strain, tying strain
    basis::template addInterpFieldsTranspose<1, 1>(pt, et.bvalue().get_data(), detn);
    basis::template addInterpFieldsTranspose<3, 3>(pt, d0.bvalue().get_data(), dd);

    basis::template addInterpFieldsGradTranspose<3, 3>(pt, d0xi.bvalue().get_data(), dd);
    basis::template addInterpFieldsGradTranspose<vars_per_node, 3>(pt, u0xi.bvalue().get_data(), res);

    basis::addInterpTyingStrainTranspose(pt, gty.bvalue().get_data(), dety);

    // reverse mode 2nd order AD for the strain energy stack
    // -----------------------------------------------------
    
    // see extractJacobian() => I think I can adapt that to get the Jacobians
    // with special calls on pvalue => hvalue to get the Jacobians (test it out)
    // also I can try to code it up myself with unit vectors if I need to as well
    // expect to get cross-hessians is you pvalue one and hvalue the other obj
    //  NOTE : Jacobian or Hessian terms just need to be setup as matrices based on hextract 
    
    // TODO : could this be made more efficient by extracting multiple Jacobians at the same time with a tuple?
    // e.g. if I pvalue() the same thing multiple times and want to get out the hvalue() for several different things
    // could make a new routine for hextract_multi() or something that takes a tuple.. (less calls then)

    // backprop the drill strain to nodal level
    A2D::Mat<TacsScalar,1,1> d2et;
    strain_energy_stack.hzero(); // hextract was failing because can't call hzero on the stack..
    strain_energy_stack.hextract(et.pvalue(), et.hvalue(), d2et);
    basis::template addInterpFieldsOuterProduct<1, 1, 1, 1>(pt, d2et.get_data(), d2etn);

    // replaces TacsShellAddDispGradHessian
    A2D::Mat<TacsScalar, 3, 3> d2d0;
    A2D::Mat<TacsScalar, 3, 6> d2d0d0xi, d2d0u0xi;
    A2D::Mat<TacsScalar, 6, 6> d2d0xi, d2d0xiu0xi, d2u0xi;

    // TODO : could make new hextract_multi(...,) with VarTuple/vector for second arg so one call for repeat pvalue inputs
    strain_energy_stack.hextract(d0.pvalue(), d0.hvalue(), d2d0);
    strain_energy_stack.hextract(d0.pvalue(), d0xi.hvalue(), d2d0d0xi);
    strain_energy_stack.hextract(d0.pvalue(), u0xi.hvalue(), d2d0u0xi);
    strain_energy_stack.hextract(d0xi.pvalue(), d0xi.hvalue(), d2d0xi);
    strain_energy_stack.hextract(d0xi.pvalue(), u0xi.hvalue(), d2d0xiu0xi);
    strain_energy_stack.hextract(u0xi.pvalue(), u0xi.hvalue(), d2u0xi);

    basis::template addInterpFieldsOuterProduct<3, 3, 3, 3>(pt, d2d0.get_data(), d2d);
    basis::template addInterpGradOuterProduct<3, 3, 3, 3>(pt, d2d0xi.get_data(), d2d);
    basis::template addInterpGradMixedOuterProduct<3, 3, 3, 3>(pt, d2d0d0xi.get_data(), d2d0d0xi.get_data(), d2d);
    basis::template addInterpGradMixedOuterProduct<3, 3, 3, 3>(pt, d2d0u0xi.get_data(), NULL, d2du);
    basis::template addInterpGradOuterProduct<3, 3, 3, 3>(pt, d2d0xiu0xi.get_data(), d2du);
    if (mat) {
      basis::template addInterpGradOuterProduct<vars_per_node, vars_per_node, 3, 3>(pt, d2u0xi.get_data(), mat);
    }

    // replaces addInterpTyingStrainHessian
    A2D::Mat<TacsScalar, 6, 6> d2gty;
    strain_energy_stack.hextract(gty.pvalue(), gty.hvalue(), d2gty);
    basis::addInterpTyingStrainHessian(pt, d2gty.get_data(), d2ety);

    // replaces TacsShellAddTypingDispCoupling
    // double check this (not quite right probably)
    A2D::Mat<TacsScalar, 6, 3> d2gtyd0;
    A2D::Mat<TacsScalar, 6, 6> d2gtyd0xi, d2gtyu0xi;
    strain_energy_stack.hextract(gty.pvalue(), d0.hvalue(), d2gtyd0); // should be able to get cross hessian A => B or B => A I believe
    strain_energy_stack.hextract(gty.pvalue(), d0xi.hvalue(), d2gtyd0xi);
    strain_energy_stack.hextract(gty.pvalue(), u0xi.hvalue(), d2gtyu0xi);

    // replace TacsShellAddTypingDispCoupling
    // TODO : add basis steps to stack so this part goes away
    TacsScalar d2gtyu[6 * usize], d2gtyd[6 * dsize];
    memset(d2gtyu, 0, 6 * usize * sizeof(TacsScalar));
    memset(d2gtyd, 0, 6 * dsize * sizeof(TacsScalar));
    const int usize = 3 * basis::NUM_NODES;
    const int dsize = 3 * basis::NUM_NODES;
    for (int k = 0; k < 6; k++) {
      // Compute the director field and the gradient of the director
      // field at the specified point
      basis::template addInterpFieldsTranspose<3, 3>(pt, d2gtyd0.get_data()[3 * k], &d2gtyd[dsize * k]);
      basis::template addInterpFieldsGradTranspose<3, 3>(pt, d2gtyd0xi.get_data()[6 * k], &d2gtyd[dsize * k]);
      basis::template addInterpFieldsGradTranspose<3, 3>(pt, d2gtyu0xi.get_data()[6 * k], &d2gtyu[usize * k]);
    }

    // Add the values into d2etyu and d2etyd
    for (int k = 0; k < usize; k++) {
      TacsScalar t1[6], t2[basis::NUM_TYING_POINTS];
      memset(t2, 0, basis::NUM_TYING_POINTS * sizeof(TacsScalar));

      for (int kk = 0; kk < 6; kk++) {
        t1[kk] = d2gtyu[usize * kk + k];
      }

      basis::addInterpTyingStrainTranspose(pt, t1, t2);

      for (int kk = 0; kk < basis::NUM_TYING_POINTS; kk++) {
        d2etyu[kk * usize + k] += t2[kk];
      }
    }

    for (int k = 0; k < dsize; k++) {
      TacsScalar t1[6], t2[basis::NUM_TYING_POINTS];
      memset(t2, 0, basis::NUM_TYING_POINTS * sizeof(TacsScalar));

      for (int kk = 0; kk < 6; kk++) {
        t1[kk] = d2gtyd[dsize * kk + k];
      }

      basis::addInterpTyingStrainTranspose(pt, t1, t2);

      for (int kk = 0; kk < basis::NUM_TYING_POINTS; kk++) {
        d2etyd[kk * dsize + k] += t2[kk];
      }
    }

    // setup before kinetic energy stack
    // ------------------------------------

    // passive variables
    A2D::Vec<TacsScalar,3> moments, u0ddot, d0ddot;

    // active variables
    A2D::A2DObj<TacsScalar> uu_term, ud_term1, ud_term2, dd_term, dTelem_dt;
    A2D::A2DObj<A2D::Vec<TacsScalar,3>> u0dot, d0dot;

    // evaluate mass moments
    con->evalMassMoments(elemIndex, pt, X.get_data(), moments.get_data());
    // interpolate first time derivatives
    basis::template interpFields<vars_per_node, 3>(pt, dvars, u0dot.value().get_data());
    basis::template interpFields<3, 3>(pt, ddot, d0dot.value().get_data());
    // interpolate second time derivatives
    basis::template interpFields<vars_per_node, 3>(pt, ddvars, u0ddot.get_data());
    basis::template interpFields<3, 3>(pt, dddot, d0ddot.get_data());

    // due to integration by parts, residual is based on dT/dt, time derivative of KE so 1st and 2nd time derivatives used
    //   double check: but Jacobian should be obtained with a cross Hessian d^2(dT/dt)/du0dot/du0ddot (and same for directors d0)
    auto kinetic_energy_stack = A2D::MakeStack(
      A2D::VecDot(u0dot, u0ddot, uu_term),
      A2D::VecDot(u0dot, d0ddot, ud_term1),
      A2D::VecDot(u0ddot, d0dot, ud_term2),
      A2D::VecDot(d0dot, d0ddot, dd_term),
      Eval(detXd * (moments[0] * uu_term + moments[1] * (ud_term1 + ud_term2) + moments[2] * dd_term), dTelem_dt)
    );

    // now reverse to from dTelem_dt => u0dot, d0dot sensitivities
    dTelem_dt.bvalue() = 1.0;
    dTelem_dt.hvalue() = 0.0;
    kinetic_energy_stack.reverse();

    // backpropagate the time derivatives to the residual
    basis::template addInterpFieldsTranspose<vars_per_node, 3>(pt, u0dot.bvalue().get_data(), res);
    basis::template addInterpFieldsTranspose<3, 3>(pt, d0dot.value().get_data(), dd);

    // now get the Hessians (but slightly different formulation here..)
    A2D::Mat<TacsScalar, 3, 3> d2u0dot, d2u0dotd0dot, d2d0dot, d2u0dot_scaled;
    kinetic_energy_stack.hextract(u0dot.pvalue(), u0dot.hvalue(), d2u0dot);
    kinetic_energy_stack.hextract(u0dot.pvalue(), d0dot.hvalue(), d2u0dotd0dot);
    kinetic_energy_stack.hextract(d0dot.pvalue(), d0dot.hvalue(), d2d0dot);
    // how to scale d2u0dot by gamma easily?
    for (int i = 0; i < 9; i++) {
      d2u0dot_scaled.get_data()[i] = gamma * d2u0dot.get_data()[i];
    }
    // A2D::MatScale<TacsScalar, 3, 3>(gamma_, d2u0dot.get_data(), d2u0dot_scaled.get_data());

    basis::template addInterpFieldsOuterProduct<vars_per_node, vars_per_node, 3, 3>(pt, d2u0dot_scaled.get_data(), mat);
    basis::template addInterpFieldsOuterProduct<3, 3, 3, 3>(pt, d2u0dotd0dot.get_data(), d2Tdotu);
    basis::template addInterpFieldsOuterProduct<3, 3, 3, 3>(pt, d2d0dot.get_data(), d2Tdotd);
  }

  // printf("Done with quad loop\n");

  // is it possible to A2D the nodal steps? => maybe too hard?

  // Add the contribution to the residual from the drill strain
  TacsShellAddDrillStrainHessian<vars_per_node, offset, basis, director, model>(
      Xdn, fn, vars, XdinvTn, Tn, u0xn, Ctn, detn, d2etn, res, mat);

  // Add the residual from the tying strain
  model::template addComputeTyingStrainTranspose<vars_per_node, basis>(
      Xpts, fn, vars, d, dety, res, dd);

  // Add the second order terms from the tying strain
  model::template addComputeTyingStrainHessian<vars_per_node, basis>(
      alpha, Xpts, fn, vars, d, dety, d2ety, d2etyu, d2etyd, mat, d2d, d2du);

  // Add the contributions to the stiffness matrix
  director::template addDirectorJacobian<vars_per_node, offset, num_nodes>(
      alpha, beta, gamma, vars, dvars, ddvars, fn, dd, d2Tdotd, d2Tdotu, d2d,
      d2du, res, mat);

  // Add the constraint associated with the rotational parametrization (if any)
  director::template addRotationConstrJacobian<vars_per_node, offset,
                                               num_nodes>(alpha, vars, res,
                                                          mat);

  printf("Done with addJacobian on elem %d\n", elemIndex);
}