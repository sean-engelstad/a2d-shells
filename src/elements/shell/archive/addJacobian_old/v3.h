/*
  Add the contributions to the residual and Jacobian matrix
*/
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
  A2D::Mat<TacsScalar,3,num_nodes> u_temp;
  for (int inode = 0; inode < num_nodes; inode++) {
    for (int ivar = 0; ivar < 3; ivar++) {
      u_temp(ivar,inode) = vars[vars_per_node * inode + ivar];
    }
  }
  A2D::A2DObj<A2D::Mat<TacsScalar,3,num_nodes>> u(u_temp);
  delete u_temp;

  // Derivative of the director field
  TacsScalar dd[dsize];
  memset(dd, 0, dsize * sizeof(TacsScalar));

  // Second derivatives required for the director
  TacsScalar d2d[dsize * dsize], d2du[usize * dsize];
  TacsScalar d2Tdotd[dsize * dsize], d2Tdotu[usize * dsize];
  memset(d2d, 0, dsize * dsize * sizeof(TacsScalar));
  memset(d2du, 0, usize * dsize * sizeof(TacsScalar));
  memset(d2Tdotd, 0, dsize * dsize * sizeof(TacsScalar));
  memset(d2Tdotu, 0, usize * dsize * sizeof(TacsScalar));

  // Zero the contributions to the tying strain derivatives
  TacsScalar dety[basis::NUM_TYING_POINTS];
  TacsScalar d2ety[basis::NUM_TYING_POINTS * basis::NUM_TYING_POINTS];
  TacsScalar d2etyu[basis::NUM_TYING_POINTS * usize];
  TacsScalar d2etyd[basis::NUM_TYING_POINTS * dsize];
  memset(dety, 0, basis::NUM_TYING_POINTS * sizeof(TacsScalar));
  memset(
      d2ety, 0,
      basis::NUM_TYING_POINTS * basis::NUM_TYING_POINTS * sizeof(TacsScalar));
  memset(d2etyu, 0, basis::NUM_TYING_POINTS * usize * sizeof(TacsScalar));
  memset(d2etyd, 0, basis::NUM_TYING_POINTS * dsize * sizeof(TacsScalar));

  // Compute the node normal directions
  TacsScalar fn[3 * num_nodes], Xdn[9 * num_nodes];
  TacsShellComputeNodeNormals<basis>(Xpts, fn, Xdn);

  // Compute the drill strain penalty at each node
  TacsScalar etn[num_nodes], detn[num_nodes];
  TacsScalar d2etn[num_nodes * num_nodes];
  memset(detn, 0, num_nodes * sizeof(TacsScalar));
  memset(d2etn, 0, num_nodes * num_nodes * sizeof(TacsScalar));

  // Store information about the transformation and derivatives at each node for
  // the drilling degrees of freedom
  TacsScalar XdinvTn[9 * num_nodes], Tn[9 * num_nodes];
  TacsScalar u0xn[9 * num_nodes], Ctn[csize];
  TacsShellComputeDrillStrain<vars_per_node, offset, basis, director, model>(
      transform, Xdn, fn, vars, XdinvTn, Tn, u0xn, Ctn, etn);

  A2D::A2DObj<A2D::Mat<TacsScalar,3,num_nodes>> d;
  TacsScalar ddot[dsize], dddot[dsize];
  director::template computeDirectorRates<vars_per_node, offset, num_nodes>(
      vars, dvars, ddvars, fn, d.value().get_data(), ddot, dddot);

  // Set the total number of tying points needed for this element
  TacsScalar ety[basis::NUM_TYING_POINTS];
  model::template computeTyingStrain<vars_per_node, basis>(Xpts, fn, vars, d, ety);

  // Loop over each quadrature point and add the residual contribution
  for (int quad_index = 0; quad_index < nquad; quad_index++) {
    
    // Get the quadrature weight
    double pt[3];
    double weight = quadrature::getQuadraturePoint(quad_index, pt);

    // printf("iquad = %d\n", quad_index);

    // interpolation and prelim quad loop section
    // ----------------------------------------

    // passive A2D Objs used in interpolation
    A2D::Vec<TacsScalar,3> X, n0;
    A2D::Mat<TacsScalar,3,2> Xxi, nxi;
    A2D::Mat<TacsScalar,3,3> T;
    
    // active A2D objs used in interpolation
    A2D::A2DObj<A2D::Vec<TacsScalar,1>> et;
    A2D::A2DObj<A2D::Vec<TacsScalar,3>> d0;
    A2D::A2DObj<A2D::Mat<TacsScalar,3,2>> d0xi, u0xi;
    A2D::A2DObj<A2D::SymMat<TacsScalar,3>> e0ty, e0ty_tmp, gty; 

    // non-disp information interpolated by basis (since less important to Jacobian and can't do all in one stack anyways)
    basis::template interpFields<3, 3>(pt, Xpts, X.get_data());
    basis::template interpFields<3, 3>(pt, fn, n0.get_data());
    basis::template interpFields<1, 1>(pt, etn, et.value().get_data());
    basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi.get_data());
    basis::template interpFieldsGrad<3, 3>(pt, fn, nxi.get_data());

    // the interpTyingStrain call is very hard to do with A2D, may need custom expression for this..
    basis::interpTyingStrain(pt, ety, gty.value().get_data());

    // get interpolation fields needed for disp interpolation
    A2D::Vec<TacsScalar, num_nodes> interp;
    A2D::Mat<TacsScalar, num_nodes, 2> interpGrad;
    basis::template getInterpFields(pt, interp.get_data()); // TODO : make routines for these
    basis::template getInterpFieldsGrad(pt, interp3Grad.get_data());

    // Compute the transformation at the quadrature point
    transform->computeTransform(Xxi.get_data(), n0.get_data(), T.get_data()); 

    // compute ABD matrix from shell theory (prospective)
    A2D::SymMat<TacsScalar,9> ABD; // normally ABD is 6x6, but this one includes transverse shear and drill strains
    con->getABDmatrix(0, pt, X.get_data(), ABD.get_data()); // TODO make this routine

    // setup before A2D strain energy stack
    // ------------------------------------

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

    // compute the strain energy from u,d nodal quantities
    auto strain_energy_stack = A2D::MakeStack(
      // part 1 - interpolate disps from nodal quantities
      A2D::MatVecMult(d, intepr3, d0),
      A2D::MatMatMult(d, interp3Grad, d0xi),
      A2D::MatMatMult(u, interp3Grad, u0xi)
      // part 2 - compute shell basis and transform matrices (passive portion)
      A2D::ShellAssembleFrame(Xxi, n0, Xd), 
      A2D::ShellAssembleFrame(nxi, zero, Xdz), 
      A2D::ShellAssembleFrame(u0xi, d0, u0xi_frame),
      A2D::ShellAssembleFrame(d0xi, zero, u1xi_frame),
      A2D::MatInv(Xd, Xdinv),
      A2D::MatDet(Xd, detXd),
      A2D::MatMatMult(Xdinv, T, XdinvT),
      // part 3 - compute u0x midplane disp gradient
      A2D::MatMatMult(u0xi_frame, Xdinv, u0x_tmp),
      A2D::MatRotateFrame(u0x_tmp, T, u0x),
      // part 4 - compute u1x director disp gradient
      A2D::MatMatMult(u1xi_frame, Xdinv, u1x_term1),
      // computes u0xi_frame * Xdinv * Xdz * Xdinv => u1x_term2 
      A2D::MatMatMult(u0xi_frame, Xdinv, u1x_tmp1), 
      A2D::MatMatMult(u1x_tmp1, Xdz, u1x_tmp2), 
      A2D::MatMatMult(u1x_tmp2, Xdinv, u1x_term2),
      // compute final u1x = T^T * (u0x * Xdinv - u1x * Xdinv * Xdz * Xdinv) * T
      A2D::MatSum(1.0, u1x_term1, -1.0, u1x_term2, u1x_sum), // for some reason this entry has no hzero?
      A2D::MatRotateFrame(u1x_sum, T, u1x),
      // part 5 - compute transformed tying strain e0ty
      A2D::MatRotateFrame(gty, XdinvT, e0ty),
      // part 6 - compute strains, stresses and then strain energy
      A2D::ShellStrain<STRAIN_TYPE>(u0x, u1x, e0ty, et, E),
      A2D::MatMatMult(ABD, E, S),
      // part 7 - compute strain energy
      A2D::VecDot(E, S, ES_dot),
      A2D::Eval(0.5 * weight * detXd * ES_dot, Uelem)
    );

    // printf("Post strain energy stack defn\n");

    // reverse mode 1st order AD for the strain energy stack
    // -------------------------------------------------
    Uelem.bvalue() = 1.0;
    Uelem.hvalue() = 0.0;
    strain_energy_stack.reverse();

    // reverse through the basis back to the director class, drill strain, tying strain
    basis::template addInterpFieldsTranspose<1, 1>(pt, et.bvalue().get_data(), detn);
    basis::template addInterpFieldsTranspose<3, 3>(pt, d0.bvalue().get_data(), dd);

    basis::template addInterpFieldsGradTranspose<3, 3>(pt, d0xi.bvalue().get_data(), dd);
    basis::template addInterpFieldsGradTranspose<vars_per_node, 3>(pt, u0xi.bvalue().get_data(), res); // u holds some of this now
    basis::addInterpTyingStrainTranspose(pt, gty.bvalue().get_data(), dety);

    // reverse mode 2nd order AD for the strain energy stack
    // ----------------------------------------------------- 

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
    basis::template addInterpGradOuterProduct<3, 3, 3, 3>(pt, d2u0xi.get_data(), d2du);
    if (mat) {
      basis::template addInterpGradOuterProduct<vars_per_node, vars_per_node, 3, 3>(pt, d2u0xi.get_data(), mat);
    }

    // replaces addInterpTyingStrainHessian and TacsShellAddTyingDispCoupling
    A2D::Mat<TacsScalar, 6, 3> d2gtyd0;
    A2D::Mat<TacsScalar, 6, 6> d2gty, d2gtyd0xi, d2gtyu0xi;
    strain_energy_stack.hextract(gty.pvalue(), gty.hvalue(), d2gty);
    strain_energy_stack.hextract(gty.pvalue(), d0.hvalue(), d2gtyd0); // should be able to get cross hessian A => B or B => A I believe
    strain_energy_stack.hextract(gty.pvalue(), d0xi.hvalue(), d2gtyd0xi);
    strain_energy_stack.hextract(gty.pvalue(), u0xi.hvalue(), d2gtyu0xi);

    basis::addInterpTyingStrainHessian(pt, d2gty.get_data(), d2ety);    

    // replace TacsShellAddTypingDispCoupling
    // would like to add basis parts so this goes away (but interpTyingStrain step is very hard to do with A2D.. + extra quad loop)
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
  }

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