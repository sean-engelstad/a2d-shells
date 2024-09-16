#ifndef TACS_DIRECTOR_H
#define TACS_DIRECTOR_H

#include "TACSElementAlgebra.h"
// #include "TACSElementVerification.h"

/*
  The director class.

  Given a reference vector, t, from the element geometry, the director
  computes the exact or approximate rate of change of the displacement
  through the thickness.
*/
class TACSLinearizedRotation {
 public:
  static const int NUM_PARAMETERS = 3;

  /**
    Compute the rotation matrices at each node

    @param vars The full variable vector
    @param C The rotation matrices at each point
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void computeRotationMat(const TacsScalar vars[], TacsScalar C[]) {
    const TacsScalar *q = &vars[offset];
    for (int i = 0; i < num_nodes; i++) {
      // C = I - q^{x}
      setMatSkew(-1.0, q, C);
      C[0] = C[4] = C[8] = 1.0;

      C += 9;
      q += vars_per_node;
    }
  }

  /**
    Compute the derivative of the rotation matrices at each node

    @param vars The full variable vector
    @param varsd The full variable vector
    @param C The rotation matrices at each point
    @param Cd The rotation matrices at each point
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void computeRotationMatDeriv(const TacsScalar vars[],
                                      const TacsScalar varsd[], TacsScalar C[],
                                      TacsScalar Cd[]) {
    const TacsScalar *q = &vars[offset];
    const TacsScalar *qd = &varsd[offset];
    for (int i = 0; i < num_nodes; i++) {
      // C = I - q^{x}
      setMatSkew(-1.0, q, C);
      C[0] = C[4] = C[8] = 1.0;

      // Cd = - qd^{x}
      setMatSkew(-1.0, qd, Cd);

      C += 9;
      Cd += 9;

      q += vars_per_node;
      qd += vars_per_node;
    }
  }

  /**
    Add the contribution to the residual from the rotation matrix

    This code adds the contribution to the residual via the derivative

    d(tr(dC^{T}C(q)))/dq_{i} = d(tr(dC^{T}*(I - q^{x})))/dq_{i}

    @param vars The full variable vector
    @param dC The derivative w.r.t. the rotation matrix
    @param res The residual array
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void addRotationMatResidual(const TacsScalar vars[],
                                     const TacsScalar dC[], TacsScalar res[]) {
    TacsScalar *r = &res[offset];

    for (int i = 0; i < num_nodes; i++) {
      r[0] += -(dC[7] - dC[5]);
      r[1] += -(dC[2] - dC[6]);
      r[2] += -(dC[3] - dC[1]);

      r += vars_per_node;
      dC += 9;
    }
  }

  /*
    Add the Jacobian of the rotation matrix to the output

    @param alpha Scalar coefficient for the Jacobian matrix
    @param vars The variable values
    @param dC The derivative of the functional w.r.t. C
    @param d2C The second derivatives of the functional w.r.t. C
    @param res The residual
    @param mat The Jacobian matrix
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void addRotationMatJacobian(const TacsScalar alpha,
                                     const TacsScalar vars[],
                                     const TacsScalar dC[],
                                     const TacsScalar d2C[], TacsScalar res[],
                                     TacsScalar mat[]) {
    const int size = vars_per_node * num_nodes;
    const int csize = 9 * num_nodes;

    TacsScalar *r = NULL;
    if (res) {
      r = &res[offset];
    }
    TacsScalar *m = &mat[offset * size + offset];

    for (int i = 0; i < num_nodes; i++) {
      if (res) {
        r[0] += -(dC[7] - dC[5]);
        r[1] += -(dC[2] - dC[6]);
        r[2] += -(dC[3] - dC[1]);
        r += vars_per_node;
      }

      for (int j = 0; j < num_nodes; j++) {
        m[vars_per_node * j] += d2C[csize * (9 * i + 5) + 9 * j + 5] -
                                d2C[csize * (9 * i + 5) + 9 * j + 7] -
                                d2C[csize * (9 * i + 7) + 9 * j + 5] +
                                d2C[csize * (9 * i + 7) + 9 * j + 7];

        m[vars_per_node * j + 1] += -d2C[csize * (9 * i + 5) + 9 * j + 2] +
                                    d2C[csize * (9 * i + 5) + 9 * j + 6] +
                                    d2C[csize * (9 * i + 7) + 9 * j + 2] -
                                    d2C[csize * (9 * i + 7) + 9 * j + 6];

        m[vars_per_node * j + 2] += d2C[csize * (9 * i + 5) + 9 * j + 1] -
                                    d2C[csize * (9 * i + 5) + 9 * j + 3] -
                                    d2C[csize * (9 * i + 7) + 9 * j + 1] +
                                    d2C[csize * (9 * i + 7) + 9 * j + 3];

        m[vars_per_node * j + size] += -d2C[csize * (9 * i + 2) + 9 * j + 5] +
                                       d2C[csize * (9 * i + 2) + 9 * j + 7] +
                                       d2C[csize * (9 * i + 6) + 9 * j + 5] -
                                       d2C[csize * (9 * i + 6) + 9 * j + 7];

        m[vars_per_node * j + 1 + size] +=
            d2C[csize * (9 * i + 2) + 9 * j + 2] -
            d2C[csize * (9 * i + 2) + 9 * j + 6] -
            d2C[csize * (9 * i + 6) + 9 * j + 2] +
            d2C[csize * (9 * i + 6) + 9 * j + 6];

        m[vars_per_node * j + 2 + size] +=
            -d2C[csize * (9 * i + 2) + 9 * j + 1] +
            d2C[csize * (9 * i + 2) + 9 * j + 3] +
            d2C[csize * (9 * i + 6) + 9 * j + 1] -
            d2C[csize * (9 * i + 6) + 9 * j + 3];

        m[vars_per_node * j + 2 * size] +=
            d2C[csize * (9 * i + 1) + 9 * j + 5] -
            d2C[csize * (9 * i + 1) + 9 * j + 7] -
            d2C[csize * (9 * i + 3) + 9 * j + 5] +
            d2C[csize * (9 * i + 3) + 9 * j + 7];

        m[vars_per_node * j + 1 + 2 * size] +=
            -d2C[csize * (9 * i + 1) + 9 * j + 2] +
            d2C[csize * (9 * i + 1) + 9 * j + 6] +
            d2C[csize * (9 * i + 3) + 9 * j + 2] -
            d2C[csize * (9 * i + 3) + 9 * j + 6];

        m[vars_per_node * j + 2 + 2 * size] +=
            d2C[csize * (9 * i + 1) + 9 * j + 1] -
            d2C[csize * (9 * i + 1) + 9 * j + 3] -
            d2C[csize * (9 * i + 3) + 9 * j + 1] +
            d2C[csize * (9 * i + 3) + 9 * j + 3];
      }

      m += vars_per_node * size;
      dC += 9;
    }
  }

  /**
    The linearized rotation class is unconstrained
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void addRotationConstraint(const TacsScalar vars[], TacsScalar res[]) {
  }

  template <int vars_per_node, int offset, int num_nodes>
  static void addRotationConstrJacobian(const TacsScalar alpha,
                                        const TacsScalar vars[],
                                        TacsScalar res[], TacsScalar mat[]) {}

  /**
    Compute the director and rates at all nodes.

    d = Q(q)*t = (C(q)^{T} - I)*t
    ddot = d/dt(Q(q))*t

    @param vars The full variable vector
    @param dvars The first time derivative of the variables
    @param ddvars The second derivatives of the variables
    @param t The reference directions
    @param d The director values
    @param ddot The first time derivative of the director
    @param dddot The second time derivative of the director
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void computeDirectorRates(const TacsScalar vars[],
                                   const TacsScalar dvars[],
                                   const TacsScalar t[], TacsScalar d[],
                                   TacsScalar ddot[]) {
    const TacsScalar *q = &vars[offset];
    const TacsScalar *qdot = &dvars[offset];
    for (int i = 0; i < num_nodes; i++) {
      crossProduct(q, t, d);
      crossProduct(qdot, t, ddot);

      t += 3;
      d += 3;
      ddot += 3;

      q += vars_per_node;
      qdot += vars_per_node;
    }
  }

  /**
    Compute the director and rates at all nodes.

    d = Q(q)*t = (C(q)^{T} - I)*t
    ddot = d/dt(Q(q))*t
    dddot = d^2/dt^2(Q(q))*t

    @param vars The full variable vector
    @param dvars The first time derivative of the variables
    @param ddvars The second derivatives of the variables
    @param t The reference directions
    @param d The director values
    @param ddot The first time derivative of the director
    @param dddot The second time derivative of the director
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void computeDirectorRates(const TacsScalar vars[],
                                   const TacsScalar dvars[],
                                   const TacsScalar ddvars[],
                                   const TacsScalar t[], TacsScalar d[],
                                   TacsScalar ddot[], TacsScalar dddot[]) {
    const TacsScalar *q = &vars[offset];
    const TacsScalar *qdot = &dvars[offset];
    const TacsScalar *qddot = &ddvars[offset];
    for (int i = 0; i < num_nodes; i++) {
      crossProduct(q, t, d);
      crossProduct(qdot, t, ddot);
      crossProduct(qddot, t, dddot);

      t += 3;
      d += 3;
      ddot += 3;
      dddot += 3;

      q += vars_per_node;
      qdot += vars_per_node;
      qddot += vars_per_node;
    }
  }

  /**
    Compute the director and rates at all nodes and the derivative.

    d = Q(q)*t = (C(q)^{T} - I)*t
    ddot = d/dt(Q(q))*t
    dddot = d^2/dt^2(Q(q))*t

    @param vars The full variable vector
    @param dvars The first time derivative of the variables
    @param ddvars The second derivatives of the variables
    @param psi The full variable vector derivative
    @param t The reference directions
    @param C The rotation matrices at each point
    @param d The director values
    @param ddot The first time derivative of the director
    @param dddot The second time derivative of the director
    @param dd The derivative of the director values
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void computeDirectorRatesDeriv(const TacsScalar vars[],
                                        const TacsScalar dvars[],
                                        const TacsScalar ddvars[],
                                        const TacsScalar psi[],
                                        const TacsScalar t[], TacsScalar d[],
                                        TacsScalar ddot[], TacsScalar dddot[],
                                        TacsScalar dpsi[]) {
    const TacsScalar *q = &vars[offset];
    const TacsScalar *qdot = &dvars[offset];
    const TacsScalar *qddot = &ddvars[offset];
    const TacsScalar *qpsi = &psi[offset];
    for (int i = 0; i < num_nodes; i++) {
      crossProduct(q, t, d);
      crossProduct(qdot, t, ddot);
      crossProduct(qddot, t, dddot);
      crossProduct(qpsi, t, dpsi);

      t += 3;
      d += 3;
      ddot += 3;
      dddot += 3;
      dpsi += 3;

      q += vars_per_node;
      qdot += vars_per_node;
      qddot += vars_per_node;
      qpsi += vars_per_node;
    }
  }

  /**
    Given the derivatives of the kinetic energy expression with
    respect to time, add the contributions to the derivative of the

    Given the partial derivatives of the Lagrangian with respect to the
    director and the time derivative of the vector, compute

    dd = d/dt(dT/d(dot{d})) - dL/dd

    In general, the residual contribution is:

    res +=
    dTdot*d(dot{d})/d(dot{q}) +
    dT/d(dot{d})*d/dt(dot{d})/d(dot{q}) +
    dd*d(d)/d(q)

    For the linearized rotation director these expressions are:

    d = q^{x} t
    dot{d} = - t^{x} dot{q}
    d(dot{d})/d(dot{q}) = - t^{x}
    d/dt(d(dot{d})/d(dot{q})) = 0

    @param vars The full variable vector
    @param dvars The first time derivative of the variables
    @param ddvars The second derivatives of the variables
    @param t The normal direction
    @param dd The contribution from the derivative of the director
    @param res The output residual
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void addDirectorResidual(const TacsScalar vars[],
                                  const TacsScalar dvars[],
                                  const TacsScalar ddvars[],
                                  const TacsScalar t[], const TacsScalar dd[],
                                  TacsScalar res[]) {
    TacsScalar *r = &res[offset];

    for (int i = 0; i < num_nodes; i++) {
      crossProductAdd(1.0, t, dd, r);

      r += vars_per_node;
      dd += 3;
      t += 3;
    }
  }

  /*
    Add terms from the Jacobian
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void addDirectorJacobian(
      TacsScalar alpha, TacsScalar beta, TacsScalar gamma,
      const TacsScalar vars[], const TacsScalar dvars[],
      const TacsScalar ddvars[], const TacsScalar t[], const TacsScalar dd[],
      const TacsScalar d2Tdotd[], const TacsScalar d2Tdotu[],
      const TacsScalar d2d[], const TacsScalar d2du[], TacsScalar res[],
      TacsScalar mat[]) {
    // Add the derivative due to and d2d
    const int dsize = 3 * num_nodes;
    const int nvars = vars_per_node * num_nodes;

    // d = crossProduct(q, t, d)
    const TacsScalar *ti = t;
    for (int i = 0; i < num_nodes; i++, ti += 3) {
      TacsScalar *jac1 = &mat[(offset + vars_per_node * i) * nvars + offset];
      TacsScalar *jac2 =
          &mat[(offset + vars_per_node * i + 1) * nvars + offset];
      TacsScalar *jac3 =
          &mat[(offset + vars_per_node * i + 2) * nvars + offset];

      const TacsScalar *tj = t;
      for (int j = 0; j < num_nodes; j++, tj += 3) {
        // Add the derivative
        TacsScalar d[9];
        d[0] = d2d[0] + gamma * d2Tdotd[0];
        d[1] = d2d[1] + gamma * d2Tdotd[1];
        d[2] = d2d[2] + gamma * d2Tdotd[2];

        d[3] = d2d[dsize] + gamma * d2Tdotd[dsize];
        d[4] = d2d[dsize + 1] + gamma * d2Tdotd[dsize + 1];
        d[5] = d2d[dsize + 2] + gamma * d2Tdotd[dsize + 2];

        d[6] = d2d[2 * dsize] + gamma * d2Tdotd[2 * dsize];
        d[7] = d2d[2 * dsize + 1] + gamma * d2Tdotd[2 * dsize + 1];
        d[8] = d2d[2 * dsize + 2] + gamma * d2Tdotd[2 * dsize + 2];

        TacsScalar tmp[9];
        mat3x3SkewMatSkewTransform(ti, d, tj, tmp);

        jac1[0] -= tmp[0];
        jac1[1] -= tmp[1];
        jac1[2] -= tmp[2];

        jac2[0] -= tmp[3];
        jac2[1] -= tmp[4];
        jac2[2] -= tmp[5];

        jac3[0] -= tmp[6];
        jac3[1] -= tmp[7];
        jac3[2] -= tmp[8];

        jac1 += vars_per_node;
        jac2 += vars_per_node;
        jac3 += vars_per_node;
        d2d += 3;
        d2Tdotd += 3;
      }

      d2d += 2 * dsize;
      d2Tdotd += 2 * dsize;
    }

    for (int i = 0; i < num_nodes; i++) {
      for (int j = 0; j < num_nodes; j++) {
        // Add the derivative
        TacsScalar d[9];
        d[0] = d2du[0] + gamma * d2Tdotu[0];
        d[1] = d2du[1] + gamma * d2Tdotu[1];
        d[2] = d2du[2] + gamma * d2Tdotu[2];

        d[3] = d2du[dsize] + gamma * d2Tdotu[dsize];
        d[4] = d2du[dsize + 1] + gamma * d2Tdotu[dsize + 1];
        d[5] = d2du[dsize + 2] + gamma * d2Tdotu[dsize + 2];

        d[6] = d2du[2 * dsize] + gamma * d2Tdotu[2 * dsize];
        d[7] = d2du[2 * dsize + 1] + gamma * d2Tdotu[2 * dsize + 1];
        d[8] = d2du[2 * dsize + 2] + gamma * d2Tdotu[2 * dsize + 2];

        TacsScalar tmp[9];
        mat3x3SkewMatTransform(&t[3 * i], d, tmp);

        for (int ii = 0; ii < 3; ii++) {
          for (int jj = 0; jj < 3; jj++) {
            int index = (vars_per_node * i + ii + offset) * nvars +
                        vars_per_node * j + jj;

            mat[index] += tmp[3 * ii + jj];
          }
        }

        for (int ii = 0; ii < 3; ii++) {
          for (int jj = 0; jj < 3; jj++) {
            int index = (vars_per_node * j + jj) * nvars + vars_per_node * i +
                        ii + offset;

            mat[index] += tmp[3 * ii + jj];
          }
        }

        d2du += 3;
        d2Tdotu += 3;
      }

      d2du += 2 * dsize;
      d2Tdotu += 2 * dsize;
    }

    // Update residual
    TacsScalar *r = &res[offset];

    for (int i = 0; i < num_nodes; i++) {
      crossProductAdd(1.0, t, dd, r);

      r += vars_per_node;
      dd += 3;
      t += 3;
    }
  }

  /*
    Add the director contributions to the derivative of the normal

    Add the adjoint sensitivity of the reference normal (dt) based on
    the adjoint sensitivity of the director (dd).

    Given that the parametrization is d = (C^{T}(q) - I) * t, compute

    dt += d(dd^{T}d)/dt = dd^{T}*C^{T}(q) = C(q) * dd

    @param vars The full variable vector
    @param t The reference directions
    @param dd The adjoint sensitivities w.r.t. the director
    @param dt The adjoint sensitivity w.r.t. the reference directions
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void addDirectorRefNormalSens(const TacsScalar vars[],
                                       const TacsScalar t[],
                                       const TacsScalar dd[], TacsScalar dt[]) {
    const TacsScalar *q = &vars[offset];

    for (int i = 0; i < num_nodes; i++) {
      crossProductAdd(1.0, dd, q, dt);

      t += 3;
      dd += 3;
      dt += 3;
      q += vars_per_node;
    }
  }

  /*
    Add the director contributions to the derivative of the normal

    Add the adjoint sensitivity of the reference normal (dt) based on
    the adjoint sensitivity of the director (dd) and the sensitivity
    of the derivative field (ddadj).

    Given that the parametrization is d = (C^{T}(q) - I) * t and the field
    dpsi = d(d)/dq^{T} * psi, compute

    dt += d(dd^{T}d)/dt = dd^{T}*C^{T}(q) = C(q) * dd

    dt += d(ddpsi^{T}*dpsi)/dt = [ d(C(q))/dq * psi ] * ddpsi

    @param vars The full variable vector
    @param psi The full variable vector derivative
    @param t The reference directions
    @param dd The adjoint sensitivities w.r.t. the director
    @param ddpsi The adjoint sensitivities w.r.t. the director derivative
    @param dt The adjoint sensitivity w.r.t. the reference directions
  */
  template <int vars_per_node, int offset, int num_nodes>
  static void addDirectorRefNormalSens(
      const TacsScalar vars[], const TacsScalar psi[], const TacsScalar t[],
      const TacsScalar dd[], const TacsScalar ddpsi[], TacsScalar dt[]) {
    const TacsScalar *q = &vars[offset];
    const TacsScalar *qpsi = &psi[offset];

    for (int i = 0; i < num_nodes; i++) {
      crossProductAdd(1.0, dd, q, dt);
      crossProductAdd(1.0, ddpsi, qpsi, dt);

      t += 3;
      dd += 3;
      ddpsi += 3;
      dt += 3;
      q += vars_per_node;
      qpsi += vars_per_node;
    }
  }

  static TacsScalar evalDrillStrain(const TacsScalar u0x[],
                                    const TacsScalar Ct[]) {
    // Compute the rotational penalty
    return 0.5 * (Ct[3] + u0x[3] - Ct[1] - u0x[1]);
  }

  static void evalDrillStrainSens(TacsScalar scale, const TacsScalar u0x[],
                                  const TacsScalar Ct[], TacsScalar du0x[],
                                  TacsScalar dCt[]) {
    dCt[0] = 0.0;
    dCt[1] = -0.5 * scale;
    dCt[2] = 0.0;
    dCt[3] = 0.5 * scale;
    dCt[4] = 0.0;
    dCt[5] = 0.0;
    dCt[6] = 0.0;
    dCt[7] = 0.0;
    dCt[8] = 0.0;

    du0x[0] = 0.0;
    du0x[1] = -0.5 * scale;
    du0x[2] = 0.0;
    du0x[3] = 0.5 * scale;
    du0x[4] = 0.0;
    du0x[5] = 0.0;
    du0x[6] = 0.0;
    du0x[7] = 0.0;
    du0x[8] = 0.0;
  }

  static TacsScalar evalDrillStrainDeriv(const TacsScalar u0x[],
                                         const TacsScalar Ct[],
                                         const TacsScalar u0xd[],
                                         const TacsScalar Ctd[],
                                         TacsScalar *ed) {
    *ed = 0.5 * (Ctd[3] + u0xd[3] - Ctd[1] - u0xd[1]);

    // Compute the rotational penalty
    return 0.5 * (Ct[3] + u0x[3] - Ct[1] - u0x[1]);
  }

  static void evalDrillStrainHessian(TacsScalar d2et, const TacsScalar u0x[],
                                     const TacsScalar Ct[], TacsScalar d2u0x[],
                                     TacsScalar d2Ct[], TacsScalar d2Ctu0x[]) {
    memset(d2u0x, 0, 81 * sizeof(TacsScalar));
    memset(d2Ct, 0, 81 * sizeof(TacsScalar));
    memset(d2Ctu0x, 0, 81 * sizeof(TacsScalar));
  }
};

#endif  // TACS_DIRECTOR_H
