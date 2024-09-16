/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2014 Georgia Tech Research Corporation

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#include "TACSElementBasis.h"

#include "TACSElementAlgebra.h"

// /*
//   Get the layout type
// */
// ElementLayout TACSElementBasis::getLayoutType() { return TACS_LAYOUT_NONE; }

TacsScalar TACSElementBasis::getFaceNormal(int face, int n,
                                           const TacsScalar Xpts[],
                                           TacsScalar X[], TacsScalar Xd[],
                                           TacsScalar normal[]) {
  const int num_params = getNumParameters();

  double pt[3];
  double tangents[2 * 3];
  getFaceQuadraturePoint(face, n, pt, tangents);

  // Compute the interpolation on the face
  interpFaceFields(face, n, pt, 3, Xpts, 1, X);
  interpFaceFieldsGrad(face, n, pt, 3, Xpts, Xd);

  if (num_params == 3) {
    // Compute the tangent directions
    TacsScalar t1[3], t2[3];
    mat3x3Mult(Xd, &tangents[0], t1);
    mat3x3Mult(Xd, &tangents[3], t2);

    // Find the cross product
    crossProduct(1.0, t1, t2, normal);

    TacsScalar A = sqrt(vec3Dot(normal, normal));
    vec3Scale(1.0 / A, normal);

    return A;
  } else if (num_params == 2) {
    // Compute the tangent direction
    TacsScalar t[2];
    mat2x2Mult(Xd, tangents, t);

    // Compute the norm of the vector
    TacsScalar A = sqrt(vec2Dot(t, t));

    normal[0] = t[1] / A;
    normal[1] = -t[0] / A;

    return A;
  }
  return 0.0;
}

TacsScalar TACSElementBasis::getJacobianTransform(int n, const double pt[],
                                                  const TacsScalar Xpts[],
                                                  TacsScalar Xd[],
                                                  TacsScalar J[]) {
  const int num_params = getNumParameters();

  interpFieldsGrad(n, pt, 3, Xpts, Xd);

  if (num_params == 3) {
    // Compute the Jacobian transformation
    TacsScalar det = inv3x3(Xd, J);

    return det;
  } else if (num_params == 2) {
    // Compute the Jacobian transformation
    TacsScalar det = inv2x2(Xd, J);

    return det;
  } else if (num_params == 1) {
    J[0] = 1.0 / Xd[0];

    return Xd[0];
  }

  return 0.0;
}

/*
  Get the field values at the specified quadrature point
*/
void TACSElementBasis::getFieldValues(int n, const double pt[],
                                      const TacsScalar Xpts[],
                                      const int vars_per_node,
                                      const TacsScalar vars[], TacsScalar X[],
                                      TacsScalar U[]) {
  // Compute the interpolation point
  interpFields(n, pt, 3, Xpts, 1, X);
  interpFields(n, pt, vars_per_node, vars, 1, U);
}

/*
  Get the gradient of the field at the quadrature point
*/
TacsScalar TACSElementBasis::getFieldGradient(
    int n, const double pt[], const TacsScalar Xpts[], const int vars_per_node,
    const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], TacsScalar X[], TacsScalar Xd[], TacsScalar J[],
    TacsScalar Ut[], TacsScalar Ud[], TacsScalar Ux[]) {
  const int num_params = getNumParameters();

  // Compute the interpolation point
  interpFields(n, pt, 3, Xpts, 1, X);

  // Compute derivative of the interpolation
  interpFieldsGrad(n, pt, 3, Xpts, Xd);

  if (vars && dvars && ddvars) {
    interpFieldsGrad(n, pt, vars_per_node, vars, Ud);
    interpFields(n, pt, vars_per_node, vars, dvars, ddvars, Ut);
  } else {
    interpFieldsGrad(n, pt, vars_per_node, vars, Ud);
    interpFields(n, pt, vars_per_node, vars, 3, &Ut[0]);

    // Zero the time derivatives
    for (int j = 0; j < vars_per_node; j++) {
      Ut[3 * j + 1] = Ut[3 * j + 2] = 0.0;
    }
  }

  if (num_params == 3) {
    // Compute the Jacobian transformation
    TacsScalar detXd = inv3x3(Xd, J);

    // U,x = U,xi * J
    TacsScalar *ux = Ux;
    for (int j = 0; j < vars_per_node; j++) {
      mat3x3MultTrans(J, &Ud[3 * j], ux);
      ux += 3;
    }

    return detXd;
  } else if (num_params == 2) {
    // Compute the Jacobian transformation
    TacsScalar detXd = inv2x2(Xd, J);

    // U,x = U,xi * J
    TacsScalar *ux = Ux;
    for (int j = 0; j < vars_per_node; j++) {
      mat2x2MultTrans(J, &Ud[2 * j], ux);
      ux += 2;
    }

    return detXd;
  } else if (num_params == 1) {
    J[0] = 1.0 / Xd[0];
    for (int j = 0; j < vars_per_node; j++) {
      Ux[j] = J[0] * Ud[j];
    }

    return Xd[0];
  }

  return 0.0;
}

TacsScalar TACSElementBasis::getFieldGradient(
    int n, const double pt[], const TacsScalar Xpts[], const int vars_per_node,
    const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], const TacsScalar psi[], TacsScalar X[],
    TacsScalar Xd[], TacsScalar J[], TacsScalar Ut[], TacsScalar Ud[],
    TacsScalar Ux[], TacsScalar Psi[], TacsScalar Psid[], TacsScalar Psix[]) {
  const int num_params = getNumParameters();

  // Compute the interpolation point
  interpFields(n, pt, 3, Xpts, 1, X);

  // Compute derivative of the interpolation
  interpFieldsGrad(n, pt, 3, Xpts, Xd);

  if (vars && dvars && ddvars && psi) {
    interpFieldsGrad(n, pt, vars_per_node, vars, Ud);
    interpFieldsGrad(n, pt, vars_per_node, psi, Psid);
    interpFields(n, pt, vars_per_node, psi, 1, Psi);
    interpFields(n, pt, vars_per_node, vars, dvars, ddvars, Ut);
  } else {
    interpFieldsGrad(n, pt, vars_per_node, vars, Ud);
    interpFieldsGrad(n, pt, vars_per_node, psi, Psid);
    interpFields(n, pt, vars_per_node, psi, 1, Psi);
    interpFields(n, pt, vars_per_node, vars, 3, &Ut[0]);

    // Zero the time derivatives
    for (int j = 0; j < vars_per_node; j++) {
      Ut[3 * j + 1] = Ut[3 * j + 2] = 0.0;
    }
  }

  if (num_params == 3) {
    // Compute the Jacobian transformation
    TacsScalar detXd = inv3x3(Xd, J);

    // U,x = U,xi * J
    TacsScalar *ux = Ux;
    for (int j = 0; j < vars_per_node; j++) {
      mat3x3MultTrans(J, &Ud[3 * j], ux);
      ux += 3;
    }

    TacsScalar *px = Psix;
    for (int j = 0; j < vars_per_node; j++) {
      mat3x3MultTrans(J, &Psid[3 * j], px);
      px += 3;
    }

    return detXd;
  } else if (num_params == 2) {
    // Compute the Jacobian transformation
    TacsScalar detXd = inv2x2(Xd, J);

    // U,x = U,xi * J
    TacsScalar *ux = Ux;
    for (int j = 0; j < vars_per_node; j++) {
      mat2x2MultTrans(J, &Ud[2 * j], ux);
      ux += 2;
    }

    TacsScalar *px = Psix;
    for (int j = 0; j < vars_per_node; j++) {
      mat2x2MultTrans(J, &Psid[2 * j], px);
      px += 2;
    }

    return detXd;
  } else if (num_params == 1) {
    J[0] = 1.0 / Xd[0];
    for (int j = 0; j < vars_per_node; j++) {
      Ux[j] = J[0] * Ud[j];
    }

    return Xd[0];
  }

  return 0.0;
}

/*
  Add the weak form of the governing equations to the residual
*/
void TACSElementBasis::addWeakResidual(int n, const double pt[],
                                       TacsScalar weight, const TacsScalar J[],
                                       const int vars_per_node,
                                       TacsScalar DUt[], TacsScalar DUx[],
                                       TacsScalar res[]) {
  const int num_params = getNumParameters();

  // Add contributions from DUt
  for (int i = 0; i < 3 * vars_per_node; i++) {
    DUt[i] *= weight;
  }
  addInterpFieldsTranspose(n, pt, 3, &DUt[0], vars_per_node, res);
  addInterpFieldsTranspose(n, pt, 3, &DUt[1], vars_per_node, res);
  addInterpFieldsTranspose(n, pt, 3, &DUt[2], vars_per_node, res);

  if (num_params == 3) {
    for (int i = 0; i < vars_per_node; i++) {
      TacsScalar dx[3];
      mat3x3Mult(J, &DUx[3 * i], dx);
      DUx[3 * i] = weight * dx[0];
      DUx[3 * i + 1] = weight * dx[1];
      DUx[3 * i + 2] = weight * dx[2];
    }
  } else if (num_params == 2) {
    for (int i = 0; i < vars_per_node; i++) {
      TacsScalar dx[2];
      mat2x2Mult(J, &DUx[2 * i], dx);
      DUx[2 * i] = weight * dx[0];
      DUx[2 * i + 1] = weight * dx[1];
    }
  } else {
    for (int i = 0; i < vars_per_node; i++) {
      DUx[i] *= J[0] * weight;
    }
  }

  addInterpFieldsGradTranspose(n, pt, vars_per_node, DUx, res);
}

void TACSElementBasis::scaleWeakMatrix(const TacsScalar weight,
                                       const TacsScalar alpha,
                                       const TacsScalar beta,
                                       const TacsScalar gamma,
                                       const int Jac_nnz, const int *Jac_pairs,
                                       TacsScalar *Jac) {
  const int num_params = getNumParameters();

  for (int ii = 0; ii < Jac_nnz; ii++) {
    int jx = Jac_pairs[2 * ii + 1];

    if (Jac[ii] != 0.0) {
      if (jx % (num_params + 3) == 0) {
        Jac[ii] *= weight * alpha;
      } else if (jx % (num_params + 3) == 1) {
        Jac[ii] *= weight * beta;
      } else if (jx % (num_params + 3) == 2) {
        Jac[ii] *= weight * gamma;
      } else {
        Jac[ii] *= weight * alpha;
      }
    }
  }
}

/*
  Add the weak form of the governing equations to the residual
*/
void TACSElementBasis::addWeakMatrix(int n, const double pt[],
                                     const TacsScalar J[],
                                     const int vars_per_node, const int Jac_nnz,
                                     const int *Jac_pairs,
                                     const TacsScalar *Jac, TacsScalar *mat) {
  const int num_nodes = getNumNodes();
  const int num_params = getNumParameters();
  const int num_vars = num_nodes * vars_per_node;
  const int row_incr = num_vars * (vars_per_node - 1);
  const int col_incr = vars_per_node;

  // Compute the transpose of the Jacobian transformation for
  // convenience
  TacsScalar JT[9];
  if (num_params == 1) {
    JT[0] = J[0];
  } else if (num_params == 2) {
    JT[0] = J[0];
    JT[1] = J[2];
    JT[2] = J[1];
    JT[3] = J[3];
  } else {
    JT[0] = J[0];
    JT[1] = J[3];
    JT[2] = J[6];

    JT[3] = J[1];
    JT[4] = J[4];
    JT[5] = J[7];

    JT[6] = J[2];
    JT[7] = J[5];
    JT[8] = J[8];
  }

  // Use a sparse representation
  for (int ii = 0; ii < Jac_nnz; ii++) {
    int ix = Jac_pairs[2 * ii];
    int jx = Jac_pairs[2 * ii + 1];

    if (Jac[ii] != 0.0) {
      // Compute the offset point in the matrix
      const int init_index =
          (num_vars * (ix / (num_params + 3)) + (jx / (num_params + 3)));
      TacsScalar *M = &mat[init_index];

      if ((ix % (num_params + 3) < 3) && (jx % (num_params + 3) < 3)) {
        addInterpOuterProduct(n, pt, Jac[ii], row_incr, col_incr, M);
      } else if (ix % (num_params + 3) < 3) {
        int j = (jx % (num_params + 3)) - 3;
        int transpose = 0;
        addInterpGradOuterProduct(n, pt, transpose, Jac[ii],
                                  &JT[num_params * j], row_incr, col_incr, M);
      } else if (jx % (num_params + 3) < 3) {
        int i = (ix % (num_params + 3)) - 3;
        int transpose = 1;
        addInterpGradOuterProduct(n, pt, transpose, Jac[ii],
                                  &JT[num_params * i], row_incr, col_incr, M);
      } else {
        int i = (ix % (num_params + 3)) - 3;
        int j = (jx % (num_params + 3)) - 3;
        addInterpGradGradOuterProduct(n, pt, Jac[ii], &JT[num_params * i],
                                      &JT[num_params * j], row_incr, col_incr,
                                      M);
      }
    }
  }
}

/*
  Add the weak form of the governing equations to the residual
*/
void TACSElementBasis::addMatVecProduct(const int vars_per_node,
                                        const int Jac_nnz, const int *Jac_pairs,
                                        const TacsScalar *data,
                                        TacsScalar *temp, const TacsScalar px[],
                                        TacsScalar py[]) {
  // Fill in the values for each entry
  interpAllFieldsGrad(vars_per_node, px, temp);

  // Get the number of parameters
  const int num_params = getNumParameters();
  const int np2 = num_params * num_params;

  // Get the number of quadrature points
  const int nquad = getNumQuadraturePoints();

  for (int n = nquad - 1; n >= 0; n--) {
    // Set the locations for the data pointers
    const TacsScalar *J = &data[n * (np2 + Jac_nnz)];
    const TacsScalar *Jac = &data[n * (np2 + Jac_nnz) + np2];

    // Set pointers into the temporary array
    TacsScalar *U = &temp[n * vars_per_node * (1 + num_params)];
    TacsScalar *Ud =
        &temp[n * vars_per_node * (1 + num_params) + vars_per_node];

    if (num_params == 3) {
      for (int i = 0; i < vars_per_node; i++) {
        TacsScalar ux[3];
        mat3x3MultTrans(J, &Ud[3 * i], ux);
        Ud[3 * i] = ux[0];
        Ud[3 * i + 1] = ux[1];
        Ud[3 * i + 2] = ux[2];
      }
    } else if (num_params == 2) {
      for (int i = 0; i < vars_per_node; i++) {
        TacsScalar ux[2];
        mat2x2MultTrans(J, &Ud[2 * i], ux);
        Ud[2 * i] = ux[0];
        Ud[2 * i + 1] = ux[1];
      }
    } else {
      for (int i = 0; i < vars_per_node; i++) {
        Ud[i] *= J[0];
      }
    }

    // For clarity, switch the names here. This is now storing U,x
    const TacsScalar *Ux = Ud;

    // Set pointers to next entry in the temporary array
    TacsScalar *DU = &temp[(n + 1) * vars_per_node * (1 + num_params)];
    TacsScalar *DUx =
        &temp[(n + 1) * vars_per_node * (1 + num_params) + vars_per_node];

    memset(DU, 0, vars_per_node * sizeof(TacsScalar));
    memset(DUx, 0, num_params * vars_per_node * sizeof(TacsScalar));

    // Compute the output based on the Jacobian input
    for (int ii = 0; ii < Jac_nnz; ii++) {
      int ix = Jac_pairs[2 * ii];
      int jx = Jac_pairs[2 * ii + 1];

      if (Jac[ii] != 0.0) {
        int i = ix / (num_params + 3);
        int j = jx / (num_params + 3);

        // Multiply by the appropriate coefficient
        if (ix % (num_params + 3) < 3) {
          if (jx % (num_params + 3) < 3) {
            DU[i] += Jac[ii] * U[j];
          } else {
            j = num_params * j + (jx % (num_params + 3)) - 3;
            DU[i] += Jac[ii] * Ux[j];
          }
        } else {
          i = num_params * i + (ix % (num_params + 3)) - 3;
          if (jx % (num_params + 3) < 3) {
            DUx[i] += Jac[ii] * U[j];
          } else {
            j = num_params * j + (jx % (num_params + 3)) - 3;
            DUx[i] += Jac[ii] * Ux[j];
          }
        }
      }
    }

    if (num_params == 3) {
      for (int i = 0; i < vars_per_node; i++) {
        TacsScalar ud[3];
        mat3x3Mult(J, &DUx[3 * i], ud);
        DUx[3 * i] = ud[0];
        DUx[3 * i + 1] = ud[1];
        DUx[3 * i + 2] = ud[2];
      }
    } else if (num_params == 2) {
      for (int i = 0; i < vars_per_node; i++) {
        TacsScalar ud[2];
        mat2x2Mult(J, &DUx[2 * i], ud);
        DUx[2 * i] = ud[0];
        DUx[2 * i + 1] = ud[1];
      }
    } else {
      for (int i = 0; i < vars_per_node; i++) {
        DUx[i] *= J[0];
      }
    }
  }

  const TacsScalar *out = &temp[vars_per_node * (1 + num_params)];

  addInterpAllFieldsGradTranspose(vars_per_node, out, py);
}

void TACSElementBasis::interpFields(const int n, const double pt[],
                                    const int vars_per_node,
                                    const TacsScalar values[], const int incr,
                                    TacsScalar field[]) {
  const int num_nodes = getNumNodes();
  double N[MAX_NUM_NODES];
  computeBasis(pt, N);
  for (int i = 0; i < vars_per_node; i++) {
    field[incr * i] = 0.0;
    for (int j = 0; j < num_nodes; j++) {
      field[incr * i] += values[vars_per_node * j + i] * N[j];
    }
  }
}

void TACSElementBasis::interpFields(const int n, const double pt[],
                                    const int vars_per_node,
                                    const TacsScalar val1[],
                                    const TacsScalar val2[],
                                    const TacsScalar val3[],
                                    TacsScalar field[]) {
  interpFields(n, pt, vars_per_node, val1, 3, &field[0]);
  interpFields(n, pt, vars_per_node, val2, 3, &field[1]);
  interpFields(n, pt, vars_per_node, val3, 3, &field[2]);
}

void TACSElementBasis::interpFaceFields(const int face, const int n,
                                        const double pt[],
                                        const int vars_per_node,
                                        const TacsScalar values[],
                                        const int incr, TacsScalar field[]) {
  interpFields(-1, pt, vars_per_node, values, incr, field);
}

void TACSElementBasis::addInterpFieldsTranspose(const int n, const double pt[],
                                                const int incr,
                                                const TacsScalar field[],
                                                const int vars_per_node,
                                                TacsScalar values[]) {
  const int num_nodes = getNumNodes();
  double N[MAX_NUM_NODES];
  computeBasis(pt, N);
  for (int i = 0; i < vars_per_node; i++) {
    for (int j = 0; j < num_nodes; j++) {
      values[vars_per_node * j + i] += field[incr * i] * N[j];
    }
  }
}

void TACSElementBasis::addInterpFaceFieldsTranspose(
    const int face, const int n, const double pt[], const int incr,
    const TacsScalar field[], const int vars_per_node, TacsScalar values[]) {
  addInterpFieldsTranspose(-1, pt, incr, field, vars_per_node, values);
}

void TACSElementBasis::interpFieldsGrad(const int n, const double pt[],
                                        const int vars_per_node,
                                        const TacsScalar values[],
                                        TacsScalar grad[]) {
  const int num_nodes = getNumNodes();
  const int num_params = getNumParameters();
  double N[MAX_NUM_NODES], Nxi[3 * MAX_NUM_NODES];
  computeBasisGradient(pt, N, Nxi);
  for (int i = 0; i < vars_per_node; i++) {
    for (int j = 0; j < num_params; j++) {
      grad[num_params * i + j] = 0.0;
    }
    for (int k = 0; k < num_nodes; k++) {
      for (int j = 0; j < num_params; j++) {
        grad[num_params * i + j] +=
            values[vars_per_node * k + i] * Nxi[num_params * k + j];
      }
    }
  }
}

void TACSElementBasis::interpFaceFieldsGrad(const int face, const int n,
                                            const double pt[],
                                            const int vars_per_node,
                                            const TacsScalar values[],
                                            TacsScalar grad[]) {
  interpFieldsGrad(-1, pt, vars_per_node, values, grad);
}

void TACSElementBasis::addInterpFieldsGradTranspose(int n, const double pt[],
                                                    const int vars_per_node,
                                                    const TacsScalar grad[],
                                                    TacsScalar values[]) {
  const int num_nodes = getNumNodes();
  const int num_params = getNumParameters();
  double N[MAX_NUM_NODES], Nxi[3 * MAX_NUM_NODES];
  computeBasisGradient(pt, N, Nxi);
  for (int i = 0; i < vars_per_node; i++) {
    for (int k = 0; k < num_nodes; k++) {
      for (int j = 0; j < num_params; j++) {
        values[vars_per_node * k + i] +=
            grad[num_params * i + j] * Nxi[num_params * k + j];
      }
    }
  }
}

void TACSElementBasis::addInterpFaceFieldsGradTranspose(const int face, int n,
                                                        const double pt[],
                                                        const int vars_per_node,
                                                        const TacsScalar grad[],
                                                        TacsScalar values[]) {
  addInterpFieldsGradTranspose(-1, pt, vars_per_node, grad, values);
}

void TACSElementBasis::addInterpOuterProduct(const int n, const double pt[],
                                             const TacsScalar weight,
                                             const int row_incr,
                                             const int col_incr,
                                             TacsScalar *mat) {
  const int num_nodes = getNumNodes();
  double N[MAX_NUM_NODES];
  computeBasis(pt, N);
  for (int i = 0; i < num_nodes; i++, mat += row_incr) {
    for (int j = 0; j < num_nodes; j++, mat += col_incr) {
      mat[0] += weight * N[i] * N[j];
    }
  }
}

void TACSElementBasis::addInterpGradOuterProduct(
    const int n, const double pt[], const int transpose,
    const TacsScalar weight, const TacsScalar scale[], const int row_incr,
    const int col_incr, TacsScalar *mat) {
  const int num_nodes = getNumNodes();
  const int num_params = getNumParameters();
  double N[MAX_NUM_NODES], Nxi[3 * MAX_NUM_NODES];
  computeBasisGradient(pt, N, Nxi);

  if (transpose) {
    if (num_params == 1) {
      for (int i = 0; i < num_nodes; i++, mat += row_incr) {
        for (int j = 0; j < num_nodes; j++, mat += col_incr) {
          mat[0] += weight * N[j] * (Nxi[i] * scale[0]);
        }
      }
    } else if (num_params == 2) {
      for (int i = 0; i < num_nodes; i++, mat += row_incr) {
        for (int j = 0; j < num_nodes; j++, mat += col_incr) {
          mat[0] += weight * N[j] *
                    (Nxi[2 * i] * scale[0] + Nxi[2 * i + 1] * scale[1]);
        }
      }
    } else if (num_params == 3) {
      for (int i = 0; i < num_nodes; i++, mat += row_incr) {
        for (int j = 0; j < num_nodes; j++, mat += col_incr) {
          mat[0] += weight * N[j] *
                    (Nxi[3 * i] * scale[0] + Nxi[3 * i + 1] * scale[1] +
                     Nxi[3 * i + 2] * scale[2]);
        }
      }
    }
  } else {
    if (num_params == 1) {
      for (int i = 0; i < num_nodes; i++, mat += row_incr) {
        for (int j = 0; j < num_nodes; j++, mat += col_incr) {
          mat[0] += weight * N[i] * (Nxi[j] * scale[0]);
        }
      }
    } else if (num_params == 2) {
      for (int i = 0; i < num_nodes; i++, mat += row_incr) {
        for (int j = 0; j < num_nodes; j++, mat += col_incr) {
          mat[0] += weight * N[i] *
                    (Nxi[2 * j] * scale[0] + Nxi[2 * j + 1] * scale[1]);
        }
      }
    } else if (num_params == 3) {
      for (int i = 0; i < num_nodes; i++, mat += row_incr) {
        for (int j = 0; j < num_nodes; j++, mat += col_incr) {
          mat[0] += weight * N[i] *
                    (Nxi[3 * j] * scale[0] + Nxi[3 * j + 1] * scale[1] +
                     Nxi[3 * j + 2] * scale[2]);
        }
      }
    }
  }
}

void TACSElementBasis::addInterpGradGradOuterProduct(
    const int n, const double pt[], const TacsScalar weight,
    const TacsScalar iscale[], const TacsScalar jscale[], const int row_incr,
    const int col_incr, TacsScalar *mat) {
  const int num_nodes = getNumNodes();
  const int num_params = getNumParameters();
  double N[MAX_NUM_NODES], Nxi[3 * MAX_NUM_NODES];
  computeBasisGradient(pt, N, Nxi);

  if (num_params == 1) {
    for (int i = 0; i < num_nodes; i++, mat += row_incr) {
      for (int j = 0; j < num_nodes; j++, mat += col_incr) {
        mat[0] += weight * (Nxi[i] * iscale[0]) * (Nxi[j] * jscale[0]);
      }
    }
  } else if (num_params == 2) {
    for (int i = 0; i < num_nodes; i++, mat += row_incr) {
      for (int j = 0; j < num_nodes; j++, mat += col_incr) {
        mat[0] += weight *
                  (Nxi[2 * i] * iscale[0] + Nxi[2 * i + 1] * iscale[1]) *
                  (Nxi[2 * j] * jscale[0] + Nxi[2 * j + 1] * jscale[1]);
      }
    }
  } else if (num_params == 3) {
    for (int i = 0; i < num_nodes; i++, mat += row_incr) {
      for (int j = 0; j < num_nodes; j++, mat += col_incr) {
        mat[0] += weight *
                  (Nxi[3 * i] * iscale[0] + Nxi[3 * i + 1] * iscale[1] +
                   Nxi[3 * i + 2] * iscale[2]) *
                  (Nxi[3 * j] * jscale[0] + Nxi[3 * j + 1] * jscale[1] +
                   Nxi[3 * j + 2] * jscale[2]);
      }
    }
  }
}

void TACSElementBasis::interpAllFieldsGrad(const int vars_per_node,
                                           const TacsScalar values[],
                                           TacsScalar out[]) {
  const int nquad = getNumQuadraturePoints();
  const int num_params = getNumParameters();

  for (int n = 0; n < nquad; n++) {
    TacsScalar *U = &out[n * vars_per_node * (1 + num_params)];
    TacsScalar *Ud = &out[n * vars_per_node * (1 + num_params) + vars_per_node];

    double pt[3];
    getQuadraturePoint(n, pt);

    interpFields(n, pt, vars_per_node, values, 1, U);
    interpFieldsGrad(n, pt, vars_per_node, values, Ud);
  }
}

void TACSElementBasis::addInterpAllFieldsGradTranspose(const int vars_per_node,
                                                       const TacsScalar in[],
                                                       TacsScalar values[]) {
  const int nquad = getNumQuadraturePoints();
  const int num_params = getNumParameters();

  for (int n = 0; n < nquad; n++) {
    const TacsScalar *U = &in[n * vars_per_node * (1 + num_params)];
    const TacsScalar *Ud =
        &in[n * vars_per_node * (1 + num_params) + vars_per_node];

    double pt[3];
    getQuadraturePoint(n, pt);

    addInterpFieldsTranspose(n, pt, 1, U, vars_per_node, values);
    addInterpFieldsGradTranspose(n, pt, vars_per_node, Ud, values);
  }
}
