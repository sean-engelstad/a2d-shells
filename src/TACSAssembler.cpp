/*
@Description : Shortened version just for linear static analysis 
               (no optimization) and for serial now
@Author : Sean Engelstad
@Date : Sep 10, 2024
*/


#include "TACSAssembler.h"

TACSAssembler::TACSAssembler(int _varsPerNode, int _numElements) : varsPerNode(_varsPerNode), numElements(_numElements) {};
TACSAssembler::~TACSAssembler() {};

/**
  Evaluates the total kinetic and potential energies of the structure

  @param Te The kinetic energy
  @param Pe The potential energy
*/
void TACSAssembler::evalEnergies(TacsScalar *Te, TacsScalar *Pe) {
  // Zero the kinetic and potential energy
  *Te = 0.0;
  *Pe = 0.0;

  // Retrieve pointers to temporary storage
  TacsScalar *vars, *dvars, *elemXpts;
  getDataPointers(elementData, &vars, &dvars, NULL, NULL, &elemXpts, NULL, NULL,
                  NULL);

  // Loop over all elements and add individual contributions to the
  // total energy
  for (int i = 0; i < numElements; i++) {
    int ptr = elementNodeIndex[i];
    int len = elementNodeIndex[i + 1] - ptr;
    const int *nodes = &elementTacsNodes[ptr];
    xptVec->getValues(len, nodes, elemXpts);
    varsVec->getValues(len, nodes, vars);
    dvarsVec->getValues(len, nodes, dvars);

    // Compute and add the element's contributions to the total
    // energy
    TacsScalar elemTe, elemPe;
    elements[i]->computeEnergies(i, time, elemXpts, vars, dvars, &elemTe,
                                 &elemPe);

    // Add up the kinetic and potential energy
    *Te += elemTe;
    *Pe += elemPe;
  }
}

/**
  Assemble the residual

  This residual includes the contributions from element tractions set
  in the auxiliary element classes. Note that the vector entries are
  zeroed first, and that the Dirichlet boundary conditions are applied
  after the assembly of the residual is complete.

  @param residual The residual vector
  @param lambda Scaling factor for the aux element contributions, by default 1
*/
void TACSAssembler::assembleRes(TACSBVec *residual, const TacsScalar lambda) {
  // Sort the list of auxiliary elements - this only performs the
  // sort if it is required (if new elements are added)

  // Zero the residual
  residual->zeroEntries();

  // Retrieve pointers to temporary storage
  TacsScalar *vars, *dvars, *ddvars, *elemRes, *elemXpts;
  getDataPointers(elementData, &vars, &dvars, &ddvars, &elemRes, &elemXpts,
                  NULL, NULL, NULL);

  // Go through and add the residuals from all the elements
  for (int i = 0; i < numElements; i++) {
    int ptr = elementNodeIndex[i];
    int len = elementNodeIndex[i + 1] - ptr;
    const int *nodes = &elementTacsNodes[ptr];
    xptVec->getValues(len, nodes, elemXpts);
    varsVec->getValues(len, nodes, vars);
    dvarsVec->getValues(len, nodes, dvars);
    ddvarsVec->getValues(len, nodes, ddvars);

    // Add the residual from the working element
    int nvars = elements[i]->getNumVariables();
    memset(elemRes, 0, nvars * sizeof(TacsScalar));
    elements[i]->addResidual(i, time, elemXpts, vars, dvars, ddvars, elemRes);

    // Add the residual values
    residual->setValues(len, nodes, elemRes, TACS_ADD_VALUES);
  }

  
  // Finish transmitting the residual
  residual->beginSetValues(TACS_ADD_VALUES);
  residual->endSetValues(TACS_ADD_VALUES);

  // Apply the boundary conditions for the residual
  residual->applyBCs(bcMap, varsVec);
}

/**
  Assemble the Jacobian matrix

  This function assembles the global Jacobian matrix and
  residual. This Jacobian includes the contributions from all
  elements. The Dirichlet boundary conditions are applied to the
  matrix by zeroing the rows of the matrix associated with a boundary
  condition, and setting the diagonal to unity. The matrix assembly
  also performs any communication required so that the matrix can be
  used immediately after assembly.

  @param alpha Coefficient for the variables
  @param beta Coefficient for the time-derivative terms
  @param gamma Coefficientfor the second time derivative term
  @param residual The residual of the governing equations
  @param A The Jacobian matrix
  @param matOr the matrix orientation NORMAL or TRANSPOSE
  @param lambda Scaling factor for the aux element contributions, by default 1
*/
void TACSAssembler::assembleJacobian(TacsScalar alpha, TacsScalar beta,
                                     TacsScalar gamma, TacsScalar *residual,
                                     TacsMat *A, MatrixOrientation matOr,
                                     const TacsScalar lambda) {
  // Zero the residual and the matrix
  if (residual) {
    residual->zeroEntries();
  }
  A->zeroEntries();

  // Retrieve pointers to temporary storage
  TacsScalar *vars, *dvars, *ddvars, *elemRes, *elemXpts;
  TacsScalar *elemWeights, *elemMat;
  getDataPointers(elementData, &vars, &dvars, &ddvars, &elemRes, &elemXpts,
                  NULL, &elemWeights, &elemMat);

  for (int i = 0; i < numElements; i++) {
    int ptr = elementNodeIndex[i];
    int len = elementNodeIndex[i + 1] - ptr;
    const int *nodes = &elementTacsNodes[ptr];
    xptVec->getValues(len, nodes, elemXpts);
    varsVec->getValues(len, nodes, vars);
    dvarsVec->getValues(len, nodes, dvars);
    ddvarsVec->getValues(len, nodes, ddvars);

    // Get the number of variables from the element
    int nvars = elements[i]->getNumVariables();

    // Compute and add the contributions to the Jacobian
    memset(elemRes, 0, nvars * sizeof(TacsScalar));
    memset(elemMat, 0, nvars * nvars * sizeof(TacsScalar));
    elements[i]->addJacobian(i, time, alpha, beta, gamma, elemXpts, vars,
                              dvars, ddvars, elemRes, elemMat);

    if (residual) {
      residual->setValues(len, nodes, elemRes, TACS_ADD_VALUES);
    }
    addMatValues(A, i, elemMat, elementIData, elemWeights, matOr);
  }

  // Do any matrix and residual assembly if required
  A->beginAssembly();
  if (residual) {
    residual->beginSetValues(TACS_ADD_VALUES);
  }

  A->endAssembly();
  if (residual) {
    residual->endSetValues(TACS_ADD_VALUES);
  }

  // Apply the boundary conditions
  if (residual) {
    residual->applyBCs(bcMap, varsVec);
  }

  // Apply the appropriate boundary conditions
  A->applyBCs(bcMap);
}
