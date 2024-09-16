/*
@Description : Shortened version just for linear static analysis 
               (no optimization) and for serial now
@Author : Sean Engelstad
@Date : Sep 10, 2024
*/


#include "TACSAssembler.h"
#include <iostream>

TACSAssembler::TACSAssembler(int _varsPerNode, int _numElements) : varsPerNode(_varsPerNode), numElements(_numElements) 
{
  meshInitializedFlag = 0;
  elements = NULL;
  time = 0.0;

  // These values will be computed later
  numMultiplierNodes = 0;
  numNodes = 0;
  maxElementNodes = 0;
  maxElementSize = 0;
  elementNodeIndex = NULL;
  elementTacsNodes = NULL;
};
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
    // xptVec->getValues(len, nodes, elemXpts);
    // varsVec->getValues(len, nodes, vars);
    // dvarsVec->getValues(len, nodes, dvars);

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
void TACSAssembler::assembleRes(TacsScalar *residual, const TacsScalar lambda) {
  // Sort the list of auxiliary elements - this only performs the
  // sort if it is required (if new elements are added)

  // Zero the residual
  // residual->zeroEntries();
  memset(residual, 0.0, getNumVariables()*sizeof(TacsScalar));

  // Retrieve pointers to temporary storage
  TacsScalar *vars, *dvars, *ddvars, *elemRes, *elemXpts;
  getDataPointers(elementData, &vars, &dvars, &ddvars, &elemRes, &elemXpts,
                  NULL, NULL, NULL);

  // Go through and add the residuals from all the elements
  for (int i = 0; i < numElements; i++) {
    int ptr = elementNodeIndex[i];
    int len = elementNodeIndex[i + 1] - ptr;
    const int *nodes = &elementTacsNodes[ptr];
    // xptVec->getValues(len, nodes, elemXpts);
    // varsVec->getValues(len, nodes, vars);
    // dvarsVec->getValues(len, nodes, dvars);
    // ddvarsVec->getValues(len, nodes, ddvars);

    // Add the residual from the working element
    int nvars = elements[i]->getNumVariables();
    memset(elemRes, 0, nvars * sizeof(TacsScalar));
    elements[i]->addResidual(i, time, elemXpts, vars, dvars, ddvars, elemRes);

    // Add the residual values
    // residual->setValues(len, nodes, elemRes, TACS_ADD_VALUES);
  }

  
  // Finish transmitting the residual
  // residual->beginSetValues(TACS_ADD_VALUES);
  // residual->endSetValues(TACS_ADD_VALUES);

  // Apply the boundary conditions for the residual
  // residual->applyBCs(bcMap, varsVec);
}

/**
  Set the element connectivity.

  Note that the number of elements that are set at this stage must be
  consistent with the number of elements passed in when TACSAssembler
  is created.  (The connectivity arrays are copied and should be freed
  by the caller.)

  @param ptr Offset into the connectivity array for this processor
  @param conn The connectivity from elements to global node index
  @return Fail flag indicating if a failure occured
*/
int TACSAssembler::setElementConnectivity(const int *ptr, const int *conn) {
  // if (meshInitializedFlag) {
  //   fprintf(stderr,
  //           "Cannot call setElementConnectivity() after initialize()\n");
  //   return 1;
  // }
  // if (tacsExtNodeNums) {
  //   fprintf(stderr,
  //           "Cannot call setElementConnectivity() after reordering\n");
  //   return 1;
  // }

  // Free the data if it has already been set before
  if (elementTacsNodes) {
    delete[] elementTacsNodes;
  }
  if (elementNodeIndex) {
    delete[] elementNodeIndex;
  }

  int size = ptr[numElements];
  elementNodeIndex = new int[numElements + 1];
  elementTacsNodes = new int[size];
  memcpy(elementNodeIndex, ptr, (numElements + 1) * sizeof(int));

  // Copy the node numbers
  memcpy(elementTacsNodes, conn, size * sizeof(int));

  // If the elements are set, check the connectivity
  if (elements) {
    for (int i = 0; i < numElements; i++) {
      int size = elementNodeIndex[i + 1] - elementNodeIndex[i];
      if (size != elements[i]->getNumNodes()) {
        fprintf(stderr,
                "TACSAssembler: Element %s does not match "
                " connectivity\n", elements[i]->getObjectName());
        return 1;
      }
    }
  }

  return 0;
}

/**
  Set the element array within TACS.

  The number of element pointers provided should be equal to the
  number of elements specified in the TACSAssembler constructor.  Note
  that the elements themselves are not copied, just the pointers to
  them.

  @param elements The array of element pointers with length = numElements.
  @return Fail flag indicating if a failure occured
*/
int TACSAssembler::setElements(TACSElement **_elements) {
  printf("Inside setElements()\n");
  std::cout << std::flush;

  if (meshInitializedFlag) {
    fprintf(stderr, "Cannot call setElements() after initialize()\n");
    return 1;
  }

  printf("Checkpt1\n");
  std::cout << std::flush;

  // Check if the number of variables per node matches
  for (int i = 0; i < numElements; i++) {
    if (_elements[i]->getVarsPerNode() != varsPerNode) {
      fprintf(stderr,
              "TACSAssembler: Element %s, num displacements %d "
              "does not match variables per node %d\n",
              _elements[i]->getObjectName(),
              _elements[i]->getVarsPerNode(), varsPerNode);
      return 1;
    }
  }

  // Copy over the element pointers into a local array
  printf("Pre copy elements\n");
  std::cout << std::flush;
  if (elements) {
    for (int i = 0; i < numElements; i++) {
      _elements[i]->incref();
      if (elements[i]) {
        elements[i]->decref();
      }
      elements[i] = _elements[i];
    }
  } else {
    elements = new TACSElement *[numElements];
    memset(elements, 0, numElements * sizeof(TACSElement *));
    for (int i = 0; i < numElements; i++) {
      _elements[i]->incref();
      elements[i] = _elements[i];
    }
  }

  // Keep track of the number of multiplier nodes
  numMultiplierNodes = 0;

  // Determine the maximum number of nodes per element
  maxElementSize = 0;
  maxElementNodes = 0;

  for (int i = 0; i < numElements; i++) {
    // Determine if the maximum number of variables and nodes needs to
    // be changed
    int elemSize = elements[i]->getNumVariables();
    if (elemSize > maxElementSize) {
      maxElementSize = elemSize;
    }

    elemSize = elements[i]->getNumNodes();
    if (elemSize > maxElementNodes) {
      maxElementNodes = elemSize;
    }

    // Count up the multiplier nodes
    int multiplier = elements[i]->getMultiplierIndex();
    if (multiplier > 0) {
      numMultiplierNodes++;
    }

  }

  // If the connectivity is set, determine if it is consistent
  if (elementNodeIndex) {
    for (int i = 0; i < numElements; i++) {
      int size = elementNodeIndex[i + 1] - elementNodeIndex[i];
      if (size != elements[i]->getNumNodes()) {
        fprintf(stderr,
                "TACSAssembler: Element %s does not match "
                "connectivity\n", elements[i]->getObjectName());
        return 1;
      }
    }
  }

  return 0;
}

/**
  Get pointers to the element data. This code provides a way to
  automatically segment an array to avoid coding mistakes.

  Note that this is coded in such a way that you can provide NULL
  arguments
*/
void TACSAssembler::getDataPointers(TacsScalar *data, TacsScalar **v1,
                                    TacsScalar **v2, TacsScalar **v3,
                                    TacsScalar **v4, TacsScalar **x1,
                                    TacsScalar **x2, TacsScalar **weights,
                                    TacsScalar **mat) {
  int s = 0;
  if (v1) {
    *v1 = &data[s];
    s += maxElementSize;
  }
  if (v2) {
    *v2 = &data[s];
    s += maxElementSize;
  }
  if (v3) {
    *v3 = &data[s];
    s += maxElementSize;
  }
  if (v4) {
    *v4 = &data[s];
    s += maxElementSize;
  }
  if (x1) {
    *x1 = &data[s];
    s += TACS_SPATIAL_DIM * maxElementNodes;
  };
  if (x2) {
    *x2 = &data[s];
    s += TACS_SPATIAL_DIM * maxElementNodes;
  };
  if (weights) {
    *weights = &data[s];
    s += maxElementIndepNodes;
  }
  if (mat) {
    *mat = &data[s];
  }
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
                                     TacsScalar *A, const TacsScalar lambda) {
  // Zero the residual and the matrix
  if (residual) {
    // residual->zeroEntries();
    memset(residual, 0.0, getNumVariables()*sizeof(TacsScalar));
  }
  // A->zeroEntries();
  memset(A, 0.0, getNumVariables()*sizeof(TacsScalar));

  // Retrieve pointers to temporary storage
  TacsScalar *vars, *dvars, *ddvars, *elemRes, *elemXpts;
  TacsScalar *elemWeights, *elemMat;
  getDataPointers(elementData, &vars, &dvars, &ddvars, &elemRes, &elemXpts,
                  NULL, &elemWeights, &elemMat);

  for (int i = 0; i < numElements; i++) {
    int ptr = elementNodeIndex[i];
    int len = elementNodeIndex[i + 1] - ptr;
    const int *nodes = &elementTacsNodes[ptr];
    // xptVec->getValues(len, nodes, elemXpts);
    // varsVec->getValues(len, nodes, vars);
    // dvarsVec->getValues(len, nodes, dvars);
    // ddvarsVec->getValues(len, nodes, ddvars);

    // Get the number of variables from the element
    int nvars = elements[i]->getNumVariables();

    // Compute and add the contributions to the Jacobian
    memset(elemRes, 0, nvars * sizeof(TacsScalar));
    memset(elemMat, 0, nvars * nvars * sizeof(TacsScalar));
    elements[i]->addJacobian(i, time, alpha, beta, gamma, elemXpts, vars,
                              dvars, ddvars, elemRes, elemMat);

    // if (residual) {
    //   residual->setValues(len, nodes, elemRes, TACS_ADD_VALUES);
    // }
    // addMatValues(A, i, elemMat, elementIData, elemWeights, matOr);
  }

  // Do any matrix and residual assembly if required
  // A->beginAssembly();
  // if (residual) {
  //   residual->beginSetValues(TACS_ADD_VALUES);
  // }

  // A->endAssembly();
  // if (residual) {
  //   residual->endSetValues(TACS_ADD_VALUES);
  // }

  // Apply the boundary conditions
  // if (residual) {
  //   residual->applyBCs(bcMap, varsVec);
  // }

  // Apply the appropriate boundary conditions
  // A->applyBCs(bcMap);
}
