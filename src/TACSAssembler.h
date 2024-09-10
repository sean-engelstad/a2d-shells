/*
@Description : Shortened version just for linear static analysis 
               (no optimization) and for serial now
@Author : Sean Engelstad
@Date : Sep 10, 2024
*/

#pragma once
#include "TACSObject.h"
#include "TACSElement.h"

class TACSAssembler : public TACSObject {
public:
    TACSAssembler(int _varsPerNode, int _numElements);
    ~TACSAssembler();

    // Evaluate kinetic and potential energy
    // -----------------------------------------
    void evalEnergies(TacsScalar *Te, TacsScalar *Pe);

    // Return elements and node numbers
    // --------------------------------
    TACSElement **getElements();
    TACSElement *getElement(int elem, TacsScalar *Xpts = NULL,
                            TacsScalar *vars = NULL, TacsScalar *dvars = NULL,
                            TacsScalar *ddvars = NULL);
    TACSElement *getElement(int elem, int *len, const int **nodes);

    // Residual and Jacobian assembly
    // ------------------------------
    void assembleRes(TacsScalar *residual, const TacsScalar lambda = 1.0);
    void assembleJacobian(TacsScalar alpha, TacsScalar beta, TacsScalar gamma,
                            TacsScalar *residual, TACSMat *A,
                            MatrixOrientation matOr = TACS_MAT_NORMAL,
                            const TacsScalar lambda = 1.0);

private:
    int varsPerNode;         // number of variables per node
    int numElements;         // number of elements
    int numNodes;            // number of nodes referenced by this process

    // The local list of elements
    TACSElement **elements;
    TACSScalar *elementData;
    int *elementIData;

    // other important data
    TacsScalar *vars0, *dvars0, *ddvars0;
    TacsScalar *res, *mat;
    TACSMat *A;
};