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

    // There are always 3 coordinates (even for 2D problems)
    static const int TACS_SPATIAL_DIM = 3;

    // Evaluate kinetic and potential energy
    // -----------------------------------------
    void evalEnergies(TacsScalar *Te, TacsScalar *Pe);

    // Return important information about the TACSAssembler object
    // -----------------------------------------------------------
    int getVarsPerNode() {return varsPerNode;};
    int getNumNodes() {return numNodes;};
    int getNumElements() {return numElements;};
    int getNumVariables() {return numNodes * varsPerNode;};
    TacsScalar *getNodeMap();
    TacsScalar *getDesignNodeMap();
    TacsScalar *getBcMap();
    TacsScalar *getInitBcMap();

    // set elements and conn
    // -----------------------------------------------
    int setElementConnectivity(const int *ptr, const int *conn);
    int setElements(TACSElement **_elements);

    // Return elements and node numbers
    // --------------------------------
    TACSElement **getElements();
    TACSElement *getElement(int elem, TacsScalar *Xpts = NULL,
                            TacsScalar *vars = NULL, TacsScalar *dvars = NULL,
                            TacsScalar *ddvars = NULL);
    TACSElement *getElement(int elem, int *len, const int **nodes);

     // Get pointers to the start-locations within the data array
    // ---------------------------------------------------------
    void getDataPointers(TacsScalar *data, TacsScalar **v1, TacsScalar **v2,
                        TacsScalar **v3, TacsScalar **v4, TacsScalar **x1,
                        TacsScalar **x2, TacsScalar **weights, TacsScalar **mat);

    // Residual and Jacobian assembly
    // ------------------------------
    void assembleRes(TacsScalar *residual, const TacsScalar lambda = 1.0);
    void assembleJacobian(TacsScalar alpha, TacsScalar beta, TacsScalar gamma,
                            TacsScalar *residual, TacsScalar *A,
                            const TacsScalar lambda = 1.0);

private:
    int varsPerNode;         // number of variables per node
    int numElements;         // number of elements
    int numNodes;            // number of nodes referenced by this process

    // The local list of elements
    TACSElement **elements;
    TacsScalar *elementData;
    int *elementIData;

    // variables/elements have been initialized
    int meshInitializedFlag;

    // Maximum element information
    int numMultiplierNodes;
    int maxElementNodes;       // maximum number of ind. and dep. element nodes
    int maxElementSize;        // maximum number of variables for any element
    int maxElementIndepNodes;  // maximum number of independent nodes


    // save the time
    double time;

    // Variables that define the CSR data structure to
    // store the element -> node information
    int *elementNodeIndex;
    int *elementTacsNodes;

    // other important data
    TacsScalar *vars0, *dvars0, *ddvars0;
    TacsScalar *res, *mat;
    TacsScalar *A;
};