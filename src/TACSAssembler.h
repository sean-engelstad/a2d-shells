/*
@Description : Shortened version just for linear static analysis 
               (no optimization) and for serial now
@Author : Sean Engelstad
@Date : Sep 10, 2024
*/

#pragma once
#include "TACSObject.h"

class TACSAssembler : public TACSObject {
public:
    TACSAssembler(int _varsPerNode, int _numElements);
    ~TACSAssembler();

private:
    int varsPerNode;         // number of variables per node
    int numElements;         // number of elements
    int numNodes;            // number of nodes referenced by this process
};