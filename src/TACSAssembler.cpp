/*
@Description : Shortened version just for linear static analysis 
               (no optimization) and for serial now
@Author : Sean Engelstad
@Date : Sep 10, 2024
*/


#include "TACSAssembler.h"

TACSAssembler::TACSAssembler(int _varsPerNode, int _numElements) : varsPerNode(_varsPerNode), numElements(_numElements) {};
TACSAssembler::~TACSAssembler() {};