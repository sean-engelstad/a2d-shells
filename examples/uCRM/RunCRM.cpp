#include "TACSMeshLoader.h"
#include "TACSAssembler.h"

int main() {
    // build the TACS mesh loader and scan the uCRM BDF file
    TACSMeshLoader myLoader{};
    myLoader.scanBDFFile("CRM_box_2nd.bdf");

    int varsPerNode = 6;
    int numElements = 1;
    TACSAssembler assembler{varsPerNode, numElements}; 
}