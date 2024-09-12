#include "TACSMeshLoader.h"
#include "TACSAssembler.h"

// this example is based off of examples/crm/crm.cpp in TACS

int main() {
    // build the TACS mesh loader and scan the uCRM BDF file
    TACSMeshLoader *mesh = new TACSMeshLoader();
    mesh->incref();
    mesh->scanBDFFile("CRM_box_2nd.bdf");

    // TODO : create elements and constitutive models

    TACSAssembler *assembler = mesh->createTACS(6);
}