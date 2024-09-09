#include "TACSMeshLoader.h"
// #include "LoadMesh.h"

int main() {
    // build the TACS mesh loader and scan the uCRM BDF file
    TACSMeshLoader myLoader{};
    // myLoader = new TACSMeshLoader();
    myLoader.scanBDFFile("CRM_box_2nd.bdf");
    // myLoader.scanBDFFile("plate.bdf");
}