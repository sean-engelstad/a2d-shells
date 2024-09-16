#include "TACSMeshLoader.h"
int main() {
    // build the TACS mesh loader and scan the uCRM BDF file
    TACSMeshLoader myLoader{};
    myLoader.scanBDFFile("CRM_box_2nd.bdf");
}