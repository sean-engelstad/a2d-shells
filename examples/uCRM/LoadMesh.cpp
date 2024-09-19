#include "TACSMeshLoader.h"

int main(int argc, char **argv) {
  // Intialize MPI and declare communicator
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;

  // build the TACS mesh loader and scan the uCRM BDF file
  TACSMeshLoader myLoader{comm};
  myLoader.scanBDFFile("CRM_box_2nd.bdf");

  MPI_Finalize();
  return 0;
}