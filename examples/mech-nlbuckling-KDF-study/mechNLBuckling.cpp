#include "getKDF.h"

// this is a nonlinear buckling example of a cylinder under mechanical loading
// with applied geometric imperfections. The load factor for nonlinear buckling is determined automatically.
// and the KDF (ratio of NL load factor / Linear load factor for buckling) or knockdown factor is computed and saved.



int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    // Get the rank
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);

    double rtVals[7] = {1000.0, 500.0, 300.0, 100.0, 50.0, 25.0, 10.0};
    
    // run each KDF simulation for mechanical nonlinear buckling
    for (int irun = 0; irun < 7; irun++) {
        double rt = rtVals[irun];
        double Lr = 2.0;
        int nelems = 8000; // 5000, 
        getNonlinearBucklingKDF(comm, irun+1, rt, 2.0, nelems);
    }   

    MPI_Finalize();

    return 0;
}