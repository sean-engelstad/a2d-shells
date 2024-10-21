#include "../mech-cylinder-include/getKDF.h"

// this is a nonlinear buckling example of a cylinder under mechanical loading
// with applied geometric imperfections. The load factor for nonlinear buckling is determined automatically.
// and the KDF (ratio of NL load factor / Linear load factor for buckling) or knockdown factor is computed and saved.


int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    // Get the rank
    MPI_Comm comm = MPI_COMM_WORLD;

    const int NRUNS = 7;
    double rtVals[NRUNS] = {1000.0, 500.0, 300.0, 100.0, 50.0, 25.0, 10.0};
    int meshSizes[NRUNS] = {40000, 20000, 10000, 10000, 10000, 10000, 10000};
    TacsScalar nasaKDF[NRUNS] = 0.0;
    TacsScalar tacsKDF[NRUNS] = 0.0;

    double t = 0.002;
    double Lr = 2.0;
    const int NUM_IMP = 3;
    // for debugging
    // TacsScalar imperfections[NUM_IMP] = {0.0 * t, 0.0, 0.0 };
    TacsScalar imperfections[NUM_IMP] = {0.5 * t, 0.0, 0.0 };
    bool useEigvals = true;
    

    // run each KDF simulation for mechanical nonlinear buckling
    for (int irun = NRUNS-1; irun >= 0; irun--) {
        double rt = rtVals[irun];
        int nelems = meshSizes[irun]; // 5000, 10000
        getNonlinearBucklingKDF(
            comm, 1, t, rt, Lr, NUM_IMP, &imperfections[0],
            useEigvals, nelems, &nasaKDF[irun], &tacsKDF[irun]
        );
    }   

    MPI_Finalize();

    return 0;
}