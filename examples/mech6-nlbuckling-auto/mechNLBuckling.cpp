#include "../mech-cylinder-include/getKDF.h"

// this is a nonlinear buckling example of a cylinder under mechanical loading
// with applied geometric imperfections. The load factor for nonlinear buckling is determined automatically.
// and the KDF (ratio of NL load factor / Linear load factor for buckling) or knockdown factor is computed and saved.


int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    // Get the rank
    MPI_Comm comm = MPI_COMM_WORLD;

    double t = 0.002;
    double rt = 100;
    double Lr = 2.0;
    const int NUM_IMP = 3;
    // for debugging
    // TacsScalar imperfections[NUM_IMP] = {0.0 * t, 0.0, 0.0 };
    TacsScalar imperfections[NUM_IMP] = {0.5 * t, 0.0, 0.0 };
    bool useEigvals = true;
    int nelems = 10000;
    std::string filePrefix = "";
    TacsScalar nasaKDF, tacsKDF;

    getNonlinearBucklingKDF(
        comm, 1, filePrefix, t, rt, Lr, NUM_IMP, &imperfections[0],
        useEigvals, nelems, &nasaKDF, &tacsKDF
    );

    MPI_Finalize();

    return 0;
}