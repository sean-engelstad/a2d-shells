#include "../therm-cylinder-include/getKDF.h"

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
    double temperature = 1.0; // K (may have to adjust depending on the)
    double E = 70e5; // 70e9 // can scale the problem (won't affect disps)
    double conv_eigval = 0.3; // 0.01, 0.1
    double conv_slope_frac = 0.1;
    // worried that if I scale down too much, won't solve as deeply though.

    // for debugging
    // TacsScalar imperfections[NUM_IMP] = {0.0 * t, 0.0, 0.0 };
    // TacsScalar imperfections[NUM_IMP] = {0.5 * t, 0.0, 0.0 };
    TacsScalar imperfections[NUM_IMP] = {0.0, 0.0, 0.5 * t };
    bool useEigvals = false; // use load-disp curve for thermal
    int nelems = 20000;
    std::string filePrefix = "";
    TacsScalar nasaKDF, tacsKDF;
    bool ringStiffened = false;
    double ringStiffenedRadiusFrac = 0.9;

    getNonlinearBucklingKDF(
        comm, 1, filePrefix, t, rt, Lr, E, temperature, 
        conv_eigval, conv_slope_frac,
        ringStiffened, ringStiffenedRadiusFrac,
        NUM_IMP, &imperfections[0],
        useEigvals, nelems, &nasaKDF, &tacsKDF
    );

    MPI_Finalize();

    return 0;
}