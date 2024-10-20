#include "TACSMeshLoader.h"
#include "TACSAssembler.h"

// dependencies to make element, constitutive objects
#include "TACSShellElementTransform.h"
#include "TACSMaterialProperties.h"
#include "TACSIsoShellConstitutive.h"
#include "TACSShellElementDefs.h"
#include "TACSBuckling.h"
#include "KSM.h"
#include "TACSContinuation.h"

#include "createCylinderDispControl.h"
#include "getKDF.h"

// this is a nonlinear buckling example of a cylinder under mechanical loading
// with applied geometric imperfections. The load factor for nonlinear buckling is determined automatically.
// and the KDF (ratio of NL load factor / Linear load factor for buckling) or knockdown factor is computed and saved.


int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    // Get the rank
    MPI_Comm comm = MPI_COMM_WORLD;

    double rtVals[6] = {500.0, 300.0, 100.0, 50.0, 25.0, 10.0};
    int meshSizes[6] = {20000, 10000, 10000, 10000, 10000, 10000};
    
    // run each KDF simulation for mechanical nonlinear buckling
    for (int irun = 0; irun < 6; irun++) {
        double rt = rtVals[irun];
        double Lr = 2.0;
        int nelems = meshSizes[irun]; // 5000, 10000
        getNonlinearBucklingKDF(comm, irun+1, rt, Lr, nelems);
    }   

    MPI_Finalize();

    return 0;
}