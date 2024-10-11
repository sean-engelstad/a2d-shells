#include "TACSMeshLoader.h"
#include "TACSAssembler.h"

// dependencies to make element, constitutive objects
#include "TACSShellElementTransform.h"
#include "TACSMaterialProperties.h"
#include "TACSIsoShellConstitutive.h"
#include "TACSShellElementDefs.h"
#include "TACSBuckling.h"
#include "KSM.h"
#include "createCylinderDispControl.h"
#include "TACSContinuation.h"

// this example is based off of examples/crm/crm.cpp in TACS

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    // Get the rank
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);

    // Parameters optionally set from the command line
    int order = 2;

    double t = 0.002; // m 
    double L = 0.4; // m
    double R = 0.2; // m
    double udisp = -1e-5; // compressive udisp

    // select nelems and it will select to retain isotropic elements (good element AR)
    // want dy = 2 * pi * R / ny the hoop elem spacing to be equal dx = L / nx the axial elem spacing
    // and want to choose # elems so that elements have good elem AR
    int nelems = 3500; // target (does round stuff)
    double pi = 3.14159265;
    double A = L / 2.0 / pi / R;
    double temp1 = sqrt(nelems * 1.0 / A);
    int ny = (int)temp1;
    double temp2 = A * ny;
    int nx = (int)temp2;
    printf("nx = %d, ny = %d\n", nx, ny);

    TacsScalar rho = 2700.0;
    TacsScalar specific_heat = 921.096;
    TacsScalar E = 70e3;
    TacsScalar nu = 0.3;
    TacsScalar ys = 270.0;
    TacsScalar cte = 24.0e-6;
    TacsScalar kappa = 230.0;
    TACSMaterialProperties *props =
        new TACSMaterialProperties(rho, specific_heat, E, nu, ys, cte, kappa);

    // TacsScalar axis[] = {1.0, 0.0, 0.0};
    // TACSShellTransform *transform = new TACSShellRefAxisTransform(axis);
    TACSShellTransform *transform = new TACSShellNaturalTransform();

    TACSShellConstitutive *con = new TACSIsoShellConstitutive(props, t);

    TACSAssembler *assembler = NULL;
    TACSCreator *creator = NULL;
    TACSElement *shell = NULL;
    // needs to be nonlinear here otherwise solve will terminate immediately
    shell = new TACSQuad4NonlinearShell(transform, con); 
    shell->incref();
    createAssembler(comm, order, nx, ny, udisp, L, R, shell, &assembler, &creator);

    // Solve the linear static analysis
    if (rank == 0) {
        printf("Solving linear system..\n");
    }
    
    // Create the design vector
    TACSBVec *x = assembler->createDesignVec();
    x->incref();

    // Get the design variable values
    assembler->getDesignVars(x);

    // Create matrix and vectors
    TACSBVec *u0 = assembler->createVec();  // displacements and rotations
    TACSBVec *f = assembler->createVec();    // loads
    u0->incref();
    f->incref();

    // create zero loads
    TacsScalar *force_vals;
    int size = f->getArray(&force_vals);
    memset(force_vals, 0.0, size * sizeof(TacsScalar));
    assembler->applyBCs(f);

    // nonlinear static
    // --------------------------------------------
    TACSSchurMat *kmat = assembler->createSchurMat();  // stiffness matrix
    kmat->incref();

    // Allocate the factorization
    int lev = 1e6;
    double fill = 10.0;
    int reorder_schur = 1;
    TACSSchurPc *pc = new TACSSchurPc(kmat, lev, fill, reorder_schur);
    pc->incref();

    int subspaceSize = 10, nrestarts = 15, isFlexible = 0;
    GMRES *solver = new GMRES(kmat, pc, subspaceSize, nrestarts, isFlexible);
    solver->incref();

    // set relative tolerances
    solver->setTolerances(1e-12, 1e-12);

    // make a KSM print object for solving buckling
    KSMPrint *ksm_print = new KSMPrintStdout("NonlinearStatic", 0, 10);
    ksm_print->incref();

    // inputs to the TACSContinuation solver    
    int max_continuation_iters = 300; // for prelim Newton solve to lambda_init
    int max_correction_iters = 100; // for the arc length method static regime
    int max_correction_restarts = 2; // restart for nonlinear static regime
    double corr_rtol = 1e-3; // this needs improvement and an atol (since already fairly low)
    double corr_dtol = 1e3; // if divergence tolerance is huge it's failing to solve and breaks out prelim newton loop
    double krylov_rtol = 1e-6; // krylov refers to the second solve where it reaches more severe loss of stability
    double krylov_atol = 1e-10; 
    double tangent_rtol = 1e-6; // for prelim newton solve section
    double tangent_atol = 1e-10;
    TACSContinuationCallback *callback = new TACSContinuationCallback(comm, "nlstatic.out");

    // make the TACSContinuation solver for nonlinear static analysis
    TACSContinuation *continuation = new TACSContinuation(
        assembler, max_continuation_iters, max_correction_iters, max_correction_restarts,
        corr_rtol, corr_dtol, krylov_rtol, krylov_atol, tangent_rtol, tangent_atol
    );
    continuation->incref();

    // try solving the nonlinear static analysis
    // looks like it doesn't displacement control here (need to adjust that and the BCs in the TACSContinuation.cpp)
    // double lambda_init = 280.0; // is this the target final load factor?
    // this is with the perfect geometry so expect buckling to occur at ~280 lambda
    // it's slowing down a lot near buckling right now => how to resolve that..
    double lambda_init = 280.0 * 0.5; // start at 50% the predicted linear buckling load (change this depending on case)
    double target_delta_lambda = 5.0;
    // right now based on 100 steps and it is canging lambda by about 4-5 each time, we reach lambd aof 
    printf("solve tangent::\n");
    continuation->solve_tangent(
        kmat, pc, solver, f, lambda_init, target_delta_lambda,
        ksm_print, callback
    );

    // TODO : need to add f5 writeouts from the nonlinear static solution
    // and/or final lambda values stop criterion
    // also the 

    // Create an TACSToFH5 object for writing output to files
    int write_flag = (TACS_OUTPUT_CONNECTIVITY | TACS_OUTPUT_NODES |
                        TACS_OUTPUT_DISPLACEMENTS | TACS_OUTPUT_STRAINS |
                        TACS_OUTPUT_STRESSES | TACS_OUTPUT_EXTRAS);
    TACSToFH5 *f5 =
        new TACSToFH5(assembler, TACS_BEAM_OR_SHELL_ELEMENT, write_flag);
    f5->incref();
    // copy solution out of nonlinear static and write to file
    // already in the assembler
    f5->writeToFile("mech-nl-cylinder.f5");

    MPI_Finalize();

    return 0;
}