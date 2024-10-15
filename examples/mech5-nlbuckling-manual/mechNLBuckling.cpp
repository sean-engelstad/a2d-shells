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
    int nelems = 5000; // prev 3500 // target (does round stuff)
    double pi = 3.14159265;
    double A = L / 2.0 / pi / R;
    double temp1 = sqrt(nelems * 1.0 / A);
    int ny = (int)temp1;
    double temp2 = A * ny;
    int nx = (int)temp2;
    printf("nx = %d, ny = %d\n", nx, ny);

    TacsScalar rho = 2700.0;
    TacsScalar specific_heat = 921.096;
    TacsScalar E = 70e9;
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
    GMRES *gmres = new GMRES(kmat, pc, subspaceSize, nrestarts, isFlexible);
    gmres->incref();
    gmres->setTolerances(1e-12, 1e-12);

    // make a KSM print object for solving buckling
    KSMPrint *ksm_print = new KSMPrintStdout("NonlinearStatic", 0, 10);
    ksm_print->incref();

    // build the TACS linear buckling analysis object
    // which we will use to check for buckling in the nonlinear load-displacement curve (u,lambda)
    // also use the linear eigenmodes as geometric imperfections in the cylinder
    // ---------------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------------
    
    // create the matrices for buckling
    TACSSchurMat *gmat = assembler->createSchurMat();  // geometric stiffness matrix
    TACSSchurMat *aux_mat = assembler->createSchurMat();  // auxillary matrix for shift and invert solver

    // optional other preconditioner settings?
    assembler->assembleMatType(TACS_STIFFNESS_MATRIX, kmat);
    assembler->assembleMatType(TACS_GEOMETRIC_STIFFNESS_MATRIX, gmat);

    subspaceSize = 10;
    nrestarts = 15; 
    isFlexible = 0;
    GMRES *lbuckle_gmres = new GMRES(aux_mat, pc, subspaceSize, nrestarts, isFlexible);
    lbuckle_gmres->incref();
    lbuckle_gmres->setTolerances(1e-12, 1e-12);

    // make the buckling solver
    TacsScalar sigma = 10.0; // need high enough num_eigvals to get it right
    int max_lanczos_vecs = 300, num_eigvals = 50; // num_eigvals = 50;
    double eig_tol = 1e-12;

    TACSLinearBuckling *buckling = new TACSLinearBuckling(assembler, sigma,
                     gmat, kmat, aux_mat, lbuckle_gmres, max_lanczos_vecs, num_eigvals, eig_tol);
    buckling->incref();

    // make a KSM print object for solving buckling
    KSMPrint *ksm_print_buckling = new KSMPrintStdout("BucklingAnalysis", 0, 10);
    ksm_print->incref();

    // solve the buckling analysis
    TACSBVec *u0_temp = NULL, *f0_temp = NULL;
    buckling->solve(f0_temp, u0_temp, ksm_print_buckling);
    // exit(0);

    // choose imperfection sizes for the cylinder based on the cylinder thickness
    TacsScalar imperfection_sizes[3] = {0.5 * t, 0.1 * t, 0.05 * t}; // t is cylinder thickness here

    // apply the first few eigenmodes as geometric imperfections to the cylinder
    TACSBVec *phi = assembler->createVec();
    TACSBVec *xpts = assembler->createNodeVec();
    TACSBVec *phi_uvw = assembler->createNodeVec();
    assembler->getNodes(xpts);
    phi->incref();
    TacsScalar error;
    for (int imode = 0; imode < 3; imode++) {
        buckling->extractEigenvector(imode, phi, &error);
        // if (imode == 0) {
        //     assembler->setVariables(phi);
        // }

        // copy the phi for all 6 shell dof into phi_uvw
        // how to copy every 3 out of 6 values from 
        TacsScalar *phi_x, *phi_uvw_x;
        int varSize = phi->getArray(&phi_x);
        int nodeSize = phi_uvw->getArray(&phi_uvw_x);
        int ixpts = 0;
        for (int iphi = 0; iphi < varSize; iphi++) {
            int idof = iphi % 6;
            if (idof > 2) { // skip rotx, roty, rotz components of eigenvector
                continue;
            }
            // printf("iphi %d, idof %d\n", iphi, idof);
            // printf("phi_x at %d / %d\n", iphi, varSize);
            // printf("phi_uvw_x at %d / %d\n", ixpts, nodeSize);
            phi_uvw_x[ixpts] = phi_x[iphi];
            ixpts++;
        }
        xpts->axpy(imperfection_sizes[imode], phi_uvw); 
        // xpts->axpy(imperfection_sizes[imode] * 100.0, phi_uvw); 
    }
    assembler->setNodes(xpts);
    // previously wrote out to f5 here in order to see the geometric imperfection

    // end of TACS linear buckling analysis for the geometric imperfections
    // ---------------------------------------------------------------------------------------
    
    // start the nonlinear static analysis, with linear buckling used to check for nonlinear buckling
    // ---------------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------------

    // first just // inputs to the TACSContinuation solver    
    int max_continuation_iters = 100; // default 300 // for prelim Newton solve to lambda_init
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

    // solve the nonlinear static analysis
    double lambda_init = 280.0 * 0.2; // start at 50% the predicted linear buckling load (change this depending on case)
    double target_delta_lambda = 5.0;
    // right now based on 100 steps and it is canging lambda by about 4-5 each time, we reach lambd aof 
    printf("solve tangent::\n");
    continuation->solve_tangent(
        kmat, pc, gmres, f, lambda_init, target_delta_lambda,
        ksm_print, callback
    );


    // end of nonlinear static analysis for nonlinear buckling
    // ---------------------------------------------------------------------------------------

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