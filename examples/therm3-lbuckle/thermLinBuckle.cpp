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

// this is a nonlinear buckling example of a cylinder under mechanical loading
// with applied geometric imperfections. The load factor for nonlinear buckling is determined automatically.
// and the KDF (ratio of NL load factor / Linear load factor for buckling) or knockdown factor is computed and saved.

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    // Get the rank
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);

    // Parameters optionally set from the command line
    int order = 2;

    double t = 0.002; // m 
    double Lr = 2.0; // default 2.0
    double rt = 100; // 100, 50, 25
    double R = t * rt; // m
    double L = R * Lr;

    double udisp = 0.0; // no compressive disp (compressive strain from heating only)

    // select nelems and it will select to retain isotropic elements (good element AR)
    // want dy = 2 * pi * R / ny the hoop elem spacing to be equal dx = L / nx the axial elem spacing
    // and want to choose # elems so that elements have good elem AR
    int nelems = 20000; // prev 3500 // target (does round stuff)
    double pi = 3.14159265;
    double A = L / 2.0 / pi / R;
    double temp1 = sqrt(nelems * 1.0 / A);
    int ny = (int)temp1;
    double temp2 = A * ny;
    int nx = (int)temp2;
    printf("nx = %d, ny = %d\n", nx, ny);

    TacsScalar rho = 2700.0;
    TacsScalar specific_heat = 921.096;
    TacsScalar E = 70e3; // 70e9
    TacsScalar nu = 0.3;
    TacsScalar ys = 270.0;
    TacsScalar cte = 23.5e-6;
    TacsScalar kappa = 230.0;
    TACSMaterialProperties *props = new TACSMaterialProperties(rho, specific_heat, E, nu, ys, cte, kappa);

    // TacsScalar axis[] = {1.0, 0.0, 0.0};
    // TACSShellTransform *transform = new TACSShellRefAxisTransform(axis);
    TACSShellTransform *transform = new TACSShellNaturalTransform();
    TACSShellConstitutive *con = new TACSIsoShellConstitutive(props, t);

    TACSAssembler *assembler = NULL;
    TACSCreator *creator = NULL;
    TACSElement *shell = NULL;

    // why do I get such different answers with the nonlinear shell than the linear one?

    // toggle on and off nonlinear vs linear element (1 vs 0)
    #define NONLINEAR_ELEM 1
    #if NONLINEAR_ELEM
        shell = new TACSQuad4NonlinearShell(transform, con); 
    #else
        shell = new TACSQuad4Shell(transform, con); 
    #endif
    shell->incref();

    createAssembler(comm, order, nx, ny, udisp, L, R, shell, &assembler, &creator);

    // set the temperatures into the structure
    TacsScalar temperature = 1.0; // default 1.0 // 1 deg K
    int numElements = assembler->getNumElements();

    #if NONLINEAR_ELEM
        TACSQuad4NonlinearShell *elem;
        for (int ielem = 0; ielem < numElements; ielem++) {
            elem = dynamic_cast<TACSQuad4NonlinearShell *>(assembler->getElement(ielem));
            elem->setTemperature(temperature);
            #ifdef TACS_USE_COMPLEX
                elem->setComplexStepGmatrix(true);
            #endif
        }
    #else
        TACSQuad4Shell *elem;
        for (int ielem = 0; ielem < numElements; ielem++) {
            elem = dynamic_cast<TACSQuad4Shell *>(assembler->getElement(ielem));
            elem->setTemperature(temperature);
            #ifdef TACS_USE_COMPLEX
                elem->setComplexStepGmatrix(true);
            #endif
        }
    #endif
    

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
    TacsScalar sigma = 5.0; // need high enough num_eigvals to get it right
    int max_lanczos_vecs = 500, num_eigvals = 100; // num_eigvals = 50;
    double eig_tol = 1e-12;

    TACSLinearBuckling *buckling = new TACSLinearBuckling(assembler, sigma,
                     gmat, kmat, aux_mat, lbuckle_gmres, max_lanczos_vecs, num_eigvals, eig_tol);
    buckling->incref();

    // make a KSM print object for solving buckling
    KSMPrint *ksm_print_buckling = new KSMPrintStdout("BucklingAnalysis", 0, 10);
    ksm_print->incref();

    // solve the buckling analysis
    buckling->setSigma(10.0);
    buckling->solve(NULL, NULL, ksm_print_buckling);
    // exit(0);

    // Create an TACSToFH5 object for writing output to files
    int write_flag = (TACS_OUTPUT_CONNECTIVITY | TACS_OUTPUT_NODES |
                        TACS_OUTPUT_DISPLACEMENTS | TACS_OUTPUT_STRAINS |
                        TACS_OUTPUT_STRESSES | TACS_OUTPUT_EXTRAS);
    TACSToFH5 *f5 =
        new TACSToFH5(assembler, TACS_BEAM_OR_SHELL_ELEMENT, write_flag);
    f5->incref();

    // write each of the buckling modes to a file
    
    TACSBVec *phi = assembler->createVec();
    phi->incref();
    TacsScalar error;
    for (int imode = 0; imode < num_eigvals; imode++) {
        buckling->extractEigenvector(imode, phi, &error);
        assembler->setVariables(phi);   
        std::string filename = "_buckling/therm-buckle" + std::to_string(imode) + ".f5";
        const char *cstr_filename = filename.c_str();
        f5->writeToFile(cstr_filename);
    }    
    

    MPI_Finalize();

    return 0;
}
