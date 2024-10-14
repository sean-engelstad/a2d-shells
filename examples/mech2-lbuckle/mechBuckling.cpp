#include "TACSMeshLoader.h"
#include "TACSAssembler.h"

// dependencies to make element, constitutive objects
#include "TACSShellElementTransform.h"
#include "TACSMaterialProperties.h"
#include "TACSIsoShellConstitutive.h"
#include "TACSShellElementDefs.h"
#include "TACSBuckling.h"
#include "KSM.h"
#include "createCylinder.h"

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
    shell = new TACSQuad4Shell(transform, con);
    shell->incref();
    createAssembler(comm, order, nx, ny, udisp, L, R, shell, &assembler, &creator);

    // set temperature into all elements for thermal buckling
    // set to 0 for mechanical buckling
    int numElements = assembler->getNumElements();
    TACSQuad4Shell *elem;
    for (int ielem = 0; ielem < numElements; ielem++) {
        elem = dynamic_cast<TACSQuad4Shell *>(assembler->getElement(ielem));
        elem->setTemperature(0.0);
    }

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
    TACSBVec *u0 = NULL;
    TACSBVec *f = NULL;
    // OR instead of solving linear static:
    //   set u0, f to NULL and it will solve it for you
    // u0 = NULL; f = NULL;

    // create the matrices for buckling
    TACSSchurMat *kmat = assembler->createSchurMat();  // stiffness matrix
    TACSSchurMat *gmat = assembler->createSchurMat();  // geometric stiffness matrix
    TACSSchurMat *aux_mat = assembler->createSchurMat();  // auxillary matrix for shift and invert solver

    // Allocate the factorization
    int lev = 1e6;
    double fill = 10.0;
    int reorder_schur = 1;
    TACSSchurPc *pc = new TACSSchurPc(kmat, lev, fill, reorder_schur);
    pc->incref();

    // optional other preconditioner settings?
    assembler->assembleMatType(TACS_STIFFNESS_MATRIX, kmat);
    assembler->assembleMatType(TACS_GEOMETRIC_STIFFNESS_MATRIX, gmat);

    int subspaceSize = 10, nrestarts = 15, isFlexible = 0;
    GMRES *solver = new GMRES(aux_mat, pc, subspaceSize, nrestarts, isFlexible);
    solver->incref();

    // set relative tolerances
    solver->setTolerances(1e-12, 1e-12);

    // solve the linear static analysis for u0 and set f to NULL
    // // f = NULL;
    // u0 = assembler->createVec();  // displacements and rotations
    // f = assembler->createVec();    // loads
    // u0->incref();
    // f->incref();
    // // Create matrix and vectors
    // TACSBVec *res = assembler->createVec();  // The residual
    // res->incref();
    // assembler->applyBCs(res);
    // assembler->assembleJacobian(1.0, 0.0, 0.0, res, kmat);
    // pc->factor();  // LU factorization of stiffness matrix
    // pc->applyFactor(res, u0);
    // u0->scale(-1.0);
    // assembler->setVariables(u0);
    // // Output for visualization
    // ElementType etype = TACS_BEAM_OR_SHELL_ELEMENT;
    // int write_flag = (TACS_OUTPUT_NODES | TACS_OUTPUT_CONNECTIVITY |
    //                     TACS_OUTPUT_DISPLACEMENTS | TACS_OUTPUT_STRAINS |
    //                     TACS_OUTPUT_STRESSES | TACS_OUTPUT_EXTRAS);
    // TACSToFH5 *f5 = new TACSToFH5(assembler, etype, write_flag);
    // f5->incref();
    // f5->writeToFile("cylinder_solution.f5");
    // // return 0;

    // make the buckling solver
    TacsScalar sigma = 10.0;
    int max_lanczos_vecs = 300, num_eigvals = 50; // num_eigvals = 50;
     // need high enough # eigvals to get it right
    double eig_tol = 1e-12;

    TACSLinearBuckling *buckling = new TACSLinearBuckling(assembler, sigma,
                     gmat, kmat, aux_mat, solver, max_lanczos_vecs, num_eigvals, eig_tol);
    buckling->incref();

    // make a KSM print object for solving buckling
    KSMPrint *ksm_print = new KSMPrintStdout("BucklingAnalysis", 0, 10);
    ksm_print->incref();

    // solve the buckling analysis
    //    if u0, f are empty (are here) => then it should do a linear static analysis first
    buckling->solve(f, u0, ksm_print);

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
        std::string filename = "_buckling/mech-buckle" + std::to_string(imode) + ".f5";
        const char *cstr_filename = filename.c_str();
        f5->writeToFile(cstr_filename);
    }    

    // decref all data
    // f5->decref();
    // x->decref();
    // aux_mat->decref();
    // kmat->decref();
    // gmat->decref();
    // pc->decref();
    // solver->decref();
    // u0->decref();
    // f->decref();
    // phi->decref();
    // buckling->decref();
    // ksm_print->decref();

    MPI_Finalize();

    return 0;
}