#include "TACSMeshLoader.h"
#include "TACSAssembler.h"

// dependencies to make element, constitutive objects
#include "TACSShellElementTransform.h"
#include "TACSMaterialProperties.h"
#include "TACSIsoShellConstitutive.h"
#include "TACSShellElementDefs.h"
#include "TACSBuckling.h"
#include "KSM.h"

// this example is based off of examples/crm/crm.cpp in TACS

int main() {
    // make the MPI communicator
    MPI_Init(NULL, NULL);
    MPI_Comm comm = MPI_COMM_WORLD;

    int rank;
    MPI_Comm_rank(comm, &rank);

    // build the TACS mesh loader and scan the uCRM BDF file
    if (rank == 0) {
        printf("Scanning BDF file\n");
    }
    
    TACSMeshLoader *mesh = new TACSMeshLoader(comm);
    mesh->incref();
    mesh->scanBDFFile("therm-cylinder.bdf");

    // get number of components
    int num_components = mesh->getNumComponents();

    // create the shell ref axis transform
    // TacsScalar refAxis[3] = {1.0, 0.0, 0.0};
    // TACSShellRefAxisTransform *transform = new TACSShellRefAxisTransform(refAxis);
    
    TACSShellTransform *transform = new TACSShellNaturalTransform();

    // set material properties for aluminum (no thermal props input this time)
    TacsScalar rho = 2718;
    TacsScalar specific_heat = 0.0;
    TacsScalar E = 72.0e9;
    TacsScalar nu = 0.33;
    TacsScalar ys = 1e11;
    TacsScalar cte = 10e-6; // double check this value
    TacsScalar kappa = 0.0;
    TACSMaterialProperties *mat = 
          new TACSMaterialProperties(rho, specific_heat, E, nu, ys, cte, kappa);
    
    // create the constitutive and element objects for each TACS component
    if (rank == 0) {
        printf("Creating constitutive and element objects for each TACS component\n");
    }
    
    for (int icomp = 0; icomp < num_components; icomp++) {
        const char *descriptor = mesh->getElementDescript(icomp);
        TacsScalar thick = 0.010; // shell thickness
        TACSIsoShellConstitutive *con = new TACSIsoShellConstitutive(mat, thick);

        // now create the shell element object
        TACSElement *shell = TacsCreateShellByName(descriptor, transform, con);

        // set the shell element into the mesh loader for that component
        mesh->setElement(icomp, shell);
    }

    mat->decref();
    transform->decref();

    if (rank == 0) {
        printf("\tdone\n");
    }
    

    // now create the TACS assembler
    if (rank == 0) {
        printf("Creating TACS assembler\n");
    }
    
    int vars_per_node = 6; // for shell elements
    TACSAssembler *assembler = mesh->createTACS(vars_per_node);
    assembler->incref();
    mesh->decref();
    if (rank == 0) {
        printf("Done\n");
    }

    // set temperature into all elements for thermal buckling
    TacsScalar temperature = 10.0; // default 1.0 // 1 deg K
    int numElements = assembler->getNumElements();
    TACSQuad4Shell *elem;
    for (int ielem = 0; ielem < numElements; ielem++) {
        elem = dynamic_cast<TACSQuad4Shell *>(assembler->getElement(ielem));
        elem->setTemperature(temperature);
        #ifdef TACS_USE_COMPLEX
            elem->setComplexStepGmatrix(true);
        #endif
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
    TACSBVec *u0 = NULL;  // displacements and rotations
    TACSBVec *f = NULL;    // loads
    // NOTE : if you want to do linearized buckling about arbitrary disps
    // then you should do linear static yourself and then give the TACSLinearBuckling your u0
    // u0->incref();
    // f->incref();

    // create the matrices for buckling
    TACSSchurMat *kmat = assembler->createSchurMat();  // stiffness matrix
    TACSSchurMat *gmat = assembler->createSchurMat();  // geometric stiffness matrix
    TACSSchurMat *aux_mat = assembler->createSchurMat();  // auxillary matrix for shift and invert solver

    // Allocate the factorization
    int lev = 1e6; // default was 10000 (4 zeros), but TACS buckling default uses 1e6 (6 zeros)
    double fill = 10.0;
    int reorder_schur = 1;
    TACSSchurPc *pc = new TACSSchurPc(kmat, lev, fill, reorder_schur);
    pc->incref();

    // int subspaceSize = 10, nrestarts = 15, isFlexible = 0; //defaults
    int subspaceSize = 100, nrestarts = 15, isFlexible = 0;
    GMRES *solver = new GMRES(aux_mat, pc, subspaceSize, nrestarts, isFlexible);
    solver->incref();

    // set relative tolerances
    solver->setTolerances(1e-12, 1e-12);

    // make the buckling solver
    TacsScalar sigma = 40.0;
    int max_lanczos_vecs = 100, num_eigvals = 50;
    double eig_tol = 1e-12; // 1e-12

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
        // Emily wants me to normalize the eigenvector here
        assembler->setVariables(phi);   
        std::string filename = "_buckling/therm-buckle" + std::to_string(imode) + ".f5";
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