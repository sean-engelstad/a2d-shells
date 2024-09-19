#include "TACSMeshLoader.h"
#include "TACSAssembler.h"

// dependencies to make element, constitutive objects
#include "TACSShellElementTransform.h"
#include "TACSMaterialProperties.h"
#include "TACSIsoShellConstitutive.h"
#include "TACSShellElementDefs.h"

// this example is based off of examples/crm/crm.cpp in TACS

int main() {
    // make the MPI communicator
    MPI_Init(NULL, NULL);
    MPI_Comm comm = MPI_COMM_WORLD;

    int rank;
    MPI_Comm_rank(comm, &rank);

    // build the TACS mesh loader and scan the uCRM BDF file
    printf("Scanning BDF file\n");
    TACSMeshLoader *mesh = new TACSMeshLoader(comm);
    mesh->incref();
    mesh->scanBDFFile("CRM_box_2nd.bdf");

    // get number of components
    int num_components = mesh->getNumComponents();

    // create the shell ref axis transform
    TacsScalar refAxis[3] = {1.0, 0.0, 0.0};
    TACSShellRefAxisTransform *transform = new TACSShellRefAxisTransform(refAxis);
    
    // set material properties for aluminum (no thermal props input this time)
    TacsScalar rho = 2718;
    TacsScalar specific_heat = 0.0;
    TacsScalar E = 72.0e9;
    TacsScalar nu = 0.33;
    TacsScalar ys = 1e11;
    TacsScalar alpha = 0.0;
    TacsScalar kappa = 0.0;
    TACSMaterialProperties *mat = 
          new TACSMaterialProperties(rho, specific_heat, E, nu, ys, alpha, kappa);
    
    // create the constitutive and element objects for each TACS component
    printf("Creating constitutive and element objects for each TACS component\n");
    for (int icomp = 0; icomp < num_components; icomp++) {
        const char *descriptor = mesh->getElementDescript(icomp);
        TacsScalar thick = 0.010; // shell thickness
        TACSIsoShellConstitutive *con = new TACSIsoShellConstitutive(mat, thick);

        // now create the shell element object
        TACSElement *shell = TacsCreateShellByName(descriptor, transform, con);

        // set the shell element into the mesh loader for that component
        mesh->setElement(icomp, shell);
    }
    printf("\tdone\n");

    // now create the TACS assembler
    printf("Creating TACS assembler\n");
    int vars_per_node = 6; // for shell elements
    TACSAssembler *assembler = mesh->createTACS(vars_per_node);
    assembler->incref();
    mesh->decref();
    printf("\tdone");

    // TODO : next part will be to do the linear static analsis.
    printf("TODO : Need to add linear static solve..\n");
    printf("Done with RunCRM!\n");
}