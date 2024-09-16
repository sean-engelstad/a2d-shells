#include "TACSShellElementTransform.h"
#include "TACSMaterialProperties.h"
#include "TACSIsoShellConstitutive.h"
#include "TACSShellElementDefs.h"

// ensure that you can build a TACS Shell element

int main() {
    // build the ACS Shell element with several prelim steps

    printf("1 - making TACS shell transform\n");
    TacsScalar refAxis[3] = {1.0, 0.0, 0.0};
    TACSShellRefAxisTransform *transform = new TACSShellRefAxisTransform(refAxis);
    printf("\tdone\n");
    
    printf("2 - make TACS material properties\n");
    TACSMaterialProperties *mat = new TACSMaterialProperties(2718, 0.0, 72e9, 0.33, 1e11, 0.0, 0.0);
    printf("\tdone\n");

    printf("3 - make TACS Isoshell constitutive\n");
    TacsScalar thick = 0.010; // thickness in (m)
    TACSIsoShellConstitutive *con = new TACSIsoShellConstitutive(mat, thick);
    printf("\tdone\n");

    printf("4 - make TACS shell element\n");
    TACSQuad4Shell *elem = new TACSQuad4Shell(transform, con);
    printf("\tdone\n");
}