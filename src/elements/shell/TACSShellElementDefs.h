#ifndef TACS_SHELL_ELEMENT_DEFS_H
#define TACS_SHELL_ELEMENT_DEFS_H

#include "TACSDirector.h"
#include "TACSShellElement.h"
#include "TACSShellElementModel.h"
#include "TACSShellElementQuadBasis.h"

/*
  Linear shell elements with appropriate quadrature schemes
*/
typedef TACSShellElement<TACSQuadLinearQuadrature, TACSShellQuadBasis<2>,
                         TACSLinearizedRotation, TACSShellLinearModel>
    TACSQuad4Shell;

/*
  Shell elements with a linearized rotation and nonlinear strain expressions
*/
typedef TACSShellElement<TACSQuadLinearQuadrature, TACSShellQuadBasis<2>,
                         TACSLinearizedRotation, TACSShellNonlinearModel>
    TACSQuad4NonlinearShell;

/**
  Create a TACS shell element based on the name of the shell.

  @param name The name of the shell element
  @param transform The transformation used for the shell
  @param con The shell constitutive object
*/
inline TACSElement *TacsCreateShellByName(const char *name,
                                          TACSShellTransform *transform,
                                          TACSShellConstitutive *con) {
  TACSElement *shell = NULL;
  if (strcmp(name, "TACSQuad4Shell") == 0 ||
             strcmp(name, "CQUAD") == 0 || strcmp(name, "CQUAD4") == 0 ||
             strcmp(name, "CQUADR") == 0) {
    shell = new TACSQuad4Shell(transform, con);
  } else if (strcmp(name, "TACSQuad4NonlinearShell") == 0) {
    shell = new TACSQuad4NonlinearShell(transform, con);
  }

  return shell;
}

#endif  // TACS_SHELL_ELEMENT_DEFS_H
