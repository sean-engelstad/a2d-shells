# TACS Shell Elements with A2D-SHELLS

## Including A2D in the Shell Elements with Director Formulation:
* Use `archive/TACSBeamElement.h` as a template since an older version of A2D was included in this element in TACS.
* Goal is to adapt the methods `computeEnergies()`, `addResidual()` and `addJacobian()` to use A2D in the `TACSShellElement.h` class.