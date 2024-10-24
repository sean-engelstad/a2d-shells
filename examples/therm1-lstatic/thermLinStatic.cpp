#include "TACSAssembler.h"
#include "TACSCreator.h"
#include "TACSElementAlgebra.h"
#include "TACSElementVerification.h"
#include "TACSIsoShellConstitutive.h"
#include "TACSShellElementDefs.h"
#include "TACSToFH5.h"
#include "tacslapack.h"

// """
// Code adapted from TACS example:
// examples/shell/cylinder.cpp
// since caps2tacs mesh generation isn't the best right now
// """

/*
  Create the TACSAssembler object and return the associated TACS
  creator object
*/
void createAssembler(MPI_Comm comm, int order, int nx, int ny, TacsScalar udisp,
                     double R, double L,
                     TACSElement *element, TACSAssembler **_assembler,
                     TACSCreator **_creator) {
  int rank;
  MPI_Comm_rank(comm, &rank);
  double defect = 0.1;

  // Set the number of nodes/elements on this proc
  int varsPerNode = element->getVarsPerNode();

  // Set up the creator object
  TACSCreator *creator = new TACSCreator(comm, varsPerNode);

  if (rank == 0) {
    // Set the number of elements
    int nnx = (order - 1) * nx + 1;
    int nny = (order - 1) * ny;
    int numNodes = nnx * nny;
    int numElements = nx * ny;

    // Allocate the input arrays into the creator object
    int *ids = new int[numElements];
    int *ptr = new int[numElements + 1];
    int *conn = new int[order * order * numElements];

    // Set the element identifiers to all zero
    memset(ids, 0, numElements * sizeof(int));

    ptr[0] = 0;
    for (int k = 0; k < numElements; k++) {
      // Back out the i, j coordinates from the corresponding
      // element number
      int i = k % nx;
      int j = k / nx;

      // Set the node connectivity
      for (int jj = 0; jj < order; jj++) {
        for (int ii = 0; ii < order; ii++) {
          if (j == ny - 1 && jj == order - 1) {
            conn[order * order * k + ii + order * jj] = ((order - 1) * i + ii);
          } else {
            conn[order * order * k + ii + order * jj] =
                ((order - 1) * i + ii) + ((order - 1) * j + jj) * nnx;
          }
        }
      }

      ptr[k + 1] = order * order * (k + 1);
    }

    // Set the connectivity
    creator->setGlobalConnectivity(numNodes, numElements, ptr, conn, ids);
    delete[] conn;
    delete[] ptr;
    delete[] ids;

    int numBcs = 2 * nny;
    int *bcNodes = new int[numBcs];
    int k = 0;

    int *bc_ptr = new int[numBcs + 1];
    int *bc_vars = new int[3 * numBcs];
    TacsScalar *bc_vals = new TacsScalar[3 * numBcs];
    bc_ptr[0] = 0;
    int i = 0; // BC counter

    for (int j = 0; j < nny; j++) { // set u, v, w, xrot to zero
      // bc at x- edge
      bc_ptr[i+1] = bc_ptr[i]; // start BC dof counter
      int node = j * nnx;
      bcNodes[k] = node;
      for (int m = 0; m < 3; m++) {
        bc_vars[bc_ptr[i+1]] = m; // DOF m is set to 0 disp
        bc_vals[bc_ptr[i+1]] = 0.0;
        bc_ptr[i+1]++;
        printf("x- node %d with DOF %d set to %.2e with total BCs %d\n", node, m, bc_vals[bc_ptr[i+1]-1], bc_ptr[i+1]);
      }
      k++; i++;

      TacsScalar vdisp = 0.0;
      TacsScalar wdisp = 0.0;

      bc_ptr[i+1] = bc_ptr[i]; // start BC dof counter
      node = nnx - 1 + j * nnx;
      bcNodes[k] = node;
      for (int m = 0; m < 3; m++) { // set x = udisp and y,z,rotx to 0
        bc_vars[bc_ptr[i+1]] = m; // DOF m is set to 0 disp
        TacsScalar disp;
        if (m == 0) {
          disp = udisp;
        } else if (m == 1) {
          disp = vdisp; // v disp to 0.01
        } else if (m == 2) {
          disp = wdisp;
        } else {
            disp = 0.0;
        }
        bc_vals[bc_ptr[i+1]] = disp;
        bc_ptr[i+1]++;
        printf("x+ node %d with DOF %d set to %.2e with total BCs %d\n", node, m, bc_vals[bc_ptr[i+1]-1], bc_ptr[i+1]);
      }
      k++; i++;
    }
    printf("udisp = %.8f", udisp);

    // Set the boundary conditions
    creator->setBoundaryConditions(numBcs, bcNodes, bc_ptr, bc_vars, bc_vals);

    delete[] bcNodes;
    delete[] bc_ptr;
    delete[] bc_vars;
    delete[] bc_vals;

    // Set the node locations
    TacsScalar *Xpts = new TacsScalar[3 * numNodes];

    for (int j = 0; j < nny; j++) {
      double v = -M_PI + (2.0 * M_PI * j) / nny;
      for (int i = 0; i < nnx; i++) {
        double u = 1.0 * i / (nnx - 1);
        // double theta =
        //     v + 0.25 * M_PI * u + defect * sin(v) * cos(2 * M_PI * u);
        // double x = L * (u + defect * cos(v) * sin(2 * M_PI * u));
        double theta = v;
        double x = L*u;

        int node = i + j * nnx;
        Xpts[3 * node] = x;
        Xpts[3 * node + 1] = R * cos(theta);
        Xpts[3 * node + 2] = -R * sin(theta);
      }
    }

    // Set the nodal locations
    creator->setNodes(Xpts);
    delete[] Xpts;
  }

  // Set the one element
  creator->setElements(1, &element);

  // Set the reordering type
  creator->setReorderingType(TACSAssembler::MULTICOLOR_ORDER,
                             TACSAssembler::GAUSS_SEIDEL);

  // Create TACS
  TACSAssembler *assembler = creator->createTACS();

  // Set the elements the node vector
  TACSBVec *X = assembler->createNodeVec();
  X->incref();
  assembler->getNodes(X);

  X->decref();

  // Set the pointers
  *_assembler = assembler;
  *_creator = creator;
}

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
  TacsScalar E = 70e9;
  TacsScalar nu = 0.3;
  TacsScalar ys = 270.0;
  TacsScalar cte = 23.5e-6;
  TacsScalar kappa = 230.0;
  TACSMaterialProperties *props =
      new TACSMaterialProperties(rho, specific_heat, E, nu, ys, cte, kappa);

  TacsScalar axis[] = {1.0, 0.0, 0.0};
  TACSShellTransform *transform = new TACSShellRefAxisTransform(axis);

  TACSShellConstitutive *con = new TACSIsoShellConstitutive(props, t);

  TACSAssembler *assembler = NULL;
  TACSCreator *creator = NULL;
  TACSElement *shell = NULL;
  shell = new TACSQuad4Shell(transform, con);
  shell->incref();
  createAssembler(comm, order, nx, ny, udisp, R, L, shell, &assembler, &creator);
  

  assembler->incref();
  creator->incref();

  // Free the creator object
  creator->decref();

  // Create matrix and vectors
  TACSBVec *ans = assembler->createVec();  // displacements and rotations
  TACSBVec *res = assembler->createVec();  // The residual
  TACSSchurMat *mat = assembler->createSchurMat();  // stiffness matrix

  // Increment reference count to the matrix/vectors
  ans->incref();
  res->incref();
  mat->incref();

  // Allocate the factorization
  int lev = 10000;
  double fill = 10.0;
  int reorder_schur = 1;
  TACSSchurPc *pc = new TACSSchurPc(mat, lev, fill, reorder_schur);
  pc->incref();

  assembler->setTemperatures(100.0); // 10 K

  // apply displacement control BCs to the residual
  assembler->applyBCs(res);

  assembler->assembleJacobian(1.0, 0.0, 0.0, res, mat);
  pc->factor();  // LU factorization of stiffness matrix
  pc->applyFactor(res, ans);

  ans->scale(-1.0);
  assembler->setVariables(ans);

  // Output for visualization
  ElementType etype = TACS_BEAM_OR_SHELL_ELEMENT;
  int write_flag = (TACS_OUTPUT_NODES | TACS_OUTPUT_CONNECTIVITY |
                    TACS_OUTPUT_DISPLACEMENTS | TACS_OUTPUT_STRAINS |
                    TACS_OUTPUT_STRESSES | TACS_OUTPUT_EXTRAS);
  TACSToFH5 *f5 = new TACSToFH5(assembler, etype, write_flag);
  f5->incref();
  f5->writeToFile("cylinder_solution.f5");

  shell->decref();
  assembler->decref();

  MPI_Finalize();
}
