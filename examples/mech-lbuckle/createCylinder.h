#pragma once

/*
  Create the TACSAssembler object and return the associated TACS
  creator object
*/
void createAssembler(MPI_Comm comm, int order, int nx, int ny, TacsScalar udisp,
                     TACSElement *element, TACSAssembler **_assembler,
                     TACSCreator **_creator) {
  int rank;
  MPI_Comm_rank(comm, &rank);

  double L = 100.0;
  double R = 100.0 / M_PI;
  double defect = 0.1;

  // Set the alpha and beta parameters
  double alpha = 4.0 / R;
  double beta = 3 * M_PI / L;

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
    int *bc_vars = new int[4 * numBcs];
    TacsScalar *bc_vals = new TacsScalar[4 * numBcs];
    bc_ptr[0] = 0;
    int i = 0; // BC counter

    for (int j = 0; j < nny; j++) { // set u, v, w, xrot to zero
      // bc at x- edge
      bc_ptr[i+1] = bc_ptr[i]; // start BC dof counter
      int node = j * nnx;
      bcNodes[k] = node;
      for (int m = 0; m < 4; m++) {
        bc_vars[bc_ptr[i+1]] = m; // DOF m is set to 0 disp
        bc_vals[bc_ptr[i+1]] = 0.0;
        bc_ptr[i+1]++;
        // printf("x- node %d with DOF %d set to %.2e with total BCs %d\n", node, m, bc_vals[bc_ptr[i+1]-1], bc_ptr[i+1]);
      }
      k++; i++;

      bc_ptr[i+1] = bc_ptr[i]; // start BC dof counter
      node = nnx - 1 + j * nnx;
      bcNodes[k] = node;
      for (int m = 0; m < 4; m++) { // set x = udisp and y,z,rotx to 0
        bc_vars[bc_ptr[i+1]] = m; // DOF m is set to 0 disp
        TacsScalar disp;
        if (m == 0) {
            disp = udisp;
        } else {
            disp = 0.0;
        }
        bc_vals[bc_ptr[i+1]] = disp;
        bc_ptr[i+1]++;
        // printf("x+ node %d with DOF %d set to %.2e with total BCs %d\n", node, m, bc_vals[bc_ptr[i+1]-1], bc_ptr[i+1]);
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