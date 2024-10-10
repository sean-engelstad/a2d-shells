#include "a2dcore.h"
#include "TACSObject.h"

int main() {
    // do a simple matrix vector operation
    A2D::Mat<TacsScalar, 3, 3> A;
    A2D::Vec<TacsScalar, 3> x, b;

    for (int i = 0; i < 3; i++) {
        x[i] = 2*i; // use A2D operator[] overloading here
        b[i] = 0.0;
        for (int j = 0; j < 3; j++) {
            A(i,j) = i+j;
        }
    }

    A2D::MatVecMult<TacsScalar, 3, 3>(A, x, b);
    for (int k = 0; k < 3; k++) {
        printf("b = %.4f", b[k]);
    }
};
