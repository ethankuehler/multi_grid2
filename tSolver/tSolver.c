#include "../multi_grid/multi_grid.h"
#include "../multi_grid/operators.h"
#include <stdlib.h>
#include <string.h>


void solver(const float* f, float* u, N_len Nlen, float w, float dxs) {
    int Ni = Nlen.i;
    int Nj = Nlen.j;
    int Nk = Nlen.k;
    int i, j, k;
#pragma omp parallel for num_threads(mul_num) schedule(auto) shared(u, f, Nlen, w, dxs) private(i, j, k) collapse(1)
    for (i = 1; i < Ni - 1; i++) {
        for (j = 1; j < Nj - 1; j++) {
            for (k = 1; k < Nk - 1; k++) {
                int n = loc(i, j, k, Nlen);
                u[n] = (1 - w)*u[n] +
                       (w/-6)*(f[n]*dxs -(u[n - Nk*Nj] + u[n - Nk] + u[n - 1] + u[n + 1] + u[n + Nk] + u[n + Nk*Nj]));
            }
        }
    }

}

void tSolve(const float* f, float* u, N_len Nlen, int iter, float w, float dx) {
    for(int _= 0; _ < iter; _++) {
        solver(f, u, Nlen, w, dx * dx);
    }
}

