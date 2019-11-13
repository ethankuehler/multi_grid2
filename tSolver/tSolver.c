#include "../multi_grid/multi_grid.h"
#include "../multi_grid/operators.h"
#include <stdlib.h>
#include <string.h>


void solver(const float* f, float* u, N_len Nlen, float w, float dxs) {
    int Ni = Nlen.i;
    int Nj = Nlen.j;
    int Nk = Nlen.k;
    int i, j, k;
#pragma omp parallel for schedule(auto) shared(u, f, Nlen, w, dxs) private(i, j, k) collapse(1)
    for (i = 1; i < Ni - 1; i++) {
        for (j = 1; j < Nj - 1; j++) {
            for (k = 1; k < Nk - 1; k++) {
                int n = loc(i, j, k, Nlen);
                int left    = loc(i, j, k + 1, Nlen);
                int right   = loc(i, j, k - 1, Nlen);
                int forward = loc(i, j + 1, k, Nlen);
                int back    = loc(i, j - 1, k, Nlen);
                int up      = loc(i + 1, j, k, Nlen);
                int down    = loc(i - 1, j, k, Nlen);
                u[n] = (1 - w)*u[n] + (w/-6)*(f[n]*dxs - (u[down] + u[back] + u[left] + u[right] + u[forward] + u[up]));
            }
        }
    }

}

void tSolve(const float* f, float* u, N_len Nlen, int iter, float w, float dx) {
    for(int _= 0; _ < iter; _++) {
        solver(f, u, Nlen, w, dx * dx);
    }
}

