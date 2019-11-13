
#include "multi_grid/operators.h"


void method1(float* u, const float* f, N_len Nlen, float w, float dxs) {
    int Ni = Nlen.i;
    int Nj = Nlen.j;
    int Nk = Nlen.k;
    int i;
    for (i = 1; i < Ni - 1; i++) {
        for (int j = 1; j < Nj - 1; j++) {
            for (int k = 1; k < Nk - 1; k++) {
                int n = loc(i, j, k, Nlen);
                u[n] = (1 - w)*u[n] +
                        (w/-6)*(f[n]*dxs -(u[n - Nk*Nj] + u[n - Nk] + u[n - 1] + u[n + 1] + u[n + Nk] + u[n + Nk*Nj]));
            }
        }
    }
}

void solve(const float* f, float* u, N_len Nlen, int iters, float w, float dx) {
    dx = dx*dx;
    for (int i = 0; i < iters; i++) {
        method1(u, f, Nlen, w, dx);
    }
}
