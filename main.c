#include <stdlib.h>
#include <stdbool.h>
#include "solve_block.h"
#include "multi_grid/multi_grid.h"
#include "tSolver/tSolver.h"
#include <omp.h>
#include "stuff.h"

void solve_top(const float* f, float* u, N_len Nlen, float dx) {
    tSolve(f, u, Nlen, 5, 1.9, dx);
}

void solve_coarse(const float* f, float* u, N_len Nlen, float dx) {
    tSolve(f, u, Nlen, 1, 1, dx);
}

void solve_base(const float* f, float* u, N_len Nlen, float dx) {
    solve(f, u, Nlen, 1, 1, dx);
}


int main() {
    for(int M = 8; M <=10; M++) {
        int Ni = pow(2,M) + 1;
        int Nj = pow(2,M) + 1;
        int Nk = pow(2,M) + 1;
        N_len Nlen = (N_len) {Ni, Nj, Nk};
        float L = 20;
        float dx = L / Ni;
        float dens = 1;
        float R = 1;
        float *u = calloc(sizeof(float), length(Nlen));
        float *u2 = calloc(sizeof(float), length(Nlen));
        float *f = calloc(sizeof(float), length(Nlen));

        //dose the initial conditions and stuff
        //its not really important what it dose
        printf("This test is on a {%i, %i, %i} grid\n", Nlen.i, Nlen.j, Nlen.k);
        inital(u, u2, f, dens, R, Nlen, L, dx, dx, M);

        char str[20];
        sprintf(str, "data_sol_%i", M);

        solve(f, u, Nlen, 100, 1.9, dx);
        save_gird(str, u, length(Nlen));

        double total = 0;
        double times = 1;

        funcs_args arg = (funcs_args) {solve_top, solve_coarse, solve_base};
        printf("starting...\n");
        for (int i = 0; i < times; i++) {
            double start = omp_get_wtime();
            multi(f, u2, Nlen, dx, arg, true); //muti call
            double end = omp_get_wtime();
            total += end - start;
            printf("time taken for multi was %f\n", end - start);
        }
        printf("avg time taken for multi was %f\n", total / times);


        //dose all the output
        sprintf(str, "data_run_%i.txt", M);
        save_gird(str, u2, length(Nlen));
        data(u, u2, f, Nlen, dx, L);
        free(u);
        free(u2);
        free(f);
    }

    return 0;
}
