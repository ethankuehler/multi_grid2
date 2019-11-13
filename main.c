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
    int Ni = 257;
    int Nj = 257;
    int Nk = 257;
    N_len Nlen = (N_len){Ni, Nj, Nk};
    float L = 20;
    float dx = L/Ni;
    float dens = 1;
    float R = 1;
    float* u = calloc(sizeof(float), length(Nlen));
    float* u2 = calloc(sizeof(float), length(Nlen));
    float* f = calloc(sizeof(float), length(Nlen));

    //dose the initial conditions and stuff
    //its not really important what it dose
    inital(u, u2, f, dens, R, Nlen, L, dx, dx);


    //printf("number of iters: %d\n", iters*6);
    tSolve(f, u, Nlen, 50, 1.9, dx);
    save_gird("data2.txt", u, length(Nlen));

    double total = 0;
    double times = 1;

    funcs_args arg = (funcs_args){solve_top, solve_coarse, solve_base};
    for(int i = 0; i < times; i++) {
        double start = omp_get_wtime();
        multi(f, u2, Nlen, dx, arg, true); //muti call
        double end = omp_get_wtime();
        total += end - start;
        printf("time taken for multi was %f\n", end - start);
    }
    printf("avg time taken for multi was %f\n", total/times);
    //tSolve(f, u2, (N_len){N, N, N}, 50, w, dx);


    //dose all the output
    save_gird("data3.txt", u2, length(Nlen));
    data(u, u2, f, Nlen, dx, L);

    return 0;
}
