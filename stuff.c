#include "multi_grid/operators.h"
#include <math.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <stdlib.h>
#include "multi_grid/multi_grid.h"
#include "stuff.h"


void save_params(const char* str, int count, ...) {
    FILE* file = fopen(str, "w");
    va_list list;
    va_start(list, count);
    for (int i = 0; i < count; i++) {
        fprintf(file, "%f\n", va_arg(list, double));
    }
    va_end(list);
    fclose(file);
}

void save_line(const char* str, float* vec,int i, int j, N_len Nlen) {
    FILE* file = fopen(str, "w");
    int n = loc(i, j, 0, Nlen);
    for (int i = n; i < n + Nlen.k; i++) {
        fprintf(file, "%f ", vec[i]);
    }
    fprintf(file, "\n");
    fclose(file);
}

void save_gird(const char* str, float* vec, int length) {
    FILE* file = fopen(str, "w");

    for (int i = 0; i < length; i++) {
        fprintf(file, "%f ", vec[i]);
    }
    fprintf(file, "\n");
    fclose(file);
}

float L2(float* f, float* u, N_len Nlen, float dx) {
    float sum = 0;
    float* r = calloc(sizeof(float*), length(Nlen));
    residual(f, u, r, Nlen, dx*dx);
    int i;
#pragma omp parallel for private(i)
    for(i = 0; i < length(Nlen); i++) {
        sum += pow(r[i],2);
    }
    free(r);
    return pow(sum, 1.0/2.0)/length(Nlen);
}

float avg_diff(const float* f1, const float* f2, N_len Nlen) {
    float sum = 0;
    int i;
#pragma omp parallel for private(i)
    for(i = 0; i < length(Nlen); i++) {
        sum += fabs(f1[i] - f2[i]);
    }
    return sum/length(Nlen);
}

location_value max_diff(const float* f1, const float* f2, N_len Nlen) {
    location_value max = (location_value){0, 0, 0, 0};
    int Ni = Nlen.i;
    int Nj = Nlen.j;
    int Nk = Nlen.k;
    int i, j, k;
#pragma omp parallel for private(i, j, k) collapse(3)
    for ( i = 0; i < Ni; i++) {
        for ( j = 0; j < Nj; j++) {
            for ( k = 0; k < Nk; k++) {
                int n = loc(i, j, k, Nlen);
                float one = f1[n];
                float two = f2[n];
                float dif = fabs(one - two);
                if (max.f < dif){
                    max.f = dif;
                    max.i = i;
                    max.j = j;
                    max.k = k;
                }
            }
        }
    }
    return max;
}

void data(float* u, float* u2, float* f, N_len Nlen, float dx, float L) {
    printf("error for multi :%lf\n", L2(f, u2, Nlen,dx));
}

void inital(float* u, float* u2, float* f, float dens, float R, N_len Nlen, float L, float dx, float shift) {
    float i_mid, j_mid, k_mid;
    i_mid = Nlen.i / 2;
    j_mid = Nlen.j / 2;
    k_mid = Nlen.k / 2;

    int Ni = Nlen.i;
    int Nj = Nlen.j;
    int Nk = Nlen.k;
    int i, j, k;

    //setting up source function
#pragma omp parallel for private(i, j, k) collapse(3)
    for ( i = 0; i < Ni; i++) {
        for ( j = 0; j < Nj; j++) {
            for ( k = 0; k < Nk; k++) {
                float x = (i - (i_mid)) * dx;
                float y = (j - (j_mid)) * dx;
                float z = (k - (k_mid)) * dx;
                float rsqrd = (x*x + y*y + z*z);
                int n = loc(i, j, k, Nlen);
                if (rsqrd < 1) {
                    f[n] = M_PI*4*dens;
                } else {
                    f[n] = 0;
                }
            }
        }
    }

    //calculating the total mass of the of the system
    double m = 0;
#pragma omp parallel for private(i, j, k) collapse(3)
    for ( i = 0; i < Ni; i++) {
        for ( j = 0; j < Nj; j++) {
            for ( k = 0; k < Nk; k++) {
                int n = loc(i,j,k, Nlen);
                if (f[n] != 0.0) {
                    m += dens*dx*dx*dx;
                }
            }
        }
    }

    //setting up th initial guess with a shift for the multi_grid to solve
    printf("The mass is %f\n", m);
#pragma omp parallel for private(i, j, k) collapse(3)
    for ( i = 0; i < Ni; i++) {
        for ( j = 0; j < Nj; j++) {
            for ( k = 0; k < Nk; k++) {
                float x = (i - (i_mid)) * dx;
                float y = (j - (j_mid)) * dx;
                float z = (k - (k_mid)) * dx + shift;
                float rsqrd = x*x + y*y + z*z;
                int n = loc(i, j, k, Nlen);
                if (rsqrd < R*R) {
                    u2[n] = -(m/(2*R*R*R))*(3*R*R - rsqrd);
                } else {
                    u2[n] = -m/sqrt(rsqrd);
                }
            }
        }
    }
    printf("dx = %f\n", dx);

    //setting up the correct solution
#pragma omp parallel for private(i, j, k) collapse(3)
    for ( i = 0; i < Ni; i++) {
        for ( j = 0; j < Nj; j++) {
            for ( k = 0; k < Nk; k++) {
                float x = (i - (i_mid)) * dx;
                float y = (j - (j_mid)) * dx;
                float z = (k - (k_mid)) * dx;
                float rsqrd = x*x + y*y + z*z;
                int n = loc(i, j, k, Nlen);
                if (rsqrd < R*R) {
                    u[n] = -(m/(2*R*R*R))*(3*R*R - rsqrd);
                } else {
                    u[n] = -m/sqrt(rsqrd);
                }
            }
        }
    }
}
