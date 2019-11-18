#ifndef C_SOR_3D_STUFF_H
#define C_SOR_3D_STUFF_H

#include "multi_grid/operators.h"
#include <math.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include "multi_grid/multi_grid.h"

typedef struct _location_value {
    int i;
    int j;
    int k;
    float f;
}location_value;

void save_params(const char* str, int count, ...);

void save_line(const char* str, float* vec,int i, int j, N_len Nlen);

void save_gird(const char* str, float* vec, int length);

float L2(float* f, float* u, N_len Nlen, float dx);

float avg_diff(const float* f1, const float* f2, N_len Nlen);

location_value max_diff(const float* f1, const float* f2, N_len Nlen);

void data(float* u, float* u2, float* f, N_len Nlen, float dx, float L);

void inital(float* u, float* u2, float* f, float dens, float R, N_len Nlen, float L, float dx, float shift, int M);

#endif //C_SOR_3D_STUFF_H
