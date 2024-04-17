#pragma once

__device__ void calculate_HUV
(int t, int pos,// 原参数
const double HM1, const double HM2, const int DT, const double QLUA, const double VMIN,
double* H_pre, double* U_pre, double* V_pre, double* NV, double* AREA, double* ZBC, double* FNC,
double** SIDE, double** SLCOS, double** SLSIN,
double* H_res, double* U_res, double* V_res, double* Z_res, double* W_res);

__device__ void calculate_WHUV
(int t, int pos, double &WH, double &WU, double &WV,// 原参数
double** SIDE, double** SLCOS, double** SLSIN);

__device__ void calculate_FLUX(int t, int pos, double (&FLUX)[4][4]);