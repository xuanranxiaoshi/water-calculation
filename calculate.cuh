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

__device__ __forceinline__ void calculate_FLUX
(int t, int pos, double (&FLUX)[4][4],
const int jt, const int NHQ, const double C0, const double C1, const double HM1, const double HM2,
/*数组长度*/ int MBQ_LEN, int MBZQ_LEN, int MBZ_LEN, int MBW_LEN, int MDI_LEN,
double* H_pre, double* U_pre, double* V_pre, double* Z_pre, double* MBQ, double* DQT, double* MBZQ, double* ZB1, double* ZBC, double* MBZ, double* DZT, double* TOPW, double* MBW, double* TOPD, double* MDI,
double** KLAS, double** QT, double** SIDE, double** ZW, double** QW, double** ZT, double** NAC, double** COSF, double** SINF);

__device__ __forceinline__ void BOUNDA
(int t, int j, int pos, double (&FLUX)[4][4],double (&QL)[3],double (&QR)[3], double &FIL, double HC, double UC, double VC, double ZC,
const int jt, const int NHQ, const double C0, const double C1,
/*数组长度*/ int MBQ_LEN, int MBZQ_LEN, int MBZ_LEN, int MBW_LEN, int MDI_LEN,
double* H_pre, double* Z_pre, double* MBQ, double* DQT, double* MBZQ, double* ZBC, double* MBZ, double* DZT, double* TOPW, double* MBW, double* TOPD, double* MDI,
double** KLAS, double** QT, double** SIDE, double** ZW, double** QW, double** ZT);

__device__ __forceinline__ void OSHER
(int t, int pos, double (&QL)[3], double (&QR)[3], double &FIL, double (&FLR_OSHER)[4],
double* H_pre);

template<int T>
__device__ __forceinline__ void QS
(int j, int pos, const double (&QL)[3], const double (&QR)[3], double &FIL, double &FIR, double (&FLR_OSHER)[4]);

__device__ __forceinline__ void QF(double h, double u, double v, double (&F)[4]);

__device__ __forceinline__ double QD(double ZL, double ZR, double ZB);

__device__ __forceinline__ void LAQP(double X, double &Y, double (&A)[5], double (&B)[5], double MS);

__device__ __forceinline__ int find_in_vec(double* list, int list_length, double target);