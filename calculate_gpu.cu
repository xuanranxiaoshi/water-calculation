#include "read_data.hh"
#include "calculate.cuh" // calculate.hh中所有函数都被改造成__device__,声明全部移到此处

using namespace DataManager;

#define FLR(x) FLUX[x][j]
#define FLUX_VAL(A, B, C, D)                                                   \
  FLR(0) = A;                                                                  \
  FLR(1) = B;                                                                  \
  FLR(2) = C;                                                                  \
  FLR(3) = D

// device 变量
int* NV_dev;
double* H_dev, *U_dev, *V_dev, *Z_dev, *W_dev, *QT_dev, *ZT_dev;
double* AREA_dev, *ZBC_dev, *ZB1_dev, *DQT_dev,  *DZT_dev,  *TOPW_dev,  *TOPD_dev,  *MBQ_dev,  *MBZQ_dev,  *MBW_dev,  *MDI_dev;
double* H_res, *U_res, *V_res, *Z_res, *W_res;

// 二维矩阵的索引指针
int **NAC_hIdx, **NAC_devIdx;
double** SIDE_hIdx, **SLCOS_hIdx, **SLSIN_hIdx, **KLAS_hIdx, **ZW_hIdx, **QW_hIdx;
double** SIDE_devIdx, **SLCOS_devIdx, **SLSIN_devIdx, **KLAS_devIdx, **ZW_devIdx, **QW_devIdx;

// 二维矩阵的实际数据指针
int *NAC_dev;
double* SIDE_dev, *SLCOS_dev, *SLSIN_dev, *KLAS_dev, *ZW_dev, *QW_dev;

size_t pitchInt, pitchDouble;
const int block_size = 512;  // CUDA maximum is 1024

void initialize_deviceVar(){
  // 显式内存分配
  cudaMalloc((void**) &H_dev, sizeof(double) * CEL);
  cudaMalloc((void**) &U_dev, sizeof(double) * CEL);
  cudaMalloc((void**) &V_dev, sizeof(double) * CEL);
  cudaMalloc((void**) &Z_dev, sizeof(double) * CEL);
  cudaMalloc((void**) &W_dev, sizeof(double) * CEL);
  cudaMalloc((void**) &QT_dev, sizeof(double) * NDAYS);
  cudaMalloc((void**) &ZT_dev, sizeof(double) * NDAYS);
  cudaMalloc((void**) &NV_dev, sizeof(int) * CEL);
  
  cudaMalloc((void**) &AREA_dev, sizeof(double)* CEL);
  cudaMalloc((void**) &ZBC_dev, sizeof(double)* CEL);
  cudaMalloc((void**) &ZB1_dev, sizeof(double)* CEL);
  cudaMalloc((void**) &DQT_dev, sizeof(double)* NQ);
  cudaMalloc((void**) &DZT_dev, sizeof(double)* NZ);
  cudaMalloc((void**) &MBQ_dev, sizeof(double)* NQ);
  cudaMalloc((void**) &TOPW_dev, sizeof(double)* CEL);
  cudaMalloc((void**) &TOPD_dev, sizeof(double)* CEL);
  cudaMalloc((void**) &MBZQ_dev, sizeof(double)* CEL);
  cudaMalloc((void**) &MBW_dev, sizeof(double)* CEL);
  cudaMalloc((void**) &MDI_dev, sizeof(double)* CEL);

  cudaMalloc((void**) &H_res, sizeof(double)* CEL);
  cudaMalloc((void**) &U_res, sizeof(double)* CEL);
  cudaMalloc((void**) &V_res, sizeof(double)* CEL);
  cudaMalloc((void**) &Z_res, sizeof(double)* CEL);
  cudaMalloc((void**) &W_res, sizeof(double)* CEL);


  // 时间无关的数据迁移
  cudaMemcpy(NV_dev, &NV[0] , sizeof(int) * CEL, cudaMemcpyHostToDevice);
  cudaMemcpy(ZBC_dev, &ZBC[0], sizeof(double)* CEL,  cudaMemcpyHostToDevice);
  cudaMemcpy(AREA_dev, &AREA[0], sizeof(double)* CEL, cudaMemcpyHostToDevice);
  cudaMemcpy(ZB1_dev, &ZB1[0], sizeof(double)* CEL,  cudaMemcpyHostToDevice);
  cudaMemcpy(MBQ_dev, &MBQ[0], sizeof(double)* NQ, cudaMemcpyHostToDevice);


  // 二维变量
  int Wd = CEL;
  int Ht = 4;

  NAC_hIdx = new int*[Ht];
  SIDE_hIdx = new double*[Ht];
  SLCOS_hIdx = new double*[Ht];
  SLSIN_hIdx = new double*[Ht];
  KLAS_hIdx = new double*[Ht];

  cudaMalloc((void**) &NAC_devIdx, sizeof(int*) * Ht);
  cudaMalloc((void**) &SIDE_devIdx, sizeof(double*) * Ht);
  cudaMalloc((void**) &SLCOS_devIdx, sizeof(double*) * Ht);
  cudaMalloc((void**) &SLSIN_devIdx, sizeof(double*) * Ht);
  cudaMalloc((void**) &KLAS_devIdx, sizeof(double*) * Ht);

  cudaMallocPitch((void**) &NAC_dev, &pitchInt, sizeof(int) * Wd, Ht);
  cudaMallocPitch((void**) &SIDE_dev, &pitchDouble, sizeof(double) * Wd, Ht);
  cudaMallocPitch((void**) &SLCOS_dev, &pitchDouble, sizeof(double) * Wd, Ht);
  cudaMallocPitch((void**) &SLSIN_dev, &pitchDouble, sizeof(double) * Wd, Ht);
  cudaMallocPitch((void**) &KLAS_dev, &pitchDouble, sizeof(double) * Wd, Ht);

  for(int row = 0; row < Ht; ++row){   
      // 1. host 端建立 device 上数据的索引
      NAC_hIdx[row] = &NAC_dev[row*(pitchInt/sizeof(int))];
      SIDE_hIdx[row] = &SIDE_dev[row*(pitchDouble/sizeof(double))];
      SLCOS_hIdx[row] = &SLCOS_dev[row*(pitchDouble/sizeof(double))];
      SLSIN_hIdx[row] = &SLSIN_dev[row*(pitchDouble/sizeof(double))];
      KLAS_hIdx[row] = &KLAS_dev[row*(pitchDouble/sizeof(double))];

      // 2. 数据从 host 拷贝到 device
		  cudaMemcpy(&NAC_dev[row*(pitchInt/sizeof(int))], NAC[row].data(), sizeof(int)*Wd, cudaMemcpyHostToDevice);
      cudaMemcpy(&SIDE_dev[row*(pitchDouble/sizeof(double))], SIDE[row].data(), sizeof(double)*Wd, cudaMemcpyHostToDevice);
      cudaMemcpy(&SLCOS_dev[row*(pitchDouble/sizeof(double))], SLCOS[row].data(), sizeof(double)*Wd, cudaMemcpyHostToDevice);
      cudaMemcpy(&SLSIN_dev[row*(pitchDouble/sizeof(double))], SLSIN[row].data(), sizeof(double)*Wd, cudaMemcpyHostToDevice);
      cudaMemcpy(&KLAS_dev[row*(pitchDouble/sizeof(double))], KLAS[row].data(), sizeof(double)*Wd, cudaMemcpyHostToDevice);
	}
  
  // 3. host 端的索引拷贝到 device 端
  cudaMemcpy(NAC_devIdx, NAC_hIdx, sizeof(int*) * Ht, cudaMemcpyHostToDevice);
  cudaMemcpy(SIDE_devIdx, SIDE_hIdx, sizeof(double*)*Ht, cudaMemcpyHostToDevice);
  cudaMemcpy(SLCOS_devIdx, SLCOS_hIdx, sizeof(double*)*Ht, cudaMemcpyHostToDevice);
  cudaMemcpy(SLSIN_devIdx, SLSIN_hIdx, sizeof(double*)*Ht, cudaMemcpyHostToDevice);
  cudaMemcpy(KLAS_devIdx, KLAS_hIdx, sizeof(double*)*Ht, cudaMemcpyHostToDevice);
  
}

__device__ __forceinline__ void calculate_FLUX
(int t, int pos, double (&FLUX)[4][4],
const int jt, const int NHQ, const double C0, const double C1, const double HM1, const double HM2, const double VMIN,
int MBQ_LEN, int MBZQ_LEN, int MBZ_LEN, int MBW_LEN, int MDI_LEN,//数组长度
double* H_pre, double* U_pre, double* V_pre, double* Z_pre, double* MBQ, double* DQT, double* MBZQ, double* ZB1, double* ZBC, double* MBZ, double* DZT, double* TOPW, double* MBW, double* TOPD, double* MDI,
double** KLAS, double** QT, double** SIDE, double** ZW, double** QW, double** ZT, double** NAC, double** COSF, double** SINF) {
  // Vec QL(3);
  // Vec QR(3);
  // QL[0] = H[TIME_PREV][pos];
  // double U1 = U[TIME_PREV][pos];
  // double V1 = V[TIME_PREV][pos];
  // double H1 = H[TIME_PREV][pos];
  // double Z1 = Z[TIME_PREV][pos];
  // double v_ZB1 = ZB1[pos];
  // Vec FLR_OSHER(4);
  double QL[3];
  double QR[3];
  QL[0] = H_pre[pos];
  double U1 = U_pre[pos];
  double V1 = V_pre[pos];
  double H1 = H_pre[pos];
  double Z1 = Z_pre[pos];
  double v_ZB1 = ZB1[pos];
  double FLR_OSHER[4];

  for (int j = 0; j < 4; j++) {
    // 局部变量声明，以及初始化运算。
    int NC = NAC[j][pos] - 1;
    double KP = KLAS[j][pos];
    double COSJ = COSF[j][pos];
    double SINJ = SINF[j][pos];
    QL[1] = U1 * COSJ + V1 * SINJ;
    QL[2] = V1 * COSJ - U1 * SINJ;
    double CL = sqrt(9.81 * H1);
    double FIL = QL[1] + 2 * CL;
    double HC, BC, ZC, UC, VC;
    double ZI = fmax(Z1, v_ZB1);

    if (NC == 0) {
      HC = 0;
      BC = 0;
      ZC = 0;
      UC = 0;
      VC = 0;
    } else {
      HC = fmax(H_pre[NC], HM1);
      BC = ZBC[NC];
      ZC = fmax(ZBC[NC], Z_pre[NC]);
      UC = U_pre[NC];
      VC = V_pre[NC];
    }

    if ((KP >= 1 && KP <= 8) || KP >= 10) {
      BOUNDA(t, j, pos, FLUX, QL, QR, FIL, HC, UC, VC, ZC,
      jt, NHQ, C0, C1,
      MBQ_LEN, MBZQ_LEN, MBZ_LEN, MBW_LEN, MDI_LEN,
      H_pre, Z_pre, MBQ, DQT, MBZQ, ZBC, MBZ, DZT, TOPW, MBW, TOPD, MDI,
      KLAS, QT, SIDE, ZW, QW, ZT, NAC, COSF, SINF);
    } else if (H1 <= HM1 && HC <= HM1) {
      FLUX_VAL(0, 0, 0, 0);
    } else if (ZI <= BC) {
      FLUX_VAL(-C1 * pow(HC, 1.5), H1 * QL[1] * fabs(QL[1]),
               0, 4.905 * H1 * H1);
    } else if (ZC <= ZBC[pos]) {
      FLUX_VAL(C1 * pow(H1, 1.5), FLR(0) * QL[1], FLR(0) * QL[2],
               0);
    } else if (H1 <= HM2) {
      if (ZC > ZI) {
        double DH = fmax(ZC - ZBC[pos], HM1);
        double UN = -C1 * sqrt(DH);
        FLUX_VAL(DH * UN, FLR(0) * UN,
                 FLR(0) * (VC * COSJ - UC * SINJ),
                 4.905 * H1 * H1);
      } else {
        FLUX_VAL(C1 * pow(H1, 1.5), 0, 0,
                 4.905 * H1 * H1);
      }
    } else if (HC <= HM2) {
      if (ZI > ZC) {
        double DH = fmax(ZI - BC, HM1);
        double UN = C1 * sqrt(DH);
        double HC1 = ZC - ZBC[pos];
        FLUX_VAL(DH * UN, FLR(0) * UN, FLR(0) * QL[2], 4.905 * HC1 * HC1);
      } else {
        FLUX_VAL(-C1 * pow(HC, 1.5),
                 H1 * QL[1] * QL[1], 0,
                 4.905 * H1 * H1);
      }
    } else {
      QR[0] = fmax(ZC - ZBC[pos], HM1);
      double UR = UC * COSJ + VC * SINJ;
      QR[1] = UR * min(HC / QR[0], 1.5);
      if (HC <= HM2 || QR[0] <= HM2) {
        QR[1] = copysign(VMIN, UR);
      }
      QR[2] = VC * COSJ - UC * SINJ;
      OSHER(t, pos, QL, QR, FIL, FLR_OSHER, H_pre);
      FLR(1) = FLR_OSHER[1] + (1 - min(HC / QR[0], 1.5)) * HC * UR * UR / 2;
      FLUX_VAL(FLR_OSHER[0], FLR(1), FLR_OSHER[2], FLR_OSHER[3]);
    }
  }
}

__device__ __forceinline__ void calculate_WHUV
(int t, int pos, double &WH, double &WU, double &WV,
const int jt, const int NHQ, const double C0, const double C1, const double HM1, const double HM2, const double VMIN,
int MBQ_LEN, int MBZQ_LEN, int MBZ_LEN, int MBW_LEN, int MDI_LEN,//数组长度
double* H_pre, double* U_pre, double* V_pre, double* Z_pre, double* MBQ, double* DQT, double* MBZQ, double* ZB1, double* ZBC, double* MBZ, double* DZT, double* TOPW, double* MBW, double* TOPD, double* MDI,
double** KLAS, double** QT, double** SIDE, double** ZW, double** QW, double** ZT, double** NAC, double** COSF, double** SINF, double** SLCOS, double** SLSIN) {
  // Vec2 FLUX(4, Vec(4));
  double FLUX[4][4];// 在kernel中不能使用vector
  calculate_FLUX(t, pos, FLUX,
  jt, NHQ, C0, C1, HM1, HM2, VMIN,
  MBQ_LEN, MBZQ_LEN, MBZ_LEN, MBW_LEN, MDI_LEN,
  H_pre, U_pre, V_pre, Z_pre, MBQ, DQT, MBZQ, ZB1, ZBC, MBZ, DZT, TOPW, MBW, TOPD, MDI,
  KLAS, QT, SIDE, ZW, QW, ZT, NAC, COSF, SINF);

  for (int j = 0; j < 4; j++) {
    FLR(1) = FLUX[1][j] + FLUX[3][j];
    FLR(2) = FLUX[2][j];
    double SL = SIDE[j][pos];
    double SLCA = SLCOS[j][pos];
    double SLSA = SLSIN[j][pos];
    WH += SL * FLUX[0][j];
    WU += SLCA * FLR(1) - SLSA * FLR(2);
    WV += SLSA * FLR(1) + SLCA * FLR(2);
  }
}

__device__ __forceinline__ void calculate_HUV
(int t, int pos,// 原参数
const int jt, const int NHQ, const double C0, const double C1, const double HM1, const double HM2, const int DT, const double QLUA, const double VMIN,
int MBQ_LEN, int MBZQ_LEN, int MBZ_LEN, int MBW_LEN, int MDI_LEN,//数组长度
int* NV, double* H_pre, double* U_pre, double* V_pre, double* Z_pre, double* MBQ, double* DQT, double* MBZQ, double* ZB1, double* ZBC, double* MBZ, double* DZT, double* TOPW, double* MBW, double* TOPD, double* MDI, double* AREA, double* FNC,
double** KLAS, double** QT, double** SIDE, double** ZW, double** QW, double** ZT, double** NAC, double** COSF, double** SINF, double** SLCOS, double** SLSIN,
double* H_res, double* U_res, double* V_res, double* Z_res, double* W_res) {
  double SIDEX, SIDES, HSIDE, DT2, DTA, WDTA, QX1, QY1, DTAU, DTAV, WSF;
  // double H1 = H[TIME_PREV][pos];
  // double U1 = U[TIME_PREV][pos];
  // double V1 = V[TIME_PREV][pos];
  double H1 = H_pre[pos];
  double U1 = U_pre[pos];
  double V1 = V_pre[pos];

  double WH = 0.0, WU = 0.0, WV = 0.0;
  calculate_WHUV(t, pos, WH, WU, WV,
  jt, NHQ, C0, C1, HM1, HM2, VMIN,
  MBQ_LEN, MBZQ_LEN, MBZ_LEN, MBW_LEN, MDI_LEN,
  H_pre, U_pre, V_pre, Z_pre, MBQ, DQT, MBZQ, ZB1, ZBC, MBZ, DZT, TOPW, MBW, TOPD, MDI,
  KLAS, QT, SIDE, ZW, QW, ZT, NAC, COSF, SINF, SLCOS, SLSIN);

  if (NV[pos] == 4) {
    SIDEX = min(0.5 * (SIDE[0][pos] + SIDE[2][pos]),
                     0.5 * (SIDE[1][pos] + SIDE[3][pos]));
  } else {
    SIDES = 0.5 * (SIDE[0][pos] + SIDE[1][pos] + SIDE[2][pos]);
    SIDEX = sqrt((SIDES - SIDE[0][pos]) * (SIDES - SIDE[1][pos]) *
                      (SIDES - SIDE[2][pos]) / SIDES);
  }
  HSIDE = max(H1, HM1);
  DT2 = SIDEX / (U1 + sqrt(9.81 * HSIDE));
  // DT2 = std::fmin(DT, DT2);
  DT2 = fmin((double)DT, DT2);
  DT2 = max(DT2, DT / 10.0);
  DTA = 1.0 * DT2 / (1.0 * AREA[pos]);
  WDTA = 1.00 * DTA;

  double H2, U2, V2, Z2, W2;
  H2 = max(H1 - WDTA * WH + QLUA, HM1);
  // Z2 = H[TIME_NOW][pos] + ZBC[pos];
  Z2 = H2 + ZBC[pos];//Problem 0: 可能是个读旧值的bug？暂时改成了读新值，存疑
  if (H2 <= HM1) {
    U2 = 0.0;
    V2 = 0.0;
  } else {
    if (H2 <= HM2) {
      U2 = copysign( min(VMIN, fabs(U1)), U1);
      V2 = copysign( min(VMIN, fabs(V1)), V1);
    } else {
      QX1 = H1 * U1;
      QY1 = H1 * V1;
      DTAU = WDTA * WU;
      DTAV = WDTA * WV;
      WSF = FNC[pos] * sqrt(U1 * U1 + V1 * V1) / pow(H1, 0.33333);
      U2 = (QX1 - DTAU - DT * WSF * U1) / H2;
      V2 = (QY1 - DTAV - DT * WSF * V1) / H2;
      U2 = copysign( min(fabs(U2), 5.0), U2);
      V2 = copysign( min(fabs(V2), 5.0), V2);
    }
  }
  W2 = sqrt(U2 * U2 + V2 * V2);
  // H[TIME_NOW][pos] = H2;
  // U[TIME_NOW][pos] = U2;
  // V[TIME_NOW][pos] = V2;
  // Z[TIME_NOW][pos] = Z2;
  // W[TIME_NOW][pos] = W2;
  H_res[pos] = H2;
  U_res[pos] = U2;
  V_res[pos] = V2;
  Z_res[pos] = Z2;
  W_res[pos] = W2;
}

__device__ __forceinline__ void BOUNDA
(int t, int j, int pos, double (&FLUX)[4][4],double (&QL)[3],double (&QR)[3], double &FIL, double HC, double UC, double VC, double ZC,
const int jt, const int NHQ, const double C0, const double C1,
int MBQ_LEN, int MBZQ_LEN, int MBZ_LEN, int MBW_LEN, int MDI_LEN,//数组长度
double* H_pre, double* Z_pre, double* MBQ, double* DQT, double* MBZQ, double* ZBC, double* MBZ, double* DZT, double* TOPW, double* MBW, double* TOPD, double* MDI,
double** KLAS, double** QT, double** SIDE, double** ZW, double** QW, double** ZT, double** NAC, double** COSF, double** SINF) {
  // Vec WZ(NHQ);
  // Vec WQ(NHQ);
  //=====================================================
  // Problem 1: __device__中不能malloc memory,但是NHQ是从文件BOUNDARY中读取的const标量，值为5，这里先写死
  // double WZ[5];
  // double WQ[5];
  //=====================================================
  double *WZ, *WQ;
  cudaMalloc((void**)WZ, sizeof(double) * NHQ);
  cudaMalloc((void**)WQ, sizeof(double) * NHQ);

  double S0 = 0.0002;
  double DX2 = 5000.0;
  double BRDTH = 100.0;
  double CL = sqrt(9.81 * H_pre[pos]);

  if (QL[1] > CL) {
    FLUX_VAL(H_pre[pos] * QL[1], FLR(0) * QL[1],
             4.905 * H_pre[pos] * H_pre[pos], FLR(0) * QL[1]);
  }
  FLR(2) = 0;
  if (QL[1] > 0) FLR(2) = H_pre[pos] * QL[1] * QL[2];
  double HB;
//
  if (KLAS[j][pos] == 10) {
    int pos_near = find_in_vec(MBQ, MBQ_LEN, pos);
    FLR(0) = -(QT[jt][pos_near] + DQT[pos_near] * t);
    FLR(0) = FLR(0) / SIDE[j][pos];
    double QB2 = FLR(0) * FLR(0);
    double HB0 = H_pre[pos];

    for (int K = 1; K <= 20; K++) {
      double W_temp = FIL - FLR(0) / HB0;
      HB = W_temp * W_temp / 39.24;
      if (fabs(HB0 - HB) <= 0.005)
        break;
      HB0 = HB0 * 0.5 + HB * 0.5;
    }
    if (HB <= 1) {
      FLR(1) = 0;
    } else {
      FLR(1) = QB2 / HB;
    }
    FLR(2) = 0;
    FLR(3) = 4.905 * HB * HB;
  } else if (KLAS[j][pos] == 3) {
    double CQ;
    double HR;
    double W_temp;
    double HR0 = H_pre[pos];
    int pos_near = find_in_vec(MBZQ, MBZQ_LEN, pos);
    for (int i = 0; i < NHQ; i++) {
      WZ[i] = ZW[i][pos_near];
      WQ[i] = QW[i][pos_near];
    }
    for (int I1 = 0; I1 < 20; I1++) {
      double ZR0 = HR0 + ZBC[pos];
      LAQP(ZR0, CQ, WZ, WQ, NHQ);
      W_temp = FIL - CQ / HR0;
      HR = W_temp * W_temp / 39.24;
      if (fabs(HR - HR0) <= 0.001)
        break;
      HR0 = HR;
    }
    FLR(0) = CQ;
    FLR(1) = CQ * CQ / HR;
    HB = (H_pre[pos] + HR) / 2;
    FLR(3) = 4.905 * HB * HB;
  } else if (KLAS[j][pos] == 1) {
    int pos_near = find_in_vec(MBZ, MBZ_LEN, pos);
    double HB1 = ZT[jt][pos_near] + DZT[pos_near] * t - ZBC[pos];
    double FIAL = QL[2] + 6.264 * sqrt(H_pre[pos]);
    double UR0 = QL[2];
    double URB = UR0;
    for (int IURB = 1; IURB <= 30; IURB++) {
      double FIAR = URB - 6.264 * sqrt(HB1);
      URB = (FIAL + FIAR) * (FIAL - FIAR) * (FIAL - FIAR) / HB1 / 313.92;
      if (fabs(URB - UR0) <= 0.0001)
        break;
      UR0 = URB;
    }
    FLR(0) = HB1 * URB;
    FLR(1) = FLR(0) * URB;
    FLR(3) = 4.905 * HB1 * HB1;
  } else if (KLAS[j][pos] == 4) {
    FLR(0) = 0;
    FLR(1) = 0;
    FLR(2) = 0;
    FLR(3) = 4.905 * H_pre[pos] * H_pre[pos];
  } else if (KLAS[j][pos] == 5) {
    QL[1] = max(QL[1], 0.0);
    FLR(0) = H_pre[pos] * QL[1];
    FLR(1) = FLR(0) * QL[1];
    FLR(3) = 4.905 * H_pre[pos] * H_pre[pos];
  } else if (KLAS[j][pos] == 6) {
    double NE = pos;
    if (NAC[j][pos] != 0) {
      // NE = std::fmin(pos, NAC[j][pos]);
      NE = fmin((double)pos, NAC[j][pos]);
    }
    int pos_near = find_in_vec(MBW, MBW_LEN, pos);
    double TOP = TOPW[pos_near];
    if (Z_pre[pos] < TOP || ZC < TOP) {
      FLR(0) = 0;
      FLR(1) = 0;
      FLR(2) = 0;
      FLR(3) = 4.905 * H_pre[pos] * H_pre[pos];
    }
    if (Z_pre[pos] > TOP && ZC < TOP) {
      FLR(0) = C0 * pow(Z_pre[pos] - TOP, 1.5);
      FLR(1) = FLR(0) * QL[1];
      FLR(2) = FLR(0) * QL[2];
      FLR(3) = 4.905 * pow(TOP - ZBC[pos], 2);
      return;
    }
    if (Z_pre[pos] < TOP && ZC > TOP) {
      FLR(0) = -C0 * pow(ZC - TOP, 1.5);
      FLR(1) = FLR(0) *
               min(UC * COSF[j][pos] + VC * SINF[j][pos], 0.0);
      FLR(2) = FLR(0) * (VC * COSF[j][pos] - UC * SINF[j][pos]);
      FLR(3) = 4.905 * pow(Z_pre[pos] - ZBC[pos], 2);
      return;
    }
    double DZ = fabs(Z_pre[pos] - ZC);
    double HD;
    double UN;
    double VT;
    if (Z_pre[pos] <= ZC) {
      HD = Z_pre[pos] - TOP;
      UN = min(UC * COSF[j][pos] + VC * SINF[j][pos], 0.0);
      VT = VC * COSF[j][pos] - UC * SINF[j][pos];
    } else {
      HD = ZC - TOP;
      UN = max(QL[1], 0.0);
      VT = QL[2];
    }
    double SH = HD + DZ;
    double CE = min(1.0, 1.05 * pow(DZ / SH, 0.33333));
    if (Z_pre[pos] < ZC && UN > 0.0) {
      UN = 0.0;
    }
    FLR(0) =
        copysign(CE * C1 * pow(SH, 1.5), Z_pre[pos] - ZC);
    FLR(1) = FLR(0) * fabs(UN);
    FLR(2) = FLR(0) * VT;
    FLR(3) = 4.905 * pow(TOP - ZBC[pos], 2);
  } else if (KLAS[j][pos] == 7) {
    int pos_near = find_in_vec(MDI, MDI_LEN, pos);
    double TOP = TOPD[pos_near];
    if (Z_pre[pos] > TOP || ZC > TOP) {
      KLAS[j][pos] = 0;
      double CQ = QD(Z_pre[pos], ZC, TOP);
      double CB = BRDTH / SIDE[j][pos];
      FLR(0) = CQ * CB;
      FLR(1) = CB * copysign(CQ * CQ / HB, CQ);
      FLR(3) = 4.905 * HB * HB;
      return;
    } else {
      FLR(0) = 0;
      FLR(1) = 0;
      FLR(2) = 0;
      FLR(3) = 4.905 * H_pre[pos] * H_pre[pos];
    }
  }
  cudaFree((void**)WZ);
  cudaFree((void**)WQ);
}

__device__ __forceinline__ void OSHER
(int t, int pos, double (&QL)[3], double (&QR)[3], double &FIL, double (&FLR_OSHER)[4],
double* H_pre) {
  double CR = sqrt(9.81 * QR[0]);
  double FIR = QR[1] - 2 * CR;
  double UA = (FIL + FIR) / 2;
  double CA = fabs((FIL - FIR) / 4);
  double CL = sqrt(9.81 * H_pre[pos]);

  FLR_OSHER[0] = 0;
  FLR_OSHER[1] = 0;
  FLR_OSHER[2] = 0;
  FLR_OSHER[3] = 0;

  int K1, K2;

  if (CA < UA) {
    K2 = 1;
  } else if (UA >= 0.0 && UA < CA) {
    K2 = 2;
  } else if (UA >= -CA && UA < 0.0) {
    K2 = 3;
  } else if (UA < -CA) {
    K2 = 4;
  }

  if (QL[1] < CL && QR[1] >= -CR) {
    K1 = 1;
  } else if (QL[1] >= CL && QR[1] >= -CR) {
    K1 = 2;
  } else if (QL[1] < CL && QR[1] < -CR) {
    K1 = 3;
  } else if (QL[1] >= CL && QR[1] < -CR) {
    K1 = 4;
  }

  switch(K1){
    case 1:
      switch(K2){
        case 1:
          QS<2>(1, pos, QL, QR, FIL, FIR, FLR_OSHER);
          break;
        case 2:
          QS<3>(1, pos, QL, QR, FIL, FIR, FLR_OSHER);
          break;
        case 3:
          QS<5>(1, pos, QL, QR, FIL, FIR, FLR_OSHER);
          break;
        case 4:
          QS<6>(1, pos, QL, QR, FIL, FIR, FLR_OSHER);
          break;
        default:
          break;
      }
      break;
    case 2:
      switch(K2){
        case 1:
          QS<1>(1, pos, QL, QR, FIL, FIR, FLR_OSHER);
          break;
        case 2:
          QS<1>(1, pos, QL, QR, FIL, FIR, FLR_OSHER);
          QS<2>(-1, pos, QL, QR, FIL, FIR, FLR_OSHER);
          QS<3>(1, pos, QL, QR, FIL, FIR, FLR_OSHER);
          break;
        case 3:
          QS<1>(1, pos, QL, QR, FIL, FIR, FLR_OSHER);
          QS<2>(-1, pos, QL, QR, FIL, FIR, FLR_OSHER);
          QS<5>(1, pos, QL, QR, FIL, FIR, FLR_OSHER);
          break;
        case 4:
          QS<1>(1, pos, QL, QR, FIL, FIR, FLR_OSHER);
          QS<2>(-1, pos, QL, QR, FIL, FIR, FLR_OSHER);
          QS<6>(1, pos, QL, QR, FIL, FIR, FLR_OSHER);
          break;
        default:
          break;
      }
      break;
    case 3:
      switch(K2){
        case 1:
          QS<2>(1, pos, QL, QR, FIL, FIR, FLR_OSHER);
          QS<6>(-1, pos, QL, QR, FIL, FIR, FLR_OSHER);
          QS<7>(1, pos, QL, QR, FIL, FIR, FLR_OSHER);
          break;
        case 2:
          QS<3>(1, pos, QL, QR, FIL, FIR, FLR_OSHER);
          QS<6>(-1, pos, QL, QR, FIL, FIR, FLR_OSHER);
          QS<7>(1, pos, QL, QR, FIL, FIR, FLR_OSHER);
          break;
        case 3:
          QS<5>(1, pos, QL, QR, FIL, FIR, FLR_OSHER);
          QS<6>(-1, pos, QL, QR, FIL, FIR, FLR_OSHER);
          QS<7>(1, pos, QL, QR, FIL, FIR, FLR_OSHER);
          break;
        case 4:
          QS<7>(1, pos, QL, QR, FIL, FIR, FLR_OSHER);
          break;
        default:
          break;
      }
      break;
    case 4:
      switch(K2){
        case 1:
          QS<1>(1, pos, QL, QR, FIL, FIR, FLR_OSHER);
          QS<6>(-1, pos, QL, QR, FIL, FIR, FLR_OSHER);
          QS<7>(1, pos, QL, QR, FIL, FIR, FLR_OSHER);
          break;
        case 2:
          QS<1>(1, pos, QL, QR, FIL, FIR, FLR_OSHER);
          QS<2>(-1, pos, QL, QR, FIL, FIR, FLR_OSHER);
          QS<3>(1, pos, QL, QR, FIL, FIR, FLR_OSHER);
          QS<6>(-1, pos, QL, QR, FIL, FIR, FLR_OSHER);
          QS<7>(1, pos, QL, QR, FIL, FIR, FLR_OSHER);
          break;
        case 3:
          QS<1>(1, pos, QL, QR, FIL, FIR, FLR_OSHER);
          QS<2>(-1, pos, QL, QR, FIL, FIR, FLR_OSHER);
          QS<5>(1, pos, QL, QR, FIL, FIR, FLR_OSHER);
          QS<6>(-1, pos, QL, QR, FIL, FIR, FLR_OSHER);
          QS<7>(1, pos, QL, QR, FIL, FIR, FLR_OSHER);
          break;
        case 4:
          QS<1>(1, pos, QL, QR, FIL, FIR, FLR_OSHER);
          QS<2>(-1, pos, QL, QR, FIL, FIR, FLR_OSHER);
          QS<7>(1, pos, QL, QR, FIL, FIR, FLR_OSHER);
          break;
        default:
          break;
      }
      break;
    default:
      break;
  }
}

#include <cstdlib>

template<int T>
__device__ __forceinline__ void QS
(int j, int pos, const double (&QL)[3], const double (&QR)[3], double &FIL, double &FIR, double (&FLUX_OSHER)[4]){
  // Vec F(4);
  double F[4];

  if constexpr (T == 1) {
    QF(QL[0], QL[1], QL[2], F);
    for (int i = 0; i < 4; i++) {
      FLUX_OSHER[i] += F[i] * j;
    }
  } else if constexpr (T == 2) {
    double US = FIL / 3;
    double HS = US * US / 9.81;
    QF(HS, US, QL[2], F);
    for (int i = 0; i < 4; i++) {
      FLUX_OSHER[i] += F[i] * j;
    }
  } else if constexpr (T == 3) {
    FIL = FIL - (FIL + FIR) / 2;
    double HA = FIL * FIL / 39.24;
    QF(HA, (FIL + FIR) / 2, QL[2], F);
    for (int i = 0; i < 4; i++) {
      FLUX_OSHER[i] += F[i] * j;
    }
  } else if constexpr (T == 4) {

  } else if constexpr (T == 5) {
    FIR = FIR - (FIL + FIR) / 2;
    double HA = FIR * FIR / 39.24;
    QF(HA, (FIL + FIR) / 2, QR[2], F);
    for (int i = 0; i < 4; i++) {
      FLUX_OSHER[i] += F[i] * j;
    }
  } else if constexpr (T == 6) {
    double US = FIR / 3;
    double HS = US * US / 9.81;
    QF(HS, US, QR[2], F);
    for (int i = 0; i < 4; i++) {
      FLUX_OSHER[i] += F[i] * j;
    }
  } else if constexpr (T == 7) {
    QF(QR[0], QR[1], QR[2], F);
    for (int i = 0; i < 4; i++) {
      FLUX_OSHER[i] += F[i] * j;
    }
  } else {
    std::exit(1);
  }
}

__device__ __forceinline__ void QF(double h, double u, double v, double (&F)[4]) {
  // COMPUTATION OF FLUX COMPONENTS
  F[0] = h * u;
  F[1] = F[0] * u;
  F[2] = F[0] * v;
  F[3] = 4.905 * h * h;
}

// int find_in_vec(const Vec& list, double target){
//   auto it = std::find(list.begin(), list.end(), target);
//   return std::distance(list.begin(), it);
// }

// Problem 2:由于__device__中不能使用vector，使用数组则需要传入长度参数，该参数需要从__global__传入
__device__ __forceinline__ int find_in_vec(double* list, int list_length, double target){
  int index = 0;
  for(; index < list_length; ++index)
    if(list[index] == target)
      break;
  return index;
}

__device__ __forceinline__ void LAQP(double X, double &Y, double *A, double *B, double MS) {
  int ILAQ = 0;
  int i;

  for (i = 0; i < MS - 3; i++) {
    if (X < A[i + 1]) {
      ILAQ = 1;
      break;
    }
  }

  if (ILAQ == 1) {
    if (i > 0 && X - A[i] < A[i + 1] - X) {
      i = i - 1;
    }

    double X0 = A[i];
    double X1 = A[i + 1];
    double X2 = A[i + 2];

    if (fabs(X0 - X1) < 0.01 || fabs(X1 - X2) < 0.01) {
      Y = B[i + 1];
    } else {
      double U_ = (X - X1) * (X - X2) / (X0 - X1) / (X0 - X2);
      double V_ = (X - X0) * (X - X2) / (X1 - X0) / (X1 - X2);
      double W_ = (X - X0) * (X - X1) / (X2 - X0) / (X2 - X1);
      Y = U_ * B[i] + V_ * B[i + 1] + W_ * B[i + 2];
    }
  } else {
    i = MS - 2;

    if (i > 0 && X - A[i] < A[i + 1] - X) {
      i = i - 1;
    }

    double X0 = A[i];
    double X1 = A[i + 1];
    double X2 = A[i + 2];

    if (fabs(X0 - X1) < 0.01 || fabs(X1 - X2) < 0.01) {
      Y = B[i + 1];
    } else {
      double U_ = (X - X1) * (X - X2) / (X0 - X1) / (X0 - X2);
      double V_ = (X - X0) * (X - X2) / (X1 - X0) / (X1 - X2);
      double W_ = (X - X0) * (X - X1) / (X2 - X0) / (X2 - X1);
      Y = U_ * B[i] + V_ * B[i + 1] + W_ * B[i + 2];
    }
  }
}

__device__ __forceinline__ double QD(double ZL, double ZR, double ZB) {
  const double CM = 0.384;
  const double SIGMA = 0.667;
  const double FI = 4.43;

  double ZU = max(ZL, ZR);
  double ZD = min(ZL, ZR);
  double H0 = ZU - ZB;
  double HS = ZD - ZB;
  double DELTA = HS / H0;

  double QD;

  if (DELTA <= SIGMA) {
    QD = copysign(CM * pow(H0, 1.5), ZL - ZR);
  } else {
    double DH = ZU - ZD;
    if (DH > 0.09) {
      QD = copysign(FI * HS * sqrt(DH), ZL - ZR);
    } else {
      QD = copysign(FI * HS * 0.3 * DH / 0.1, ZL - ZR);
    }
  }

  return QD;
}

// Kernel函数，具体可以放到其他核上跑
// 0 <= pos <= CEL
// void time_step(int t, int pos) {
//   calculate_HUV(t, pos);
// }

// #include <omp.h>

void closeFile(){
  using namespace DataManager;
  ZUV_file.close();
  H2U2V2_file.close();
  XY_TEC_file.close();
}

__global__ void translate_step(
    // int t,  const double QLUA, const double VMIN, const double C0, const double C1,
    // int MBQ_LEN, int MBZQ_LEN, int MBZ_LEN, int MBW_LEN, int MDI_LEN,//数组长度
    const int CEL, const int DT, const int jt, const int NHQ, const double HM1, const double HM2,   // 常量
    double* H_pre, double* U_pre, double* V_pre, double* Z_pre, double* W_pre,  // 前一时刻结果
    int* NV, double* AREA, double* ZBC, double* ZB1, double* DQT, double* DZT, double* TOPW, double* TOPD, double* MBQ, double* MBZQ, double* MBW, double* MDI, //double* FNC, double* MBZ, 
    double* d_QT, double* d_ZT,
    double** SIDE, double** SLCOS, double** SLSIN, double** KLAS, double** ZW, double** QW, int** NAC, //double** QT, double** ZT, double** COSF, double** SINF,
    double* H_res, double* U_res, double* V_res, double* Z_res, double* W_res){   // 更新结果
    
  int pos = threadIdx.x+blockDim.x*blockIdx.x;// 线程pos计算HUVWZ[pos]位置的格子
  // calculate_HUV
//   calculate_HUV(
//     t, pos,
//     jt, NHQ, C0, C1, HM1, HM2, DT, QLUA, VMIN,
//     MBQ_LEN, MBZQ_LEN, MBZ_LEN, MBW_LEN, MDI_LEN,
//     NV, H_pre, U_pre, V_pre, Z_pre, MBQ, DQT, MBZQ, ZB1, ZBC, MBZ, DZT, TOPW, MBW, TOPD, MDI, AREA, FNC,
//     KLAS, QT, SIDE, ZW, QW, ZT, NAC, COSF, SINF, SLCOS, SLSIN,
//     H_res, U_res, V_res, Z_res, W_res);
    // 数据读取测试
    if(pos == CEL - 1){
        printf("SIDE[0][idx] =  %f \n", SIDE[0][pos]);
      }
  }


int main() {
  //1. 数据初始化
  READ_DATA("TIME", MDT, NDAYS)
  READ_DATA("GIRD", NOD, CEL)
  READ_DATA("DEPTH", HM1, HM2)
  READ_DATA("BOUNDARY", NZ, NQ, dummy, NHQ, dummy, NDI)
  READ_DATA("CALTIME", dummy, DT)
  DZT.resize(NZ, 0);
  DQT.resize(NQ, 0);
  W.resize(2, Vec(CEL, 0));
  pre2();
  take_boundary_for_two_d();

  ZUV_file.open("ZUV.OUT", std::ios::out | std::ofstream::trunc);
  H2U2V2_file.open("H2U2V2.OUT", std::ios::out | std::ofstream::trunc);
  XY_TEC_file.open("XY-TEC.DAT", std::ios::out | std::ofstream::trunc);
  atexit(closeFile);

  int K0 = MDT / DT;
  int pos;

  // 跟时间无关的数据迁移
  initialize_deviceVar();

  printf("grid: %d, block: %d\n", (CEL + block_size -1)/block_size, block_size);
  for (jt = 0; jt < 100; jt++) {

    for (int l = 0; l < NZ; l++) {
      if (jt != NDAYS) {
        DZT[l] = (ZT[jt + 1][l] - ZT[jt][l]) / K0;
      }
    }

    for (int l = 1; l < NQ; l++) {
      if (jt != NDAYS) {
        DQT[l] = (QT[jt + 1][l] - QT[jt][l]) / K0;
      }
    }

    // 跟时间 jt 相关的数据迁移
    if(NQ != 0){
      cudaMemcpy(DQT_dev, DQT.data(), sizeof(double)* NQ, cudaMemcpyHostToDevice);
      cudaMemcpy(QT_dev, QT[jt].data(), sizeof(double) * NDAYS, cudaMemcpyHostToDevice);
    }
    if(NZ != 0){
      cudaMemcpy(DZT_dev, DZT.data(), sizeof(double)* NZ, cudaMemcpyHostToDevice);
      cudaMemcpy(ZT_dev, ZT[jt].data(), sizeof(double) * NDAYS, cudaMemcpyHostToDevice);
    }

    // start_time = omp_get_wtime();
    clock_t start_time = clock();
    for (kt = 1; kt <= K0; kt++) {

        cudaMemcpy(H_dev, &H[TIME_PREV][0], sizeof(double) * CEL, cudaMemcpyHostToDevice);
        cudaMemcpy(U_dev, &U[TIME_PREV][0], sizeof(double) * CEL, cudaMemcpyHostToDevice);
        cudaMemcpy(V_dev, &V[TIME_PREV][0], sizeof(double) * CEL, cudaMemcpyHostToDevice);
        cudaMemcpy(Z_dev, &Z[TIME_PREV][0], sizeof(double) * CEL, cudaMemcpyHostToDevice);
        cudaMemcpy(W_dev, &W[TIME_PREV][0], sizeof(double) * CEL, cudaMemcpyHostToDevice);

        // 2. kernel 调用
        

        translate_step<<<(CEL + block_size -1)/block_size, block_size>>>( 
            // todo: 补充前两行参数
            CEL, DT, jt, NHQ, HM1, HM2,
            H_dev, U_dev, V_dev, Z_dev, W_dev,
            NV_dev, AREA_dev, ZBC_dev, ZB1_dev, DQT_dev, DZT_dev, TOPW_dev, TOPD_dev, MBQ_dev, MBZQ_dev, MBW_dev, MDI_dev,
            QT_dev, ZT_dev,
            SIDE_devIdx, SLCOS_devIdx, SLSIN_devIdx, KLAS_devIdx, ZW_devIdx, QW_devIdx, NAC_devIdx,
            H_res, U_res, V_res, Z_res, W_res
        );
        // // 3. 结果处理
        cudaDeviceSynchronize(); // 等待所有的核函数执行完毕

        cudaMemcpy(H[TIME_NOW].data(), H_res, sizeof(double) * CEL, cudaMemcpyDeviceToHost);
        cudaMemcpy(U[TIME_NOW].data(), U_res, sizeof(double) * CEL, cudaMemcpyDeviceToHost);
        cudaMemcpy(V[TIME_NOW].data(), V_res, sizeof(double) * CEL, cudaMemcpyDeviceToHost);
        cudaMemcpy(Z[TIME_NOW].data(), Z_res, sizeof(double) * CEL, cudaMemcpyDeviceToHost);
        cudaMemcpy(W[TIME_NOW].data(), W_res, sizeof(double) * CEL, cudaMemcpyDeviceToHost);

    }
    // end_time = omp_get_wtime();
    // cout << end_time - start_time << " seconds" << endl;
    clock_t end_time = clock();
    double duration = ((double)(end_time-start_time))/CLOCKS_PER_SEC;
    printf("Done. Compute took %f seconds\n", duration);
    
    output_data(kt);
  }
  return 0;
}