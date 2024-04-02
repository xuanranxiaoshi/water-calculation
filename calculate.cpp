#include "read_data.hh"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <pthread.h>
#include <sstream>

// TIME -> Y -> X
// const value
const double C0 = 1.33;
const double C1 = 1.7;
const double VMIN = 0.001;
const double QLUA = 0.0;
const string data_root_path = "./origin/changjiangkou/";

// Time scalar
int jt;
int kt;
int CEL;
int NDAYS;
int NDI;
int NOD;
int MDT;
int DT;
int NNZ0;
int NNQ0;

// scalar
double HM1, HM2;
double BI;
double NHQ;
int NZ;
double NQ;
double STIME;

// no-state matrix
Vec WH;
Vec WU;
Vec WV;
Vec ZB1;
Vec ZBC;
Vec QL;
Vec QR;
Vec MBQ;
Vec MBZ;
Vec MBW;
Vec MDI;
vector<int> NV;
Vec NNZ;
Vec NNQ;
Vec AREA;
Vec FNC;
Vec DQT;
Vec MBZQ;
Vec DZT;
Vec TOPW;
Vec TOPD;
Vec XP;
Vec XP0;
Vec YP;
Vec YP0;
Vec FNC0;
Vec QZSTIME1;
Vec QZSTEMP1;
Vec QZSTIME2;
Vec QZSTEMP2;

// Origin scalar, but change to Matrix
Vec CL;
Vec FIL;
Vec FIR;
Vec ZC;
Vec HC;
Vec UC;
Vec BC;
Vec VC;

// one-dimention time-wise matrix
Vec2 H;
Vec2 U;
Vec2 V;
Vec2 Z;
Vec2 KLAS;
Vec2 NAC;
Vec2 W;
Vec2 QT;
Vec2 ZW;
Vec2 ZT;
Vec2 QW;
Vec2 NAP;

Vec2 COSF;
Vec2 SINF;

// Origin one-dimention, but change to Vec2
Vec2 FLR_OSHER;
// one-dimention j-wise matrix
Vec2 SIDE;
Vec2 SLCOS;
Vec2 SLSIN;

// other
Vec3 FLUX;

// DUMMY VALUE.
int dummy;

void calculate_FLUX(int t, int pos);
void calculate_WHUV(int t, int pos);
void calculate_HUV(int t, int pos);
void BOUNDA(int t, int j, int pos);
void OSHER(int t, int pos);
void CHOICE(Vec list, double target, int &index);
void LAQP(double X, double &Y, Vec A, Vec B, double MS);
double QD(double ZL, double ZR, double ZB);
void QS(int k, int j, int pos);
void QF(double H, double U, double V, Vec &F);
double BOUNDRYinterp(double THOURS, int NZQSTEMP, Vec ZQSTIME, Vec ZQSTEMP);

void calculate_FLUX(int t, int pos) {
  for (int j = 0; j < 4; j++) {
    // 局部变量初始化运算。
    double KP = KLAS[j][pos];
    double NC = NAC[j][pos];
    QL[0] = H[TIME_PREV][pos];
    QL[1] = U[TIME_PREV][pos] * COSF[j][pos] + V[TIME_PREV][pos] * SINF[j][pos];
    QL[2] = V[TIME_PREV][pos] * COSF[j][pos] - U[TIME_PREV][pos] * SINF[j][pos];
    CL[pos] = std::sqrt(9.81 * H[TIME_PREV][pos]);
    FIL[pos] = QL[1] + 2 * CL[pos];
    double ZI = std::fmax(Z[TIME_PREV][pos], ZB1[pos]);
    if (NC == 0) {
      HC[pos] = 0;
      BC[pos] = 0;
      ZC[pos] = 0;
      UC[pos] = 0;
      VC[pos] = 0;
    } else {
      HC[pos] = std::fmax(H[TIME_PREV][NC], HM1);
      BC[pos] = ZBC[NC];
      ZC[pos] = std::fmax(ZBC[NC], Z[TIME_PREV][NC]);
      UC[pos] = U[TIME_PREV][NC];
      VC[pos] = V[TIME_PREV][NC];
    }

    if (KP >= 1 && KP <= 8 || KP >= 10) {
      BOUNDA(t, j, pos);
    } else if (H[TIME_PREV][pos] <= HM1 && HC[pos] <= HM1) {
      FLUX_VAL(0, 0, 0, 0);
    } else if (ZI <= BC[pos]) {
      FLUX_VAL(-C1 * pow(HC[pos], 1.5), H[TIME_PREV][pos] * QL[1] * fabs(QL[1]),
               0, 4.905 * H[TIME_PREV][pos] * H[TIME_PREV][pos]);
    } else if (ZC[pos] <= BI) {
      FLUX_VAL(C1 * pow(H[TIME_PREV][pos], 1.5), FLR(0) * QL[1], FLR(0) * QL[2],
               0);
    } else if (H[TIME_PREV][pos] <= HM2) {
      if (ZC[pos] > ZI) {
        double DH = std::fmax(ZC[pos] - BI, HM1);
        double UN = -C1 * std::sqrt(DH);
        FLUX_VAL(DH * UN, FLR(0) * UN,
                 FLR(0) * (VC[pos] * COSF[j][pos] - UC[pos] * SINF[j][pos]),
                 4.905 * H[TIME_PREV][pos] * H[TIME_PREV][pos]);
      } else {
        FLUX_VAL(-C1 * pow(HC[pos], 1.5), 0, 0,
                 4.905 * H[TIME_PREV][pos] * H[TIME_PREV][pos]);
      }
    } else if (HC[pos] <= HM2) {
      if (ZI > ZC[pos]) {
        double DH = std::fmax(ZC[pos] - BI, HM1);
        double UN = -C1 * std::sqrt(DH);
        double HC1 = ZC[pos] - BI;
        FLUX_VAL(DH * UN, FLR(0) * UN, FLR(0) * QL[2], 4.905 * HC1 * HC1);
      } else {
        FLUX_VAL(-HC[pos] * C1 * sqrt(HC[pos]),
                 H[TIME_PREV][pos] * QL[1] * QL[1], 0,
                 4.905 * H[TIME_PREV][pos] * H[TIME_PREV][pos]);
      }
    } else {
      QR[0] = std::max(ZC[pos] - BI, HM1);
      double UR = UC[pos] * COSF[j][pos] + VC[pos] * SINF[j][pos];
      QR[1] = UR * std::min(HC[pos] / QR[0], 1.5);
      if (HC[pos] <= HM2 || QR[0] <= HM2) {
        QR[1] = std::copysign(VMIN, UR);
      }
      QR[2] = VC[pos] * COSF[j][pos] - UC[pos] * SINF[j][pos];
      OSHER(t, pos);
      FLR(1) = FLR_OSHER[1][pos] + (1 - std::min(HC[pos] / QR[0], 1.5)) * HC[pos] * UR * UR / 2;
      FLUX_VAL(FLR_OSHER[0][pos], FLR(1), FLR_OSHER[2][pos], FLR_OSHER[3][pos]);
    }
  }
}

void calculate_WHUV(int t, int pos) {
  WH[pos] = 0;
  WU[pos] = 0;
  WV[pos] = 0;

  for (int j = 0; j < 4; j++) {
    FLR(1) = FLUX[1][j][pos] + FLUX[3][j][pos];
    FLR(2) = FLUX[2][j][pos];
    double SL = SIDE[j][pos];
    double SLCA = SLCOS[j][pos];
    double SLSA = SLSIN[j][pos];
    WH[pos] += SL * FLUX[0][j][pos];
    WU[pos] += SLCA * FLR(1) - SLSA * FLR(2);
    WV[pos] += SLSA * FLR(1) - SLCA * FLR(2);
  }
}

void calculate_HUV(int t, int pos) {
  double SIDEX, SIDES, HSIDE, DT2, DTA, WDTA;
  double QX1, QY1, DTAU, DTAV, FNCC, WSF;

  if (NV[pos] == 4) {
    SIDEX = std::min(0.5 * (SIDE[0][pos] + SIDE[2][pos]),
                     0.5 * (SIDE[1][pos] + SIDE[3][pos]));
  } else {
    SIDES = 0.5 * (SIDE[0][pos] + SIDE[1][pos] + SIDE[2][pos]);
    SIDEX = std::sqrt((SIDES - SIDE[0][pos]) * (SIDES - SIDE[1][pos]) *
                      (SIDES - SIDE[2][pos]) / SIDES);
  }
  HSIDE = std::max(H[TIME_PREV][pos], HM1);
  DT2 = SIDEX / (U[TIME_PREV][pos] + std::sqrt(9.81 * HSIDE));
  DT2 = std::fmin(DT, DT2);
  DT2 = std::max(DT2, DT / 10.0);
  DTA = 1.0 * DT2 / (1.0 * AREA[pos]);
  WDTA = 1.00 * DTA;

  H[TIME_NOW][pos] = std::max(H[TIME_PREV][pos] - WDTA * WH[pos] + QLUA, HM1);
  Z[TIME_NOW][pos] = H[TIME_NOW][pos] + ZBC[pos];
  if (H[TIME_NOW][pos] <= HM1) {
    U[TIME_NOW][pos] = 0.0;
    V[TIME_NOW][pos] = 0.0;
  } else {
    if (H[TIME_NOW][pos] <= HM2) {
      U[TIME_NOW][pos] = std::copysign(
          std::min(VMIN, std::abs(U[TIME_PREV][pos])), U[TIME_PREV][pos]);
      V[TIME_NOW][pos] = std::copysign(
          std::min(VMIN, std::abs(V[TIME_PREV][pos])), V[TIME_PREV][pos]);
    } else {
      QX1 = H[TIME_PREV][pos] * U[TIME_PREV][pos];
      QY1 = H[TIME_PREV][pos] * V[TIME_PREV][pos];
      DTAU = WDTA * WU[pos];
      DTAV = WDTA * WV[pos];

      FNCC = FNC[pos];
      WSF = FNCC *
            std::sqrt(U[TIME_PREV][pos] * U[TIME_PREV][pos] +
                      V[TIME_PREV][pos] * V[TIME_PREV][pos]) /
            std::pow(H[TIME_PREV][pos], 0.33333);

      U[TIME_NOW][pos] =
          (QX1 - DTAU - DT * WSF * U[TIME_PREV][pos]) / H[TIME_NOW][pos];
      V[TIME_NOW][pos] =
          (QY1 - DTAV - DT * WSF * V[TIME_PREV][pos]) / H[TIME_NOW][pos];

      U[TIME_NOW][pos] = std::copysign(
          std::min(std::abs(U[TIME_NOW][pos]), 5.0), U[TIME_NOW][pos]);
      V[TIME_NOW][pos] = std::copysign(
          std::min(std::abs(V[TIME_NOW][pos]), 5.0), V[TIME_NOW][pos]);
    }
  }
  W[TIME_NOW][pos] = sqrt(U[TIME_NOW][pos] * U[TIME_NOW][pos] +
                          V[TIME_NOW][pos] * V[TIME_NOW][pos]);
}

void BOUNDA(int t, int j, int pos) {
  Vec WZ;
  Vec WQ;
  Vec QB;

  double S0 = 0.0002;
  double DX2 = 5000.0;
  double BRDTH = 100.0;
  int KP;

  if (QL[1] > CL[pos]) {
    FLUX_VAL(H[TIME_PREV][pos] * QL[1], FLR(0) * QL[1],
             4.905 * H[TIME_PREV][pos] * H[TIME_PREV][pos], FLR(0) * QL[1]);
  }
  FLR(2) = 0;
  if (QL[1] > 0) FLR(2) = H[TIME_PREV][pos] * QL[1] * QL[2]; 
  int II;
  double HB;
  if (KP == 10) {
    CHOICE(MBQ, pos, II);
    FLR(0) = -(QT[jt][II] + DQT[II] * t);
    FLR(0) = FLR(0) / SIDE[j][pos];
    double QB2 = FLR(0) * FLR(0);
    double HB0 = H[TIME_PREV][pos];

    for (int K = 1; K <= 20; K++) {
      double W_temp = FIL[pos] - FLR(0) / HB0;
      HB = W_temp * W_temp / 39.24;
      if (std::abs(HB0 - HB) <= 0.005)
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
  } else if (KP == 3) {
    double CQ;
    double HR;
    double W;
    double HR0 = H[TIME_PREV][pos];
    CHOICE(MBZQ, pos, II);
    for (int i = 0; i < NHQ; i++) {
      WZ[i] = ZW[i][II];
      WQ[i] = QW[i][II];
    }
    for (int I1 = 0; I1 < 20; I1++) {
      double ZR0 = HR0 + BI;
      LAQP(ZR0, CQ, WZ, WQ, NHQ);
      W = FIL[pos] - CQ / HR0;
      HR = W * W / 39.24;
      if (std::abs(HR - HR0) <= 0.001)
        break;
      HR0 = HR;
    }
    FLR(0) = CQ;
    FLR(1) = CQ * CQ / HR;
    HB = (H[TIME_PREV][pos] + HR) / 2;
    FLR(3) = 4.905 * HB * HB;
  } else if (KP == 1) {
    CHOICE(MBZ, pos, II);
    double HB1 = ZT[jt][II] + DZT[II] * t - BI;
    double FIAL = QL[2] + 6.264 * sqrt(H[TIME_PREV][pos]);
    double UR0 = QL[2];
    double URB = UR0;
    for (int IURB = 1; IURB <= 30; IURB++) {
      double FIAR = URB - 6.264 * sqrt(HB1);
      URB = (FIAL + FIAR) * (FIAL - FIAR) * (FIAL - FIAR) / HB1 / 313.92;
      if (std::abs(URB - UR0) <= 0.0001)
        break;
      UR0 = URB;
    }
    FLR(0) = HB1 * URB;
    FLR(1) = FLR(0) * URB;
    FLR(3) = 4.905 * HB1 * HB1;
  } else if (KP == 4) {
    FLR(0) = 0;
    FLR(1) = 0;
    FLR(2) = 0;
    FLR(3) = 4.905 * H[TIME_PREV][pos] * H[TIME_PREV][pos];
  } else if (KP == 5) {
    QL[1] = std::max(QL[1], 0.0);
    FLR(0) = H[TIME_PREV][pos] * QL[1];
    FLR(1) = FLR(0) * QL[1];
    FLR(3) = 4.905 * H[TIME_PREV][pos] * H[TIME_PREV][pos];
  } else if (KP == 6) {
    double NC = NAC[j][pos];
    double NE;
    if (NC != 0) {
      NE = std::fmin(pos, NC);
    }
    CHOICE(MBW, NE, II);
    double TOP = TOPW[II];
    if (Z[TIME_PREV][pos] < TOP || ZC[pos] < TOP) {
      KP = 4;
      FLR(0) = 0;
      FLR(1) = 0;
      FLR(2) = 0;
      FLR(3) = 4.905 * H[TIME_PREV][pos] * H[TIME_PREV][pos];
    }
    if (Z[TIME_PREV][pos] > TOP && ZC[pos] < TOP) {
      FLR(0) = C0 * std::pow(Z[TIME_PREV][pos] - TOP, 1.5);
      FLR(1) = FLR(0) * QL[1];
      FLR(2) = FLR(0) * QL[2];
      FLR(3) = 4.905 * std::pow(TOP - BI, 2);
      return;
    }
    if (Z[TIME_PREV][pos] < TOP && ZC[pos] > TOP) {
      FLR(0) = -C0 * std::pow(ZC[pos] - TOP, 1.5);
      FLR(1) = FLR(0) *
               std::min(UC[pos] * COSF[j][pos] + VC[pos] * SINF[j][pos], 0.0);
      FLR(2) = FLR(0) * (VC[pos] * COSF[j][pos] - UC[pos] * SINF[j][pos]);
      FLR(3) = 4.905 * std::pow(Z[TIME_PREV][pos] - BI, 2);
      return;
    }
    double DZ = std::abs(Z[TIME_PREV][pos] - ZC[pos]);
    double HD;
    double UN;
    double VT;
    if (Z[TIME_PREV][pos] <= ZC[pos]) {
      HD = Z[TIME_PREV][pos] - TOP;
      UN = std::min(UC[pos] * COSF[j][pos] + VC[pos] * SINF[j][pos], 0.0);
      VT = VC[pos] * COSF[j][pos] - UC[pos] * SINF[j][pos];
    } else {
      HD = ZC[pos] - TOP;
      UN = std::max(QL[1], 0.0);
      VT = QL[2];
    }
    double SH = HD + DZ;
    double CE = std::min(1.0, 1.05 * std::pow(DZ / SH, 0.33333));
    if (Z[TIME_PREV][pos] < ZC[pos] && UN > 0.0) {
      UN = 0.0;
    }
    FLR(0) =
        std::copysign(CE * C1 * std::pow(SH, 1.5), Z[TIME_PREV][pos] - ZC[pos]);
    FLR(1) = FLR(0) * std::abs(UN);
    FLR(2) = FLR(0) * VT;
    FLR(3) = 4.905 * std::pow(TOP - BI, 2);
  } else if (KP == 7) {
    CHOICE(MDI, pos, II);
    double TOP = TOPD[II];
    if (Z[TIME_PREV][pos] > TOP || ZC[pos] > TOP) {
      KP = 0;
      KLAS[j][pos] = 0;
      double CQ = QD(Z[TIME_PREV][pos], ZC[pos], TOP);
      double CB = BRDTH / SIDE[j][pos];
      FLR(0) = CQ * CB;
      FLR(1) = CB * std::copysign(CQ * CQ / HB, CQ);
      FLR(3) = 4.905 * HB * HB;
      return;
    } else {
      KP = 4;
      FLR(0) = 0;
      FLR(1) = 0;
      FLR(2) = 0;
      FLR(3) = 4.905 * H[TIME_PREV][pos] * H[TIME_PREV][pos];
    }
  }
}

void OSHER(int t, int pos) {
  double CR = sqrt(9.81 * QR[0]);
  FIR[pos] = QR[1] - 2 * CR;
  double UA = (FIL[pos] + FIR[pos]) / 2;
  double CA = fabs((FIL[pos] - FIR[pos]) / 4);

  int K1, K2;
  if (QL[1] < CL[pos] && QR[1] >= -CR) {
    K1 = 1;
  } else {
    if (CA < UA) {
      K2 = 10;
    }
    if (UA >= 0.0 && UA < CA) {
      K2 = 20;
    }
    if (UA >= -CA && UA < 0.0) {
      K2 = 30;
    }
    if (UA < -CA) {
      K2 = 40;
    }
  }

  if (QL[1] >= CL[pos] && QR[1] >= -CR) {
    K1 = 2;
  } else {
    if (CA < UA) {
      K2 = 10;
    }
    if (UA >= 0.0 && UA < CA) {
      K2 = 20;
    }
    if (UA >= -CA && UA < 0.0) {
      K2 = 30;
    }
    if (UA < -CA) {
      K2 = 40;
    }
  }

  if (QL[1] < CL[pos] && QR[1] < -CR) {
    K1 = 3;
  } else {
    if (CA < UA) {
      K2 = 10;
    }
    if (UA >= 0.0 && UA < CA) {
      K2 = 20;
    }
    if (UA >= -CA && UA < 0.0) {
      K2 = 30;
    }
    if (UA < -CA) {
      K2 = 40;
    }
  }

  if (QL[1] >= CL[pos] && QR[1] < -CR) {
    K1 = 4;
  } else {
    if (CA < UA) {
      K2 = 10;
    }
    if (UA >= 0.0 && UA < CA) {
      K2 = 20;
    }
    if (UA >= -CA && UA < 0.0) {
      K2 = 30;
    }
    if (UA < -CA) {
      K2 = 40;
    }
  }

  int K12 = K1 + K2;

  if (K12 == 11) {
    QS(2, 1, pos);
    return;
  }

  if (K12 == 21) {
    QS(3, 1, pos);
    return;
  }

  if (K12 == 31) {
    QS(5, 1, pos);
    return;
  }

  if (K12 == 41) {
    QS(6, 1, pos);
    return;
  }

  if (K12 == 12) {
    QS(1, 1, pos);
    return;
  }

  if (K12 == 22) {
    QS(1, 1, pos);
    QS(2, -1, pos);
    QS(3, 1, pos);
    return;
  }

  if (K12 == 32) {
    QS(1, 1, pos);
    QS(2, -1, pos);
    QS(5, 1, pos);
    return;
  }

  if (K12 == 42) {
    QS(1, 1, pos);
    QS(2, -1, pos);
    QS(6, 1, pos);
    return;
  }

  if (K12 == 13) {
    QS(2, 1, pos);
    QS(6, -1, pos);
    QS(7, 1, pos);
    return;
  }

  if (K12 == 23) {
    QS(3, 1, pos);
    QS(6, -1, pos);
    QS(7, 1, pos);
    return;
  }

  if (K12 == 33) {
    QS(5, 1, pos);
    QS(6, -1, pos);
    QS(7, 1, pos);
    return;
  }

  if (K12 == 43) {
    QS(7, 1, pos);
    return;
  }

  if (K12 == 14) {
    QS(1, 1, pos);
    QS(6, -1, pos);
    QS(7, 1, pos);
    return;
  }

  if (K12 == 24) {
    QS(1, 1, pos);
    QS(2, -1, pos);
    QS(3, 1, pos);
    QS(6, -1, pos);
    QS(7, 1, pos);
    return;
  }

  if (K12 == 34) {
    QS(1, 1, pos);
    QS(2, -1, pos);
    QS(5, 1, pos);
    QS(6, -1, pos);
    QS(7, 1, pos);
    return;
  }

  if (K12 == 44) {
    QS(1, 1, pos);
    QS(2, -1, pos);
    QS(7, 1, pos);
  }
}

void QS(int k, int j, int pos) {
  // COMPUTATION OF FLUX BASED ON LEFT, RIGHT OR INTERMEDIATE STATE
  Vec F(4);

  if (k == 1) {
    QF(QL[0], QL[1], QL[2], F);
    for (int i = 0; i < 4; i++) {
      FLR_OSHER[i][pos] += F[i] * j;
    }
  }

  if (k == 2) {
    double US = FIL[pos] / 3;
    double HS = US * US / 9.81;
    QF(HS, US, QL[2], F);
    for (int i = 0; i < 4; i++) {
      FLR_OSHER[i][pos] += F[i] * j;
    }
  }

  if (k == 3) {
    FIL[pos] = FIL[pos] - (FIL[pos] + FIR[pos]) / 2;
    double HA = FIL[pos] * FIL[pos] / 39.24;
    QF(HA, (FIL[pos] + FIR[pos]) / 2, QL[2], F);
    for (int i = 0; i < 4; i++) {
      FLR_OSHER[i][pos] += F[i] * j;
    }
  }

  if (k == 5) {
    FIR[pos] = FIR[pos] - (FIL[pos] + FIR[pos]) / 2;
    double HA = FIR[pos] * FIR[pos] / 39.24;
    QF(HA, (FIL[pos] + FIR[pos]) / 2, QR[2], F);
    for (int i = 0; i < 4; i++) {
      FLR_OSHER[i][pos] += F[i] * j;
    }
  }

  if (k == 6) {
    double US = FIR[pos] / 3;
    double HS = US * US / 9.81;
    QF(HS, US, QR[2], F);
    for (int i = 0; i < 4; i++) {
      FLR_OSHER[i][pos] += F[i] * j;
    }
  }

  if (k == 7) {
    QF(QR[0], QR[1], QR[2], F);
    for (int i = 0; i < 4; i++) {
      FLR_OSHER[i][pos] += F[i] * j;
    }
  }
}

void QF(double H, double U, double V, Vec &F) {
  // COMPUTATION OF FLUX COMPONENTS
  F[0] = H * U;
  F[1] = F[0] * U;
  F[2] = F[0] * V;
  F[3] = 4.905 * H * H;
}

void CHOICE(Vec list, double target, int &index) {
  auto it = std::find(list.begin(), list.end(), target);
  index = std::distance(list.begin(), it);
}

void LAQP(double X, double &Y, Vec A, Vec B, double MS) {
  const int NHQ = 5;
  int ILAQ = 0;
  int I;

  for (I = 0; I < MS - 3; I++) {
    if (X < A[I + 1]) {
      ILAQ = 1;
      break;
    }
  }

  if (ILAQ == 1) {
    if (I > 0 && X - A[I] < A[I + 1] - X) {
      I = I - 1;
    }

    double X0 = A[I];
    double X1 = A[I + 1];
    double X2 = A[I + 2];

    if (std::abs(X0 - X1) < 0.01 || std::abs(X1 - X2) < 0.01) {
      Y = B[I + 1];
    } else {
      double U = (X - X1) * (X - X2) / (X0 - X1) / (X0 - X2);
      double V = (X - X0) * (X - X2) / (X1 - X0) / (X1 - X2);
      double W = (X - X0) * (X - X1) / (X2 - X0) / (X2 - X1);
      Y = U * B[I] + V * B[I + 1] + W * B[I + 2];
    }
  } else {
    I = MS - 2;

    if (I > 0 && X - A[I] < A[I + 1] - X) {
      I = I - 1;
    }

    double X0 = A[I];
    double X1 = A[I + 1];
    double X2 = A[I + 2];

    if (std::abs(X0 - X1) < 0.01 || std::abs(X1 - X2) < 0.01) {
      Y = B[I + 1];
    } else {
      double U = (X - X1) * (X - X2) / (X0 - X1) / (X0 - X2);
      double V = (X - X0) * (X - X2) / (X1 - X0) / (X1 - X2);
      double W = (X - X0) * (X - X1) / (X2 - X0) / (X2 - X1);
      Y = U * B[I] + V * B[I + 1] + W * B[I + 2];
    }
  }
}

double QD(double ZL, double ZR, double ZB) {
  const double CM = 0.384;
  const double SIGMA = 0.667;
  const double FI = 4.43;

  double ZU = std::max(ZL, ZR);
  double ZD = std::min(ZL, ZR);
  double H0 = ZU - ZB;
  double HS = ZD - ZB;
  double DELTA = HS / H0;

  double QD;

  if (DELTA <= SIGMA) {
    QD = std::copysign(CM * std::pow(H0, 1.5), ZL - ZR);
  } else {
    double DH = ZU - ZD;
    if (DH > 0.09) {
      QD = std::copysign(FI * HS * std::sqrt(DH), ZL - ZR);
    } else {
      QD = std::copysign(FI * HS * 0.3 * DH / 0.1, ZL - ZR);
    }
  }

  return QD;
}

// Kernel函数，具体可以放到其他核上跑
// 0 <= pos <= CEL
void time_step(int t, int pos) {
  calculate_FLUX(t, pos);
  calculate_WHUV(t, pos);
  calculate_HUV(t, pos);
}

template <typename T> T readFromLine(const std::string &line) {
  std::istringstream iss(line);
  T result;
  iss >> result;
  return result;
}

void loadFromFilePNAC(const std::string p) {
  std::ifstream file(data_root_path + p);
  ASSERT_READ(file);
  std::string line;
  std::getline(file, line);
  int NO;
  NAC.resize(4, std::vector<double>(CEL));
  for (int i = 0; i < CEL; i++) {
    std::getline(file, line);
    std::istringstream iss(line);
    iss >> NO;
    for (int j = 0; j < 4; j++) {
      iss >> NAC[j][i];
    }
  }
  file.close();
}

void loadFromFilePNAP(const std::string p) {
  std::ifstream file(data_root_path + p);
  ASSERT_READ(file);
  std::string line;
  std::getline(file, line);
  int NO;
  NAP.resize(4, std::vector<double>(CEL));
  for (int i = 0; i < CEL; i++) {
    std::getline(file, line);
    std::istringstream iss(line);
    iss >> NO;
    for (int j = 0; j < 4; j++) {
      iss >> NAP[j][i];
    }
  }
  file.close();
}

void loadFromFilePKLAS(const std::string p) {
  std::ifstream file(data_root_path + p);
  ASSERT_READ(file);
  std::string line;
  std::getline(file, line);
  int NO;
  KLAS.resize(4, std::vector<double>(CEL));
  for (int i = 0; i < CEL; i++) {
    std::getline(file, line);
    std::istringstream iss(line);
    iss >> NO;
    for (int j = 0; j < 4; j++) {
      iss >> KLAS[j][i];
    }
  }
  file.close();
}

void loadFromFilePZBC(const std::string p) {
  std::ifstream file(data_root_path + p);
  ASSERT_READ(file);
  std::string line;
  std::getline(file, line);
  int NO;
  ZBC.resize(CEL, 0);
  for (int i = 0; i < CEL; i++) {
    std::getline(file, line);
    std::istringstream iss(line);
    iss >> ZBC[i];
  }
  file.close();
}

void loadFromFileMBZ(const std::string p) {
  std::ifstream file(data_root_path + p);
  ASSERT_READ(file);
  std::string line;
  std::getline(file, line);
  std::istringstream iss(line);
  iss >> NNZ0;
  int NO;
  MBZ.resize(CEL, 0);
  NNZ.resize(CEL, 0);
  for (int i = 0; i < CEL; i++) {
    std::getline(file, line);
    std::istringstream iss(line);
    iss >> NO;
    iss >> MBZ[i];
    iss >> NNZ[i];
  }
  file.close();
}

void loadFromFileMBQ(const std::string p) {
  std::ifstream file(data_root_path + p);
  ASSERT_READ(file);
  std::string line;
  std::getline(file, line);
  std::istringstream iss(line);
  iss >> NNQ0;
  int NO;
  MBQ.resize(CEL, 0);
  NNQ.resize(CEL, 0);
  for (int i = 0; i < CEL; i++) {
    std::getline(file, line);
    std::istringstream iss(line);
    iss >> NO;
    iss >> MBQ[i];
    iss >> NNQ[i];
  }
  file.close();
}

void loadFromFilePXY(const std::string p) {
  std::ifstream file(data_root_path + p);
  ASSERT_READ(file);
  std::string line;
  std::getline(file, line);
  int NO;
  XP.resize(CEL, 0);
  YP.resize(CEL, 0);
  for (int i = 0; i < CEL; i++) {
    std::getline(file, line);
    std::istringstream iss(line);
    iss >> NO;
    iss >> XP[i];
    iss >> YP[i];
  }
  file.close();
}

void loadFromFileInitLevel(const std::string p) {
  std::ifstream file(data_root_path + p);
  ASSERT_READ(file);
  std::string line;
  std::getline(file, line);
  int NO;
  Z.resize(2, vector<double>(CEL));
  for (int i = 0; i < CEL; i++) {
    std::getline(file, line);
    std::istringstream iss(line);
    iss >> Z[0][i];
  }
  file.close();
}

void loadFromFileU1(const std::string p) {
  std::ifstream file(data_root_path + p);
  ASSERT_READ(file);
  std::string line;
  std::getline(file, line);
  int NO;
  U.resize(2, vector<double>(CEL));
  for (int i = 0; i < CEL; i++) {
    std::getline(file, line);
    std::istringstream iss(line);
    iss >> U[0][i];
  }
  file.close();
}

void loadFromFileV1(const std::string p) {
  std::ifstream file(data_root_path + p);
  ASSERT_READ(file);
  std::string line;
  std::getline(file, line);
  int NO;
  V.resize(2, vector<double>(CEL));
  for (int i = 0; i < CEL; i++) {
    std::getline(file, line);
    std::istringstream iss(line);
    iss >> V[0][i];
  }
  file.close();
}

void loadFromFileCV(const std::string p) {
  std::ifstream file(data_root_path + p);
  ASSERT_READ(file);
  std::string line;
  std::getline(file, line);
  int NO;
  FNC0.resize(CEL, 0);
  for (int i = 0; i < CEL; i++) {
    std::getline(file, line);
    std::istringstream iss(line);
    iss >> FNC0[i];
  }
  file.close();
}

void pre2() {
  vector<int> NW;
  NW.resize(4, 0);
  loadFromFilePNAC("/SOURCES/PNAC.DAT");
  loadFromFilePNAP("/SOURCES/PNAP.DAT");
  loadFromFilePKLAS("/SOURCES/PKLAS.DAT");
  loadFromFilePZBC("/SOURCES/PZBC.DAT");
  loadFromFilePXY("/SOURCES/PXY.DAT");
  loadFromFileMBZ("/SOURCES/MBZ.DAT");
  loadFromFileMBQ("/SOURCES/MBQ.DAT");
  loadFromFileInitLevel("/INITIALLEVEL.DAT");
  loadFromFileU1("/INITIALU1.DAT");
  loadFromFileV1("/INITIALV1.DAT");
  loadFromFileCV("/CV.DAT");

  std::sort(MBZ.begin(), MBZ.end());
  std::sort(NNZ.begin(), NNZ.end());

  double XIMIN = *std::min_element(XP.begin(), XP.end());
  double YIMIN = *std::min_element(YP.begin(), YP.end());
  for (int i = 0; i < NOD; i++) {
    XP[i] -= XIMIN;
    YP[i] -= YIMIN;
  }

  ZB1.resize(CEL, 0);
  NV.resize(CEL, 0);
  NW.resize(CEL, 0);
  SIDE.resize(4, vector<double>(CEL));
  SINF.resize(4, vector<double>(CEL));
  COSF.resize(4, vector<double>(CEL));
  AREA.resize(CEL, 0);
  for (int i = 0; i < CEL; i++) {
    if (NAP[0][i] == 0) {
      continue;
    }
    ZB1[i] = ZBC[i] + HM1;
    NV[i] = 4;
    int NA = NAP[3][i];
    if (NA == 0 || NA == NAP[0][i]) {
      NV[i] = 3;
    }
    for (int j = 0; j < NV[i]; j++) {
      NW[j] = NAP[j][i];
    }
    double XP1 = XP[NW[0]];
    double XP2 = XP[NW[1]];
    double XP3 = XP[NW[2]];
    double YP1 = YP[NW[0]];
    double YP2 = YP[NW[1]];
    double YP3 = YP[NW[2]];
    AREA[i] = (XP1 * YP2 - YP1 * XP2 - XP1 * YP3 + YP1 * XP3 + XP2 * YP3 -
               YP2 * XP3) /
              2.0;
    if (NV[i] == 4) {
      double XP4 = XP[NW[3]];
      double YP4 = YP[NW[3]];
      AREA[i] = (XP1 * YP3 - YP1 * XP3 - XP1 * YP4 + YP1 * XP4 + XP3 * YP4 -
                 YP3 * XP4) /
                2.0;
    }
    for (int j = 0; j < NV[i]; j++) {
      int N1 = NW[j];
      int N2 = NW[(j + 1) % NV[i]];
      double DX = XP[N1] - XP[N2];
      double DY = YP[N2] - YP[N1];
      SIDE[j][i] = std::sqrt(DX * DX + DY * DY);
      if (SIDE[j][i] > 0.0) {
        SINF[j][i] = DX / SIDE[j][i];
        COSF[j][i] = DY / SIDE[j][i];
      }
    }
  }

  for (int i = 0; i < CEL; i++) {
    ZB1[i] = ZBC[i] + HM1;
  }

  // if (JT0 > 1)
  // {
  //   std::ifstream file("HUVZ.DAT", std::ios::binary);
  //   file.read(reinterpret_cast<char *>(&H1[0]), CEL * sizeof(int));
  //   file.read(reinterpret_cast<char *>(&U1[0]), CEL * sizeof(int));
  //   file.read(reinterpret_cast<char *>(&V1[0]), CEL * sizeof(int));
  //   file.read(reinterpret_cast<char *>(&Z1[0]), CEL * sizeof(int));
  //   file.read(reinterpret_cast<char *>(&KLAS[0]), CEL * sizeof(int));
  //   file.close();
  // }
  // else
  // {

  // }
  H.resize(2, Vec(CEL));
  SLCOS.resize(4, vector<double>(CEL));
  SLSIN.resize(4, vector<double>(CEL));
  FNC.resize(CEL, 0);
  for (int i = 0; i < CEL; i++) {
    if (NAP[0][i] == 0) {
      continue;
    }
    if (Z[0][i] <= ZBC[i]) {
      H[0][i] = HM1;
      Z[0][i] = ZB1[i];
    } else {
      H[0][i] = Z[0][i] - ZBC[i];
    }
  }
  for (int i = 0; i < CEL; i++) {
    if (NAP[0][i] == 0) {
      continue;
    }
    FNC[i] = 9.81 * FNC0[i] * FNC0[i];
    Z[1][i] = Z[0][i];
    H[1][i] = H[0][i];
    U[1][i] = U[0][i];
    V[1][i] = V[0][i];
    for (int j = 0; j < NV[i]; j++) {
      SLCOS[j][i] = SIDE[j][i] * COSF[j][i];
      SLSIN[j][i] = SIDE[j][i] * SINF[j][i];
    }
  }
}

// BUG:待修复
void take_boundary_for_two_d() {
  // INCLUDE 相关的文件或数据结构的声明和定义
  std::string pointname;
  double STIME1;

  int NZTEMP, NQTEMP;
  ZT.resize(NZ, vector<double>(NDAYS));
  for (int k = 1; k <= NNZ0; k++) {
    std::stringstream ss;
    ss << std::setfill('0') << std::setw(4) << k;
    pointname = data_root_path + "BOUNDE/NZ/NZ" + ss.str() + ".DAT";
    cout << pointname << endl;
    std::fstream NZ_file(pointname);
    ASSERT_READ(NZ_file)
    std::getline(NZ_file, current_line);
    NZTEMP = readFromLine<int>(current_line);
    QZSTIME1.resize(NZTEMP, 0);
    QZSTEMP1.resize(NZTEMP, 0);
    for (int i = 0; i < NZTEMP; i++) {
      std::getline(NZ_file, current_line);
      std::istringstream iss(current_line);
      iss >> QZSTIME1[i];
      iss >> QZSTEMP1[i];
    }
    for (int i = 0; i < NDAYS; i++) {
      STIME1 = STIME + (i - 1) / (24.0 * 3600.0 / MDT);
      for (int j = 0; j < NZ; j++) {
        if (NNZ[j] == k)
          ZT[j][i] = BOUNDRYinterp(STIME1, NZTEMP, QZSTIME1, QZSTEMP1);
      }
    }
  }

  QZSTIME1.resize(NQTEMP, 0);
  QZSTEMP1.resize(NQTEMP, 0);
  QT.resize(NQ, vector<double>(NDAYS));
  for (int k = 1; k <= NNQ0; k++) {
    std::stringstream ss;
    ss << std::setfill('0') << std::setw(4) << k;
    pointname = data_root_path + "BOUNDE/NQ/NQ" + ss.str() + ".DAT";
    cout << pointname << endl;
    std::fstream NQ_file(pointname);
    ASSERT_READ(NQ_file)
    std::getline(NQ_file, current_line);
    NQTEMP = readFromLine<int>(current_line);
    for (int i = 0; i < NQTEMP; i++) {
      std::getline(NQ_file, current_line);
      std::istringstream iss(current_line);
      iss >> QZSTIME2[i];
      iss >> QZSTEMP2[i];
    }
    for (int i = 0; i < NDAYS; i++) {
      STIME1 = STIME + (i - 1) / (24.0 * 3600.0 / MDT);
      for (int j = 0; j < NQ; j++) {
        if (NNQ[j] == k)
          QT[j][i] = BOUNDRYinterp(STIME1, NQTEMP, QZSTIME2, QZSTEMP2);
      }
    }
  }
  // for (K = 1; K <= NNQ0; K++) {
  //   sprintf(pointname, "%04d", K);
  //   FILE *file = fopen(
  //       ("../BOUNDE/NQ/NQ" + std::string(pointname) + ".dat").c_str(), "r");
  //   if (file == NULL) {
  //     // 错误处理
  //     exit(EXIT_FAILURE);
  //   }
  //   fscanf(file, "%d", &NQTEMP);
  //   for (I = 1; I <= NQTEMP; I++) {
  //     fscanf(file, "%f %f", &QZSTIME2[I], &QZSTEMP2[I]);
  //   }
  //   for (I = 1; I <= NDAYS; I++) {
  //     STIME1 = STIME + (I - 1) / (24.0 * 3600.0 / MDT);
  //     double QTTEMP;
  //     BOUNDRYinterp(STIME1, &QTTEMP, NQTEMP, QZSTIME2, QZSTEMP2);
  //     for (J = 1; J <= NQ; J++) {
  //       if (NNQ[J] == K) {
  //         QT[I][J] = QTTEMP;
  //       }
  //     }
  //   }
  //   fclose(file);
  // }
}

double BOUNDRYinterp(double THOURS, int NZQSTEMP, Vec ZQSTIME, Vec ZQSTEMP) {
  double ZQSTEMP1;
  for (int i = 0; i < NZQSTEMP; i++) {
    if (THOURS >= ZQSTIME[i] && THOURS <= ZQSTIME[i + 1]) {
      ZQSTEMP1 = ZQSTEMP[i] + (ZQSTEMP[i + 1] - ZQSTEMP[i]) /
                                  (ZQSTIME[i + 1] - ZQSTIME[i]) *
                                  (THOURS - ZQSTIME[i]);
    }
  }
  return ZQSTEMP1;
}

int main() {
  READ_DATA("TIME", MDT, NDAYS)
  READ_DATA("GIRD", NOD, CEL)
  READ_DATA("DEPTH", HM1, HM2)
  READ_DATA("BOUNDARY", NZ, NQ, dummy, NHQ, dummy, NDI)
  READ_DATA("CALTIME", dummy, DT)
  pre2();
  take_boundary_for_two_d();
  DZT.resize(NZ, 0);
  DQT.resize(NQ, 0);
  HC.resize(CEL, 0);
  BC.resize(CEL, 0);
  ZC.resize(CEL, 0);
  UC.resize(CEL, 0);
  VC.resize(CEL, 0);
  QR.resize(3, 0);
  QL.resize(3, 0);
  CL.resize(CEL, 0);
  FIL.resize(CEL, 0);
  FIR.resize(CEL, 0);
  FLUX.resize(4, Vec2(4, Vec(CEL)));
  WH.resize(CEL, 0);
  WU.resize(CEL, 0);
  WV.resize(CEL, 0);
  W.resize(2, Vec(CEL, 0));
  FLR_OSHER.resize(4, Vec(CEL, 0));
  // 必须串行执行的部分：
  // K0 = 2000
  double TAL = 0; // 将TAL（总面积）初始化为0
  double K0 = (double)MDT / DT;
  for (int II = 1; II <= CEL; II++) { // 开始一个循环，II从1到CEL
    TAL += AREA[II];
  }
  for (jt = 0; jt < NDAYS; jt++) {
    // 由于需要区分天数和小时数，故每一天都至少需要做一次同步。
    // TODO: 补充边界插值
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

    for (kt = 1; kt < 2000; kt++) {
      // 可以考虑放到核上跑
      for (int pos = 0; pos < 24000; pos++) {
        time_step(kt, pos);
      }
    }

    cout << jt << endl;
  }
  return 0;
}