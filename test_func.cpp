#include <algorithm>
#include <cstdlib>
#include <iterator>
#include <vector>
#include <cmath>
#include <iostream>

#define TIME_NOW (t%2)
#define TIME_PREV ((t-1)%2)
#define FLR(x) FLUX[x][j][pos]
#define FLUX_VAL(A, B, C, D)  FLR(0) = A; \
                              FLR(1) = B; \
                              FLR(2) = C; \
                              FLR(3) = D 

using std::vector;
using Vec = vector<double>;
using Vec2 = vector<vector<double>>;
using Vec3 = vector<vector<vector<double>>>;

// TIME -> Y -> X
// const value
const double C1 = 1.7;

// Time scalar
int jt;
int kt;

// scalar
double HM1, HM2;
double BI;
double DT;
double NHQ;

// no-state matrix
Vec WH;
Vec WU;
Vec WV;
Vec ZB1;
Vec ZBC;
Vec QL;
Vec MBQ;
Vec NV;
Vec AREA;
Vec FNC;
Vec DQT;
Vec MBZQ;

// Origin scalar, but change to Matrix
Vec CL;
Vec FIL;

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
Vec2 QW;

// one-dimention j-wise matrix
Vec2 SIDE;
Vec2 SLCOS;
Vec2 SLSIN;

// other
Vec3 FLUX;


void calculate_FLUX(int t, int pos);
void calculate_WHUV(int t, int pos);
void calculate_HUV(int t, int pos);
void BOUNDA(int t, int j, int pos);
void OSHER(int t, int pos);
void CHOICE(Vec list, double target, int &index);
void LAQP(double X, double& Y, Vec A, Vec B, double MS);

void calculate_FLUX(int t, int pos) {
  for (int j = 0; j < 4; j++){
    // 局部变量初始化运算。
    double KP = KLAS[j][pos];
    double NC = NAC[j][pos];
    double ZI = fmax(Z[TIME_PREV][pos], ZB1[pos]);
    int HC, BC, ZC, UC, VC;
    if (NC == 0) {
      HC = 0;
      BC = 0;
      ZC = 0;
      UC = 0;
      VC = 0;
    } else {
      HC = std::fmax(H[TIME_PREV][HM1], HM1); 
      BC = ZBC[NC];
      ZC = std::fmax(ZBC[NC], Z[TIME_PREV][NC]);
      UC = U[TIME_PREV][NC];
      VC = V[TIME_PREV][NC];
    }

    if(KP >= 1 && KP <= 8 || KP >= 10) {
      BOUNDA(t, j, pos);
    } else if (H[TIME_PREV][pos] <= HM1 && HC <= HM1) {
      FLUX_VAL(0, 0, 0, 0);
    } else if (ZI <= BC) {
      // BUG: C1 undefined.
      FLUX_VAL(-C1 * pow(HC, 1.5),
               H[TIME_PREV][pos] * QL[1] * fabs(QL[1]),
               0,
               4.905 * H[TIME_PREV][pos] * H[TIME_PREV][pos]);

    } else if (ZC <= BI) {
      FLUX_VAL(C1 * pow(H[TIME_PREV][pos], 1.5),
               FLR(0) * QL[1],
               FLR(0) * QL[2],
               0);

    } else if (H[TIME_PREV][pos] <= HM2) {
      if (ZC > ZI) {
        double DH = fmax(ZC - BI, HM1);
        double UN = -C1 * std::sqrt(DH);
        FLUX_VAL(DH * UN,
                 FLR(0) * UN,
                 FLR(0) * (VC * COSA - UC * SINA),
                 4.905 *  H[TIME_PREV][pos] * H[TIME_PREV][pos]);
      } else {
        FLUX_VAL(-C1 * pow(HC, 1.5),
                 0, 
                 0, 
                 4.905 *  H[TIME_PREV][pos] * H[TIME_PREV][pos]);
      }
    } else if (HC <= HM2) {
      if (ZI > ZC) {
        double DH = fmax(ZC - BI, HM1);
        double UN = -C1 * std::sqrt(DH);
        double HC1 = ZC - BI;
        FLUX_VAL(DH * UN,
                 FLR(0) * UN,
                 FLR(0) * QL[2],
                 4.905 * HC1 * HC1);
      } else {
        FLUX_VAL(-HC * C1 * sqrt(HC),
                 H[TIME_PREV][pos] * QL[1] * QL[1],
                 0,
                 4.905 *  H[TIME_PREV][pos] * H[TIME_PREV][pos]);
      }
    } else {
      OSHER(t, pos);
    }
  }
}

void calculate_WHUV(int t, int pos) {
  WH[pos] = 0;
  WU[pos] = 0;
  WV[pos] = 0;

  for (int j = 0; j < 4; j++){
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
      SIDEX = std::min(0.5 * (SIDE[0][pos] + SIDE[2][pos]), 0.5 * (SIDE[1][pos] + SIDE[3][pos]));
  } else {
      SIDES = 0.5 * (SIDE[0][pos] + SIDE[1][pos] + SIDE[2][pos]);
      SIDEX = std::sqrt((SIDES - SIDE[0][pos]) * (SIDES - SIDE[1][pos]) * (SIDES - SIDE[2][pos]) / SIDES);
  }
  HSIDE = std::max(H[TIME_PREV][pos], HM1);
  DT2 = SIDEX / (U[TIME_PREV][pos] + std::sqrt(9.81 * HSIDE));
  DT2 = std::min(DT, DT2);
  DT2 = std::max(DT2, DT / 10.0);
  DTA = 1.0 * DT2 / (1.0 * AREA[pos]);
  WDTA = 1.00 * DTA;

  H[TIME_NOW][pos] = std::max(H[TIME_PREV][pos] - WDTA * WH[i] + QLUA, HM1);
  Z[TIME_NOW][pos] = H[TIME_NOW][pos] + ZBC[TIME_NOW][pos];
    if (H[TIME_NOW][pos] <= HM1) {
        U[TIME_NOW][pos] = 0.0;
        V[TIME_NOW][pos] = 0.0;
    } else {
        if (H[TIME_NOW][pos] <= HM2) {
            U[TIME_NOW][pos] = std::copysign(std::min(VMIN, std::abs(U[TIME_PREV][pos])), U[TIME_PREV][pos]);
            V[TIME_NOW][pos] = std::copysign(std::min(VMIN, std::abs(V[TIME_PREV][pos])), V[TIME_PREV][pos]);
        } else {
            QX1 = H[TIME_PREV][pos] * U[TIME_PREV][pos];
            QY1 = H[TIME_PREV][pos] * V[TIME_PREV][pos];
            DTAU = WDTA * WU[pos];
            DTAV = WDTA * WV[pos];
            
            FNCC = FNC[pos];
            WSF = FNCC * std::sqrt(U[TIME_PREV][pos] * U[TIME_PREV][pos] + V[TIME_PREV][pos] * V[TIME_PREV][pos]) / std::pow(H[TIME_PREV][pos], 0.33333);
            
            U[TIME_NOW][pos] = (QX1 - DTAU - DT * WSF * U[TIME_PREV][pos]) / H[TIME_NOW][pos];
            V[TIME_NOW][pos] = (QY1 - DTAV - DT * WSF * V[TIME_PREV][pos]) / H[TIME_NOW][pos];
            
            U[TIME_NOW][pos] = std::copysign(std::min(std::abs(U[TIME_NOW][pos]), 5.0), U[TIME_NOW][pos]);
            V[TIME_NOW][pos] = std::copysign(std::min(std::abs(V[TIME_NOW][pos]), 5.0), V[TIME_NOW][pos]);
        }
    }
  W[TIME_NOW][pos] = sqrt(U[TIME_NOW][pos] * U[TIME_NOW][pos] + V[TIME_NOW][pos] * V[TIME_NOW][pos]);
}

void BOUNDA(int t, int j, int pos) {
  Vec WZ;
  Vec WQ;
  Vec QB;


  double S0 = 0.0002;
  double DX2 = 5000.0;
  double BRDTH = 100.0;
  int KP;

  if(QL[1] > CL[pos]){
    FLUX_VAL( H[TIME_PREV][pos] * QL[1],
              FLR(0) * QL[1], 
              4.905 * H[TIME_PREV][pos] * H[TIME_PREV][pos], 
              FLR(0) * QL[1]);
  }

  int II;
  double HB;
  switch (KP) {
    case 10:
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
      break;

    case 3:
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
    break;

    case 1:
    choice(NZ, MBZ, I, II);
    double HB1 = ZT[JT][II] + DZT[II] * KT - BI;
    double FIAL = QL[2] + 6.264 * sqrt(HI);
    double UR0 = QL[2];
    double URB = UR0;
    for (int IURB = 1; IURB <= 30; IURB++) {
        double FIAR = URB - 6.264 * sqrt(HB1);
        URB = (FIAL + FIAR) * (FIAL - FIAR) * (FIAL - FIAR) / HB1 / 313.92;
        if (std::abs(URB - UR0) <= 0.0001)
            break;
        UR0 = URB;
    }
    FLR[0] = HB1 * URB;
    FLR[1] = FLR[0] * URB;
    FLR[3] = 4.905 * HB1 * HB1;
    return;
    break;
    
    case 4:
    FLR[0] = 0;
FLR[1] = 0;
FLR[2] = 0;
FLR[3] = 4.905 * HI * HI;
    break;

    case 5:
    QL[1] = std::max(QL[1], 0.0);
FLR[0] = HI * QL[1];
FLR[1] = FLR[0] * QL[1];
FLR[3] = 4.905 * HI * HI;
    break;

    case 6:
    if (KP == 6) {
    NE = I;
    if (NC != 0) {
        NE = std::min(I, NC);
    }
    CHOICE(NWE, MBW, NE, I1);
    TOP = TOPW(I1);
    if (ZI < TOP || ZC < TOP) {
        KP = 4;
        FLR[0] = 0;
        FLR[1] = 0;
        FLR[2] = 0;
        FLR[3] = 4.905 * HI * HI;
    }
    if (ZI > TOP && ZC < TOP) {
        FLR[0] = C0 * std::pow(ZI - TOP, 1.5);
        FLR[1] = FLR[0] * QL[1];
        FLR[2] = FLR[0] * QL[2];
        FLR[3] = 4.905 * std::pow(TOP - BI, 2);
        return;
    }
    if (ZI < TOP && ZC > TOP) {
        FLR[0] = -C0 * std::pow(ZC - TOP, 1.5);
        FLR[1] = FLR[0] * std::min(UC * COSA + VC * SINA, 0.0);
        FLR[2] = FLR[0] * (VC * COSA - UC * SINA);
        FLR[3] = 4.905 * std::pow(ZI - BI, 2);
        return;
    }
    DZ = std::abs(ZI - ZC);
    if (ZI <= ZC) {
        HD = ZI - TOP;
        UN = std::min(UC * COSA + VC * SINA, 0.0);
        VT = VC * COSA - UC * SINA;
    } else {
        HD = ZC - TOP;
        UN = std::max(QL[1], 0.0);
        VT = QL[2];
    }
    SH = HD + DZ;
    CE = std::min(1.0, 1.05 * std::pow(DZ / SH, 0.33333));
    if (ZI < ZC && UN > 0.0) {
        UN = 0.0;
    }
    FLR[0] = std::copysign(CE * C1 * std::pow(SH, 1.5), ZI - ZC);
    FLR[1] = FLR[0] * std::abs(UN);
    FLR[2] = FLR[0] * VT;
    FLR[3] = 4.905 * std::pow(TOP - BI, 2);
    return;
}
    break;

    case 7:
        CHOICE(NDI, MDI, I, I1);
    TOP = TOPD(I1);
    if (ZI > TOP || ZC > TOP) {
        KP = 0;
        KLAS[J][I] = 0;
        CQ = QD(ZI, ZC, TOP);
        CB = BRDTH / SL;
        FLR[0] = CQ * CB;
        FLR[1] = CB * std::copysign(CQ * CQ / HB, CQ);
        FLR[3] = 4.905 * HB * HB;
        return;
    } else {
        KP = 4;
        FLR[0] = 0;
        FLR[1] = 0;
        FLR[2] = 0;
        FLR[3] = 4.905 * HI * HI;
    }
    break;
  }

}

void OSHER(int t, int pos){

}

void CHOICE(Vec list, double target, int &index){
  auto it = std::find(list.begin(), list.end(), target);
  index = std::distance(list.begin(), it);
}

void LAQP(double X, double& Y, Vec A, Vec B, double MS) {
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

// Kernel函数，具体可以放到其他核上跑
// 0 <= pos <= CEL
void time_step(int t, int pos){
  calculate_FLUX(t, pos);
  calculate_WHUV(t, pos);
  calculate_HUV(t, pos);
}

int main() {
  // 必须串行执行的部分：
  // K0 = 2000
  for(jt = 0; jt < 365; jt ++){
    // 由于需要区分天数和小时数，故每一天都至少需要做一次同步。
    for(kt = 0; kt < 2000; kt ++){
      // 可以考虑放到核上跑
      for (int pos = 0; pos < 24000; pos++){
        time_step(kt, pos);
      }
    }
  }
  return 0;  
}