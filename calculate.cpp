#include "read_data.hh"
#include "common.hh"

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
int NHQ;
int NZ;
int NQ;
double STIME;

// no-state matrix
Vec ZB1;
Vec ZBC;
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
Vec YP;
Vec FNC0;
Vec QZSTIME1;
Vec QZSTEMP1;
Vec QZSTIME2;
Vec QZSTEMP2;


// one-dimention time-wise matrix
Vec2 H;
Vec2 U;
Vec2 V;
Vec2 Z;
Vec2 KLAS;
vector<vector<int>> NAC;
Vec2 W;
Vec2 QT;
Vec2 ZW;
Vec2 ZT;
Vec2 QW;
vector<vector<int>> NAP;

Vec2 COSF;
Vec2 SINF;

// one-dimention j-wise matrix
Vec2 SIDE;
Vec2 SLCOS;
Vec2 SLSIN;

// other

// DUMMY VALUE.
int dummy;

void calculate_FLUX(int t, int pos, Vec2 &FLUX);
void calculate_WHUV(int t, int pos, double &WH, double &WU, double &WV);
void calculate_HUV(int t, int pos);
void BOUNDA(int t, int j, int pos, Vec2 &FLUX, Vec &QL, Vec &QR, double &FIL
        , double HC,  double UC,  double VC,  double ZC);
void OSHER(int t, int pos, Vec &QL, Vec &QR, double &FIL, Vec &FLR_OSHER);
void CHOICE(Vec list, double target, int &index);
void LAQP(double X, double &Y, Vec A, Vec B, double MS);
double QD(double ZL, double ZR, double ZB);
template <int T> void QS(int j, int pos, const Vec &QL, const Vec &QR, double &FIL, double &FIR);
void QF(double H, double U, double V, Vec &F);
double BOUNDRYinterp(double THOURS, int NZQSTEMP, Vec ZQSTIME, Vec ZQSTEMP);

void calculate_FLUX(int t, int pos, Vec2 &FLUX) {
  Vec QL(3);
  Vec QR(3);
  QL[0] = H[TIME_PREV][pos];
  double U1 = U[TIME_PREV][pos];
  double V1 = V[TIME_PREV][pos];
  double H1 = H[TIME_PREV][pos];
  double Z1 = Z[TIME_PREV][pos];
  double v_ZB1 = ZB1[pos];
  Vec FLR_OSHER(4);

  for (int j = 0; j < 4; j++) {
    // 局部变量声明，以及初始化运算。
    int NC = NAC[j][pos] - 1;
    double KP = KLAS[j][pos];
    double COSJ = COSF[j][pos];
    double SINJ = SINF[j][pos];
    QL[1] = U1 * COSJ + V1 * SINJ;
    QL[2] = V1 * COSJ - U1 * SINJ;
    double CL = std::sqrt(9.81 * U1);
    double FIL = QL[1] + 2 * CL;
    double HC, BC, ZC, UC, VC;
    double ZI = std::fmax(Z1, v_ZB1);

    if (NC == 0) {
      HC = 0;
      BC = 0;
      ZC = 0;
      UC = 0;
      VC = 0;
    } else {
      HC = std::fmax(H1, HM1);
      BC = ZBC[NC];
      ZC = std::fmax(ZBC[NC], Z[TIME_PREV][NC]);
      UC = U[TIME_PREV][NC];
      VC = V[TIME_PREV][NC];
    }

    if ((KP >= 1 && KP <= 8) || KP >= 10) {
      BOUNDA(t, j, pos, FLUX, QL, QR, FIL, HC, UC, VC, ZC);
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
        double DH = std::fmax(ZC - ZBC[pos], HM1);
        double UN = -C1 * std::sqrt(DH);
        FLUX_VAL(DH * UN, FLR(0) * UN,
                 FLR(0) * (VC * COSJ - UC * SINJ),
                 4.905 * H1 * H1);
      } else {
        FLUX_VAL(C1 * pow(HC, 1.5), 0, 0,
                 4.905 * H1 * H1);
      }
    } else if (HC <= HM2) {
      if (ZI > ZC) {
        double DH = std::fmax(ZI - BC, HM1);
        double UN = C1 * std::sqrt(DH);
        double HC1 = ZC - ZBC[pos];
        FLUX_VAL(DH * UN, FLR(0) * UN, FLR(0) * QL[2], 4.905 * HC1 * HC1);
      } else {
        FLUX_VAL(-C1 * pow(HC, 1.5),
                 H1 * QL[1] * QL[1], 0,
                 4.905 * H1 * H1);
      }
    } else {
      QR[0] = std::fmax(ZC - ZBC[pos], HM1);
      double UR = UC * COSJ + VC * SINJ;
      QR[1] = UR * std::min(HC / QR[0], 1.5);
      if (HC <= HM2 || QR[0] <= HM2) {
        QR[1] = std::copysign(VMIN, UR);
      }
      QR[2] = VC * COSJ - UC * SINJ;
      OSHER(t, pos, QL, QR, FIL, FLR_OSHER);
      FLR(1) = FLR_OSHER[1] + (1 - std::min(HC / QR[0], 1.5)) * HC * UR * UR / 2;
      FLUX_VAL(FLR_OSHER[0], FLR(1), FLR_OSHER[2], FLR_OSHER[3]);
    }
  }
}

void calculate_WHUV(int t, int pos, double &WH, double &WU, double &WV) {
  Vec2 FLUX(4, Vec(4));
  calculate_FLUX(t, pos, FLUX);

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

void calculate_HUV(int t, int pos) {
  double SIDEX, SIDES, HSIDE, DT2, DTA, WDTA, QX1, QY1, DTAU, DTAV, WSF;
  double H1 = H[TIME_PREV][pos];
  double U1 = U[TIME_PREV][pos];
  double V1 = V[TIME_PREV][pos];
  double WH = 0.0, WU = 0.0, WV = 0.0;
  calculate_WHUV(t, pos, WH, WU, WV);

  if (NV[pos] == 4) {
    SIDEX = std::min(0.5 * (SIDE[0][pos] + SIDE[2][pos]),
                     0.5 * (SIDE[1][pos] + SIDE[3][pos]));
  } else {
    SIDES = 0.5 * (SIDE[0][pos] + SIDE[1][pos] + SIDE[2][pos]);
    SIDEX = std::sqrt((SIDES - SIDE[0][pos]) * (SIDES - SIDE[1][pos]) *
                      (SIDES - SIDE[2][pos]) / SIDES);
  }
  HSIDE = std::max(H1, HM1);
  DT2 = SIDEX / (U1 + std::sqrt(9.81 * HSIDE));
  DT2 = std::fmin(DT, DT2);
  DT2 = std::max(DT2, DT / 10.0);
  DTA = 1.0 * DT2 / (1.0 * AREA[pos]);
  WDTA = 1.00 * DTA;

  double H2, U2, V2, Z2, W2;
  H2 = std::max(H1 - WDTA * WH + QLUA, HM1);
  Z2 = H[TIME_NOW][pos] + ZBC[pos];
  if (H2 <= HM1) {
    U2 = 0.0;
    V2 = 0.0;
  } else {
    if (H2 <= HM2) {
      U2 = std::copysign( std::min(VMIN, std::abs(U1)), U1);
      V2 = std::copysign( std::min(VMIN, std::abs(V1)), V1);
    } else {
      QX1 = H1 * U1;
      QY1 = H1 * V1;
      DTAU = WDTA * WU;
      DTAV = WDTA * WV;
      WSF = FNC[pos] * std::sqrt(U1 * U1 + V1 * V1) / std::pow(H1, 0.33333);
      U2 = (QX1 - DTAU - DT * WSF * U1) / H2;
      V2 = (QY1 - DTAV - DT * WSF * V1) / H2;
      U2 = std::copysign( std::min(std::abs(U2), 5.0), U2);
      V2 = std::copysign( std::min(std::abs(V2), 5.0), V2);
    }
  }
  W2 = std::sqrt(U2 * U2 + V2 * V2);
  H[TIME_NOW][pos] = H2;
  U[TIME_NOW][pos] = U2;
  V[TIME_NOW][pos] = V2;
  W[TIME_NOW][pos] = W2;
}

void BOUNDA(int t, int j, int pos, Vec2 &FLUX, Vec &QL, Vec &QR, double &FIL,
            double HC, double UC, double VC, double ZC) {
  Vec WZ(NHQ);
  Vec WQ(NHQ);

  double S0 = 0.0002;
  double DX2 = 5000.0;
  double BRDTH = 100.0;
  double CL = std::sqrt(9.81 * H[TIME_PREV][pos]);

  if (QL[1] > CL) {
    FLUX_VAL(H[TIME_PREV][pos] * QL[1], FLR(0) * QL[1],
             4.905 * H[TIME_PREV][pos] * H[TIME_PREV][pos], FLR(0) * QL[1]);
  }
  FLR(2) = 0;
  if (QL[1] > 0) FLR(2) = H[TIME_PREV][pos] * QL[1] * QL[2];
  int II;
  double HB;
//
  if (KLAS[j][pos] == 10) {
    CHOICE(MBQ, pos, II);
    FLR(0) = -(QT[jt][II] + DQT[II] * t);
    FLR(0) = FLR(0) / SIDE[j][pos];
    double QB2 = FLR(0) * FLR(0);
    double HB0 = H[TIME_PREV][pos];

    for (int K = 1; K <= 20; K++) {
      double W_temp = FIL - FLR(0) / HB0;
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
  } else if (KLAS[j][pos] == 3) {
    double CQ;
    double HR;
    double W_temp;
    double HR0 = H[TIME_PREV][pos];
    CHOICE(MBZQ, pos, II);
    for (int i = 0; i < NHQ; i++) {
      WZ[i] = ZW[i][II];
      WQ[i] = QW[i][II];
    }
    for (int I1 = 0; I1 < 20; I1++) {
      double ZR0 = HR0 + ZBC[pos];
      LAQP(ZR0, CQ, WZ, WQ, NHQ);
      W_temp = FIL - CQ / HR0;
      HR = W_temp * W_temp / 39.24;
      if (std::abs(HR - HR0) <= 0.001)
        break;
      HR0 = HR;
    }
    FLR(0) = CQ;
    FLR(1) = CQ * CQ / HR;
    HB = (H[TIME_PREV][pos] + HR) / 2;
    FLR(3) = 4.905 * HB * HB;
  } else if (KLAS[j][pos] == 1) {
    CHOICE(MBZ, pos, II);
    double HB1 = ZT[jt][II] + DZT[II] * t - ZBC[pos];
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
  } else if (KLAS[j][pos] == 4) {
    FLR(0) = 0;
    FLR(1) = 0;
    FLR(2) = 0;
    FLR(3) = 4.905 * H[TIME_PREV][pos] * H[TIME_PREV][pos];
  } else if (KLAS[j][pos] == 5) {
    QL[1] = std::max(QL[1], 0.0);
    FLR(0) = H[TIME_PREV][pos] * QL[1];
    FLR(1) = FLR(0) * QL[1];
    FLR(3) = 4.905 * H[TIME_PREV][pos] * H[TIME_PREV][pos];
  } else if (KLAS[j][pos] == 6) {
    double NE = pos;
    if (NAC[j][pos] != 0) {
      NE = std::fmin(pos, NAC[j][pos]);
    }
    CHOICE(MBW, NE, II);
    double TOP = TOPW[II];
    if (Z[TIME_PREV][pos] < TOP || ZC < TOP) {
      FLR(0) = 0;
      FLR(1) = 0;
      FLR(2) = 0;
      FLR(3) = 4.905 * H[TIME_PREV][pos] * H[TIME_PREV][pos];
    }
    if (Z[TIME_PREV][pos] > TOP && ZC < TOP) {
      FLR(0) = C0 * std::pow(Z[TIME_PREV][pos] - TOP, 1.5);
      FLR(1) = FLR(0) * QL[1];
      FLR(2) = FLR(0) * QL[2];
      FLR(3) = 4.905 * std::pow(TOP - ZBC[pos], 2);
      return;
    }
    if (Z[TIME_PREV][pos] < TOP && ZC > TOP) {
      FLR(0) = -C0 * std::pow(ZC - TOP, 1.5);
      FLR(1) = FLR(0) *
               std::min(UC * COSF[j][pos] + VC * SINF[j][pos], 0.0);
      FLR(2) = FLR(0) * (VC * COSF[j][pos] - UC * SINF[j][pos]);
      FLR(3) = 4.905 * std::pow(Z[TIME_PREV][pos] - ZBC[pos], 2);
      return;
    }
    double DZ = std::abs(Z[TIME_PREV][pos] - ZC);
    double HD;
    double UN;
    double VT;
    if (Z[TIME_PREV][pos] <= ZC) {
      HD = Z[TIME_PREV][pos] - TOP;
      UN = std::min(UC * COSF[j][pos] + VC * SINF[j][pos], 0.0);
      VT = VC * COSF[j][pos] - UC * SINF[j][pos];
    } else {
      HD = ZC - TOP;
      UN = std::max(QL[1], 0.0);
      VT = QL[2];
    }
    double SH = HD + DZ;
    double CE = std::min(1.0, 1.05 * std::pow(DZ / SH, 0.33333));
    if (Z[TIME_PREV][pos] < ZC && UN > 0.0) {
      UN = 0.0;
    }
    FLR(0) =
        std::copysign(CE * C1 * std::pow(SH, 1.5), Z[TIME_PREV][pos] - ZC);
    FLR(1) = FLR(0) * std::abs(UN);
    FLR(2) = FLR(0) * VT;
    FLR(3) = 4.905 * std::pow(TOP - ZBC[pos], 2);
  } else if (KLAS[j][pos] == 7) {
    CHOICE(MDI, pos, II);
    double TOP = TOPD[II];
    if (Z[TIME_PREV][pos] > TOP || ZC > TOP) {
      KLAS[j][pos] = 0;
      double CQ = QD(Z[TIME_PREV][pos], ZC, TOP);
      double CB = BRDTH / SIDE[j][pos];
      FLR(0) = CQ * CB;
      FLR(1) = CB * std::copysign(CQ * CQ / HB, CQ);
      FLR(3) = 4.905 * HB * HB;
      return;
    } else {
      FLR(0) = 0;
      FLR(1) = 0;
      FLR(2) = 0;
      FLR(3) = 4.905 * H[TIME_PREV][pos] * H[TIME_PREV][pos];
    }
  }
}

void OSHER(int t, int pos, Vec &QL, Vec &QR, double &FIL, Vec &FLR_OSHER) {
  double CR = sqrt(9.81 * QR[0]);
  double FIR = QR[1] - 2 * CR;
  double UA = (FIL + FIR) / 2;
  double CA = fabs((FIL - FIR) / 4);
  double CL = std::sqrt(9.81 * H[TIME_PREV][pos]);

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
          QS<2>(1, pos, QL, QR, FIL, FIR);
          break;
        case 2:
          QS<3>(1, pos, QL, QR, FIL, FIR);
          break;
        case 3:
          QS<5>(1, pos, QL, QR, FIL, FIR);
          break;
        case 4:
          QS<6>(1, pos, QL, QR, FIL, FIR);
          break;
        default:
          break;
      }
      break;
    case 2:
      switch(K2){
        case 1:
          QS<1>(1, pos, QL, QR, FIL, FIR);
          break;
        case 2:
          QS<1>(1, pos, QL, QR, FIL, FIR);
          QS<2>(-1, pos, QL, QR, FIL, FIR);
          QS<3>(1, pos, QL, QR, FIL, FIR);
          break;
        case 3:
          QS<1>(1, pos, QL, QR, FIL, FIR);
          QS<2>(-1, pos, QL, QR, FIL, FIR);
          QS<5>(1, pos, QL, QR, FIL, FIR);
          break;
        case 4:
          QS<1>(1, pos, QL, QR, FIL, FIR);
          QS<2>(-1, pos, QL, QR, FIL, FIR);
          QS<6>(1, pos, QL, QR, FIL, FIR);
          break;
        default:
          break;
      }
      break;
    case 3:
      switch(K2){
        case 1:
          QS<2>(1, pos, QL, QR, FIL, FIR);
          QS<6>(-1, pos, QL, QR, FIL, FIR);
          QS<7>(1, pos, QL, QR, FIL, FIR);
          break;
        case 2:
          QS<3>(1, pos, QL, QR, FIL, FIR);
          QS<6>(-1, pos, QL, QR, FIL, FIR);
          QS<7>(1, pos, QL, QR, FIL, FIR);
          break;
        case 3:
          QS<5>(1, pos, QL, QR, FIL, FIR);
          QS<6>(-1, pos, QL, QR, FIL, FIR);
          QS<7>(1, pos, QL, QR, FIL, FIR);
          break;
        case 4:
          QS<7>(1, pos, QL, QR, FIL, FIR);
          break;
        default:
          break;
      }
      break;
    case 4:
      switch(K2){
        case 1:
          QS<1>(1, pos, QL, QR, FIL, FIR);
          QS<6>(-1, pos, QL, QR, FIL, FIR);
          QS<7>(1, pos, QL, QR, FIL, FIR);
          break;
        case 2:
          QS<1>(1, pos, QL, QR, FIL, FIR);
          QS<2>(-1, pos, QL, QR, FIL, FIR);
          QS<3>(1, pos, QL, QR, FIL, FIR);
          QS<6>(-1, pos, QL, QR, FIL, FIR);
          QS<7>(1, pos, QL, QR, FIL, FIR);
          break;
        case 3:
          QS<1>(1, pos, QL, QR, FIL, FIR);
          QS<2>(-1, pos, QL, QR, FIL, FIR);
          QS<5>(1, pos, QL, QR, FIL, FIR);
          QS<6>(-1, pos, QL, QR, FIL, FIR);
          QS<7>(1, pos, QL, QR, FIL, FIR);
          break;
        case 4:
          QS<1>(1, pos, QL, QR, FIL, FIR);
          QS<2>(-1, pos, QL, QR, FIL, FIR);
          QS<7>(1, pos, QL, QR, FIL, FIR);
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
void QS(int j, int pos, const Vec &QL, const Vec &QR, double &FIL, double &FIR){
  Vec result(4);
  Vec F(4);

  if constexpr (T == 1) {
    QF(QL[0], QL[1], QL[2], F);
    for (int i = 0; i < 4; i++) {
      result[i] += F[i] * j;
    }
  } else if constexpr (T == 2) {
    double US = FIL / 3;
    double HS = US * US / 9.81;
    QF(HS, US, QL[2], F);
    for (int i = 0; i < 4; i++) {
      result[i] += F[i] * j;
    }
  } else if constexpr (T == 3) {
    FIL = FIL - (FIL + FIR) / 2;
    double HA = FIL * FIL / 39.24;
    QF(HA, (FIL + FIR) / 2, QL[2], F);
    for (int i = 0; i < 4; i++) {
      result[i] += F[i] * j;
    }
  } else if constexpr (T == 4) {

  } else if constexpr (T == 5) {
    FIR = FIR - (FIL + FIR) / 2;
    double HA = FIR * FIR / 39.24;
    QF(HA, (FIL + FIR) / 2, QR[2], F);
    for (int i = 0; i < 4; i++) {
      result[i] += F[i] * j;
    }
  } else if constexpr (T == 6) {
    double US = FIR / 3;
    double HS = US * US / 9.81;
    QF(HS, US, QR[2], F);
    for (int i = 0; i < 4; i++) {
      result[i] += F[i] * j;
    }
  } else if constexpr (T == 7) {
    QF(QR[0], QR[1], QR[2], F);
    for (int i = 0; i < 4; i++) {
      result[i] += F[i] * j;
    }
  } else {
    std::exit(1);
  }
}

void QF(double h, double u, double v, Vec &F) {
  // COMPUTATION OF FLUX COMPONENTS
  F[0] = h * u;
  F[1] = F[0] * u;
  F[2] = F[0] * v;
  F[3] = 4.905 * h * h;
}

void CHOICE(Vec list, double target, int &index) {
  auto it = std::find(list.begin(), list.end(), target);
  index = std::distance(list.begin(), it);
}

void LAQP(double X, double &Y, Vec A, Vec B, double MS) {
  int NHQ = 5;
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

template <typename T> T readFromLine(const std::string &line) {
  std::istringstream iss(line);
  T result;
  iss >> result;
  return result;
}

void loadFromFilePNAC(const std::string& p) {
  std::ifstream file(data_root_path + p);
  ASSERT_READ(file);
  std::string line;
  std::getline(file, line);
  int NO;
  NAC.resize(4, std::vector<int>(CEL));
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

void loadFromFilePNAP(const std::string& p) {
  std::ifstream file(data_root_path + p);
  ASSERT_READ(file);
  std::string line;
  std::getline(file, line);
  int NO;
  NAP.resize(4, std::vector<int>(CEL));
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

void loadFromFilePKLAS(const std::string& p) {
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

void loadFromFilePZBC(const std::string& p) {
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

void loadFromFileMBZ(const std::string& p) {
  std::ifstream file(data_root_path + p);
  ASSERT_READ(file);
  std::string line;
  std::getline(file, line);
  std::istringstream iss(line);
  iss >> NNZ0;
  int NO;
  MBZ.resize(NZ, 0);
  NNZ.resize(NZ, 0);
  for (int i = 0; i < NZ; i++) {
    std::getline(file, line);
    std::istringstream iss(line);
    iss >> NO;
    iss >> MBZ[i];
    iss >> NNZ[i];
  }
  file.close();
}

void loadFromFileMBQ(const std::string& p) {
  std::ifstream file(data_root_path + p);
  ASSERT_READ(file);
  std::string line;
  std::getline(file, line);
  std::istringstream iss(line);
  iss >> NNQ0;
  int NO;
  MBQ.resize(NQ, 0);
  NNQ.resize(NQ, 0);
  for (int i = 0; i < NQ; i++) {
    std::getline(file, line);
    std::istringstream iss(line);
    iss >> NO;
    iss >> MBQ[i];
    iss >> NNQ[i];
  }
  file.close();
}

void loadFromFilePXY(const std::string& p) {
  std::ifstream file(data_root_path + p);
  ASSERT_READ(file);
  std::string line;
  std::getline(file, line);
  int NO;
  XP.resize(NOD, 0);
  YP.resize(NOD, 0);
  for (int i = 0; i < NOD; i++) {
    std::getline(file, line);
    std::istringstream iss(line);
    iss >> NO;
    iss >> XP[i];
    iss >> YP[i];
  }
  file.close();
}

void loadFromFileInitLevel(const std::string& p) {
  std::ifstream file(data_root_path + p);
  ASSERT_READ(file);
  std::string line;
  std::getline(file, line);
  int NO;
  Z.resize(2, vector<double>(CEL));
  for (int i = 0; i < CEL; i++) {
    std::getline(file, line);
    std::istringstream iss(line);
    iss >> Z[1][i];
  }
  file.close();
}

void loadFromFileU1(const std::string& p) {
  std::ifstream file(data_root_path + p);
  ASSERT_READ(file);
  std::string line;
  std::getline(file, line);
  int NO;
  U.resize(2, vector<double>(CEL));
  for (int i = 0; i < CEL; i++) {
    std::getline(file, line);
    std::istringstream iss(line);
    iss >> U[1][i];
  }
  file.close();
}

void loadFromFileV1(const std::string& p) {
  std::ifstream file(data_root_path + p);
  ASSERT_READ(file);
  std::string line;
  std::getline(file, line);
  int NO;
  V.resize(2, vector<double>(CEL));
  for (int i = 0; i < CEL; i++) {
    std::getline(file, line);
    std::istringstream iss(line);
    iss >> V[1][i];
  }
  file.close();
}

void loadFromFileCV(const std::string& p) {
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
    double XP1 = XP[NW[0] - 1];
    double XP2 = XP[NW[1] - 1];
    double XP3 = XP[NW[2] - 1];
    double YP1 = YP[NW[0] - 1];
    double YP2 = YP[NW[1] - 1];
    double YP3 = YP[NW[2] - 1];
    AREA[i] = (XP1 * YP2 - YP1 * XP2 - XP1 * YP3 + YP1 * XP3 + XP2 * YP3 -
               YP2 * XP3) /
              2.0;
    if (NV[i] == 4) {
      double XP4 = XP[NW[3] - 1];
      double YP4 = YP[NW[3] - 1];
      AREA[i] += (XP1 * YP3 - YP1 * XP3 - XP1 * YP4 + YP1 * XP4 + XP3 * YP4 -
                 YP3 * XP4) /
                2.0;
    }
    for (int j = 0; j < NV[i]; j++) {
      int N1 = NW[j] - 1;
      int N2 = NW[(j + 1) % NV[i]] - 1;
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

  H.resize(2, Vec(CEL));
  SLCOS.resize(4, vector<double>(CEL));
  SLSIN.resize(4, vector<double>(CEL));
  FNC.resize(CEL, 0);
  for (int i = 0; i < CEL; i++) {
    if (NAP[0][i] == 0) {
      continue;
    }
    if (Z[1][i] <= ZBC[i]) {
      H[1][i] = HM1;
      Z[1][i] = ZB1[i];
    } else {
      H[1][i] = Z[1][i] - ZBC[i];
    }
  }
  for (int i = 0; i < CEL; i++) {
    FNC[i] = 9.81 * FNC0[i] * FNC0[i];
    Z[0][i] = Z[1][i];
    H[0][i] = H[1][i];
    U[0][i] = U[1][i];
    V[0][i] = V[1][i];
    for (int j = 0; j < NV[i]; j++) {
      SLCOS[j][i] = SIDE[j][i] * COSF[j][i];
      SLSIN[j][i] = SIDE[j][i] * SINF[j][i];
    }
  }
}

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
      STIME1 = STIME + i / (24.0 * 3600.0 / MDT);
      for (int j = 0; j < NZ; j++) {
        if (NNZ[j] == k)
          ZT[j][i] = BOUNDRYinterp(STIME1, NZTEMP, QZSTIME1, QZSTEMP1);
      }
    }
  }

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
    QZSTIME1.resize(NQTEMP, 0);
    QZSTEMP1.resize(NQTEMP, 0);
    for (int i = 0; i < NQTEMP; i++) {
      std::getline(NQ_file, current_line);
      std::istringstream iss(current_line);
      iss >> QZSTIME2[i];
      iss >> QZSTEMP2[i];
    }
    for (int i = 0; i < NDAYS; i++) {
      STIME1 = STIME + i / (24.0 * 3600.0 / MDT);
      for (int j = 0; j < NQ; j++) {
        if (NNQ[j] == k)
          QT[j][i] = BOUNDRYinterp(STIME1, NQTEMP, QZSTIME2, QZSTEMP2);
      }
    }
  }
}

double BOUNDRYinterp(double THOURS, int NZQSTEMP, Vec ZQSTIME, Vec ZQSTEMP) {
  double result = 0;
  for (int i = 0; i < NZQSTEMP; i++) {
    if (THOURS >= ZQSTIME[i] && THOURS <= ZQSTIME[i + 1]) {
      result = ZQSTEMP[i] + (ZQSTEMP[i + 1] - ZQSTEMP[i]) /
                                  (ZQSTIME[i + 1] - ZQSTIME[i]) *
                                  (THOURS - ZQSTIME[i]);
    }
  }
  return result;
}

// Kernel函数，具体可以放到其他核上跑
// 0 <= pos <= CEL
void time_step(int t, int pos) {
  calculate_HUV(t, pos);
}

#include <omp.h>

int main() {
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

  double end_time;
  double start_time;
  omp_set_num_threads(16);
  // 必须串行执行的部分：
  // K0 = 2000
  double K0 = (double) MDT / DT;
  int pos;

  for (jt = 0; jt < 100; jt++) {
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

    start_time = omp_get_wtime();
    for (kt = 1; kt <= K0; kt++) {
      // 可以考虑放到核上跑
      #pragma omp parallel for
      for (pos = 0; pos < CEL; pos++) {
        time_step(kt, pos);
      }
    }
    end_time = omp_get_wtime();
    cout << end_time - start_time << " seconds" << endl;
  }
  return 0;
}