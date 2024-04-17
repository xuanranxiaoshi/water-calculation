#include "read_data.hh"
#include "common.hh"
#include "calculate.hh"

using namespace DataManager;

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
    double CL = std::sqrt(9.81 * H1);
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
      HC = std::fmax(H[TIME_PREV][NC], HM1);
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
        FLUX_VAL(C1 * pow(H1, 1.5), 0, 0,
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
  Z[TIME_NOW][pos] = Z2;
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
  double HB;
//
  if (KLAS[j][pos] == 10) {
    int pos_near = find_in_vec(MBQ, pos);
    FLR(0) = -(QT[jt][pos_near] + DQT[pos_near] * t);
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
    int pos_near = find_in_vec(MBZQ, pos);
    for (int i = 0; i < NHQ; i++) {
      WZ[i] = ZW[i][pos_near];
      WQ[i] = QW[i][pos_near];
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
    int pos_near = find_in_vec(MBZ, pos);
    double HB1 = ZT[jt][pos_near] + DZT[pos_near] * t - ZBC[pos];
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
    int pos_near = find_in_vec(MBW, pos);
    double TOP = TOPW[pos_near];
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
    int pos_near = find_in_vec(MDI, pos);
    double TOP = TOPD[pos_near];
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
void QS(int j, int pos, const Vec &QL, const Vec &QR, double &FIL, double &FIR, Vec& FLUX_OSHER){
  Vec F(4);

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

void QF(double h, double u, double v, Vec &F) {
  // COMPUTATION OF FLUX COMPONENTS
  F[0] = h * u;
  F[1] = F[0] * u;
  F[2] = F[0] * v;
  F[3] = 4.905 * h * h;
}

int find_in_vec(const Vec& list, double target){
  auto it = std::find(list.begin(), list.end(), target);
  return std::distance(list.begin(), it);
}

void LAQP(double X, double &Y, Vec A, Vec B, double MS) {
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

    if (std::abs(X0 - X1) < 0.01 || std::abs(X1 - X2) < 0.01) {
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

    if (std::abs(X0 - X1) < 0.01 || std::abs(X1 - X2) < 0.01) {
      Y = B[i + 1];
    } else {
      double U_ = (X - X1) * (X - X2) / (X0 - X1) / (X0 - X2);
      double V_ = (X - X0) * (X - X2) / (X1 - X0) / (X1 - X2);
      double W_ = (X - X0) * (X - X1) / (X2 - X0) / (X2 - X1);
      Y = U_ * B[i] + V_ * B[i + 1] + W_ * B[i + 2];
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
  calculate_HUV(t, pos);
}

#include <omp.h>

void closeFile(){
  ZUV_file.close();
  H2U2V2_file.close();
  XY_TEC_file.close();
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

  double end_time;
  double start_time;
  int K0 = MDT / DT;
  int pos;

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

    start_time = omp_get_wtime();
    for (kt = 1; kt <= K0; kt++) {

      // for (pos = 0; pos < CEL; pos++) {
      //   time_step(kt, pos);
      // }
      // todo: 调用 kernel 函数

    }
    end_time = omp_get_wtime();
    cout << end_time - start_time << " seconds" << endl;
    
    output_data(kt);
  }
  return 0;
}


__global__ void translate_step(const int CEL, const int DT, const int jt, const int NHQ, const double HM1, const double HM2,
      double* H_pre, double* U_pre, double* V_pre, double* Z_pre, double* W_pre,
      int* NV, double* AREA, double* ZBC, double* ZB1, double* DQT, double* DZT, double* TOPW, double* TOPD, double* MBQ, double* MBZQ, double* MBW, double* MDI,
      double* d_QT, double* d_ZT,
      double** SIDE, double** SLCOS, double** SLSIN, double** KLAS, double** NAC, double** ZW, double** QW,
      double* H_res, double* U_res, double* V_res, double* Z_res, double* W_res){
        

      // calculate_HUV


      // calculate_WHUV


      // calculate_FLUX: 


      // BOUNDA


      // OSHER

      }