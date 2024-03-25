#include <vector>
#include <cmath>

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
// scalar
double HM1, HM2;
double BI;

// no-state matrix
Vec WH;
Vec WU;
Vec WV;
Vec ZB1;
Vec ZBC;
Vec QL;

// one-dimention time-wise matrix
Vec2 H;
Vec2 U;
Vec2 V;
Vec2 Z;
Vec2 KLAS;
Vec2 NAC;
Vec2 W;

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
  H[TIME_NOW][pos] = std::max(H[TIME_PREV][pos] - WDTA * WH[i] + QLUA, HM1);
  Z[TIME_NOW][pos] = H[TIME_NOW][pos] + ZBC[TIME_NOW][pos];
  if (H[TIME_NOW][pos] <= HM1) {
    U[TIME_NOW][pos] = 0;
    V[TIME_NOW][pos] = 0;
  } else {
    if (H[TIME_NOW][i] <= HM2) {
      U[TIME_NOW][pos] = std::copysign(std::min(VMIN, std::abs(U[TIME_PREV][pos])), U[TIME_PREV][pos]);
      V[TIME_NOW][pos] = std::copysign(std::min(VMIN, std::abs(V[TIME_PREV][pos])), V[TIME_PREV][pos]);
    } else {
      U[TIME_NOW][pos] = (QX1 - DTAU - DT * WSF * U[TIME_PREV][pos]) / H[TIME_NOW][pos];
      V[TIME_NOW][pos] = (QY1 - DTAV - DT * WSF * V[TIME_PREV][pos]) / H[TIME_NOW][pos];      
    }
  }
  W[TIME_NOW][pos] = sqrt(U[TIME_NOW][pos] * U[TIME_NOW][pos] + V[TIME_NOW][pos] * V[TIME_NOW][pos]);
}

void BOUNDA(int t, int j, int pos) {
  double S0 = 0.0002;
  double DX2 = 5000.0;
  double BRDTH = 100.0;
  int KP;
  if(QL[1] > CL){
    FLUX_VAL( HI * QL[1],
              FLR(0) * QL[1], 
              4.905 * HI * HI, 
              FLR(0) * QL[1]);
  }

  switch (KP) {
    case 10:
    
    break;

    case 3:
    break;

    case 1:
    break;
    
    case 4:
    break;

    case 5:
    break;

    case 6:
    break;

    case 7:
    break;
  }

}

void OSHER(int t, int pos){

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
  for(int t = 0; t <= 2000; t ++){
    // 可以考虑放到核上跑
    for (int pos = 0; pos <= 24000; pos++){
      time_step(t, pos);
    }
  }
  return 0;  
}