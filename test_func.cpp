#include <cmath>
#include<vector>

#define TIME_NOW (t%2)
#define TIME_PREV ((t-1)%2)
using std::vector;
using Vec = vector<double>;
using Vec2 = vector<vector<double>>;
using Vec3 = vector<vector<vector<double>>>;

// TIME -> Y -> X
Vec WH;
Vec WU;
Vec WV;
Vec2 H;
Vec2 U;
Vec2 V;
Vec2 Z;
Vec2 ZBC;
Vec2 W;
Vec3 FLUX;

void calculate_FLUX(int t, int pos);
void calculate_WHUV(int t, int pos);
void calculate_HUV(int t, int pos);
void BOUNDA(int t, int j, int pos);
void OSHER(int t, int pos);

#define FLR(x) FLUX[x][j][pos]
#define FLUX_VAL(A, B, C, D)  FLR(0) = A; \
                              FLR(1) = B; \
                              FLR(2) = C; \
                              FLR(3) = D 

void calculate_FLUX(int t, int pos) {
  for (int j = 0; j < 4; j++){
      if(KP >= 1 && KP <= k8 || KP >= 10) {
      BOUNDA(t, j, pos);
    } else if (HI <= HM1 && HC <= HM1) {
      FLUX_VAL(0, 0, 0, 0);
    } else if (ZI <= BC) {
      FLUX_VAL(0, 0, 0, 0);
    } else if (ZC <= BI) {
      FLUX_VAL(0, 0, 0, 0);
    } else if (HI <= HM2) {
      FLUX_VAL(0, 0, 0, 0);
    } else if (HC <= HM2) {
      FLUX_VAL(0, 0, 0, 0);
    } else {
      OSHER();
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
    WH[pos] += SL * FLUX[0][j][pos];
    WU[pos] += SLCA * FLR[1] - SLSA * FLR[2];
    WV[pos] += SLSA * FLR[1] - SLCA * FLR[2];
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