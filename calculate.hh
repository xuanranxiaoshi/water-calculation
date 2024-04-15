#pragma once
void calculate_FLUX(int t, int pos, Vec2 &FLUX);
void calculate_WHUV(int t, int pos, double &WH, double &WU, double &WV);
void calculate_HUV(int t, int pos);
void BOUNDA(int t, int j, int pos, Vec2 &FLUX, Vec &QL, Vec &QR, double &FIL
        , double HC,  double UC,  double VC,  double ZC);
int find_in_vec(const Vec& list, double target);
void OSHER(int t, int pos, Vec &QL, Vec &QR, double &FIL, Vec &FLR_OSHER);
void LAQP(double X, double &Y, Vec A, Vec B, double MS);
double QD(double ZL, double ZR, double ZB);
template <int T> void QS(int j, int pos, const Vec &QL, const Vec &QR, double &FIL, double &FIR, Vec& FLUX_OSHER);
void QF(double H, double U, double V, Vec &F);