#pragma once
#include "common.hh"
using std::cout;
using std::endl;

namespace DataManager{
    // TIME -> Y -> X
// const value
  const double C0 = 1.33;
  const double C1 = 1.7;
  const double VMIN = 0.001;
  const double QLUA = 0.0;
  const string data_root_path = "./origin/changjiangkou/";

// 时间层面的各个变量。
  int jt;
  int kt;

// 从文件中读取的标量
  int CEL;
  int NDAYS;
  int NDI;
  int NOD;
  int MDT;
  int DT;
  int NNZ0;
  int NNQ0;
  double HM1, HM2;
  int NHQ;
  int NZ;
  int NQ;
  double STIME;

  // 从文件夹直接读入的变量
  Vec ZBC; // SOURCES/PZBC.DAT
  Vec2i NAC; // SOURCES/PNAC.DAT
  Vec2i NAP; // SOURCES/PNAP.DAT
  Vec2 KLAS; // SOURCES/PKLAS.DAT
  Vec XP;    // SOURCES/PXY.DAT
  Vec YP;    // SOURCES/PXY.DAT
  Vec MBQ;   // SOURCES/MBQ.DAT
  Vec NNQ;   // SOURCES/MBQ.DAT
  Vec MBZ;   // SOURCES/MBZ.DAT
  Vec NNZ;   // SOURCES/MBZ.DAT
  Vec FNC0;  // CV.DAT

  // 经过预处理函数pre与take_boundary_for_two_d这两个函数处理后的只读数据
  Vec2 ZT;
  Vec2 QT;
  Vec2 SIDE;
  Vec2 COSF;
  Vec2 SINF;
  Vec2 SLCOS;
  Vec2 SLSIN;
  Vec ZB1;
  Veci NV; // 初始化的时候貌似是全部填上4
  Vec FNC;
  Vec AREA;
  Vec DZT;
  Vec DQT;

  // 暂时不起作用的变量
  Vec MBW; // ??
  Vec MDI; // ??
  Vec TOPW; // ??
  Vec TOPD; // ??
  Vec MBZQ; // ??


// one-dimention time-wise matrix
// 时间层面需要存若干份的数据。其中，HUV是需要计算的三个变量。
  Vec2 H;
  Vec2 U;
  Vec2 V;
  Vec2 Z;
  Vec2 W;
  Vec2 ZW;
  Vec2 QW;

  // 管理输出文件:
  std::ofstream ZUV_file;
  std::ofstream SHA_file;
  std::ofstream H2U2V2_file;
  std::ofstream XY_TEC_file;
  std::ofstream SHA_TEC_file;


    // 用来填充。
  int dummy;

  static std::fstream file_to_read;
  static std::string current_line;
  static std::string current_path;

  #define TIME_NOW (kt % 2)
  #define TIME_PREV ((kt - 1) % 2)
  #define FLR(x) FLUX[x][j]
  #define FLUX_VAL(A, B, C, D)                                                   \
    FLR(0) = A;                                                                  \
    FLR(1) = B;                                                                  \
    FLR(2) = C;                                                                  \
    FLR(3) = D

  #define READ(X)                                                                \
    std::getline(file_to_read, current_line);                                    \
    X = readFromLine<decltype(X)>(current_line);
  #define ASSERT_READ(X)                                                         \
    if (!(X.is_open())) {                                                        \
      std::cerr << "ERROR FILE." << std::endl;                                   \
      exit(1);                                                                   \
    }

  // 用于获取可变参数个数
  #define ArgNums(...) _ArgNums(, ##__VA_ARGS__, 6, 5, 4, 3, 2, 1, 0)
  #define _ArgNums(_0, _1, _2, _3, _4, _5, _6, COUNT, ...) COUNT

  // #和##运算符直接使用会阻止宏展开，所以封装一层
  #define concat(a, b) _concat(a, b)
  #define _concat(a, b) a##b

  #define READ_DATA(FILE, ...)                                                   \
    do {                                                                         \
      current_path = string("./origin/changjiangkou/") + FILE + ".DAT";          \
      file_to_read.open(current_path);                                           \
      ASSERT_READ(file_to_read);                                                 \
      concat(READ_DATA_, ArgNums(__VA_ARGS__))(__VA_ARGS__);                     \
      file_to_read.close();                                                      \
    } while (0);

  #define READ_DATA_1(arg) READ(arg)
  #define READ_DATA_2(_1, _2) READ_DATA_1(_1) READ_DATA_1(_2)
  #define READ_DATA_3(_1, ...) READ_DATA_1(_1) READ_DATA_2(__VA_ARGS__)
  #define READ_DATA_4(_1, ...) READ_DATA_1(_1) READ_DATA_3(__VA_ARGS__)
  #define READ_DATA_5(_1, ...) READ_DATA_1(_1) READ_DATA_4(__VA_ARGS__)
  #define READ_DATA_6(_1, ...) READ_DATA_1(_1) READ_DATA_5(__VA_ARGS__)

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
    Vec QZSTIME1(NZTEMP, 0);
    Vec QZSTEMP1(NZTEMP, 0);
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
    Vec QZSTIME2(NQTEMP, 0);
    Vec QZSTEMP2(NQTEMP, 0);
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

  // 文件声明：
  // 28: /OUTPUT/H2U2V2.OUT
  // 20: /OUTPUT/ZUV.OUT
  // 24: /OUTPUT/SHA.OUT
  // 901: /OUTPUT/XY-TEC.DAT
  // 902: /OUTPUT/SHA-TEC.DAT
  void output_data(int t) {
    if (!ZUV_file.is_open() || !H2U2V2_file.is_open() || !XY_TEC_file.is_open()) {
      std::cerr << "BAD FILE OUTPUT" << endl;
      exit(1);
    }
    H2U2V2_file << std::setw(1) << " " << std::setw(3) << "JT=" << std::setw(5) << jt
                << std::setw(2) << " " << "KT=" << std::setw(5) << kt
                << std::setw(2) << " " << "DT=" << std::setw(3) << DT << std::setw(1) << " "
                << "SEC" << std::setw(2) << " " << "T=" << std::setw(2) << 0 << std::setw(1) << "H"
                << std::setw(2) << " " << "NSF=" << std::setw(2) << 0 << "/" << 0 << std::setw(2)
                << " " << "WEC=" << std::fixed << std::setprecision(2) << 0 << "/" << std::setw(2)
                << " " << "CQL=" << std::fixed << std::setprecision(2) << 0
                << std::setw(2) << " " << "INE=" << std::setw(1) << 0 << std::endl;
    ZUV_file << std::setw(1) << " " << std::setw(3) << "JT=" << std::setw(5) << jt
             << std::setw(2) << " " << "KT=" << std::setw(5) << kt
             << std::setw(2) << " " << "DT=" << std::setw(3) << DT << std::setw(1) << " "
             << "SEC" << std::setw(2) << " " << "T=" << std::setw(2) << 0 << std::setw(1) << "H"
             << std::setw(2) << " " << "NSF=" << std::setw(2) << 0 << "/" << 0 << std::setw(2)
             << " " << "WEC=" << std::fixed << std::setprecision(2) << 0 << "/" << std::setw(2)
             << " " << "CQL=" << std::fixed << std::setprecision(2) << 0
             << std::setw(2) << " " << "INE=" << std::setw(1) << 0 << std::endl;
    H2U2V2_file << std::setw(5) << " "; // 插入一个宽度为5的空格
    H2U2V2_file << std::endl; // 换行
    H2U2V2_file << std::setw(5) << " "; // 插入一个宽度为5的空格
    H2U2V2_file << std::setw(3) << "H2="; // 输出"H2="字符串
    for (int i = 0; i < CEL; i++) {
      if (i % 10 == 0) {
        H2U2V2_file << std::endl; // 每输出10个数换行
        H2U2V2_file << std::setw(5) << " "; // 每行开头插入5个空格
      }
      H2U2V2_file << std::setw(10) << std::fixed << std::setprecision(2) << H[TIME_NOW][i]; // 输出浮点数，总宽度为10，保留2位小数
    }

    H2U2V2_file << std::setw(5) << " "; // 插入一个宽度为5的空格
    H2U2V2_file << std::endl; // 换行
    H2U2V2_file << std::setw(5) << " "; // 插入一个宽度为5的空格
    H2U2V2_file << std::setw(3) << "U2="; // 输出"H2="字符串
    for (int i = 0; i < CEL; i++) {
      if (i % 10 == 0) {
        H2U2V2_file << std::endl; // 每输出10个数换行
        H2U2V2_file << std::setw(5) << " "; // 每行开头插入5个空格
      }
      H2U2V2_file << std::setw(10) << std::fixed << std::setprecision(2) << U[TIME_NOW][i]; // 输出浮点数，总宽度为10，保留2位小数
    }

    H2U2V2_file << std::setw(5) << " "; // 插入一个宽度为5的空格
    H2U2V2_file << std::endl; // 换行
    H2U2V2_file << std::setw(5) << " "; // 插入一个宽度为5的空格
    H2U2V2_file << std::setw(3) << "V2="; // 输出"H2="字符串
    for (int i = 0; i < CEL; i++) {
      if (i % 10 == 0) {
        H2U2V2_file << std::endl; // 每输出10个数换行
        H2U2V2_file << std::setw(5) << " "; // 每行开头插入5个空格
      }
      H2U2V2_file << std::setw(10) << std::fixed << std::setprecision(2) << V[TIME_NOW][i]; // 输出浮点数，总宽度为10，保留2位小数
    }

    ZUV_file << std::setw(5) << " "; // 插入一个宽度为5的空格
    ZUV_file << std::endl; // 换行
    ZUV_file << std::setw(5) << " "; // 插入一个宽度为5的空格
    ZUV_file << std::setw(3) << "Z2="; // 输出"H2="字符串
    for (int i = 0; i < CEL; i++) {
      if (i % 10 == 0) {
        ZUV_file << std::endl; // 每输出10个数换行
        ZUV_file << std::setw(5) << " "; // 每行开头插入5个空格
      }
      ZUV_file << std::setw(10) << std::fixed << std::setprecision(2) << Z[TIME_NOW][i];
    }

    ZUV_file << std::setw(5) << " "; // 插入一个宽度为5的空格
    ZUV_file << std::endl; // 换行
    ZUV_file << std::setw(5) << " "; // 插入一个宽度为5的空格
    ZUV_file << std::setw(3) << "W2="; // 输出"H2="字符串
    for (int i = 0; i < CEL; i++) {
      if (i % 10 == 0) {
        ZUV_file << std::endl; // 每输出10个数换行
        ZUV_file << std::setw(5) << " "; // 每行开头插入5个空格
      }
      ZUV_file << std::setw(10) << std::fixed << std::setprecision(2) << W[TIME_NOW][i];
    }

    cout <<"JT=" << jt << " KT=" << kt << endl;
  }
}
