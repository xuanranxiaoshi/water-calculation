#pragma once

#include "common.hh"
using std::cout;
using std::endl;

static std::fstream file_to_read;
static std::string current_line;
static std::string current_path;

#define TIME_NOW (t % 2)
#define TIME_PREV ((t - 1) % 2)
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
