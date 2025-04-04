#include <algorithm>
#include <iostream>
#include "FakeDPX.hpp"

using std::cout; using std::endl;
using std::max;
using std::min;

/* 3 PARAMETERS FUNCTIONS */

int FakeDPX::__vimax3_s32(const int a, const int b, const int c) {
  return max(max(a, b), c);
}

unsigned int FakeDPX::__vimax3_s16x2(const unsigned int a, const unsigned int b, const unsigned int c) {

  unsigned int ret;

  short a_high = a >> 16;
  short b_high = b >> 16;
  short c_high = c >> 16;
  ret = max(max(a_high, b_high), c_high) << 16;

  short a_low = a & 0xFFFF;
  short b_low = b & 0xFFFF;
  short c_low = c & 0xFFFF;

  return ret | max(max(a_low, b_low), c_low);
}

unsigned int FakeDPX::__vimax3_u32(const unsigned int a, const unsigned int b, const unsigned int c) {
  return max(max(a, b), c);
}

unsigned int FakeDPX::__vimax3_u16x2(const unsigned int a, const unsigned int b, const unsigned int c) {

  unsigned int ret;

  unsigned short a_high = a >> 16;
  unsigned short b_high = b >> 16;
  unsigned short c_high = c >> 16;
  ret = max(max(a_high, b_high), c_high) << 16;

  
  unsigned short a_low = a & 0xFFFF;
  unsigned short b_low = b & 0xFFFF;
  unsigned short c_low = c & 0xFFFF;

  return ret | (max(max(a_low, b_low), c_low) & 0xFFFF);
}

int FakeDPX::__vimin3_s32(const int a, const int b, const int c) {
  return min(min(a, b), c);
}

unsigned int FakeDPX::__vimin3_s16x2(const unsigned int a, const unsigned int b, const unsigned int c) {

  unsigned int ret;

  short a_high = a >> 16;
  short b_high = b >> 16;
  short c_high = c >> 16;
  ret = min(min(a_high, b_high), c_high) << 16;

  short a_low = a & 0xFFFF;
  short b_low = b & 0xFFFF;
  short c_low = c & 0xFFFF;

  return ret | (min(min(a_low, b_low), c_low) & 0xFFFF);
}

unsigned int FakeDPX::__vimin3_u32(const unsigned int a, const unsigned int b, const unsigned int c) {
  return min(min(a, b), c);
}

unsigned int FakeDPX::__vimin3_u16x2(const unsigned int a, const unsigned int b, const unsigned int c) {

  unsigned int ret;

  unsigned short a_high = a >> 16;
  unsigned short b_high = b >> 16;
  unsigned short c_high = c >> 16;
  ret = min(min(a_high, b_high), c_high) << 16;

  unsigned short a_low = a & 0xFFFF;
  unsigned short b_low = b & 0xFFFF;
  unsigned short c_low = c & 0xFFFF;

  return ret | (min(min(a_low, b_low), c_low) & 0xFFFF);
}


/* 2 PARAMETERS + ReLU FUNCTIONS */

int FakeDPX::__vimax_s32_relu(const int a, const int b) {
  return FakeDPX::__vimax3_s32(a, b, 0);
}

unsigned int FakeDPX::__vimax_s16x2_relu(const unsigned int a, const unsigned int b) {
  return FakeDPX::__vimax3_s16x2(a, b, 0);
}

int FakeDPX::__vimin_s32_relu(const int a, const int b) {
  return max(min(a, b), 0);
}

unsigned int FakeDPX::__vimin_s16x2_relu(const unsigned int a, const unsigned int b) {

  unsigned int ret;

  short a_high = a >> 16;
  short b_high = b >> 16;
  ret = max(min(a_high, b_high), (short)0) << 16;

  short a_low = a & 0xFFFF;
  short b_low = b & 0xFFFF;

  return ret | (max(min(a_low, b_low), (short)0) & 0xFFFF);
}


/* 3 PARAMETERS + ReLU FUNCTIONS */

int FakeDPX::__vimax3_s32_relu(const int a, const int b, const int c) {
  return max(FakeDPX::__vimax3_s32(a, b, c), 0);
}

// I know I don't need to do 2 ReLUs but didn't want to write it all out again
unsigned int FakeDPX::__vimax3_s16x2_relu(const unsigned int a, const unsigned int b, const unsigned int c) {
  return FakeDPX::__vimax_s16x2_relu(FakeDPX::__vimax_s16x2_relu(a, b), c);
}

int FakeDPX::__vimin3_s32_relu(const int a, const int b, const int c) {
  return max(FakeDPX::__vimin3_s32(a, b, c), 0);
}

// I know I don't need to do 2 ReLUs but didn't want to write it all out again
unsigned int FakeDPX::__vimin3_s16x2_relu(const unsigned int a, const unsigned int b, const unsigned int c) {
  return FakeDPX::__vimin_s16x2_relu(FakeDPX::__vimin_s16x2_relu(a, b), c);
}


/* 2 PARAMETERS + RETURNING SMALLER/LARGER PARAMETER FUNCTIONS */

int FakeDPX::__vibmax_s32(const int a, const int b, bool *const pred) {
  if (a >= b) {
    *pred = true;
    return a;
  }

  *pred = false;
  return b;
}

unsigned int FakeDPX::__vibmax_u32(const unsigned int a, const unsigned int b, bool *const pred) {
  if (a >= b) {
    *pred = true;
    return a;
  }

  *pred = false;
  return b;
}

int FakeDPX::__vibmin_s32(const int a, const int b, bool *const pred) {
  if (a <= b) {
    *pred = true;
    return a;
  }

  *pred = false;
  return b;
}

unsigned int FakeDPX::__vibmin_u32(const unsigned int a, const unsigned int b, bool *const pred) {
  if (a <= b) {
    *pred = true;
    return a;
  }

  *pred = false;
  return b;
}

unsigned int FakeDPX::__vibmax_s16x2(const unsigned int a, const unsigned int b, bool *const pred_hi, bool *const pred_low) {

  unsigned int ret;

  short a_high = a >> 16;
  short b_high = b >> 16;

  if (a_high >= b_high) {
    *pred_hi = true;
    ret = a_high << 16;
  } else {
    *pred_hi = false;
    ret = b_high << 16;
  }

  short a_low = a & 0xFFFF;
  short b_low = b & 0xFFFF;

  if (a_low >= b_low) {
    *pred_low = true;
    return ret | (a_low & 0xFFFF);
  }

  *pred_low = false;
  return ret | (b_low & 0xFFFF);
}

unsigned int FakeDPX::__vibmax_u16x2(const unsigned int a, const unsigned int b, bool *const pred_hi, bool *const pred_low) {

  unsigned int ret;

  unsigned short a_high = a >> 16;
  unsigned short b_high = b >> 16;

  if (a_high >= b_high) {
    *pred_hi = true;
    ret = a_high << 16;
  } else {
    *pred_hi = false;
    ret = b_high << 16;
  }

  unsigned short a_low = a & 0xFFFF;
  unsigned short b_low = b & 0xFFFF;

  if (a_low >= b_low) {
    *pred_low = true;
    return ret | (a_low & 0xFFFF);
  }

  *pred_low = false;
  return ret | (b_low & 0xFFFF);
}

unsigned int FakeDPX::__vibmin_s16x2(const unsigned int a, const unsigned int b, bool *const pred_hi, bool *const pred_low) {

  unsigned int ret;

  short a_high = a >> 16;
  short b_high = b >> 16;

  if (a_high <= b_high) {
    *pred_hi = true;
    ret = a_high << 16;
  } else {
    *pred_hi = false;
    ret = b_high << 16;
  }

  short a_low = a & 0xFFFF;
  short b_low = b & 0xFFFF;

  if (a_low <= b_low) {
    *pred_low =true;
    return ret | (a_low & 0xFFFF);
  }

  *pred_low =false;
  return ret | (b_low & 0xFFFF);
}

unsigned int FakeDPX::__vibmin_u16x2(const unsigned int a, const unsigned int b, bool *const pred_hi, bool *const pred_low) {

  unsigned int ret;

  unsigned short a_high = a >> 16;
  unsigned short b_high = b >> 16;

  if (a_high <= b_high) {
    *pred_hi = true;
    ret = a_high << 16;
  } else {
    *pred_hi = false;
    ret = b_high << 16;
  }

  unsigned short a_low = a & 0xFFFF;
  unsigned short b_low = b & 0xFFFF;

  if (a_low <= b_low) {
    *pred_low =true;
    return ret | (a_low & 0xFFFF);
  }

  *pred_low =false;
  return ret | (b_low & 0xFFFF);
}


/* 3 PARAMETERS, COMPARING (FIRST + SECOND) WITH THIRD FUNCTIONS */

int FakeDPX::__viaddmax_s32(const int a, const int b, const int c) {
  return max(a+b, c);
}

unsigned int FakeDPX::__viaddmax_u32(const unsigned int a, const unsigned int b, const unsigned int c) {
  return max(a+b, c);
}

unsigned int FakeDPX::__viaddmax_s16x2(const unsigned int a, const unsigned int b, const unsigned int c) {

  unsigned int ret;

  short ab_high = (short)(a >> 16) + (short)(b >> 16);
  short c_high = c >> 16;
  ret = max(ab_high, c_high) << 16;

  short ab_low = (short)(a & 0xFFFF) + (short)(b & 0xFFFF);
  short c_low = c & 0xFFFF;

  return ret | (max(ab_low, c_low) & 0xFFFF);
}

unsigned int FakeDPX::__viaddmax_u16x2(const unsigned int a, const unsigned int b, const unsigned int c) {

  unsigned int ret;

  unsigned short ab_high = (unsigned short)(a >> 16) + (unsigned short)(b >> 16);
  unsigned short c_high = c >> 16;
  ret = max(ab_high, c_high) << 16;

  unsigned short ab_low = (unsigned short)(a & 0xFFFF) + (unsigned short)(b & 0xFFFF);
  unsigned short c_low = c & 0xFFFF;

  return ret | (max(ab_low, c_low) & 0xFFFF);
}

int FakeDPX::__viaddmin_s32(const int a, const int b, const int c) {
  return min(a+b, c);
}

unsigned int FakeDPX::__viaddmin_u32(const unsigned int a, const unsigned int b, const unsigned int c) {
  return min(a+b, c);
}

unsigned int FakeDPX::__viaddmin_s16x2(const unsigned int a, const unsigned int b, const unsigned int c) {

  unsigned int ret;

  short ab_high = (short)(a >> 16) + (short)(b >> 16);
  short c_high = c >> 16;
  ret = min(ab_high, c_high) << 16;

  short ab_low = (short)(a & 0xFFFF) + (short)(b & 0xFFFF);
  short c_low = c & 0xFFFF;

  return ret | (min(ab_low, c_low) & 0xFFFF);
}

unsigned int FakeDPX::__viaddmin_u16x2(const unsigned int a, const unsigned int b, const unsigned int c) {

  unsigned int ret;

  unsigned short ab_high = (unsigned short)(a >> 16) + (unsigned short)(b >> 16);
  unsigned short c_high = c >> 16;
  ret = min(ab_high, c_high) << 16;

  unsigned short ab_low = (unsigned short)(a & 0xFFFF) + (unsigned short)(b & 0xFFFF);
  unsigned short c_low = c & 0xFFFF;

  return ret | (min(ab_low, c_low) & 0xFFFF);
}
  

/* 3 PARAMETERS, COMPARING (FIRST + SECOND) WITH THIRD AND ReLU FUNCTIONS */

int FakeDPX::__viaddmax_s32_relu(const int a, const int b, const int c) {
  return max(max(a + b, c), 0);
}

unsigned int FakeDPX::__viaddmax_s16x2_relu(const unsigned int a, const unsigned int b, const unsigned int c) {

  unsigned int ret;

  short ab_high = (short)(a >> 16) + (short)(b >> 16);
  short c_high = c >> 16;
  ret = max(max(ab_high, c_high), (short)0) << 16;

  short ab_low = (short)(a & 0xFFFF) + (short)(b & 0xFFFF);
  short c_low = c & 0xFFFF;

  return ret | (max(max(ab_low, c_low), (short)0) & 0xFFFF);
}

int FakeDPX::__viaddmin_s32_relu(const int a, const int b, const int c) {
  return max(min(a + b, c), 0);
}

unsigned int FakeDPX::__viaddmin_s16x2_relu(const unsigned int a, const unsigned int b, const unsigned int c) {

  unsigned int ret;

  short ab_high = (short)(a >> 16) + (short)(b >> 16);
  short c_high = c >> 16;
  ret = max(min(ab_high, c_high), (short)0) << 16;

  short ab_low = (short)(a & 0xFFFF) + (short)(b & 0xFFFF);
  short c_low = c & 0xFFFF;

  return ret | (max(min(ab_low, c_low), (short)0) & 0xFFFF);
}