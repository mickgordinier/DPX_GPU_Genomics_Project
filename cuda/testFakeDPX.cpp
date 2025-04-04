#include <cassert>
#include <iostream>
#include "FakeDPX.hpp"

using std::cout;

int main() {
  cout << "TESTING ALL FakeDPX INSTRUCTIONS...\n";

  assert(FakeDPX::__vimax3_s32(1, 2, 3) == 3);
  assert(FakeDPX::__vimax3_s32(2, 3, 1) == 3);
  assert(FakeDPX::__vimax3_s32(-5, -10, -30) == -5);

  // Ensuring looking at both lower and upper
  // 0xFF00 is a negative number , 0x00FF is a positive number
  // 0xFFFF > 0xFFFE > 0xFFFD  (Assuming 2's complement)
  assert(FakeDPX::__vimax3_s16x2(0, 0x00FF00FF, 0xFF00FF00) == 0x00FF00FF);
  assert(FakeDPX::__vimax3_s16x2(0, 0xFFFF00FF, 0xFFFFFF00) == 0x000000FF);
  assert(FakeDPX::__vimax3_s16x2(0xFFFD00FF, 0xFFFE00FF, 0xFFFFFF00) == 0xFFFF00FF);

  assert(FakeDPX::__vimax3_u32(1, 2, 3) == 3);
  assert(FakeDPX::__vimax3_u32(2, 3, 0) == 3);

  assert(FakeDPX::__vimax3_u16x2(0, 0x00FF00FF, 0xFF00FF00) == 0xFF00FF00);
  assert(FakeDPX::__vimax3_u16x2(0, 0xFFFF00FF, 0xFFFFFF00) == 0xFFFFFF00);
  assert(FakeDPX::__vimax3_u16x2(0xFFFD00FF, 0xFFFE00FF, 0xFFFFFF00) == 0xFFFFFF00);

  assert(FakeDPX::__vimin3_s32(1, 2, 3) == 1);
  assert(FakeDPX::__vimin3_s32(2, 3, 1) == 1);
  assert(FakeDPX::__vimin3_s32(-5, -10, -30) == -30);

  assert(FakeDPX::__vimin3_s16x2(0, 0x00FF00FF, 0xFF00FF00) == 0xFF00FF00);
  assert(FakeDPX::__vimin3_s16x2(0, 0xFFFF00FF, 0xFFFFFF00) == 0xFFFFFF00);
  assert(FakeDPX::__vimin3_s16x2(0xFFFD00FF, 0xFFFE00FF, 0xFFFFFF00) == 0xFFFDFF00);

  assert(FakeDPX::__vimin3_u32(1, 2, 3) == 1);
  assert(FakeDPX::__vimin3_u32(2, 3, 0) == 0);

  assert(FakeDPX::__vimin3_u16x2(0, 0x00FF00FF, 0xFF00FF00) == 0);
  assert(FakeDPX::__vimin3_u16x2(0, 0xFFFF00FF, 0xFFFFFF00) == 0);
  assert(FakeDPX::__vimin3_u16x2(0xFFFD00FF, 0xFFFE00FF, 0xFFFFFF00) == 0XFFFD00FF);
  
  assert(FakeDPX::__vimax_s32_relu(1, 2) == 2);
  assert(FakeDPX::__vimax_s32_relu(2, 3) == 3);
  assert(FakeDPX::__vimax_s32_relu(-10, -30) == 0);

  assert(FakeDPX::__vimax_s16x2_relu(0x00FF00FF, 0xFF00FF00) == 0x00FF00FF);
  assert(FakeDPX::__vimax_s16x2_relu(0xFFFF00FF, 0xFFFFFF00) == 0x000000FF);
  assert(FakeDPX::__vimax_s16x2_relu(0xFFFD00FF, 0xFFFFFF00) == 0x000000FF);
  
  assert(FakeDPX::__vimin_s32_relu(1, 2) == 1);
  assert(FakeDPX::__vimin_s32_relu(2, 3) == 2);
  assert(FakeDPX::__vimin_s32_relu(-10, -30) == 0);

  assert(FakeDPX::__vimin_s16x2_relu(0x00FF00FF, 0xFF00FF00) == 0);
  assert(FakeDPX::__vimin_s16x2_relu(0xFFFF00FF, 0xFFFFFF00) == 0);
  assert(FakeDPX::__vimin_s16x2_relu(0xFFFD00FF, 0xFFFF0001) == 0x00000001);



  cout << "PASSED ALL ASSERTIONS FOR INSTRUCTION CHECKING!!\n";

}