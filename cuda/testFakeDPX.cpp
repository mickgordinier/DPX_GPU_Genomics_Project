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

  cout << "PASSED ALL ASSERTIONS FOR INSTRUCTION CHECKING!!\n";

}