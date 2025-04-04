#include <cassert>
#include <iostream>
#include "FakeDPX.hpp"

using std::cout;

int main() {
  cout << "TESTING ALL FakeDPX INSTRUCTIONS...\n";
  assert(FakeDPX::__vimax3_s32(1, 2, 3) == 3);
  assert(FakeDPX::__vimax3_s32(2, 3, 1) == 3);
  assert(FakeDPX::__vimax3_s32(-5, -10, -30) == -5);

  cout << "PASSED ALL ASSERTIONS FOR INSTRUCTION CHECKING!!\n";

}