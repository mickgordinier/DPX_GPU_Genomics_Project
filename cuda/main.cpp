#include "LinearSmithWaterman.h"

int main(){
    LinearSmithWaterman LSW("AATCG", "AACG", 3, -1, -2);
    LSW.align();
}