#include <stdio.h>
#include <math.h>

#include "softposit.h"
#include "softposit_ext.h"
int main() {
  posit_t x;posit_t y;
  x = rapl_convertDoubleToP32 (1.0E+16);;
  posit_t tmp1 = rapl_convertDoubleToP32 (1);
  posit_t tmp2 = rapl_p32_add(x,tmp1);
  posit_t tmp3 = rapl_p32_sub(rapl_p32_sqrt(tmp2),rapl_p32_sqrt(x));
  y = tmp3;
  double tmp4 = rapl_convertP32ToDouble (y);
  printf("%e\n", tmp4);
  return 0;
}
