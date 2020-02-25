#ifndef __SOFTPOSIT_EXT_H__
#define __SOFTPOSIT_EXT_H__

#include <stdlib.h>
#include <math.h>
#include "softposit.h"

typedef struct{
  posit32_t pos;
  size_t arr[2];
}posit_t;

posit_t rapl_convertDoubleToP32(double a);
double rapl_convertP32ToDouble(posit_t a);
int rapl_p32_to_i32(posit_t a);
posit_t rapl_p32_sqrt(posit_t a);
posit_t rapl_p32_add(posit_t a, posit_t b);
posit_t rapl_p32_sub(posit_t a, posit_t b);
posit_t rapl_p32_mul(posit_t a, posit_t b);
posit_t rapl_p32_div(posit_t a, posit_t b);
posit_t p32_fabs(posit_t val);
posit_t p32_cos(posit_t val);
posit_t p32_acos(posit_t val);
posit_t p32_cos(posit_t val);
posit_t p32_sin(posit_t val);
posit_t p32_exp(posit_t val);
posit_t p32_tan(posit_t val);
posit_t p32_atan(posit_t val);
posit_t p32_log(posit_t val1);
posit_t p32_log10(posit_t val1);
posit_t p32_pow(posit_t val1, posit_t val2);
posit_t p32_atan2(posit_t val1, posit_t val2);
posit_t p32_floor(posit_t val);
bool rapl_p32_lt(posit_t a, posit_t b);
bool rapl_p32_le(posit_t a, posit_t b);
bool rapl_p32_eq(posit_t a, posit_t b);

#endif
