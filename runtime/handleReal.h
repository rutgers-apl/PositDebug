#ifndef REAL_HEADER_H
#define REAL_HEADER_H 
#include <iostream>
#include <map>
#include <queue>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <vector>
#include <stack>
#include <list>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#define MMAP_FLAGS (MAP_PRIVATE| MAP_ANONYMOUS| MAP_NORESERVE)
#define MAX_STACK_SIZE 500
#define MAX_SIZE 1000

/* Assumption: float types are 4-byte aligned. */
const size_t SS_PRIMARY_TABLE_ENTRIES = ((size_t) 4194304);//2^22
const size_t SS_SEC_TABLE_ENTRIES = ((size_t) 16*(size_t) 1024 * (size_t) 1024); // 2^24
const size_t PRIMARY_INDEX_BITS = 20;
const size_t SECONDARY_INDEX_BITS = 24;
const size_t SECONDARY_MASK = 0xffffff;

enum fp_op{FADD, FSUB, FMUL, FDIV, CONSTANT, SQRT, FLOOR, TAN, SIN, COS, ATAN, ABS, LOG, ASIN, EXP, POW, UNKNOWN};

FILE* m_errfile;
FILE* m_brfile;

struct error_info {
  double error;
  unsigned int cbad;
  unsigned int regime;
  unsigned int lineno;
  unsigned int colnumber;
  bool debugInfoAvail;
};

struct smem_entry{
  mpfr_t val;
  int error;
  posit32_t computed;
  unsigned int lineno;
  enum fp_op opcode;
  bool is_init;
  bool initV;
  size_t timestamp;
  struct smem_entry* lhs;
  struct smem_entry* rhs;
};

smem_entry *m_shadow_stack;
smem_entry **m_shadow_memory;

# if !defined(PREC_128) && !defined(PREC_256) && !defined(PREC_512) && !defined(PREC_1024) 
#  define PREC_512
# endif

# if !defined(PRECISION) 
# ifdef PREC_128
#   define PRECISION 128
# endif
# ifdef PREC_256
#   define PRECISION 256
# endif
# ifdef PREC_512
#   define PRECISION 512
# endif
# ifdef PREC_1024
#   define PRECISION 1024
# endif
#endif


size_t m_stack_top = 0;
size_t timestampCounter = 0;
bool m_init_flag = false;

size_t infCount = 0;
size_t nanCount = 0;
size_t ccCount = 0;
size_t errorCount63 = 0;
size_t errorCount55 = 0;
size_t errorCount45 = 0;
size_t errorCount35 = 0;
size_t flipsCount = 0;
size_t convCount = 0;
size_t countMinPosMaxPos = 0;


std::list <smem_entry*> m_expr;

extern "C" void pd_init();
unsigned long m_ulpd(double, double);
int m_update_error(smem_entry *, double, unsigned int); 
void m_print_real(mpfr_t*);
int m_isnan(mpfr_t);
extern "C" unsigned int pd_check_error(void *, posit32_t);
unsigned int m_count_p32_regime_bits(posit32_t);

/* posit_t type */

#include <stdlib.h>
#include <math.h>
#include "softposit.h"

typedef struct{
  posit32_t pos;
  size_t arr[2];
}posit_t;

extern "C" posit_t rapl_convertDoubleToP32(double a);
extern "C" double rapl_convertP32ToDouble(posit_t a);
extern "C" int rapl_p32_to_i32(posit_t a);
extern "C" posit_t rapl_p32_sqrt(posit_t a);
extern "C" posit_t rapl_p32_add(posit_t a, posit_t b);
extern "C" posit_t rapl_p32_sub(posit_t a, posit_t b);
extern "C" posit_t rapl_p32_mul(posit_t a, posit_t b);
extern "C" posit_t rapl_p32_div(posit_t a, posit_t b);
extern "C" posit_t p32_fabs(posit_t val);
extern "C" posit_t p32_cos(posit_t val);
extern "C" posit_t p32_acos(posit_t val);
extern "C" posit_t p32_cos(posit_t val);
extern "C" posit_t p32_sin(posit_t val);
extern "C" posit_t p32_exp(posit_t val);
extern "C" posit_t p32_tan(posit_t val);
extern "C" posit_t p32_atan(posit_t val);
extern "C" posit_t p32_log(posit_t val1);
extern "C" posit_t p32_log10(posit_t val1);
extern "C"  posit_t p32_pow(posit_t val1, posit_t val2);
extern "C"  posit_t p32_atan2(posit_t val1, posit_t val2);
extern "C"  posit_t p32_floor(posit_t val);
extern "C"  bool rapl_p32_lt(posit_t a, posit_t b);
extern "C"  bool rapl_p32_le(posit_t a, posit_t b);
extern "C"  bool rapl_p32_eq(posit_t a, posit_t b);


#endif
