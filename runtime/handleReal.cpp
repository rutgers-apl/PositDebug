#include <string.h>
#include <fstream>
#include <queue>
#include <iostream>
#include <stdlib.h>
#include <execinfo.h>
#include <limits.h>
#include <sys/time.h>
#include <sys/ioctl.h>
#include <linux/perf_event.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <asm/unistd.h>
#include "softposit.h"
#include "handleReal.h"

#define debug 0
#define debugtrace 1
#define debugerror 0
#define ERRORTHRESHOLD 35
#define ERRORTHRESHOLD1 35
#define ERRORTHRESHOLD2 45
#define ERRORTHRESHOLD3 55
#define ERRORTHRESHOLD4 63

mpfr_t op1_mpfr, op2_mpfr, res_mpfr;
mpfr_t computed, temp_diff;




/***MPFR functions***/
int m_isnan(mpfr_t real){
    return mpfr_nan_p(real);
}

void m_set_mpfr(mpfr_t *val1, mpfr_t *val2) {
  mpfr_set(*val1, *val2, MPFR_RNDN);
}

void m_set_shadow(smem_entry *op, double d, unsigned int lineno) {
  mpfr_set_d(op->val, d, MPFR_RNDN);
  op->error = 0;
  op->lineno = lineno;
  op->opcode = CONSTANT;
  op->computed = convertDoubleToP32(d);
  op->initV = true;
  op->error = 0;
  if(debug){
    std::cout<<"m_set_shadow: setting:"<<d<<"\n";
    std::cout<<"\n";
  }
}

double m_get_double(mpfr_t mpfr_val) { return mpfr_get_d(mpfr_val, MPFR_RNDN); }

long double m_get_longdouble(smem_entry *real) {
  return mpfr_get_ld(real->val, MPFR_RNDN);
}

smem_entry* m_get_shadowaddress(size_t address){

  size_t addrInt = (address/sizeof(posit_t));
  size_t primary_index = (addrInt >> SECONDARY_INDEX_BITS );
  smem_entry* primary_ptr = m_shadow_memory[primary_index];
  if (primary_ptr == NULL) {
    size_t sec_length = (SS_SEC_TABLE_ENTRIES) * sizeof(smem_entry);
    primary_ptr = (smem_entry*)mmap(0, sec_length, PROT_READ| PROT_WRITE,
        MAP_PRIVATE|MAP_ANONYMOUS|MAP_NORESERVE, -1, 0);
    m_shadow_memory[primary_index] = primary_ptr;
  }
  size_t offset = (addrInt) & SECONDARY_MASK;
  smem_entry* realAddr = primary_ptr + offset;
  if(!realAddr->is_init){
    realAddr->is_init = true;
    mpfr_init2(realAddr->val, PRECISION);
  }
  return realAddr;
}

void m_print_real(mpfr_t mpfr_val){
  mpfr_out_str (stdout, 10, 15, mpfr_val, MPFR_RNDN);
} 

void handle_math_d(fp_op opCode, smem_entry *op,
          smem_entry *res,  posit_t* computedRes,
                  unsigned int lineno) {

  switch(opCode){
    case SQRT:
      mpfr_sqrt(res->val, op->val, MPFR_RNDN);
      break;
    case FLOOR:
      mpfr_floor(res->val, op->val);
      break;
    case TAN:
      mpfr_tan(res->val, op->val, MPFR_RNDN);
      break;
    case SIN:
      mpfr_sin(res->val, op->val, MPFR_RNDN);
      break;
    case COS:
      mpfr_cos(res->val, op->val, MPFR_RNDN);
      break;
    case ATAN:
      mpfr_atan(res->val, op->val, MPFR_RNDN);
      break;
    case ABS:
      mpfr_abs(res->val, op->val, MPFR_RNDN);
      break;
    case LOG:
      mpfr_log(res->val, op->val, MPFR_RNDN);
      break;
    case ASIN: 
      mpfr_asin(res->val, op->val, MPFR_RNDN);
      break;
    case EXP: 
      mpfr_exp(res->val, op->val, MPFR_RNDN);
      break;
    default:
      break;
  }
  unsigned int regime = m_count_p32_regime_bits(computedRes->pos);
  int bitsError = m_update_error(res, convertP32ToDouble(computedRes->pos), regime);
  res->lineno = lineno;
  res->timestamp = timestampCounter++;
  res->error = bitsError;
  res->opcode = opCode;
  res->computed = computedRes->pos;
  res->lhs = op;
  res->rhs = nullptr;
}
//Count preicison bits in <32,2>
unsigned int m_count_prec_bits(unsigned ps, unsigned nbit) {
  unsigned bitsAvail = nbit;
  unsigned countPrec = 0;
  unsigned countExp = 0;

  // (1) Get Sign bit: Logical shift right by (nbit - 1)
  unsigned sign = ps >> (nbit - 1);
  bitsAvail--;

  // (1.5) If Sign bit is 1, perform two's complement
  if (sign == 1) {
    ps = -ps;
    ps = (ps << (32 - nbit)) >> (32 - nbit);
  }

  // (2) Get regime bits:
  // Move regime+exponent+fraction bits all the way to the left.
  unsigned num_regime = 0;
  unsigned temp = ps << (32 - bitsAvail);
  unsigned regime_sign = temp >> 31;

  while (bitsAvail > 0 && (temp >> 31) == regime_sign) {
    temp = temp << 1;
    num_regime++;
    bitsAvail--;
  }

  while (bitsAvail > 0 &&  countExp < 2) {
    temp = temp << 1;
    countExp++;
    bitsAvail--;
  }

  while (bitsAvail > 0) {
    temp = temp << 1;
    countPrec++;
    bitsAvail--;
  }
  return countPrec;
}
// Given an nbit posit, count the number of regime bits.
// Take note that 0 and NaR will report (n - 1) bits of regime.
// nbit can be up to 32 bits.
unsigned int m_count_regime_bits(unsigned ps, unsigned nbit) {
  unsigned bitsAvail = nbit;

  // (1) Get Sign bit: Logical shift right by (nbit - 1)
  unsigned sign = ps >> (nbit - 1);
  bitsAvail--;

  // (1.5) If Sign bit is 1, perform two's complement
  if (sign == 1) {
    ps = -ps;
    ps = (ps << (32 - nbit)) >> (32 - nbit);
  }

  // (2) Get regime bits:
  // Move regime+exponent+fraction bits all the way to the left.
  unsigned num_regime = 0;
  unsigned temp = ps << (32 - bitsAvail);
  unsigned regime_sign = temp >> 31;

  while (bitsAvail > 0 && (temp >> 31) == regime_sign) {
    temp = temp << 1;
    num_regime++;
    bitsAvail--;
  }
  return num_regime;
}

unsigned int m_count_p32_regime_bits(posit32_t val) {
      return m_count_regime_bits(val.v, 32);
}

unsigned int m_count_p32_prec_bits(posit32_t val) {
      return m_count_prec_bits(val.v, 32);
}

unsigned int m_get_exact_bits(
    posit32_t op,
    smem_entry *shadow){
  double positVal_d = convertP32ToDouble(op);

  mpfr_set_d(computed, positVal_d, MPFR_RNDN);

  mpfr_sub(temp_diff, shadow->val, computed, MPFR_RNDN);

  mpfr_exp_t exp_real = mpfr_get_exp(shadow->val);
  mpfr_exp_t exp_computed = mpfr_get_exp(computed);
  mpfr_exp_t exp_diff = mpfr_get_exp(temp_diff);

  int precBits = m_count_p32_prec_bits(op);
  if(mpfr_cmp(computed, shadow->val) == 0){
    return precBits;
  }
  else if(exp_real != exp_computed){

    return 0;
  }
  else{
    if(mpfr_cmp_ui(temp_diff, 0) != 0) {
      if(precBits < abs(exp_real -  exp_diff)){
        return precBits;
      }
      else{
        return abs(exp_real -  exp_diff);   
      }
    }
    else{
      return 0;
    }
  }
}

mpfr_exp_t m_get_cancelled_bits(posit32_t op1, posit32_t op2, posit32_t res){
  mpfr_set_d(op1_mpfr, convertP32ToDouble(op1), MPFR_RNDN);

  mpfr_set_d(op2_mpfr, convertP32ToDouble(op2), MPFR_RNDN);

  mpfr_set_d(res_mpfr, convertP32ToDouble(res), MPFR_RNDN);

  mpfr_exp_t exp_op1 = mpfr_get_exp(op1_mpfr);
  mpfr_exp_t exp_op2 = mpfr_get_exp(op2_mpfr);
  mpfr_exp_t exp_res = mpfr_get_exp(res_mpfr);

  mpfr_exp_t max_exp;
  if( mpfr_regular_p(op1_mpfr) == 0 ||
      mpfr_regular_p(op2_mpfr) == 0 ||
      mpfr_regular_p(res_mpfr) == 0)
    return 0;

  if(exp_op1 > exp_op2)
    max_exp = exp_op1;
  else
    max_exp = exp_op2;

  if(max_exp > exp_res)
    return abs(max_exp - exp_res);
  else
    return 0;
}

unsigned int m_get_cbad(mpfr_exp_t cbits, 
                              unsigned int ebitsOp1, 
                              unsigned int ebitsOp2){
  unsigned int min_ebits;
  if (ebitsOp1 > ebitsOp2)
    min_ebits = ebitsOp2;
  else
    min_ebits = ebitsOp1;
  int badness = 1 + cbits - min_ebits;
  if(badness > 0)
    return badness;
  else
    return 0;
}
bool m_is_min_pos_or_max_pos(posit32_t val){
  if(p32_eq(val, castP32(0x7FFFFFFF))){
    return true; //maxpos || -minpos
  }
  if(p32_eq(val, castP32(0x80000001))){ //-maxpos || minpos
    return true;
  }
  return false;
}

void m_compute(fp_op opCode, 
	          smem_entry *op1, 
                  smem_entry *op2, 
 		  smem_entry *res,
		  posit_t* computedRes,
		  unsigned int lineno) {

  switch(opCode) {                                                                                            
    case FADD: 
      mpfr_add (res->val, op1->val, op2->val, MPFR_RNDN);

      break;
    case FSUB: 
      mpfr_sub (res->val, op1->val, op2->val, MPFR_RNDN);
      break;
    case FMUL: 
      mpfr_mul (res->val, op1->val, op2->val, MPFR_RNDN);
      break;
    case FDIV: 
      mpfr_div (res->val, op1->val, op2->val, MPFR_RNDN);
      break;
    default:
      // do nothing
      break;
  } 
  unsigned int regime = m_count_p32_regime_bits(computedRes->pos);
  int bitsError = m_update_error(res, convertP32ToDouble(computedRes->pos), regime);
  res->lineno = lineno;
  res->error = bitsError;
  res->opcode = opCode;
  res->timestamp = timestampCounter++;
  res->computed = computedRes->pos;
  res->lhs = op1;
  res->rhs = op2;
}

long m_get_mpfr_long(smem_entry *op) {
  return mpfr_get_si(op->val, MPFR_RNDN);
}

extern "C" void pd_init_mpfr(smem_entry *op) {
  if (!m_init_flag) 
    pd_init();
  mpfr_init2(op->val, PRECISION);
}


extern "C" void pd_clear_mpfr(smem_entry *op) {
  mpfr_clear(op->val);
}

extern "C" void pd_set_const(void *toAddr, double d, unsigned int lineno) {
  if (!m_init_flag) 
    pd_init();
  size_t toAddrInt = (size_t) toAddr;
  smem_entry* dst = m_get_shadowaddress(toAddrInt);

  mpfr_set_d(dst->val, d, MPFR_RNDN);
  dst->lineno = lineno;
  dst->opcode = CONSTANT;
  dst->lhs = nullptr;
  dst->rhs = nullptr;
  dst->timestamp = timestampCounter++;
  dst->computed = convertDoubleToP32(d);
  dst->initV = true;
  dst->error = 0;
  if(debug){
    std::cout<<"pd_set_const: setting"<<d<<"\n";
    std::cout<<"\n";
  }
}

extern "C" void pd_copy_phi(smem_entry* src, smem_entry* dst){
  if(src != NULL){
    m_set_mpfr(&(dst->val), &(src->val));
  }
}

extern "C" void pd_store_real(void* toAddr, smem_entry* src){
  size_t toAddrInt = (size_t) toAddr;
  smem_entry* dst = m_get_shadowaddress(toAddrInt);
  if(src != NULL){
    //copy val
    m_set_mpfr(&(dst->val), &(src->val));
    //copy everything else except res key and opcode
    dst->error = src->error;
    dst->lineno = src->lineno;
    dst->computed = src->computed;
  }
  else{
    std::cout<<"Error !!! pd_store_real trying to read invalid memory\n";
  }
}

//copy metadata from adress to shadow address
//TODO:I think we don't need this, but I am leaving it as it is for now
extern "C" void pd_load_shadow(smem_entry *src, void *Addr){
  if(src != NULL){
    size_t AddrInt = (size_t) Addr;
    smem_entry* dst = m_get_shadowaddress(AddrInt);

    m_set_mpfr(&(src->val), &(dst->val));

    src->error = dst->error;
    src->lineno = dst->lineno;
    src->computed = dst->computed;

  }
  else{
    if(debug){
      printf("__load_p32:Error !!! __load_p32 trying to load from invalid stack\n");
    }
  }
}

unsigned int m_check_cc(posit32_t op1,
                posit32_t op2, 
                posit32_t res, 
                smem_entry *shadowOp1,
                smem_entry *shadowOp2,
                smem_entry *shadowVal){

  // If op1 or op2 is NaR, then it is not catastrophic cancellation
  if (isNaRP32UI(op1.v) || isNaRP32UI(op2.v)) return 0;
  // If result is 0 and it has error, then it is catastrophic cancellation
  if (p32_eq(res, convertDoubleToP32(0)) && shadowVal->error != 0) return 1;
  
  unsigned int ebitsOp1 = m_get_exact_bits(op1, shadowOp1);
  unsigned int ebitsOp2 = m_get_exact_bits(op2, shadowOp2);
  mpfr_exp_t cbits = m_get_cancelled_bits(op1, op2, res);
  unsigned int cbad = m_get_cbad(cbits, ebitsOp1, ebitsOp2);
  return cbad;
}

extern "C" void pd_handle_memcpy(void* dstAddr, void *srcAddr){
  size_t srcInt = (size_t) srcAddr;
  smem_entry* src = m_get_shadowaddress(srcInt);

  size_t dstInt = (size_t) dstAddr;
  smem_entry* dst = m_get_shadowaddress(dstInt);

  m_set_mpfr(&(dst->val), &(src->val));
  dst->timestamp = timestampCounter++;
  dst->computed = src->computed;
  dst->error = src->error;
  dst->lineno = src->lineno;
  dst->lhs = src->lhs;
  dst->rhs = src->rhs;
  dst->opcode = src->opcode;
  dst->is_init = true;
}

extern "C" void pd_mpfr_rapl_p32_add( void* op1Addr, 
                                void* op2Addr, 
                                void* resAddr, 
                                posit_t* posOp1, 
                                posit_t* posOp2, 
                                posit_t* computedRes, 
                                unsigned long long int instId, 
                                bool debugInfoAvail, 
                                unsigned int lineno, 
                                unsigned int colnumber) {

  size_t Op1Int = (size_t) op1Addr;
  smem_entry* op1Idx = m_get_shadowaddress(Op1Int);

  size_t Op2Int = (size_t) op2Addr;
  smem_entry* op2Idx = m_get_shadowaddress(Op2Int);

  size_t ResInt = (size_t) resAddr;
  smem_entry* res = m_get_shadowaddress(ResInt);

  m_compute(FADD, op1Idx, op2Idx, res, computedRes, lineno);
 
  if(m_is_min_pos_or_max_pos(computedRes->pos)){
    countMinPosMaxPos++;
  }
  posit32_t zero = convertDoubleToP32(0);
  unsigned int cbad = 0;
  //check for subtraction
  if((p32_lt(posOp1->pos, zero) && 
        !p32_le(posOp2->pos, zero)) || 
        (!p32_le(posOp1->pos, zero) && 
         p32_lt(posOp2->pos, zero))){
    cbad = m_check_cc(posOp1->pos, posOp2->pos, computedRes->pos, op1Idx, op2Idx, res);
    if(cbad > 0)
      ccCount++;
  }

  if (isNaRP32UI(computedRes->pos.v)) {
    nanCount++;
  }
}

extern "C" void pd_mpfr_rapl_p32_sub( void* op1Addr, 
                                void* op2Addr, 
                                void* resAddr, 
                                posit_t* posOp1, 
                                posit_t* posOp2, 
                                posit_t* computedRes, 
                                unsigned long long int instId, 
                                bool debugInfoAvail, 
                                unsigned int lineno, 
                                unsigned int colnumber) {
  size_t Op1Int = (size_t) op1Addr;
  smem_entry* op1Idx = m_get_shadowaddress(Op1Int);

  size_t Op2Int = (size_t) op2Addr;
  smem_entry* op2Idx = m_get_shadowaddress(Op2Int);

  //TODO: I could not find a way  to get type for memcpy operands
  //call void @llvm.memcpy.p0i8.p0i8.i64(i8* align 16 bitcast 
  //([10 x [3 x %struct.posit_t]]* @r_kj to i8*), i8* nonnull align 8 %0, i64 24, i1 false), !tbaa.struct !37, !psan_inst_id !43
  //To handle such cases: Adding a way to set value on demand
  /*Update: I check either operand of memcpy, it seems to work but will work in future cases not sure*/

  size_t ResInt = (size_t) resAddr;
  smem_entry* res = m_get_shadowaddress(ResInt);

  m_compute(FSUB, op1Idx, op2Idx, res, computedRes, lineno);

  if(m_is_min_pos_or_max_pos(computedRes->pos)){
    countMinPosMaxPos++;
  }
  //check for subtraction
  posit32_t zero = convertDoubleToP32(0);
  unsigned int cbad = 0;
  if((!p32_le(posOp1->pos, zero) && 
        !p32_le(posOp2->pos, zero)) || 
        (p32_lt(posOp1->pos, zero) && 
         p32_lt(posOp2->pos, zero))){
    cbad = m_check_cc(posOp1->pos, posOp2->pos, computedRes->pos, op1Idx, op2Idx, res);
    if(cbad > 0)
      ccCount++;
  }
}

extern "C" void pd_mpfr_rapl_p32_mul( void* op1Addr, 
                                void* op2Addr, 
                                void* resAddr, 
                                posit_t* posOp1, 
                                posit_t* posOp2, 
                                posit_t* computedRes, 
                                unsigned long long int instId, 
                                bool debugInfoAvail, 
                                unsigned int lineno, 
                                unsigned int colnumber) {
  size_t Op1Int = (size_t) op1Addr;
  smem_entry* op1Idx = m_get_shadowaddress(Op1Int);

  size_t Op2Int = (size_t) op2Addr;
  smem_entry* op2Idx = m_get_shadowaddress(Op2Int);

  size_t ResInt = (size_t) resAddr;
  smem_entry* res = m_get_shadowaddress(ResInt);

  m_compute(FMUL, op1Idx, op2Idx, res, computedRes, lineno);
  if(m_is_min_pos_or_max_pos(computedRes->pos)){
    countMinPosMaxPos++;
  }
}

extern "C" void pd_mpfr_rapl_p32_div( void* op1Addr, 
                                void* op2Addr, 
                                void* resAddr, 
                                posit_t* posOp1, 
                                posit_t* posOp2, 
                                posit_t* computedRes, 
                                unsigned long long int instId, 
                                bool debugInfoAvail, 
                                unsigned int lineno, 
                                unsigned int colnumber) {
  size_t Op1Int = (size_t) op1Addr;
  smem_entry* op1Idx = m_get_shadowaddress(Op1Int);

  size_t Op2Int = (size_t) op2Addr;
  smem_entry* op2Idx = m_get_shadowaddress(Op2Int);

  size_t ResInt = (size_t) resAddr;
  smem_entry* res = m_get_shadowaddress(ResInt);

  m_compute(FDIV, op1Idx, op2Idx, res ,computedRes, lineno);
  if(m_is_min_pos_or_max_pos(computedRes->pos)){
    countMinPosMaxPos++;
  }
}

extern "C" bool pd_mpfr_eq(void* op1Addr, void* op2Addr,
			  bool computedRes, 
        unsigned long long int instId,
        bool debugInfoAvail, 
        unsigned int lineno, 
        unsigned int colnumber){
  size_t Op1Int = (size_t) op1Addr;
  smem_entry* op1 = m_get_shadowaddress(Op1Int);

  size_t Op2Int = (size_t) op2Addr;
  smem_entry* op2 = m_get_shadowaddress(Op2Int);

  bool realRes = false;
  int ret = mpfr_cmp(op1->val, op2->val);

  if(ret == 0) 
    realRes = true;
  if(realRes != computedRes){
    flipsCount++;
  }
  return realRes;
}

extern "C" bool pd_mpfr_le(void* op1Addr, void* op2Addr,
			  bool computedRes, 
        unsigned long long int instId,
        bool debugInfoAvail,  
        unsigned int lineno, 
        unsigned int colnumber){

  size_t Op1Int = (size_t) op1Addr;
  smem_entry* op1 = m_get_shadowaddress(Op1Int);

  size_t Op2Int = (size_t) op2Addr;
  smem_entry* op2 = m_get_shadowaddress(Op2Int);
  bool realRes = false;
  int ret = mpfr_cmp(op1->val, op2->val);

  if(ret <= 0) 
    realRes = true;

  if(realRes != computedRes){
    flipsCount++;
  }
  return realRes;
}


extern "C" bool pd_mpfr_lt(void* op1Addr, void* op2Addr,
			  bool computedRes, 
        unsigned long long int instId, 
        bool debugInfoAvail, 
        unsigned int lineno, 
        unsigned int colnumber){
  size_t Op1Int = (size_t) op1Addr;
  smem_entry* op1 = m_get_shadowaddress(Op1Int);

  size_t Op2Int = (size_t) op2Addr;
  smem_entry* op2 = m_get_shadowaddress(Op2Int);

  bool realRes = false;
  int ret = mpfr_cmp(op1->val, op2->val);

  if(ret < 0) 
    realRes = true;

  if(realRes != computedRes){
    flipsCount++;
  }
  return realRes;
}

extern "C" void pd_mpfr_sqrt(void* op1Addr, void* resAddr, posit_t* computedRes,  
    unsigned long long int instId, bool debugInfoAvail, unsigned int lineno, unsigned int colnumber){
  size_t Op1Int = (size_t) op1Addr;
  smem_entry* op1Idx = m_get_shadowaddress(Op1Int);

  size_t ResInt = (size_t) resAddr;
  smem_entry* res = m_get_shadowaddress(ResInt);

  handle_math_d(SQRT, op1Idx, res, computedRes, lineno);

}

extern "C" void pd_mpfr_fabs(void* op1Addr, void* resAddr, posit_t* computedRes, 
    unsigned long long int instId, bool debugInfoAvail, unsigned int lineno, unsigned int colnumber){
  size_t Op1Int = (size_t) op1Addr;
  smem_entry* op1Idx = m_get_shadowaddress(Op1Int);

  size_t ResInt = (size_t) resAddr;
  smem_entry* res = m_get_shadowaddress(ResInt);

  handle_math_d(ABS, op1Idx, res, computedRes, lineno);
}

extern "C" void pd_mpfr_pow(void* op1Addr, void* op2Addr, void* resAddr, posit_t* computedRes, 
    unsigned long long int instId, bool debugInfoAvail, unsigned int lineno, unsigned int colnumber){
  size_t Op1Int = (size_t) op1Addr;
  smem_entry* op1Idx = m_get_shadowaddress(Op1Int);

  size_t Op2Int = (size_t) op2Addr;
  smem_entry* op2Idx = m_get_shadowaddress(Op2Int);

  size_t ResInt = (size_t) resAddr;
  smem_entry* res = m_get_shadowaddress(ResInt);

  mpfr_pow(res->val, op1Idx->val, op2Idx->val, MPFR_RNDN);

  unsigned int regime = m_count_p32_regime_bits(computedRes->pos);
  int bitsError = m_update_error(res, convertP32ToDouble(computedRes->pos), regime);
  res->lineno = lineno;
  res->timestamp = timestampCounter++;
  res->error = bitsError;
  res->opcode = POW;
  res->computed = computedRes->pos;
  res->lhs = op1Idx;
  res->rhs = op2Idx;
}

extern "C" void pd_mpfr_exp(void* op1Addr, void* resAddr, posit_t* computedRes, 
    unsigned long long int instId, bool debugInfoAvail, unsigned int lineno, unsigned int colnumber){
  size_t Op1Int = (size_t) op1Addr;
  smem_entry* op1Idx = m_get_shadowaddress(Op1Int);

  size_t ResInt = (size_t) resAddr;
  smem_entry* res = m_get_shadowaddress(ResInt);

  if(debug){
    std::cout<<"handle_math_d: exp operand:"<<"\n";
    m_print_real(op1Idx->val);
  }
  handle_math_d(EXP, op1Idx, res, computedRes, lineno);

}

extern "C" void pd_mpfr_floor(void* op1Addr, void* resAddr, posit_t* computedRes, 
    unsigned long long int instId, bool debugInfoAvail, unsigned int lineno, unsigned int colnumber){
  size_t Op1Int = (size_t) op1Addr;
  smem_entry* op1Idx = m_get_shadowaddress(Op1Int);

  size_t ResInt = (size_t) resAddr;
  smem_entry* res = m_get_shadowaddress(ResInt);

  handle_math_d(FLOOR, op1Idx, res, computedRes, lineno);
}

extern "C" void pd_mpfr_atan(void* op1Addr, void* resAddr, posit_t* computedRes, 
    unsigned long long int instId, bool debugInfoAvail, unsigned int lineno, unsigned int colnumber){
  size_t Op1Int = (size_t) op1Addr;
  smem_entry* op1Idx = m_get_shadowaddress(Op1Int);

  size_t ResInt = (size_t) resAddr;
  smem_entry* res = m_get_shadowaddress(ResInt);

  handle_math_d(ATAN, op1Idx, res, computedRes, lineno);
}

extern "C" void pd_mpfr_tan(void* op1Addr, void* resAddr, posit_t* computedRes, 
    unsigned long long int instId, bool debugInfoAvail, unsigned int lineno, unsigned int colnumber){
  size_t Op1Int = (size_t) op1Addr;
  smem_entry* op1Idx = m_get_shadowaddress(Op1Int);

  size_t ResInt = (size_t) resAddr;
  smem_entry* res = m_get_shadowaddress(ResInt);

  handle_math_d(TAN, op1Idx, res, computedRes, lineno);
}

extern "C" void pd_mpfr_cos(void* op1Addr, void* resAddr, posit_t* computedRes, 
    unsigned long long int instId, bool debugInfoAvail, unsigned int lineno, unsigned int colnumber){
  size_t Op1Int = (size_t) op1Addr;
  smem_entry* op1Idx = m_get_shadowaddress(Op1Int);

  size_t ResInt = (size_t) resAddr;
  smem_entry* res = m_get_shadowaddress(ResInt);

  handle_math_d(COS, op1Idx, res, computedRes, lineno);
}

extern "C" void pd_mpfr_sin(void* op1Addr, void* resAddr, posit_t* computedRes, 
    unsigned long long int instId, bool debugInfoAvail, unsigned int lineno, unsigned int colnumber){
  size_t Op1Int = (size_t) op1Addr;
  smem_entry* op1Idx = m_get_shadowaddress(Op1Int);

  size_t ResInt = (size_t) resAddr;
  smem_entry* res = m_get_shadowaddress(ResInt);

  handle_math_d(SIN, op1Idx, res, computedRes, lineno);
}

extern "C" long pd_mpfr_rapl_p32_to_i32(void* op1Addr, int computed, unsigned long long int instId, unsigned int lineno, unsigned int colnumber){
  size_t Op1Int = (size_t) op1Addr;
  smem_entry* op1 = m_get_shadowaddress(Op1Int);

  long ret = m_get_mpfr_long(op1);
  if(ret != computed){
    convCount++;
  }
  return ret;
}

std::string m_get_string_opcode(size_t opcode){
  switch(opcode){
    case FADD:
      return "P32_ADD";
    case FMUL:
      return "P32_MUL";
    case FSUB:
      return "P32_SUB";
    case FDIV:
      return "P32_DIV";
    case CONSTANT:
      return "CONSTANT";
    case SQRT:  
      return "SQRT";
    case FLOOR:  
      return "FLOOR";
    case TAN:  
      return "TAN";
    case SIN:  
      return "SIN";
    case COS:  
      return "COS";
    case ATAN:  
      return "ATAN";
    case ABS:  
      return "ABS";
    case LOG:  
      return "LOG";
    case ASIN:  
      return "ASIN";
    case EXP:  
      return "EXP";
    case POW:  
      return "POW";
    default:
      return "Unknown";
  }
}

int m_get_depth(smem_entry *current){
  int depth = 0;
  m_expr.push_back(current);
  int level;
  while(!m_expr.empty()){
    level = m_expr.size();
    std::cout<<"\n";
    while(level > 0){
      smem_entry *cur = m_expr.front();
      if(cur == NULL){
        return depth;
      }
      if(cur->lhs != NULL){
        if(cur->lhs->timestamp < cur->timestamp){
          m_expr.push_back(cur->lhs);
        }
      }
      if(cur->rhs != NULL){
        if(cur->rhs->timestamp < cur->timestamp){
          m_expr.push_back(cur->rhs);
        }
      }
      m_expr.pop_front();
      level--;
    }
    depth++;
  }
  return depth;
}

extern "C" void pd_trace(smem_entry *current){
  m_expr.push_back(current);
  int level;
  while(!m_expr.empty()){
    level = m_expr.size();
    smem_entry *cur = m_expr.front();
    std::cout<<"\n";
    if(cur == NULL){
      return;
    }
    std::cout<<" "<<cur->lineno<<" "<<m_get_string_opcode(cur->opcode)<<" ";
    fflush(stdout);
    if(cur->lhs != NULL){
      std::cout<<" "<<cur->lhs->lineno<<" ";
      if(cur->lhs->timestamp < cur->timestamp){
        m_expr.push_back(cur->lhs);
        fflush(stdout);
      }
    }
    if(cur->rhs != NULL){
      std::cout<<" "<<cur->rhs->lineno<<" ";
      if(cur->rhs->timestamp < cur->timestamp){
        m_expr.push_back(cur->rhs);
        fflush(stdout);
      }
    }
    std::cout<<"(real:";
    m_print_real(cur->val);
    printf(" computed: %e", convertP32ToDouble(cur->computed));
    std::cout<<", error:"<<cur->error<<" "<<")";
    fflush(stdout);
    m_expr.pop_front();
    level--;
  }
  int depth = m_get_depth(current);
  std::cout<<"depth:"<<depth<<"\n";
}

extern "C" void pd_func_init(long totalArgs) {
  m_stack_top = m_stack_top + totalArgs;
}

extern "C" void pd_func_exit(long totalArgs) {
  m_stack_top = m_stack_top - totalArgs;
}

//copy local var to shadow stack for return
extern "C" void pd_set_return(smem_entry* src, size_t totalArgs, posit32_t op) {
  smem_entry *dst = &(m_shadow_stack[m_stack_top - totalArgs]); //save return m_stack_top - totalArgs 
  if(src != NULL){
    m_set_mpfr(&(dst->val), &(src->val));
    dst->computed = src->computed;
    dst->error = src->error;
    dst->lineno = src->lineno;
    dst->opcode = src->opcode;

  }
  else{
    std::cout<<"__set_return copying src is null:"<<"\n";
    m_set_shadow(dst, convertP32ToDouble(op), 0); //when we return tmp 
    //we don't need to set metadata as it is null
  }
}

//copy return from shadow stack to local alloca var
extern "C" void pd_get_return(smem_entry* dst) {
  smem_entry *src = &(m_shadow_stack[m_stack_top]); //save return m_stack_top - totalArgs 
  m_set_mpfr(&(dst->val), &(src->val));
  dst->computed = src->computed;
  dst->error = src->error;
  dst->lineno = src->lineno;
  dst->timestamp = timestampCounter++;
  dst->opcode = src->opcode;
}

extern "C" void pd_get_argument(size_t argIdx, void *dstAddr) {
  size_t dstInt = (size_t) dstAddr;
  smem_entry* dst = m_get_shadowaddress(dstInt);
  smem_entry *src = &(m_shadow_stack[m_stack_top-argIdx]);

  m_set_mpfr(&(dst->val), &(src->val));
  dst->timestamp = timestampCounter++;
  dst->computed = src->computed;
  dst->error = src->error;
  dst->lineno = src->lineno;
  dst->lhs = src->lhs;
  dst->rhs = src->rhs;
  dst->opcode = src->opcode;

  if(!src->is_init){ //caller maybe is not instrumented for set_arg
    std::cout<<"Error !!!! __get_argument trying to get argument which was not set\n";
  }
}

extern "C" void pd_set_argument(size_t argIdx, void* srcAddr) {
  size_t srcInt = (size_t) srcAddr;
  smem_entry* src = m_get_shadowaddress(srcInt);
  smem_entry *dst = &(m_shadow_stack[argIdx+m_stack_top]);
  assert(argIdx < MAX_STACK_SIZE && "Arg index is more than MAX_STACK_SIZE");

  if(src != NULL){
    m_set_mpfr(&(dst->val), &(src->val));
    dst->timestamp = timestampCounter++;
    dst->computed = src->computed;
    dst->error = src->error;
    dst->lineno = src->lineno;
    dst->lhs = src->lhs;
    dst->rhs = src->rhs;
    dst->opcode = src->opcode;
    dst->is_init = true;
  }
  else{
    std::cout<<"__set_argument Error!!! trying to set null argument\n";
  }
}

extern "C" void pd_init() {
  if (!m_init_flag) {

    m_errfile = fopen ("error.log","w");

    //    printf("sizeof posit_t %lu", sizeof(posit_t));
    m_init_flag = true;
    size_t length = MAX_STACK_SIZE * sizeof(smem_entry);
    size_t memLen = SS_PRIMARY_TABLE_ENTRIES * sizeof(smem_entry);
    m_shadow_stack =
      (smem_entry *)mmap(0, length, PROT_READ | PROT_WRITE, MMAP_FLAGS, -1, 0);
    m_shadow_memory =
      (smem_entry **)mmap(0, memLen, PROT_READ | PROT_WRITE, MMAP_FLAGS, -1, 0);
    assert(m_shadow_stack != (void *)-1);
    assert(m_shadow_memory != (void *)-1);

    for(int i =0; i<MAX_STACK_SIZE; i++){
      mpfr_init2(m_shadow_stack[i].val, PRECISION);
    }
    m_stack_top = 0;

    mpfr_init2(op1_mpfr, PRECISION);
    mpfr_init2(op2_mpfr, PRECISION);
    mpfr_init2(res_mpfr, PRECISION);
    mpfr_init2(temp_diff, PRECISION);
    mpfr_init2(computed, PRECISION);
  }
}


unsigned long m_ulpd(double x, double y) {
  if (x == 0)
    x = 0; // -0 == 0
  if (y == 0)
    y = 0; // -0 == 0

  /* if (x != x && y != y) return 0; */
  if (x != x)
    return ULLONG_MAX - 1; // Maximum error
  if (y != y)
    return ULLONG_MAX - 1; // Maximum error

  long long xx = *((long long *)&x);
  xx = xx < 0 ? LLONG_MIN - xx : xx;

  long long yy = *((long long *)&y);
  yy = yy < 0 ? LLONG_MIN - yy : yy;
  return xx >= yy ? xx - yy : yy - xx;
}

extern "C" void pd_finish() {

  fprintf(m_errfile, "Error above %d bits found %zd\n", ERRORTHRESHOLD4, errorCount63);
  fprintf(m_errfile, "Error above %d bits found %zd\n", ERRORTHRESHOLD3, errorCount55);
  fprintf(m_errfile, "Error above %d bits found %zd\n", ERRORTHRESHOLD2, errorCount45);
  fprintf(m_errfile, "Error above %d bits found %zd\n", ERRORTHRESHOLD1, errorCount35);
  fprintf(m_errfile, "Total CC found %zd\n", ccCount);
  fprintf(m_errfile, "Total NaR found %zd\n", nanCount);
  fprintf(m_errfile, "Total branch flips found %zd\n", flipsCount);
  fprintf(m_errfile, "Total conversion errors %zd\n", convCount);
  fprintf(m_errfile, "Total minpos or maxpos found %zd\n", countMinPosMaxPos);
  mpfr_clear(computed);
  mpfr_clear(temp_diff);
  mpfr_clear(op1_mpfr);
  mpfr_clear(op2_mpfr);
  mpfr_clear(res_mpfr);
  fclose(m_errfile);
}

extern "C" unsigned int  pd_check_branch(bool realBr, 
			bool computedBr, 
			smem_entry *realRes1, 
			smem_entry *realRes2){
  if(realBr != computedBr)
    return 1;
  return 0;
}

extern "C" unsigned int  pd_check_conversion(long real, 
					     long computed, 
					     smem_entry *realRes){
  if(real != computed){
    return 1;
  }
  return 0;
}

extern "C" unsigned int  pd_check_error(void* realAddr, 
					posit32_t computedRes){
			
  size_t realInt = (size_t) realAddr;
  smem_entry* realRes = m_get_shadowaddress(realInt);
  if(realRes->error >= ERRORTHRESHOLD4){
    errorCount63++;
    return 1; 
  }
  else if(realRes->error >= ERRORTHRESHOLD3){
    errorCount55++;
    return 2; 
  }
  else if(realRes->error >= ERRORTHRESHOLD2){
    errorCount45++; 
    return 3; 
  }
  else if(realRes->error >= ERRORTHRESHOLD1){
    errorCount35++; 
    return 4; 
  }
  return 0;
}

//I need a way to compare with herbgrind, hence storing error
//in a similar way as herbgrind
int m_update_error(smem_entry *real, double computedVal, unsigned int regime){
   //compare
  double shadowRounded = m_get_double(real->val);
  unsigned long ulpsError = m_ulpd(shadowRounded, computedVal);

  double bitsError = log2(ulpsError + 1);
 
  if(m_isnan(real->val)){
    if(debug)
      std::cout<<"nan computed\n";
  } 
  if (debugerror){
    std::cout<<"\nThe shadow value is ";
    m_print_real(real->val);
    if (computedVal != computedVal){
      std::cout<<", but NaN was computed.\n";
    } else {
      std::cout<<", but ";
      printf("%e", computedVal);
      std::cout<<" was computed.\n";
      std::cout<<"m_update_error: computedVal:"<<computedVal<<"\n";
    }
    printf("%f bits error (%lu ulps)\n",
        bitsError, ulpsError);
    std::cout<<"m_update_error: regime:"<<regime<<"\n";
    std::cout<<"****************\n\n"; 
  }
  return bitsError;
}

/***posit_t functions */

extern "C" posit_t rapl_p32_sqrt(posit_t a){
  posit_t c;
  c.pos = p32_sqrt(a.pos);
  return c;
}

extern "C" posit_t rapl_p32_add(posit_t a, posit_t b){
  posit_t c;
  c.pos = p32_add(a.pos, b.pos);
  return c;
}

extern "C" posit_t rapl_p32_sub(posit_t a, posit_t b){
  posit_t c;
  c.pos = p32_sub(a.pos, b.pos);
  return c;
}

extern "C" posit_t rapl_p32_mul(posit_t a, posit_t b){
  posit_t c;
  c.pos = p32_mul(a.pos, b.pos);
  return c;
}

extern "C"  posit_t rapl_p32_div(posit_t a, posit_t b){
  posit_t c;
  c.pos = p32_div(a.pos, b.pos);
  return c;
}

extern "C" posit_t rapl_convertDoubleToP32(double a){
  posit_t c;
  c.pos = convertDoubleToP32(a);
  return c;
}

extern "C" double rapl_convertP32ToDouble(posit_t a){
    return convertP32ToDouble(a.pos);
}

extern "C"  int rapl_p32_to_i32(posit_t a){
  return p32_to_i32(a.pos);
}

extern "C" bool rapl_p32_lt(posit_t a, posit_t b){
  return p32_lt(a.pos, b.pos);
}

extern "C" bool rapl_p32_le(posit_t a, posit_t b){
  return p32_le(a.pos, b.pos);
}
extern "C" bool rapl_p32_eq(posit_t a, posit_t b){
  return p32_eq(a.pos, b.pos);
}
extern "C" posit_t p32_strtod(const char* str, char** endptr){
  double d = strtod(str, endptr);
  return rapl_convertDoubleToP32(d);
}
extern "C"  posit_t p32_cos(posit_t val){
  double d = cos(rapl_convertP32ToDouble(val));
  return rapl_convertDoubleToP32(d);
}

extern "C" posit_t p32_sin(posit_t val){
  double d = sin(rapl_convertP32ToDouble(val));
  return rapl_convertDoubleToP32(d);
}

extern "C" posit_t p32_atan2(posit_t val1, posit_t val2){
	double d = atan2(rapl_convertP32ToDouble(val1), rapl_convertP32ToDouble(val2));
	return rapl_convertDoubleToP32(d);
}
extern "C"  posit_t p32_tan(posit_t val){
	double d = tan(rapl_convertP32ToDouble(val));
	return rapl_convertDoubleToP32(d);
}

extern "C" posit_t p32_acos(posit_t val){
	double d = acos(rapl_convertP32ToDouble(val));
	return rapl_convertDoubleToP32(d);
}

extern "C" posit_t p32_fabsf(posit_t val){
	double d = fabs(rapl_convertP32ToDouble(val));
	return rapl_convertDoubleToP32(d);
}

extern "C" posit_t p32_fabs(posit_t val){
	double d = fabs(rapl_convertP32ToDouble(val));
	return rapl_convertDoubleToP32(d);
}

extern "C" posit_t p32_exp(posit_t val){
	double d = exp(rapl_convertP32ToDouble(val));
	return rapl_convertDoubleToP32(d);
}
extern "C"  posit_t p32_expf(posit_t val){
	double d = expf(rapl_convertP32ToDouble(val));
	return rapl_convertDoubleToP32(d);
}

extern "C"  posit_t p32_floor(posit_t val){
	double d = floor(rapl_convertP32ToDouble(val));
	return rapl_convertDoubleToP32(d);
}

extern "C" posit_t p32_cosf(posit_t val){
	double d = cosf(rapl_convertP32ToDouble(val));
	return rapl_convertDoubleToP32(d);
}

extern "C"  posit_t p32_sinf(posit_t val){
	double d = sinf(rapl_convertP32ToDouble(val));
	return rapl_convertDoubleToP32(d);
}

extern "C"  posit_t p32_pow(posit_t val1, posit_t val2){
	double d = pow(rapl_convertP32ToDouble(val1), rapl_convertP32ToDouble(val2));
	return rapl_convertDoubleToP32(d);
}
extern "C" posit_t p32_atan2f(posit_t val1, posit_t val2){
	double d = atan2f(rapl_convertP32ToDouble(val1), rapl_convertP32ToDouble(val2));
	return rapl_convertDoubleToP32(d);
}

extern "C"  posit_t p32_atan(posit_t val1){
	double d = atan(rapl_convertP32ToDouble(val1));
	return rapl_convertDoubleToP32(d);
}
extern "C"  posit_t p32_log10(posit_t val1){
	double d = log10(rapl_convertP32ToDouble(val1));
	return rapl_convertDoubleToP32(d);
}
extern "C"  posit_t p32_log(posit_t val1){
	double d = log(rapl_convertP32ToDouble(val1));
	return rapl_convertDoubleToP32(d);
}
extern "C" bool isNaN(posit_t val){
  posit32_t nar = castP32(0x80000000);
  if(p32_eq(val.pos, nar))
    return true;
  else
    return false;
}

