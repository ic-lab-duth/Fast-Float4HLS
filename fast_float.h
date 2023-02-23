// Copyright 2022 Democritus University of Thrace
// Integrated Circuits Lab
// 
// Fast-Float4HLS v1.0
// -------------------
// Floating Point Library of fast operators for High Level Synthesis
#ifndef __FASTFLOAT_H__
#define __FASTFLOAT_H__

#include <ac_fixed.h>
#include <ac_std_float.h>

template <int W>
bool lzc_reduce (ac_int<W,false> a) {
  ac_int<W/2,false> odd, even;


  // split odd even bits of a
  #pragma unroll yes
  SPLIT:for (int i=0; i < W/2; i++) {
    even[i] = a[2*i];
    odd[i] = a[2*i+1];
  }

  // prefix AND from MSB-to-LSB of inverted even bits
  even = even.bit_complement();
  ac_int<W/2,false> even_prefix;
  ac_int<1,false> t=true;
  
  #pragma unroll yes
  PREFIXAND:for (int i=W/2-1; i >=0; i--) {
    if (i == W/2-1) t = even[i];
    else t = t & even[i];
    even_prefix[i] = t;
  }

  // fix alignment of prefixed and terms
  even_prefix = even_prefix >> 1;
  even_prefix[W/2-1] = 1;

  // prepare terms for each bit position
  ac_int<W/2,false> tmp = even_prefix & odd;

  // return the wide OR of those terms
  return tmp.or_reduce();
}

// Version 1: try to do it
template<int N>
struct lzc_s {
  template<typename T>
  static void lzc(ac_int<N,false> a, T &out) {
    ac_int<N/2, false> a0;
    out[ac::log2_ceil<N>::val] = lzc_reduce<N>(a);

    #pragma unroll yes
    for (int i = 0; i < N / 2; i++) {
      a0[i] = a[2*i] | a[2*i + 1];
    }
    lzc_s<N/2>::lzc(a0,out);
  }
};

template<>
struct lzc_s<1> {
  template<typename T>
  static void lzc(ac_int<1,false> a,  T &out) {
    out[0] = a[0];
  }
};


template<int N=8>
ac_int<ac::log2_ceil<N>::val+1,false> lzcount(ac_int<N,false> x) {
  ac_int<ac::log2_ceil<N>::val+1,false> b,res;
  lzc_s<N>::lzc(x,b);

  // reverse bits
  #pragma unroll yes
  for (int i=0; i< ac::log2_ceil<N>::val+1; i++) {
    res[i] = b[ac::log2_ceil<N>::val - i]; 
  }

  // complement them and send them out
  return (res.or_reduce() == 0) ? (ac_int<ac::log2_ceil<N>::val+1,false>)0 :  res.bit_complement();
}


template<int N>
struct max_s {
  template<typename T>
  static T max(T *a) {
    T m0 = max_s<N/2>::max(a);
    T m1 = max_s<N-N/2>::max(a+N/2);

    return m0 > m1 ? m0 : m1;
  }
};

template<> 
struct max_s<1> {
  template<typename T>
  static T max(T *a) {
    return a[0];
  }
};

template<int N, typename T>
T max(T *a) {
  return max_s<N>::max(a);
};

// Fast Float class
template<int M, int E>
class fast_float {
public:
  typedef ac_int<M,false> man_t;
  typedef ac_int<E,false> exp_t;
  typedef ac_int<1,false> sgn_t;
  
  man_t mantissa;
  exp_t exponent;
  sgn_t sign;
   
  static const int width = M + E + 1;
  static const int man_width = M;
  static const int exp_width = E;
  static const int e_bias = (1 << (E-1)) - 1;
  static const int frac_bits = M + 1;

  enum RND_ENUM {EVEN, ODD, INF};
  
  // Constructors
  fast_float() {};
  fast_float(const ac_int<M+E+1,false> &in) {
    this->operator=(in);
  }
  fast_float(const fast_float<M,E> &in) {
    this->operator=(in);
  }
  ~fast_float() {};

  void set(sgn_t s, exp_t e, man_t m) {
    sign = s;
    exponent = e;
    mantissa = m;
  }

  void setNan()    { this->set(0, (1<<E)-1, 1); }
  void setPosINF() { this->set(0, (1<<E)-1, 0); }
  void setNegINF() { this->set(1, (1<<E)-1, 0); }
  
  /*
  * .to_float() is a non-synthesizable function that
  * is used for converting the fast_float into a c++ float
  */
  float to_float() {
    
    ac_int<32,false> data;
    for (int i=0; i<M; i++)
      data[22-i] = mantissa[M-1-i];

    for (int i=M; i<23; i++)
      data[22-i] = 0;
    
    ac_int<8,false> tmpExp = (exponent==0) ? ( ac_int<8,false>)0 : ( ac_int<8,false>)(exponent - e_bias + 127);
    for (int i=0; i<8; i++)
      data[23+i] = (exponent == (1<<E)+-1) ? 1 : tmpExp[i];

    data[31] = sign;

    ac_std_float<32,8> ieee_data;
    ieee_data.set_data(data);
    return  ieee_data.to_float();
  }

  /** FLOATING POINT OPERATIONS **/

  /* Floating Point Addition
  *  adds the value of 'a' to this fast_float and
  *  return the result to the 'output'. 
  *  The addition is implemented using dual path (Far-Near).
  */
  template<RND_ENUM RND_MODE=EVEN, bool DENORMALS=false>
  void fpa_dual(const fast_float<M,E> a, fast_float<M,E> &output) {
  
    // Set create needed lengths.
    const int signif_uniform_W = M + 4;

    // Calculate effective sign.
    const ac_int<1,false> sign_eff =  sign ^ a.sign; 

    // Infinity check:
    ac_int<1,false> inf_sign =  sign;
    bool Is_inf = false;
    if (( exponent).and_reduce()) {
      Is_inf = true;
    } else if ((a.exponent).and_reduce()) {
      inf_sign = a.sign;
      Is_inf = true;
    }
  

    // End of infinity check.

    // Create normalized significant with guard bits.

    ac_int<signif_uniform_W,true> norm_significant_1;
    ac_int<signif_uniform_W,true> norm_significant_2;

    norm_significant_1 =  mantissa;
    norm_significant_1 <<= 2;
    
    norm_significant_2 = a.mantissa;
    norm_significant_2 <<= 2;
    
    const ac_int<1,false> in_1_subnorm = ( exponent == 0);
    const ac_int<1,false> in_2_subnorm = (a.exponent == 0);
    const ac_int<1,false> both_ins_subnorm = in_1_subnorm & in_2_subnorm;
    const ac_int<1,false> only_one_in_subnorm = in_1_subnorm ^ in_2_subnorm;

    norm_significant_1[signif_uniform_W - 2] = (DENORMALS) ? ((in_1_subnorm) ? 0 : 1) : 1;
    norm_significant_2[signif_uniform_W - 2] = (DENORMALS) ? ((in_2_subnorm) ? 0 : 1) : 1;

    
    // Start of R-path:

    // Compute exponent difference.

    const ac_int<E + 1,true> R_exp_diff_1 =  exponent - a.exponent;
    const ac_int<E + 1,true> R_exp_diff_2 = a.exponent -  exponent;
    ac_int<E + 1,true> R_abs_exp_diff;

    // Initialize swap variables.

    ac_int<signif_uniform_W,true> R_nsignif_small;
    ac_int<signif_uniform_W,true> R_nsignif_large;
    ac_int<E,false> R_exp_large;
    ac_int<1,false> R_sign_large = 0;
    
    const int R_align_sft_extra_W = M + 3;
    
    ac_int<signif_uniform_W + R_align_sft_extra_W,true> R_small_signif_aligner;
    
    // Swap and allign.

    if (!R_exp_diff_1[E]) {

		R_nsignif_large = norm_significant_1;
        R_exp_large =  exponent;
        R_sign_large =  sign;
        R_nsignif_small = norm_significant_2;
        R_abs_exp_diff = R_exp_diff_1;

    } else if (!R_exp_diff_2[E]) {

        R_nsignif_large = norm_significant_2;
        R_exp_large = a.exponent;
        R_sign_large = a.sign;
        R_nsignif_small = norm_significant_1;
        R_abs_exp_diff = R_exp_diff_2;

    }

    R_small_signif_aligner = R_nsignif_small;
    
    if (R_abs_exp_diff < R_align_sft_extra_W) {
      R_small_signif_aligner <<= R_align_sft_extra_W;
      if (DENORMALS) R_small_signif_aligner <<= only_one_in_subnorm;
      R_small_signif_aligner >>= R_abs_exp_diff;
    }
    
    R_nsignif_small = R_small_signif_aligner.template slc<signif_uniform_W>(R_align_sft_extra_W);
    
    const ac_int<1,false> R_sticky = (R_small_signif_aligner.template slc<R_align_sft_extra_W>(0)).or_reduce();
    const ac_int<1,false> R_add_1 = !R_sticky;



    // Significant Addition. Compute sum and both rounded sums.

    ac_int<signif_uniform_W,false> R_sum;
    ac_int<signif_uniform_W,false> R_rounded_sum; // rounded sum in case no normalization left shift is needed.
    ac_int<signif_uniform_W,false> R_unnormal_rounded_sum; // rounded sum requires left shift normalition.
    
    if (sign_eff) {
    	
    	R_sum = R_nsignif_large + (~R_nsignif_small) + R_add_1;
    	R_rounded_sum = R_nsignif_large + (~R_nsignif_small) + R_add_1 + 4;
    	R_unnormal_rounded_sum = R_nsignif_large + (~R_nsignif_small) + R_add_1 + 2;
    	
    } else {
    	
    	R_sum = R_nsignif_large + R_nsignif_small;
    	R_rounded_sum = R_nsignif_large + R_nsignif_small + 4;
    	R_unnormal_rounded_sum = R_nsignif_large + R_nsignif_small + 2;
    	
	  }

    bool Ris_zero = !(R_sum.or_reduce());

    // Normalize and overflow check R_sum, R_rounded_sum
    // and R_unnormal_rounded_sum.

    const ac_int<E,false> R_exp_large_plus_1 = R_exp_large + 1;
    const ac_int<E,false> R_exp_large_min_1 = R_exp_large - 1;
	
    bool R_sum_overf = false;
    bool R_sum_is_unnormal = false;
    bool R_rounded_sum_overf = false;
    bool R_is_unnormal_rounded = false;
    
    const ac_int<1,false> R_sum_sticky_hold = R_sum[0];
    
    if (DENORMALS) {
      if (R_sum[signif_uniform_W - 1]) {

        R_sum >>= 1;
        R_sum[0] = R_sum[0] | R_sum_sticky_hold;
        R_sum_overf = true;

      } else if (both_ins_subnorm && R_sum[signif_uniform_W - 2]) {

          R_sum_overf = true;

      } else if (!both_ins_subnorm && !R_sum[signif_uniform_W - 2]) {

          R_sum <<= 1;
          R_sum_is_unnormal = true;

      }
    } else {
      if (R_sum[signif_uniform_W - 1]) {
          R_sum >>= 1;
          R_sum[0] = R_sum[0] | R_sum_sticky_hold;
          R_sum_overf = true;

      } else if (!R_sum[signif_uniform_W - 2]) {

          R_sum <<= 1;
          R_sum_is_unnormal = true;

      }
    }
    
    if (DENORMALS) {
      if (R_rounded_sum[signif_uniform_W - 1]) {

        R_rounded_sum >>= 1;
        R_rounded_sum_overf = true;

      } else if (both_ins_subnorm && R_rounded_sum[signif_uniform_W - 2]) {

          R_rounded_sum_overf = true;

      }
    } else {
      if (R_rounded_sum[signif_uniform_W - 1]) {

        R_rounded_sum >>= 1;
        R_rounded_sum_overf = true;

      }
    }
    
    if (DENORMALS) {
      if (!both_ins_subnorm && !R_unnormal_rounded_sum.template slc<2>(signif_uniform_W - 2)) {

          R_unnormal_rounded_sum <<= 1;
          R_is_unnormal_rounded = true;
      }
    } else {
      if (!R_unnormal_rounded_sum.template slc<2>(signif_uniform_W - 2)) {

          R_unnormal_rounded_sum <<= 1;
          R_is_unnormal_rounded = true;
      }
    }

    // Rounding Decision.
	
    bool rnd = false;
    
    switch (RND_MODE) {
      
      case EVEN:
        // Round to nearest tie to even.
          rnd = R_sum[1] & (R_sum[0] | R_sum[2] | R_sticky);
          break;
      case ODD:
          // Round to nearest tie to odd.
          rnd = R_sum[1] & (R_sum[0] | (!R_sum[2]) | R_sticky);
          break;
      case INF:
          // Round infinity.
          rnd = R_sum[1];
          break;

      }

    ac_int<M,false> R_magn_result;
    ac_int<E,false> R_exp_result;
    
    if (rnd) {

      if (R_is_unnormal_rounded) {
        
        R_magn_result = R_unnormal_rounded_sum.template slc<M>(2);
        R_exp_result = R_exp_large_min_1;

      } else {

        R_magn_result = R_rounded_sum.template slc<M>(2);
        R_exp_result = (R_rounded_sum_overf) ? R_exp_large_plus_1 : R_exp_large;
         
      }

    } else {

      R_magn_result = R_sum.template slc<M>(2);
        
      if (R_sum_overf) {

        R_exp_result = R_exp_large_plus_1;

      } else if (R_sum_is_unnormal) {

        R_exp_result = R_exp_large_min_1;

      } else {
        
        R_exp_result = R_exp_large;

      }

    }
    

    // End of R-path.


    // Start of N-path:

    // 2-bit exponent difference.

    const ac_int<3,true> N_exp_diff =  exponent.template slc<2>(0) - a.exponent.template slc<2>(0);

    // Swap and 1-bit allign.

    ac_int<signif_uniform_W,true> N_nsignif_small;
    ac_int<signif_uniform_W,true> N_nsignif_large;
    ac_int<E,false> N_exp_large;
    ac_int<1,false> N_sign_large = 0;

    if (N_exp_diff[1]) {

        N_nsignif_large = norm_significant_2;
        N_exp_large = a.exponent;
        N_sign_large = a.sign;
        N_nsignif_small = norm_significant_1;

    } else {

        N_nsignif_large = norm_significant_1;
        N_exp_large =  exponent;
        N_sign_large =  sign;
        N_nsignif_small = norm_significant_2;

    }
	  
    if (DENORMALS) {
      if (N_exp_diff.or_reduce() && (!only_one_in_subnorm)) {
    	  N_nsignif_small >>= 1;

      }
    } else {
      if (N_exp_diff.or_reduce()) {
        
        N_nsignif_small >>= 1;

      }
    }
    
	  // Significant substraction.
	
	  const ac_int<signif_uniform_W,true> N_large_sub_small = N_nsignif_large - N_nsignif_small;
    const ac_int<signif_uniform_W,true> N_small_sub_large = N_nsignif_small - N_nsignif_large;

    ac_int<signif_uniform_W - 1,false> N_chosen_sub;

    if (!N_small_sub_large[signif_uniform_W - 1]) {

        N_sign_large = !N_sign_large;
        N_chosen_sub = N_small_sub_large;

    } else if (!N_large_sub_small[signif_uniform_W - 1]) {
		
        N_chosen_sub = N_large_sub_small;
        
    }

    // Normalize chosen sub.

    const bool N_no_normalize = N_chosen_sub[signif_uniform_W - 2];

	  bool N_zero = false;

    // const ac_int<ac::nbits<signif_uniform_W-1>::val, false> N_lz_cnt = (N_chosen_sub).leading_sign(N_zero);
    // M+3 = 23 + 3 = 26 
    // M+3 = 7  + 3 = 10 
    // M+3 = 3  + 3 = 6 
    const int ww = (M+3 <= 8)  ? 8  :
                   (M+3 <= 16) ? 16 :
                   (M+3 <= 32) ? 32 : 64;
    const int uu = ww - M - 3;
    ac_int<ww,false> tcz = 0;
    #pragma hls_unroll
    for (int ii=0; ii<M+3; ii++)
      tcz[ii+uu]=N_chosen_sub[ii];
    const ac_int<ac::nbits<signif_uniform_W-1>::val, false> N_lz_cnt = lzcount<ww>(tcz);
    // REMOVE LINE -3
    // ADD LINES -2 AND -1
    // CHANGES FOR LEADING ZERO COUNT

    if (N_zero) {
		
		N_exp_large = 0;
		N_sign_large = 0;
		
    } else {

      N_chosen_sub <<= N_lz_cnt;
      N_exp_large -= N_lz_cnt;
        
    }
  
    if (DENORMALS)
      N_chosen_sub >>= ((N_exp_large == 0) && !both_ins_subnorm) ? 1 : 0;


    // End of N-path.

    // Result Decision.

    const bool Is_R1 = (R_abs_exp_diff >= 2);
    const bool Is_R2 = sign_eff & !Is_R1 & N_no_normalize;
    const bool Is_R = !sign_eff | Is_R1 | Is_R2;

    
    
    if (R_exp_result.and_reduce() & Is_R) {
		Is_inf = true;	
    }

    sgn_t tmp_s = (Is_inf) ? inf_sign                      : ((Is_R) ? R_sign_large  : N_sign_large);
    man_t tmp_m = (Is_inf) ? ac_int<M,false>(0)            : ((Is_R) ? R_magn_result : N_chosen_sub.template slc<M>(2));
    exp_t tmp_e = (Is_inf) ? ac_int<E,false>((1 << E) - 1) : ((Is_R) ? R_exp_result  : N_exp_large);


    bool zero_res = (Is_R) ? Ris_zero : !(N_chosen_sub.or_reduce());

    sgn_t tmp_s_A = (in_1_subnorm) ? a.sign     : (sgn_t)0;
    man_t tmp_m_A = (in_1_subnorm) ? a.mantissa : (man_t)0;
    exp_t tmp_e_A = (in_1_subnorm) ? a.exponent : (exp_t)0;

    sgn_t tmp_s_B = (in_2_subnorm) ? sign     : (sgn_t)0;
    man_t tmp_m_B = (in_2_subnorm) ? mantissa : (man_t)0;
    exp_t tmp_e_B = (in_2_subnorm) ? exponent : (exp_t)0;
 
    output.sign     = (DENORMALS) ? tmp_s : ((in_1_subnorm) ? tmp_s_A : ((in_2_subnorm) ? tmp_s_B : ((zero_res) ? (sgn_t)0 : tmp_s)));
    output.exponent = (DENORMALS) ? tmp_e : ((in_1_subnorm) ? tmp_e_A : ((in_2_subnorm) ? tmp_e_B : ((zero_res) ? (exp_t)0 : tmp_e)));
    output.mantissa = (DENORMALS) ? tmp_m : ((in_1_subnorm) ? tmp_m_A : ((in_2_subnorm) ? tmp_m_B : tmp_m));
  }
  

  template<RND_ENUM RND_MODE=EVEN, bool DENORMALS=false>
  void fpma_dual(const fast_float<M,E> a, const fast_float<M,E> b, fast_float<M,E> &output) {
                        
    const bool t_in_a_is_inf = (exponent).and_reduce();
    const bool t_in_c_is_inf = (a.exponent).and_reduce();
    const bool t_in_b_is_inf = (b.exponent).and_reduce();

    // Start by calculating the multiplication
    
    //Infinite check
    const bool mul_is_inf = t_in_a_is_inf | t_in_c_is_inf;		 	
        
    // Find mul sign.
    const ac_int<1,false> mul_sign = sign ^ a.sign;
    
    // Create significants of mul.
    ac_int<M+1,false> mul_signif_a = mantissa;
    ac_int<M+1,false> mul_signif_b = a.mantissa;

    const ac_int<1,false> a_denorm = (exponent == 0);
    const ac_int<1,false> b_denorm = (a.exponent == 0);
    const ac_int<1,false> c_denorm = (b.exponent == 0);

    const ac_int<1,false> a_and_b_denorm = a_denorm & b_denorm;
    const ac_int<1,false> a_or_b_denorm = a_denorm ^ b_denorm;
    
    mul_signif_a[M] = (DENORMALS) ? ((a_denorm) ? 0 : 1) : 1;
    mul_signif_b[M] = (DENORMALS) ? ((b_denorm) ? 0 : 1) : 1;

    if (DENORMALS) {
      mul_signif_a <<= a_denorm;
      mul_signif_b <<= b_denorm;
    }
    
    // Calculate exponent result 
    /*static*/ const ac_int<E,false> mul_exp_base = (1 << (E-1)) - 1;
    /*static*/ const ac_int<E,false> mul_exp_base_min_1 = (1 << (E-1)) - 2;
    const ac_int<E+2,true> mul_exp_overf = exponent + a.exponent - mul_exp_base_min_1;
    const ac_int<E+2,true> mul_exp = exponent + a.exponent - mul_exp_base;

    // Do the multiplication.
    /*static*/ const int mul_input_W = M + 1;
    /*static*/ const int mul_W = (mul_input_W) << 1;
    
    ac_int<mul_W,false> mul_prod1, mul_prod;
    
    mul_prod1 = mul_signif_a*mul_signif_b;

    mul_prod = (DENORMALS) ? mul_prod1 : ((a_denorm || b_denorm) ? (ac_int<mul_W,false>)0 : mul_prod1);
    
    ac_int<1,false> mul_overf = mul_prod[mul_W-1];
    mul_prod <<= !mul_overf;
    ac_int<E+3,true> mul_exp_result1 = (mul_overf) ? mul_exp_overf : mul_exp;
    ac_int<E+3,true> mul_exp_result  = (DENORMALS) ? mul_exp_result1 : ((a_denorm || b_denorm) ? (ac_int<E+3,true>)0 : mul_exp_result1);
    
    // ONLY CHANGE FIX IT
    // ONLY CHANGE FIX IT
    // ONLY CHANGE FIX IT
    // ONLY CHANGE FIX IT
    // ONLY CHANGE FIX IT
    bool mul_zero = (DENORMALS) ? false : (a_denorm || b_denorm);
    int mul_prod_lz_cnt = 0;

    if (DENORMALS) {
      // int mul_prod_lz_cnt = (mul_prod).leading_sign(mul_zero);
      

      // 2*M+2 = 2*(23 + 1) = 48 
      // 2*M+2 = 2*(7  + 1) = 16 
      // 2*M+2 = 2*(3  + 1) = 8 
      const int ww = (2*M+2 <= 8)  ? 8  :
                     (2*M+2 <= 16) ? 16 :
                     (2*M+2 <= 32) ? 32 : 64;
      const int uu = ww - mul_W;
      ac_int<ww,false> tcz = 0;
      #pragma hls_unroll
      for (int ii=0; ii<mul_W; ii++)
        tcz[ii+uu]=mul_prod[ii];
      int mul_prod_lz_cnt = lzcount<ww>(tcz);//.leading_sign(mul_zero);

         
      if (mul_zero) {
        mul_exp_result = 0;  
      } else {
        mul_prod <<= mul_prod_lz_cnt;
        mul_exp_result -= mul_prod_lz_cnt; 
      } 
    }
    // ONLY CHANGE FIX IT
    // ONLY CHANGE FIX IT
    // ONLY CHANGE FIX IT
    // ONLY CHANGE FIX IT
    // ONLY CHANGE FIX IT


    // Perform the addition
    
    //Infinite check
    bool add_is_inf = mul_is_inf | t_in_b_is_inf;
    
    // Calculate add effective operation.
    const ac_int<1,false> add_eff_sign = mul_sign ^ b.sign;
    
    // Create add significant.
    ac_int<mul_W,false> add_signif_c =  (DENORMALS) ? b.mantissa : ((c_denorm) ? (man_t)0 : b.mantissa);
    add_signif_c[M] = (DENORMALS) ? ((c_denorm) ? 0 : 1) : 1;
    
    if (DENORMALS) add_signif_c <<= c_denorm;

    add_signif_c <<= M + 1;


    // Start of Far Path:

    // Compute exponent difference.

    const ac_int<E + 3,true> Far_exp_diff_1 = mul_exp_result - b.exponent;
    const ac_int<E + 3,true> Far_exp_diff_2 = b.exponent - mul_exp_result;

    ac_int<E+2,false> Far_abs_exp_diff;
      
    // Initialize swap variables.

    ac_int<mul_W+1,true> Far_nsignif_small;
    ac_int<mul_W,false> Far_nsignif_large;
    ac_int<E+1,false> Far_exp_large;
    ac_int<1,false> Far_sign_large;
    ac_int<1,false> Far_large_is_mul;
    
    const int Far_align_extra_W = M + 3;

    ac_int<mul_W + Far_align_extra_W,false> Far_small_signif_aligner;
    
    // Swap and allign.

    if (!Far_exp_diff_1[E+2]) {

        Far_nsignif_large = mul_prod;
        Far_exp_large = mul_exp_result;
        Far_sign_large = mul_sign;
        Far_nsignif_small = add_signif_c;
        Far_abs_exp_diff = Far_exp_diff_1.template slc<E+1>(0);
        Far_large_is_mul = 1;

    } else if (!Far_exp_diff_2[E+2]) {

        Far_nsignif_large = add_signif_c;
        Far_exp_large = b.exponent;
        Far_sign_large = b.sign;
        Far_nsignif_small = mul_prod;
        Far_abs_exp_diff = Far_exp_diff_2.template slc<E+1>(0);
        Far_large_is_mul = 0;

    }
    
      
    Far_small_signif_aligner = Far_nsignif_small;
    
    
    const int Far_align_max_sft = (Far_large_is_mul) ? mul_W : Far_align_extra_W;
    
    Far_small_signif_aligner <<= Far_align_extra_W;
    
    if (Far_abs_exp_diff < Far_align_max_sft) {
      Far_small_signif_aligner >>= Far_abs_exp_diff;
    } else { 
      Far_small_signif_aligner >>= Far_align_max_sft; 
    }
    

    Far_nsignif_small = Far_small_signif_aligner.template slc<mul_W>(Far_align_extra_W);
    
    
    const ac_int<1,false> Far_sticky = (Far_small_signif_aligner.template slc<Far_align_extra_W>(0)).or_reduce();
    const ac_int<1,false> Far_add_1 = (add_eff_sign) ? (!Far_sticky) : 0;
    
    const ac_int<E+1,false> Far_exp_large_overf = Far_exp_large + 1;
    const ac_int<E+1,false> Far_exp_large_underf = Far_exp_large - 1;
    ac_int<E+1,false> Far_exp_result;
    
    
    const ac_int<mul_W+1,true> Far_nsignif_small_compl2 = (add_eff_sign) ? (~Far_nsignif_small) : Far_nsignif_small;
    
    const ac_int<mul_W+2,true> Far_sum_1 = Far_nsignif_large + Far_nsignif_small_compl2 + Far_add_1;
    const ac_int<mul_W+2,true> Far_sum_2 = Far_nsignif_small - Far_nsignif_large;
    ac_int<1,false> Far_sign_result;
    
    ac_int<mul_W+1,false> Far_sum;

    if (!Far_sum_2[mul_W+1]) {
      Far_sign_result = !Far_sign_large;
      Far_sum = Far_sum_2;
    }  
    if (!Far_sum_1[mul_W+1]) {
      Far_sign_result = Far_sign_large;
      Far_sum = Far_sum_1;        
    }
    
    
    const ac_int<1,false> Far_sticky_hold = Far_sum[0];
    
    
    if (Far_sum[mul_W]) {
      
      Far_sum >>= 1;
      Far_sum[0] = Far_sum[0] | Far_sticky_hold;
      Far_exp_result = Far_exp_large_overf;

    } else if (!Far_sum[mul_W-1]) {
      
      Far_sum <<= 1;
      Far_exp_result = Far_exp_large_underf;	

    } else {
      
      Far_exp_result = Far_exp_large;

    }

    
    // Start of Near Path:
    
    // 2-bit exponent difference.
    
    const ac_int<2,false> mul_exp_result_slc = mul_exp_result.template slc<2>(0);
    
    const ac_int<3,true> Near_exp_diff = mul_exp_result_slc - b.exponent.template slc<2>(0);

    // Swap and 1-bit allign.
    ac_int<mul_W+1,true> Near_nsignif_large = (Near_exp_diff[1]) ? add_signif_c : mul_prod;
    ac_int<E+1,false>    Near_exp_large     = (Near_exp_diff[1]) ? (ac_int<E+1,false>)b.exponent : (ac_int<E+1,false>)mul_exp_result;
    ac_int<1,false>      Near_sign_large    = (Near_exp_diff[1]) ? b.sign : mul_sign;
    ac_int<mul_W+1,true> Near_nsignif_small = (Near_exp_diff[1]) ? mul_prod : add_signif_c;

    
    ac_int<1,false> Near_sticky_hold = 0;
    

    if (Near_exp_diff.or_reduce()) {
      Near_sticky_hold = Near_nsignif_small[0];
      Near_nsignif_small >>= 1;
    }
    
    const ac_int<1,false> Near_add_1 = !Near_sticky_hold;
    
    
    // Significant substraction.

    const ac_int<mul_W+1,true> Near_large_sub_small = Near_nsignif_large + (~Near_nsignif_small) + Near_add_1;
    const ac_int<mul_W+1,true> Near_small_sub_large = Near_nsignif_small - Near_nsignif_large;
    ac_int<E+1,false> Near_exp_result;


    ac_int<mul_W,false> Near_chosen_sub;
    ac_int<1,false> Near_sign_result;
    
    if (!Near_large_sub_small[mul_W]) {	
      Near_sign_result = Near_sign_large;
      Near_chosen_sub = Near_large_sub_small;   
    } else if (!Near_small_sub_large[mul_W]) {	
      Near_sign_result = !Near_sign_large;
      Near_chosen_sub = Near_small_sub_large;  
    }
    
    
    bool Near_zero = false;
    //const int Near_lz_cnt = (Near_chosen_sub).leading_sign(Near_zero);
    

    // 2*M+2 = 2*(23 + 1) = 48 
    // 2*M+2 = 2*(7  + 1) = 16 
    // 2*M+2 = 2*(3  + 1) = 8 
    const int ww = (2*M+2 <= 8)  ? 8  :
                   (2*M+2 <= 16) ? 16 :
                   (2*M+2 <= 32) ? 32 : 64;
    const int uu = ww - mul_W;
    ac_int<ww,false> tcz = 0;
    #pragma hls_unroll
    for (int ii=0; ii<mul_W; ii++)
      tcz[ii+uu]=Near_chosen_sub[ii];
    const int Near_lz_cnt = lzcount<ww>(tcz);
    
    

    if (Near_zero) {
      Near_exp_result = 0;
    } else {
      Near_chosen_sub <<= Near_lz_cnt;
      Near_exp_result = Near_exp_large - Near_lz_cnt;	
    }

    const bool Is_Near = (Far_abs_exp_diff < 2) & add_eff_sign;

    ac_int<mul_W,false> Chosen_signif_result = (Is_Near) ? (ac_int<mul_W+1, false>)Near_chosen_sub : Far_sum;
    ac_int<E+1,false> Chosen_exp_result = (Is_Near) ? Near_exp_result : Far_exp_result;
    ac_int<1,false> Chosen_sign_result = (Is_Near) ? Near_sign_result : Far_sign_result;
    ac_int<1,false> Chosen_sticky = (Is_Near) ? Near_sticky_hold : Far_sticky;
    
    // Round.
    
    /*static*/ const int inj_point = M + 1;
    
    ac_int<1,false> round;
    switch (RND_MODE) {
      
      case EVEN:
        // Round to nearest tie to even.
          round = Chosen_signif_result[inj_point-1] & (Chosen_signif_result[inj_point] | Chosen_signif_result[inj_point-2] | Chosen_sticky);
          break;
      case ODD:
          // Round to nearest tie to odd.
          round = Chosen_signif_result[inj_point-1] & (Chosen_signif_result[inj_point] | (!Chosen_signif_result[inj_point-2]) | Chosen_sticky);
          break;
      case INF:
          // Round infinity.
          round = Chosen_signif_result[inj_point-1];
          break;

      }
      

    ac_int<M+1,false> Final_signif = Chosen_signif_result.template slc<M+1>(inj_point);
    
    bool zero_res = !(Final_signif.or_reduce() | round);
    Chosen_sign_result = (zero_res) ? mul_sign : Chosen_sign_result;
              
    ac_int<M+2,false> Final_signif_result = Final_signif + round;
    
    ac_int<E+1,false> Chosen_exp_result_overf = Chosen_exp_result + 1;
    ac_int<E+1,false> Final_exp_result;
    

    if (Final_signif_result[M+1]) {
      Final_signif_result >>= 1;
      Final_exp_result = Chosen_exp_result_overf;
    } else { 
      Final_exp_result = Chosen_exp_result; 
    }

    add_is_inf = (Final_exp_result >= ((1 << E) - 1)) ? true : add_is_inf;
    
    const ac_int<1,false> inf_sign = (t_in_b_is_inf) ? b.sign : mul_sign;

    
    output.sign = (add_is_inf) ? inf_sign :((zero_res) ? (sgn_t)0 : Chosen_sign_result);
    output.exponent = (add_is_inf) ? (exp_t)((1 << E) - 1) : ((zero_res) ? (exp_t)0 : Final_exp_result.template slc<E>(0));
    output.mantissa = (add_is_inf) ? (man_t)0 : Final_signif_result.template slc<M>(0);
  }

 
  template<RND_ENUM RND_MODE=EVEN, bool DENORMALS=false>
  void fpmul(const fast_float<M,E> a, fast_float<M,E> &output) {

    static const ac_int<E,false> mul_exp_base = (1 << (E-1)) - 1;
    static const int signif_W = M + 1;
    static const int mul_W = (signif_W) << 1;

    ac_int<1,false> a_is_zero, c_is_zero;

    // Sign Calculation
    const ac_int<1,false> sign_res = sign ^ a.sign;

    // Create significants of mul.
    ac_int<M+1,false> mul_signif_a = mantissa;
    ac_int<M+1,false> mul_signif_c = a.mantissa;
    
    ac_int<1,false> is_a_subnormal = (exponent == 0);
    ac_int<1,false> is_c_subnormal = (a.exponent == 0);
    
    mul_signif_a[M] = (DENORMALS) ? ((is_a_subnormal) ? 0 : 1) : 1;
    mul_signif_c[M] = (DENORMALS) ? ((is_c_subnormal) ? 0 : 1) : 1;
    
    if (DENORMALS) {
      mul_signif_a <<= is_a_subnormal;
      mul_signif_c <<= is_c_subnormal;

      a_is_zero = ~(mul_signif_a.or_reduce());
      c_is_zero = ~(mul_signif_c.or_reduce());
    }
    
    // Infinity Check.
    
    bool Is_inf = ((exponent).and_reduce() | (a.exponent).and_reduce());
    
    // Calculate exponent result and overf exponent result.
    ac_int<E+2,true> mul_exp =  (DENORMALS) ? 
                                ((a_is_zero | c_is_zero) ? (ac_int<E+2,true>)0 : (ac_int<E+2,true>)(exponent + a.exponent - mul_exp_base)) :
                                ((is_a_subnormal | is_c_subnormal)) ? (ac_int<E+2,true>)0 : (ac_int<E+2,true>)(exponent + a.exponent - mul_exp_base);
    if (mul_exp < 0) 
       mul_exp = 0;
    else if (mul_exp >= ((1<<E)-1)) 
       mul_exp = ((1<<E)-1);
    
    // Overflow Exponent Calculation
    ac_int<E+2,true> mul_exp_overf = (DENORMALS) ? 
                                     ((a_is_zero | c_is_zero) ? (ac_int<E+2,true>)0 : (ac_int<E+2,true>)(exponent + a.exponent - mul_exp_base + 1)) :
                                     ((is_a_subnormal | is_c_subnormal)) ? (ac_int<E+2,true>)0 : (ac_int<E+2,true>)(exponent + a.exponent - mul_exp_base + 1);
    if (mul_exp_overf < 0) 
       mul_exp_overf = 0;
    else if (mul_exp_overf >= ((1<<E)-1)) 
       mul_exp_overf = ((1<<E)-1);

    ac_int<E+2,true> right_shift;
    if (DENORMALS) {
      right_shift = (a_is_zero | c_is_zero) ? (ac_int<E+2,true>)0 : (ac_int<E+2,true>)(mul_exp_base - exponent - a.exponent);
      if (right_shift < 0) 
        right_shift = 0;
      else if (right_shift > M+2) 
        right_shift = M+2;

    }

    ac_int<1,false> mul_zero = (DENORMALS) ? (a_is_zero | c_is_zero) : (is_a_subnormal | is_c_subnormal);
    
    // Do the multiplication.
    ac_int<mul_W,false> mul_res = mul_signif_a*mul_signif_c;
    ac_int<mul_W,false> mul_prod = (DENORMALS) ? ((a_is_zero | c_is_zero) ? (ac_int<mul_W,false>)0 : mul_res) : ((is_a_subnormal | is_c_subnormal) ? (ac_int<mul_W,false>)0 : mul_res);
    

    // Overflow Check.
    ac_int<1,false> mul_overflow = mul_prod[mul_W-1];
    ac_int<E+2,true> exp_res = (mul_overflow) ? mul_exp_overf : mul_exp;
  
    // Normalization.
    if (DENORMALS) {
      // const ac_int<ac::nbits<mul_W>::val,false> mul_lzc = mul_prod.leading_sign();
      

      // 2*M+2 = 2*(23 + 1) = 48 
      // 2*M+2 = 2*(7  + 1) = 16 
      // 2*M+2 = 2*(3  + 1) = 8 
      const int ww = (2*M+2 <= 8)  ? 8  :
                     (2*M+2 <= 16) ? 16 :
                     (2*M+2 <= 32) ? 32 : 64;
      const int uu = ww - mul_W;
      ac_int<ww,false> tcz = 0;
      #pragma hls_unroll
      for (int ii=0; ii<mul_W; ii++)
        tcz[ii+uu]=mul_prod[ii];
      const ac_int<ac::nbits<mul_W>::val,false> mul_lzc = lzcount<ww>(tcz);

      if (exp_res==0) {     
        mul_prod >>= right_shift;     
      } else {
        //mul_prod <<= mul_lzc;
        if (mul_lzc >= exp_res) {	 
          mul_prod <<= exp_res;
	        exp_res = 0;  
      	} else {
          mul_prod <<= (mul_lzc <= 1) ? (ac_int<ac::nbits<mul_W>::val,false>)0 : (ac_int<ac::nbits<mul_W>::val,false>)(mul_lzc-1);
	        exp_res -= (mul_lzc <= 1) ? (ac_int<ac::nbits<mul_W>::val,false>)0 : (ac_int<ac::nbits<mul_W>::val,false>)(mul_lzc-1);	
        }
      }
      
      // Subnormal Result Check.
      mul_prod >>= ((exp_res == 0) && (!mul_overflow));
    }
    
   
    
    // Injection Rounding.
    const int inj_point = M;
    
    ac_int<1,false> mul_round;
    switch (RND_MODE) {
      
      case EVEN:
        // Round to nearest tie to even.
          mul_round = (mul_overflow) ? mul_prod[inj_point] & (mul_prod[inj_point+1] | mul_prod[inj_point-1] | (mul_prod.template slc<inj_point-1>(0)).or_reduce()) :
                                       mul_prod[inj_point-1] & (mul_prod[inj_point] | mul_prod[inj_point-2] | (mul_prod.template slc<inj_point-2>(0)).or_reduce());
          break;
      case ODD:
          // Round to nearest tie to odd.
          mul_round = (mul_overflow) ? mul_prod[inj_point] & (mul_prod[inj_point+1] | !mul_prod[inj_point-1] | (mul_prod.template slc<inj_point-1>(0)).or_reduce()) :
                                       mul_prod[inj_point-1] & (mul_prod[inj_point] | !mul_prod[inj_point-2] | (mul_prod.template slc<inj_point-2>(0)).or_reduce());
          break;
      case INF:
          // Round infinity.
          mul_round = (mul_overflow) ? mul_prod[inj_point] : mul_prod[inj_point-1];
          break;

      }

    // TODO: check if moving addition after saves area              
    ac_int<M+2,false> mul_rounded = (mul_overflow) ? mul_prod.template slc<M+1>(inj_point+1) + mul_round : mul_prod.template slc<M+1>(inj_point) + mul_round;

    ac_int<1,false> extra_shift = 0;
    ac_int<1,false> extra_exp = 0;
    // Rounding Overflow Check.
    if (DENORMALS) {
      if ((exp_res == 0) && mul_rounded[M]) {
        extra_shift = 0;
        extra_exp = 1;
      }
    } else {
      if (mul_rounded[M+1]) {
        extra_shift = 1;
        extra_exp = 1;
      }
    } 
    
    exp_t tmp_exp = exp_res.template slc<E>(0);
    exp_t mul_exp_result = (tmp_exp == (1<<E)-1) ? tmp_exp : (exp_t)(tmp_exp+extra_exp);
    
    // Infinity Result Check.
    Is_inf = mul_exp_result.and_reduce() |  Is_inf; 

    // Return result.
    output.sign = sign_res;
    output.exponent = (Is_inf) ? (exp_t)((1 << E) - 1) : mul_exp_result;
    output.mantissa = (Is_inf) ? (man_t)0 : mul_rounded.template slc<M>(extra_shift);
  }
  
  // TODO : Support for denormals
  template<int N, RND_ENUM RND_MODE=EVEN, bool DENORMALS=false>
  void dotProd(const fast_float<M,E> A[N], const fast_float<M,E> B[N]) {
    
    static const ac_int<E,false> mul_exp_base = e_bias;
    static const ac_int<E,false> mul_exp_base_min_1 = e_bias-1;
    static const int mul_input_W = M + 1;
    static const int mul_W = (mul_input_W) << 1;

    bool AisInf[N], BisInf[N], AisZero[N], BisZero[N];
    
    sgn_t infSign;
    sgn_t mulSign[N];

    ac_int<M+1,false> fracA[N], fracB[N];

    
    ac_int<N,false> mulisInf;
    ac_int<E+2,true> mul_exp_overf[N]; 
    ac_int<E+2,true> mul_exp[N];
    ac_int<E+2,true> mul_exp_result[N];
    
    #pragma hls_unroll
    INP_DEC: for (int i=0; i< N; i++){
      // Check for infinite values
      AisInf[i] = A[i].exponent.and_reduce();
      BisInf[i] = B[i].exponent.and_reduce();
      
      // Check for zero values (DENORMALS==0)
      AisZero[i] = A[i].exponent == 0;
      BisZero[i] = B[i].exponent == 0;

      mulSign[i] = A[i].sign ^ B[i].sign;

      mulisInf[i] = AisInf[i] | BisInf[i];
      infSign =mulSign[i]; // TODO : check sign for infinity

      fracA[i]    = A[i].mantissa;
      fracA[i][M] = 1;

      fracB[i]    = B[i].mantissa;
      fracB[i][M] = 1;
    
      // Calculate exponent result 

      
      mul_exp[i] = (AisZero[i] || BisZero[i]) ? (ac_int<E+2,true> )0 :  A[i].exponent + B[i].exponent - mul_exp_base;
      if (mul_exp[i]<0) 
        mul_exp[i] = 0;
      else if (mul_exp[i] >= ((1<<E)-1)) 
        mul_exp[i] = ((1<<E)-1);
      
      mul_exp_overf[i] = (AisZero[i] || BisZero[i]) ? (ac_int<E+2,true> )0 : A[i].exponent + B[i].exponent - mul_exp_base_min_1;
      if (mul_exp_overf[i]<0) 
        mul_exp_overf[i] = 0;
      else if (mul_exp_overf[i] >= ((1<<E)-1)) 
        mul_exp_overf[i] = ((1<<E)-1);
    }
    
    
    // Do the multiplication.
    ac_int<mul_W,false> mul_prod[N];
    ac_int<1,false> mul_overf[N];

    #pragma hls_unroll
    DO_MUL: for (int i=0; i<N; i++) {
      
      mul_prod[i] = fracA[i]*fracB[i];

      if (AisZero[i] || BisZero[i]) {
        mul_prod[i] = 0;
        mul_overf[i] =0; 
      } else {
        mul_overf[i] = mul_prod[i][mul_W-1];
        mul_prod[i] <<= !mul_overf[i];
      }
        
      mul_exp_result[i] = (mul_overf[i]) ? mul_exp_overf[i] : mul_exp[i];
      
      if (mul_exp_result[i] >= ((1<<E)-1)) {
        mulisInf[i] = 1;
      }
    } // DO_MUL


    // Start the addition

    ac_int<1,false> addisInf = mulisInf.or_reduce();  // check infinity

    ac_int<E+2,true> maxExp = max<N>(mul_exp_result); // find max exponent
    
    ac_int<E+2,true> shifted_exp_res[N];
    ac_int<mul_W,false> shifted_prod[N];
    ac_int<E+1,true> diffE[N];
    #pragma hls_unroll
    ALLIGN: for (int i=0; i<N; i++) {
      diffE[i] = maxExp-mul_exp_result[i];
      shifted_prod[i] = mul_prod[i] >> diffE[i];
      shifted_exp_res[i] = mul_exp_result[i] + diffE[i];
    }
        
    ac_int<mul_W+1,true> add_op[N];
    ac_int<mul_W+1+N,true> acc=0;

    #pragma hls_unroll
    TWOs_COMPL_ADD: for (int i=0; i<N; i++) {
      if (mulSign[i]) {
        add_op[i] = -shifted_prod[i];
      } else {
        add_op[i] = shifted_prod[i];
      }   
      acc += add_op[i];
    }

    
    ac_int<1,false> res_sign = acc[mul_W+N];
    
    ac_int<mul_W+1+N,false> res1 = acc;
    ac_int<mul_W+1+N,false> res2 = acc.bit_complement();
    
    ac_int<mul_W+N,false> res;
    res = (res_sign) ? res2.template slc<mul_W+N>(0) : res1.template slc<mul_W+N>(0);
    
    // const int lead_zer = res.leading_sign();
    // mul_W+N = (23 + 1)*2 + 4 = 52 (or 56 for dot<8>)
    // mul_W+N = (7  + 1)*2 + 4 = 20 (or 24 for dot<8>)
    // mul_W+N = (3  + 1)*2 + 4 = 12 (or 16 for dot<8>)
    const int ww = (2*M+2+N <= 8)  ? 8  :
                   (2*M+2+N <= 16) ? 16 :
                   (2*M+2+N <= 32) ? 32 : 64;
    const int uu = ww - mul_W - N;
    ac_int<ww,false> tcz = 0;
    #pragma hls_unroll
    for (int ii=0; ii<mul_W+N; ii++)
      tcz[ii+uu]=res[ii];
    const int lead_zer = lzcount<ww>(tcz);

    res <<= lead_zer;
    maxExp -= (lead_zer);
    maxExp += (lead_zer==0) ? N-1 : N;
    
    // res width = 2M+2+N, we need to keep 1+M bits, so we start 
    //  from the M+N+2-1 index and keep M+1 
    static const int ind = M+N+1;
    ac_int<M+2,false> rounded_res = res.template slc<M+1>(ind);
    
    ac_int<1,false> mul_round;
    switch (RND_MODE) {
      
      case EVEN:
        // Round to nearest tie to even.
          mul_round = res[ind-1] & ((res[ind] | res[ind-2] | (res.template slc<ind-2>(0)).or_reduce()));
          break;
      case ODD:
          // Round to nearest tie to odd.
          mul_round = res[ind-1] & ((res[ind] | res[ind-2] | ~(res.template slc<ind-2>(0)).or_reduce()));
          break;
      case INF:
          // Round infinity.
          mul_round = res[ind-1];
          break;

      }                  
    rounded_res += (mul_round+res_sign);

    exp_t tmp_exp = maxExp.template slc<E>(0);
    exp_t over_exp = (tmp_exp == (1<<E)-1) ? tmp_exp : (exp_t)(tmp_exp+1);

    mantissa = (addisInf) ? (man_t)(0) : (rounded_res[M+1]) ? rounded_res.template slc<M>(1) : rounded_res.template slc<M>(0); 
    exponent = (addisInf) ? (exp_t)((1<<E) -1) : (rounded_res[M+1]) ?  over_exp : tmp_exp;
    sign = (addisInf) ? infSign : res_sign;

  }

  
  // OVERLOAD OPERATORS
  void operator = (const float &inFP) {
    
    #ifndef __SYNTHESIS__
    bool iszer = (inFP == 0) ? true : false;

    unsigned int_part = (inFP < 0) ? (-inFP) : inFP;
    float    fra_part = (inFP < 0) ? (-inFP - int_part) : (inFP - int_part);

    ac_int<32,false>     int_bin = int_part;
    ac_fixed<32,0,false> fra_bin = fra_part;

    int right_shift = -1;
    if (int_bin > 0) {
      for (int i = 0; i < 32; i++) {
        if (int_bin[i] == 1) {
          right_shift = i;
        }
      }
    }
    int left_shift = -1;
    for (int i = 0; i < 32; i++) {
      if (fra_bin[i] == 1) {
        left_shift = 32- i;
      }
    }

    ac_int<64,false> bin;

    for(int i = 0; i < 32; i++) {
      bin[i] = fra_bin[i];
    }
    for(int i = 0; i < 32; i++) {
      bin[32+i] = int_bin[i];
    }

    int e = 0;
    if (right_shift >= 0) {
      bin = bin >> right_shift;
      e = right_shift;
    } else if (left_shift >= 0) {
      bin = bin << left_shift;
      e = -left_shift;
    }

    e += e_bias;
    iszer = (iszer || (e<0)) ? true : false;
    
    int maxE = ((1 << (E)) - 1);
    bool isinf = (e >= maxE) ? true : false;
    
    sign     = (iszer) ? (sgn_t)0 : ((inFP < 0) ? (sgn_t)1                : (sgn_t)0);
    exponent = (iszer) ? (exp_t)0 : ((isinf)    ? (exp_t)((1 << (E)) - 1) : (exp_t)e);
    mantissa = (iszer) ? (man_t)0 : ((isinf)    ? (man_t)0                : (man_t)(bin.template slc<M>(32-M)));
    
    #endif
  }
  
  void operator = (const ac_int<M+E+1,false> &in) {
    sign = in[M+E];
    exponent = in.template slc<E>(M);
    mantissa = in.template slc<M>(0);
  }
  
  void operator = (const fast_float<M,E> &in) {
    mantissa = in.mantissa;
    exponent = in.exponent;
    sign = in.sign;
  }

  fast_float<M, E> operator + (const fast_float<M, E> &b) {
    fast_float<M, E> r;
    fpa_dual(b,r);
    return r;
  }
  
  fast_float<M, E> operator - (const fast_float<M, E> &b) {
    fast_float<M, E> r, tmp;
    tmp = b;
    tmp.sign = ~tmp.sign;
    fpa_dual(tmp,r);
    return r;
  }
  
  fast_float<M, E> operator * (const fast_float<M, E> &b) {
      fast_float<M, E> r;
      fpmul(b,r);
      return r;
  }

  fast_float<M, E> &operator += (const fast_float<M, E> &b) {
    *this = this->operator+(b);
    return *this;
  }

  fast_float<M, E> &operator -= (const fast_float<M, E> &b) {
    *this = this->operator-(b);
    return *this;
  }

  fast_float<M, E> &operator *= (const fast_float<M, E> &b) {
    *this = this->operator*(b);
    return *this;
  }
  
  bool operator == (const fast_float<M,E> b) {
    return ((sign == b.sign) && (mantissa == b.mantissa) && (exponent == b.exponent));
  }

  bool operator != (const fast_float<M,E> b) {
    return ((sign ^ b.sign) || (mantissa != b.mantissa) || (exponent != b.exponent));
  }

  bool operator >  (const fast_float<M,E> b) {
    bool greaterA = (exponent>b.exponent || (exponent==b.exponent && mantissa>b.mantissa));
    bool greaterB = (exponent<b.exponent || (exponent==b.exponent && mantissa<b.mantissa));
    bool greater  = (sign == 0) ? greaterA : greaterB;
    return (sign < b.sign || (sign==b.sign && greater));
  }

  bool operator <  (const fast_float<M,E> b) {
    bool lessA = (exponent<b.exponent || (exponent==b.exponent && mantissa<b.mantissa));
    bool lessB = (exponent>b.exponent || (exponent==b.exponent && mantissa>b.mantissa));
    bool less  = (sign == 0) ? lessA : lessB;
    return (sign > b.sign || (sign==b.sign && less));
  }

  bool operator >= (const fast_float<M,E> b) {
    bool grOReqA = exponent>b.exponent || (exponent==b.exponent && mantissa>=b.mantissa);
    bool grOReqB = exponent<b.exponent || (exponent==b.exponent && mantissa<=b.mantissa);
    bool grOReq  = (sign==0) ? grOReqA : grOReqB;
    return (sign<b.sign || (sign==b.sign && grOReq));
  }

  bool operator <= (const fast_float<M,E> b) {
    bool lsOReqA = exponent<b.exponent || (exponent==b.exponent && mantissa<=b.mantissa);
    bool lsOReqB = exponent>b.exponent || (exponent==b.exponent && mantissa>=b.mantissa);
    bool lsOReq  = (sign==0) ? lsOReqA : lsOReqB;
    return (sign>b.sign || (sign==b.sign && lsOReq));
  }

};

// double-precission IEEE754 
typedef fast_float<52,11> ffp64;  
// single-precission IEEE754 
typedef fast_float<23,8>  ffp32;  
// half-precission IEEE754 
typedef fast_float<10,5>  ffp16;  
// bfloat16
typedef fast_float<7, 8>  ffp16b; 

#endif
