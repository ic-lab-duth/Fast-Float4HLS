#ifndef __FASTFLOAT_H__
#define __FASTFLOAT_H__

#include <ac_int.h>
#include <ac_std_float.h>

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

private:
  template<int N=M, bool S=false>
  void full_adder(ac_int<2*(M+1),false> a,
                  ac_int<2*(M+1),false> b,
                  ac_int<2*(M+1),false> c, 
                  ac_int<2*(M+1),false> &sum_bits,
                  ac_int<2*(M+1),false> &carry_bits) {
    sum_bits= a^b^c;
    carry_bits= ((a & b)|(a & c)|(b & c))<<1;
  }

  template<int N=M, bool S=false>
  void multiply(ac_int<N,S> operand_1, ac_int<N,S> operand_2, ac_int<2*N,S> &result) {
    
    ac_int<N/2,false> fh_operand_1, sh_operand_1;
    ac_int<N/2,false> fh_operand_2, sh_operand_2;
    ac_int<N,false> p1, p2, p3, p4;

    ac_int<2*N,false> partial_product_1=0;
    ac_int<2*N,false> partial_product_2=0;
    ac_int<2*N,false> partial_product_3=0;
    ac_int<2*N,false> partial_product_4=0;
    ac_int<2*N,false> sum_first=0;
    ac_int<2*N,false> carry_first=0;
    ac_int<2*N,false> sum=0;
    ac_int<2*N,false> carry=0;
    ac_int<2*N,false> result_mid;

    const int half=N/2;
    ac_int<N,false> op1 = operand_1;
    ac_int<N,false> op2 = operand_2;

    fh_operand_1=op1.template slc<half>(0);
    sh_operand_1=op1.template slc< half>(half);
    fh_operand_2=op2.template slc<half>(0);
    sh_operand_2=op2.template slc<half>(half);

    p1=fh_operand_1*fh_operand_2;
    p2=fh_operand_1*sh_operand_2;
    p3=sh_operand_1*fh_operand_2;
    p4=sh_operand_1*sh_operand_2;
    
    partial_product_1.set_slc(0,p1);
    partial_product_2.set_slc(N/2,p2);
    partial_product_3.set_slc(N/2,p3);
    partial_product_4.set_slc(N,p4);

    full_adder(partial_product_1,partial_product_2,partial_product_3,sum_first,carry_first);

    full_adder(partial_product_4,sum_first,carry_first,sum,carry);

    result_mid=sum+carry;
    result=result_mid;
  }

public:
  
  // Constructors
  fast_float() {};
  fast_float(const float &in) {
    this->operator=(in);
  };
  
// #ifndef __SYNTHESIS__
  /*
  * .to_float() is a non-synthesizable function that
  * is used for converting the fast_float into a c++ float
  */
  float to_float() {

    // create significand
    ac_fixed<M+2,2,true> fix;
    fix[M+1] = 0;
    fix[M] = ((mantissa== 0)&& (exponent==0)) ? 0 : 1;
    for (int i=0; i<M; i++)
      fix[i] = mantissa[i];

    // create 2^exp
    float exp= (exponent < 127) ? 0.5 : 2;
    int iter = (exponent < 127) ? e_bias - exponent : exponent - e_bias;
    float e = 1;
    for (int i=0; i< iter; i++) 
      e *= exp;

    // define sign  
    float sg = (sign == 0) ? (float)1.0 : (float)-1.0;
    
    // compute float number
    float fo = sg*((float)fix.to_double())*e;

    if (exponent==255) {
      if (mantissa ==0) {
        return (sign == 0) ? INFINITY : -INFINITY;
      } else {
        return NAN;
      }
    } else {
      return fo;
    }
  }
// #endif

  /** FLOATING POINT OPERATIONS **/

  /* Floating Point Addition
  *  adds the value of 'a' to this fast_float and
  *  return the result to the 'output'. 
  *  The addition is implemented using dual path (Far-Near).
  */
  template<RND_ENUM RND_MODE=EVEN, bool DENORMALS=false>
  void fpa_dual(fast_float<M,E> a, fast_float<M,E> &output) {
  
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
    const ac_int<ac::nbits<signif_uniform_W-1>::val, false> N_lz_cnt = (N_chosen_sub).leading_sign(N_zero);
    
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

    sgn_t tmp_s_A = (in_1_subnorm) ? a.sign     : (sgn_t)0;
    man_t tmp_m_A = (in_1_subnorm) ? a.mantissa : (man_t)0;
    exp_t tmp_e_A = (in_1_subnorm) ? a.exponent : (exp_t)0;

    sgn_t tmp_s_B = (in_2_subnorm) ? sign     : (sgn_t)0;
    man_t tmp_m_B = (in_2_subnorm) ? mantissa : (man_t)0;
    exp_t tmp_e_B = (in_2_subnorm) ? exponent : (exp_t)0;
 
    output.sign     = (DENORMALS) ? tmp_s : ((in_1_subnorm) ? tmp_s_A : ((in_2_subnorm) ? tmp_s_B : tmp_s));
    output.exponent = (DENORMALS) ? tmp_e : ((in_1_subnorm) ? tmp_e_A : ((in_2_subnorm) ? tmp_e_B : tmp_e));
    output.mantissa = (DENORMALS) ? tmp_m : ((in_1_subnorm) ? tmp_m_A : ((in_2_subnorm) ? tmp_m_B : tmp_m));
  }
  

  template<RND_ENUM RND_MODE=EVEN, bool DENORMALS=false>
  void fpma_dual(fast_float<M,E> a, fast_float<M,E> b, fast_float<M,E> &output) {
                        
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

    //const ac_int<1,false> mul_and_c_denorm = in_1_subnorm & c_denorm;
    //const ac_int<1,false> mul_or_c_denorm = in_1_subnorm ^ c_denorm;

    
    mul_signif_a[M] = (DENORMALS) ? ((a_denorm) ? 0 : 1) : 1;
    mul_signif_b[M] = (DENORMALS) ? ((b_denorm) ? 0 : 1) : 1;

    if (DENORMALS) {
      mul_signif_a <<= a_denorm;
      mul_signif_b <<= b_denorm;
    }
    
    // Calculate exponent result 
    static const ac_int<E,false> mul_exp_base = (1 << (E-1)) - 1;
    static const ac_int<E,false> mul_exp_base_min_1 = (1 << (E-1)) - 2;
    const ac_int<E+2,true> mul_exp_overf = exponent + a.exponent - mul_exp_base_min_1;
    const ac_int<E+2,true> mul_exp = exponent + a.exponent - mul_exp_base;

    // Do the multiplication.
    static const int mul_input_W = M + 1;
    static const int mul_W = (mul_input_W) << 1;
    
    ac_int<mul_W,false> mul_prod1, mul_prod;
    
    multiply(mul_signif_a, mul_signif_b, mul_prod1);

    mul_prod = (DENORMALS) ? mul_prod1 : ((a_denorm || b_denorm) ? (ac_int<mul_W,false>)0 : mul_prod1);
    
    ac_int<1,false> mul_overf = mul_prod[mul_W-1];
    mul_prod <<= !mul_overf;
    ac_int<E+3,true> mul_exp_result1 = (mul_overf) ? mul_exp_overf : mul_exp;
    ac_int<E+3,true> mul_exp_result  = (DENORMALS) ? mul_exp_result1 : ((a_denorm || b_denorm) ? (ac_int<E+3,true>)0 : mul_exp_result1);
    
    bool mul_zero = false;
    const int mul_prod_lz_cnt = (mul_prod).leading_sign(mul_zero);

    if (mul_zero) {
      mul_exp_result = 0;  
    } else {
      mul_prod <<= mul_prod_lz_cnt;
      mul_exp_result -= mul_prod_lz_cnt; 
    }


    // Perform the addition
    
    //Infinite check
    bool add_is_inf = mul_is_inf | t_in_b_is_inf;
    
    // Calculate add effective operation.
    const ac_int<1,false> add_eff_sign = mul_sign ^ b.sign;
    
    // Create add significant.
    ac_int<mul_W,false> add_signif_c = (DENORMALS) ? b.mantissa : ((c_denorm) ? (man_t)0 : b.mantissa);
    
    add_signif_c[M] = (DENORMALS) ? ((c_denorm) ? 0 : 1) : 1;
    
    if (DENORMALS) add_signif_c <<= c_denorm;

    add_signif_c <<= M + 1;


    // Start of Far Path:

    // Compute exponent difference.

    const ac_int<E + 3,true> Far_exp_diff_1 = mul_exp_result - b.exponent;
    const ac_int<E + 3,true> Far_exp_diff_2 = b.exponent - mul_exp_result;

    ac_int<E+2,false> Far_abs_exp_diff;// = (sel) ? Far_exp_diff_1.template slc<E+1>(0) : Far_exp_diff_2.template slc<E+1>(0);
      
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
    const int Near_lz_cnt = (Near_chosen_sub).leading_sign(Near_zero);
    

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
    
    static const int inj_point = M + 1;
    
    ac_int<1,false> round;
    // ac_int<1,false> round = Chosen_signif_result[inj_point-1] & 
            // (Chosen_signif_result[inj_point] | (Chosen_signif_result.template slc<inj_point-1>(0)).or_reduce() /*Chosen_signif_result[inj_point-2]*/ | Chosen_sticky);
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

    
    output.sign = (add_is_inf) ? inf_sign : Chosen_sign_result;
    output.exponent = (add_is_inf) ? (ac_int<E+1, false>)((1 << E) - 1) : Final_exp_result;
    output.mantissa = (add_is_inf) ? (ac_int<M+2, false>)0 : Final_signif_result;
  }

 
  template<RND_ENUM RND_MODE=EVEN, bool DENORMALS=false>
  void fpmul(fast_float<M,E> a, fast_float<M,E> &output, bool &invalid) {

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
    }
    
    // Infinity Check.
    
    bool Is_inf = ((exponent).and_reduce() | (a.exponent).and_reduce());
    
    // Calculate exponent result and overf exponent result.
    
    static const ac_int<E,false> mul_exp_base = (1 << (E-1)) - 1;
    ac_int<E+2,true> mul_pos_exp_result = exponent + a.exponent - mul_exp_base;
    const ac_int<E+2,true> mul_neg_exp_sft_num = mul_exp_base - exponent - a.exponent;
    
    // Do the multiplication.
    
    static const int signif_W = M + 1;
    static const int mul_W = (signif_W) << 1;
    
    ac_int<mul_W,false> mul_res;
    multiply(mul_signif_a,mul_signif_c,mul_res);

    ac_int<mul_W+1,false> mul_prod = (DENORMALS) ? mul_res : ((is_a_subnormal || is_c_subnormal) ? (ac_int<mul_W,false>)0 : mul_res);

    mul_prod <<= 1;
    
    // Overflow Check.
    
    bool mul_overflow = false;
    
    if (mul_prod[mul_W]) {
      mul_prod >>= 1;
      mul_overflow = true;
    }
    
    // Create mul_prod_normalizer witch has M more LSBs 
    // to keep shifted out sticky bits in case of subnormal result.
    
    // ac_int<mul_W+M+2,false> mul_prod_normalizer = mul_prod.template slc<mul_W>(0);
    // mul_prod_normalizer <<= M+2;
    
    ac_int<mul_W,false> mul_prod_normalizer = mul_prod.template slc<mul_W>(0);

    //TODO: MUL NEEDS LEADING ZEROS?????

    // Normalization.
    
    bool mul_zero = false;
    const bool pos_exp = (mul_pos_exp_result >= 0);
    // const ac_int<ac::nbits<mul_W>::val,false> mul_lzc = (mul_prod_normalizer.template slc<mul_W>(M+2)).leading_sign(mul_zero);
    const ac_int<ac::nbits<mul_W>::val,false> mul_lzc = mul_prod_normalizer.leading_sign(mul_zero);
    
    if (mul_zero) {
      
      mul_pos_exp_result = 0;
    
    } else if (!pos_exp) {

      mul_pos_exp_result = 0;
      mul_prod_normalizer >>= mul_neg_exp_sft_num;
      
    } else {

      if (mul_lzc >= mul_pos_exp_result) {
        
        mul_prod_normalizer <<= mul_pos_exp_result;
        mul_pos_exp_result = 0;
          
      } else {
        
        mul_prod_normalizer <<= mul_lzc;
        mul_pos_exp_result -= mul_lzc;
        
      }
    }
    
    // Subnormal Result Check.
    
    mul_prod_normalizer >>= ((mul_pos_exp_result == 0) && (!mul_overflow));
    
    // Overflow Exponent Calculation
    
    const ac_int<E+2,true> mul_pos_overf_exp_res = mul_pos_exp_result + 1;
    
    // Sign Calculation
    
    const ac_int<1,false> sign_res = sign ^ a.sign;
    
    // Injection Rounding.
    
    static const int mul_inj_point = M + 1;
    // ac_int<1,false> mul_round = mul_prod_normalizer[mul_inj_point-1] & ((mul_prod_normalizer[mul_inj_point]) 
    //               | (mul_prod_normalizer.template slc<mul_inj_point-1>(0)).or_reduce());
    
    ac_int<1,false> mul_round = mul_prod_normalizer[mul_inj_point-1] & ((mul_prod_normalizer[mul_inj_point]) 
                  | (mul_prod_normalizer.template slc<mul_inj_point-1>(0)).or_reduce());
                  
    ac_int<M+2,false> mul_signif_result = mul_prod_normalizer.template slc<M+1>(mul_inj_point) + mul_round;
    
    // Rounding Overflow Check.
    
    if (mul_signif_result[M+1]) {
      
      mul_signif_result >>= 1;
      mul_overflow = true;
      
    } else if ((mul_pos_exp_result == 0) && mul_signif_result[M]) {
      
      mul_overflow = true;
      
    }
    
    const ac_int<E+2,true> mul_exp_result = (mul_overflow && pos_exp) ? mul_pos_overf_exp_res : mul_pos_exp_result;
    
    // Infinity Result Check.
    
    Is_inf = ((mul_exp_result.template slc<E>(0)).and_reduce() | mul_exp_result[E]) ? true : Is_inf;
          
    // Return result.
    
    output.sign = sign_res;
    
    if (Is_inf) {
      
      output.exponent = (1 << E) - 1;
      output.mantissa = 0;
      
    } else {
      
      output.exponent = mul_exp_result;				
      output.mantissa = mul_signif_result.template slc<M>(0);
      
    }


  }

  template<int N>
  void dotProd(fast_float<M,E> A[N], fast_float<M,E> B[N]) {

     bool AisInf[N], BisInf[N], AisZero[N], BisZero[N];
     sgn_t mulSign[N];
    
    #pragma hls_unroll
    INP_DEC: for (int i=0; i< N; i++){
      AisInf[i] = A[i].exponent.and_reduce();
      BisInf[i] = B[i].exponent.and_reduce();

      AisZero[i] = A[i].exponent == 0;
      BisZero[i] = B[i].exponent == 0;

      mulSign[i] = A[i].sign ^ B[i].sign;
    }
    
    ac_int<1,false> infSign;
    bool mulisInf[N];
    #pragma hls_unroll
    M_INFcheck: for (int i=0; i< N; i++) {
      mulisInf[i] = AisInf[i] | BisInf[i];
      infSign =mulSign[i];
    }

    ac_int<M+1,false> fracA[N], fracB[N];
    #pragma hls_unroll
    GET_FRAC: for (int i=0; i<N; i++) {
      fracA[i]    = A[i].mantissa;
      fracA[i][M] = 1;

      fracB[i]    = B[i].mantissa;
      fracB[i][M] = 1;
    }

    // Calculate exponent result 
    static const ac_int<E,false> mul_exp_base = (1 << (E-1)) - 1;
    static const ac_int<E,false> mul_exp_base_min_1 = (1 << (E-1)) - 2;

    static const int mul_input_W = M + 1;
    static const int mul_W = (mul_input_W) << 1;
     
    
    ac_int<E+2,true> mul_exp_overf[N]; 
    ac_int<E+2,true> mul_exp[N];
    #pragma hls_unroll
    M_EXPcalc: for (int i=0; i<N; i++) {
      mul_exp_overf[i] = A[i].exponent + B[i].exponent - mul_exp_base_min_1;
      mul_exp[i] = A[i].exponent + B[i].exponent - mul_exp_base;
    }
    
    // Do the multiplication.
    ac_int<1,false> mul_overf[N];
    ac_int<E+2,true> mul_exp_result[N];

    ac_int<mul_W,false> mul_prod[N];
    #pragma hls_unroll
    DO_MUL: for (int i=0; i<N; i++) {
      multiply(fracA[i], fracB[i], mul_prod[i]);

      if (AisZero[i] || BisZero[i])  {
        mul_prod[i] = 0;
        mul_exp_result[i] = 0;
      }
          

      mul_overf[i] = mul_prod[i][mul_W-1];
      mul_prod[i] <<= !mul_overf[i];

      mul_exp_result[i] = (mul_overf[i]) ? mul_exp_overf[i] : mul_exp[i];

      if (mul_exp_result[i] > 255 ) {
        mul_exp_result[i] = 255;
        mulisInf[i] = true;
      }
      if (mul_exp_result[i]<0) {
        mul_exp_result[i] = 0;
      }
      
      // bool mul_zero;
      // const int mul_prod_lz_cnt = (mul_prod[i]).leading_sign(mul_zero);
      // if (mul_zero) 
      //   mul_exp_result[i] = 0;  

    }

    // Start the addition
    bool addisInf = false;
    #pragma hls_unroll
    A_INFcheck: for (int i=0; i<N; i++)
      addisInf |= mulisInf[i];

    
    ac_int<E+1,true> diffE[N][N];
    #pragma hls_unroll
    MAXE_H: for (int i=0; i<N; i++) {
      #pragma hls_unroll
      MAXE_V: for (int j=0; j<N; j++) {
        diffE[i][j] = (mul_exp_result[i]) - (mul_exp_result[j]);
      }
    }
    
    bool maxExp_OH[N];
    #pragma hls_unroll
    MAX_E: for (int i=0; i<N; i++) {
      maxExp_OH[i] = true;
      #pragma hls_unroll
      AND_MSB: for (int j=0; j<N; j++) {
        maxExp_OH[i] &= ~diffE[i][j][E];
      }
    }

    bool only1 = true;
    #pragma hls_unroll
    A_ALLIGN: for (int i=0; i<N; i++) {
      #pragma hls_unroll
      for (int j=0; j<N; j++){     
        if (maxExp_OH[i] && only1) { 
          mul_prod[j] >>= diffE[i][j];
          mul_exp_result[j] +=diffE[i][j];
          if (j==N-1) 
            only1 = false;
        }
      }
    }
    
    ac_int<mul_W+1,true> add_op[N];
    #pragma hls_unroll
    TWOs_COMPL: for (int i=0; i<N; i++) {
      if (mulSign[i]) {
        add_op[i] = -mul_prod[i];
      } else {
        add_op[i] = mul_prod[i];
      }
        
    }

    ac_int<mul_W+1+N,true> res1=0;
    #pragma hls_unroll
    MANT_ADD: for (int i=0; i<N; i++) {
      res1 += add_op[i];
    }

    ac_int<1,false> res_sign = res1[mul_W+N];

    ac_int<mul_W+N,false> res;
    res = (res_sign) ? (-res1).template slc<mul_W+N>(0) : res1.template slc<mul_W+N>(0);

    
    const int lead_zer = res.leading_sign();

    res <<= lead_zer;
    mul_exp_result[0] -= (lead_zer);
    mul_exp_result[0] += (lead_zer==0) ? N-1 : N;
   

    mantissa = (addisInf) ? (ac_int<M, false>)(0) : res.template slc<M>(M+N+1); 
    exponent = (addisInf) ? (ac_int<E, false>)((1<<E) -1) : (ac_int<E, false>)(mul_exp_result[0].template slc<E>(0));
    sign = (addisInf) ? infSign : res_sign;
    
// #ifndef __SYNTHESIS__
//     std::cout << "sig= " << sign << std::endl;
//     std::cout << "exp= " << exponent << std::endl;
//     std::cout << "man= " << mantissa << std::endl;
// #endif
    
  }

  
  // OVERLOAD OPERATORS

// #ifndef __SYNTHESIS__
  void operator = (const float &inFP) {
    ac_int<32,false> in= ((ac_ieee_float<binary32>)inFP).data_ac_int();
    exponent = in.template slc<E>(23);
    mantissa = (exponent==0) ? (man_t)0 : in.template slc<M>(23-M);
    sign = (exponent==0) ? (sgn_t)0 : in[31];
  }

  void operator = (const double &inFP) {
    ac_int<64,false> in= ((ac_ieee_float<binary64>)inFP).data_ac_int();
    mantissa = in.template slc<M>(52-M);
    exponent = in.template slc<E>(52);
    sign = in[63];
  }
// #endif

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

    bool x;
    fpa_dual(tmp,r,x);
    return r;
  }
  
  fast_float<M, E> operator * (const fast_float<M, E> &b) {
      fast_float<M, E> r;
      bool x;
      fpmul(b,r,x);
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
    if (sign <= b.sign){
      if (sign < b.sign) {
        return true;
      } else {
        if (exponent >= b.exponent) {
          if (exponent == b.exponent) {
            if (mantissa > b.mantissa) {
              return true;
            } else {
              return false;
            }
          } else {
            return true;
          }
        } else {
          return false;
        }
      }
    } else {
      return false;
    }
  }

  bool operator <  (const fast_float<M,E> b) {
    if (sign >= b.sign){
      if (sign > b.sign) {
        return true;
      } else {
        if (exponent <= b.exponent) {
          if (exponent == b.exponent) {
            if (mantissa < b.mantissa) {
              return true;
            } else {
              return false;
            }
          } else {
            return true;
          }
        } else {
          return false;
        }
      }
    } else {
      return false;
    }
  }

  bool operator >= (const fast_float<M,E> b) {
    if (sign <= b.sign){
      if (exponent >= b.exponent) {
        if (exponent == b.exponent) {
          if (mantissa >= b.mantissa) {
            return true;
          } else {
            return false;
          }
        } else {
          return true;
        }
      } else {
        return false;
      }
    } else {
      return false;
    }
  }

  bool operator <= (const fast_float<M,E> b) {
    if (sign >= b.sign){
      if (exponent <= b.exponent) {
        if (exponent == b.exponent) {
          if (mantissa <= b.mantissa) {
            return true;
          } else {
            return false;
          }
        } else {
          return true;
        }
      } else {
        return false;
      }
    } else {
      return false;
    }
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
