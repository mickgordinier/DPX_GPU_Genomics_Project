// Defining the SIMD intrinsic functions that are only supported in device code.
// These instructions specifically are accelerated by the DPX hardware

// Instruction Definition Can be found: https://docs.nvidia.com/cuda/cuda-math-api/cuda_math_api/group__CUDA__MATH__INTRINSIC__SIMD.html
// DPX instruction (Section 7.25): https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html?highlight=dpx#dpx

// NOTE: I believe the initial "vi_" stands for Vector (Vector SIMD processing operating) integer

// NOTE: These functions consist of MIN/MAX instructions for 16/32 bit signed/unsigned integers
// NOTE: Additional features include: ReLU, Pointer returning, and addition maxing

class FakeDPX {

  public:

  /* 3 PARAMETERS FUNCTIONS */
  
  // Computes max(max(a, b), c)
  static int __vimax3_s32(const int a, const int b, const int c);

  // Performs per-halfword max(max(a, b), c)
  // NOTE: INTENTIALLY LEFT AS UNSIGNED INT, WE ARE 2 DOING MAX FUNCTIONS!!!!!!
  // NOTE: MAX done on upper 16-bits, and lower 16-bit
  /* 
    Splits 4 bytes of each argument into 2 parts, each consisting of 2 bytes. 
    These 2 byte parts are interpreted as signed shorts. 
    For corresponding parts function performs a 3-way max ( = max(max(a_part, b_part), c_part) ). 
    Partial results are recombined and returned as unsigned int. 
  */
  static unsigned int __vimax3_s16x2(const unsigned int a, const unsigned int b, const unsigned int c);
  
  // Computes max(max(a, b), c)
  static unsigned int __vimax3_u32(const unsigned int a, const unsigned int b, const unsigned int c);

  // Performs per-halfword max(max(a, b), c)
  static unsigned int __vimax3_u16x2(const unsigned int a, const unsigned int b, const unsigned int c);

  // Computes min(min(a, b), c)
  static int __vimin3_s32(const int a, const int b, const int c);

  // Performs per-halfword min(min(a, b), c)
  static unsigned int __vimin3_s16x2(const unsigned int a, const unsigned int b, const unsigned int c);

  // Computes min(min(a, b), c)
  static unsigned int __vimin3_u32(const unsigned int a, const unsigned int b, const unsigned int c);

  // Performs per-halfword min(min(a, b), c)
  static unsigned int __vimin3_u16x2(const unsigned int a, const unsigned int b, const unsigned int c);


  /* 2 PARAMETERS + ReLU FUNCTIONS */
  // NOTE: Only exist for signed numbers as ReLU is not needed for unsigned numbers

  // Computes max(max(a, b), 0)
  static int __vimax_s32_relu(const int a, const int b);

  // Performs per-halfword max(max(a, b), 0)
  static unsigned int __vimax_s16x2_relu(const unsigned int a, const unsigned int b);

  // Computes max(min(a, b), 0)
  // NOTE: WE STILL MAX ON 0 AS WE DON'T WANT ANY NEGATIVE NUMBERS
  static int __vimin_s32_relu(const int a, const int b);

  // Performs per-halfword max(min(a, b), 0)
  static unsigned int __vimin_s16x2_relu(const unsigned int a, const unsigned int b);


  /* 3 PARAMETERS + ReLU FUNCTIONS */
  static int __vimax3_s32_relu(const int a, const int b, const int c);
  static unsigned int __vimax3_s16x2_relu(const unsigned int a, const unsigned int b, const unsigned int c);
  static int __vimin3_s32_relu(const int a, const int b, const int c);
  static unsigned int __vimin3_s16x2_relu(const unsigned int a, const unsigned int b, const unsigned int c);


  /* 2 PARAMETERS + RETURNING SMALLER/LARGER PARAMETER FUNCTIONS */

  // Computes max(a, b), also sets the value pointed to by pred to (a >= b).
  // NOTE: "bool *const" means a constant POINTER to a non-constant boolean value
  static int __vibmax_s32(const int a, const int b, bool *const pred);
  static unsigned int __vibmax_u32(const unsigned int a, const unsigned int b, bool *const pred);

  // Computes min(a, b), also sets the value pointed to by pred to (a <= b).
  static int __vibmin_s32(const int a, const int b, bool *const pred);
  static unsigned int __vibmin_u32(const unsigned int a, const unsigned int b, bool *const pred);

  // Performs per-halfword max(a, b), also sets the value pointed to by pred_hi and pred_lo to the per-halfword result of (a >= b).
  // NOTE: As we are doing 2 MAX functions, we have to keep track of 2 pointers for the upper/lower 16-bits
  static unsigned int __vibmax_s16x2(const unsigned int a, const unsigned int b, bool *const pred_hi, bool *const pred_lo);
  static unsigned int __vibmax_u16x2(const unsigned int a, const unsigned int b, bool *const pred_hi, bool *const pred_lo);

  // Performs per-halfword min(a, b), also sets the value pointed to by pred_hi and pred_lo to the per-halfword result of (a <= b).
  // NOTE: As we are doing 2 MIN functions, we have to keep track of 2 pointers for the upper/lower 16-bits
  static unsigned int __vibmin_s16x2(const unsigned int a, const unsigned int b, bool *const pred_hi, bool *const pred_lo);
  static unsigned int __vibmin_u16x2(const unsigned int a, const unsigned int b, bool *const pred_hi, bool *const pred_lo);


  /* 3 PARAMETERS, COMPARING (FIRST + SECOND) WITH THIRD FUNCTIONS */

  // Computes max(a + b, c)
  static int __viaddmax_s32(const int a, const int b, const int c);
  static unsigned int __viaddmax_u32(const unsigned int a, const unsigned int b, const unsigned int c);

  // Performs per-halfword max(a + b, c)
  static unsigned int __viaddmax_s16x2(const unsigned int a, const unsigned int b, const unsigned int c);
  static unsigned int __viaddmax_u16x2(const unsigned int a, const unsigned int b, const unsigned int c);

  // min calculations
  static int __viaddmin_s32(const int a, const int b, const int c);
  static unsigned int __viaddmin_u32(const unsigned int a, const unsigned int b, const unsigned int c);
  static unsigned int __viaddmin_s16x2(const unsigned int a, const unsigned int b, const unsigned int c);
  static unsigned int __viaddmin_u16x2(const unsigned int a, const unsigned int b, const unsigned int c);

  
  /* 3 PARAMETERS, COMPARING (FIRST + SECOND) WITH THIRD AND ReLU FUNCTIONS */

  // Computes max(max(a + b, c), 0)
  static int __viaddmax_s32_relu(const int a, const int b, const int c);

  // Performs per-halfword max(max(a + b, c), 0)
  static unsigned int __viaddmax_s16x2_relu(const unsigned int a, const unsigned int b, const unsigned int c);

  // Computes max(min(a + b, c), 0)
  static int __viaddmin_s32_relu(const int a, const int b, const int c);

  // Performs per-halfword max(min(a + b, c), 0)
  static unsigned int __viaddmin_s16x2_relu(const unsigned int a, const unsigned int b, const unsigned int c);
};