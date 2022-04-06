# Fast-Float4HLS
Fast-Float is a C++ header only library for floating point arithmetic operations developed to be used for High Level Synthesis implemetations. 

Inside the header file <fast-float.h> is described a templated datatype 'fast_float<M,E>' together with a set of operations. The datatype supports different floating point representations depending on the definitions of the mantissa M and exponent E widths through the template.
For example:
* fast_float<23,8>  is used for the single precision representation.
* fast_float<52,11> is used for the double precision representation.
* fast_float<7,8>  is used for the representation of Brain Floating Point, bfloat16, developed by Google Brain.

This library utilizes features from the ac_int and ac_std_float libraries that are available in HLSLibs (''https://github.com/hlslibs'').

# Supported Operators

* Addition
* Multiplication
* Multiply-Add
* Dot Product

Addition and multiplication are overloaded also to + and * operators. The same holds for comparisons and assignment operators

## Denormals
The operations of addition, multiplication and MAC support denormalized values through a template parameter ''DENORMALS''. By default, the operators do not support denormalized values, while converting them into zero values.

The operators +, -, *, +=, -=, *= compute the corresponding operation without support for denormalized values.

# Pending Features

* Increase the output precision of dot product.
* support rounding on dot product
* Optimize the multiply operation

# Examples

* FIR filter

