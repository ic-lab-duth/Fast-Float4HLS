# Fast-Float
Fast-Float is a C++ header only library for floating point arithmetic operations developed to be used for High Level Synthesis implemetations. 

Inside the header file <fast-float.h> is described a templated datatype 'fast_float<M,E>' together with a set of operations. The datatype supports different floating point representations depending on the definitions of the mantissa M and exponent E widths through the template.
For example:
* fast_float<23,8>  is used for the single precision representation.
* fast_float<52,11> is used for the double precision representation.
* fast_float<7,8>  is used for the representation of Brain Floating Point, bfloat16, developed by Google Brain.

To use the fast-float datatype include the header file.
>#include "fast_float.h"



The library includes features from the <ac_int.h> and <ac_std_float.h> libraries from HLSLibs (''https://github.com/hlslibs'').

# Supported Features

* Addition
:  Dual Path
* Multiplication
* MAC
* Dot Product
* a set of overloaded operators for arithmetic operations, comparisons and assignments

## Denormals
The operations of addition, multiplication and MAC support denormalized values through a template parameter ''DENORMALS''. By default, the operators do not support denormalized values, while converting them into zero values.

The operators +, -, *, +=, -=, *= compute the corresponding operation without support for denormalized values.

# Pending Features

* Increase the output precision of dot product.
* support rounding on dot product

