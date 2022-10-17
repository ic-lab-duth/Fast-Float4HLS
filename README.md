# Fast-Float4HLS
Fast-Float4HLS is a C++ header only library for floating point arithmetic operations developed to be used for High Level Synthesis implemetations. 

Inside the header file <fast-float.h> is described a templated datatype 'fast_float<M,E>' together with a set of operations. The datatype supports different floating point representations depending on the definitions of the mantissa M and exponent E widths through the template.
For example:
* fast_float<23,8>  is used for the single precision representation.
* fast_float<52,11> is used for the double precision representation.
* fast_float<7,8>  is used for the representation of Brain Floating Point, bfloat16, developed by Google Brain.

Fast-Float4HLS depends only on ac_fixed library that is available in [HLSLibs](https://github.com/hlslibs/ac_types).

Also the post-synthesis RTL co-simultion of the given examples require the sc_verify flow of Catapult HLS. The necessary header 
(mc_scverify) is publicly available in [ac_simutils](https://github.com/hlslibs/ac_simutils/tree/master/include).

Downloading the needed header-only libraries can be done by running the provided script ```set_libs.sh```

# Supported Operators

* Addition (A+B)
* Multiplication (A*B)
* Multiply-Add (A*B + C)
* Vector Dot Product (A[0]*B[0]+A[1]*B[1]+...+A[N-1]*B[N-1])

Addition and multiplication are overloaded also to + and * operators. The same holds for comparisons and assignment operators

## Denormals
The operations of addition, multiplication and MAC support denormalized values through a template parameter ''DENORMALS''. By default, the operators do not support denormalized values, while converting them into zero values.

Currently the overloaded operators +, -, *, +=, -=, *= compute the corresponding operation without support for denormalized values.

# Pending Features

* Allow for possible increase the output precision of dot product.
* Implement Division

# Reference

The architecture and performance of the fused vector dot product unit, implemented as part of the FastFloat4HLS library, was published in the MDPI Journal of Low Power Electronics and Applications on Oct. 2022. You can find the paper [here](https://www.mdpi.com/2079-9268/12/4/56/pdf). To cite this work, please use:

```
@article{fused-fp-dot,
author = {Filippas, Dionysios and Nicopoulos, Chrysostomos and Dimitrakopoulos, Giorgos},
title = {Templatized Fused Vector Floating-Point Dot Product for High-Level Synthesis},
journal = {Journal of Low Power Electronics and Applications},
volume = {12},
year = {2022},
number = {4},
article-number = {56},
```

# Contributors

Currently active: [Dionysios Filippas](https://github.com/dionisisfil), [Christodoulos Peltekis](https://github.com/chrispelt)
and [Giorgos Dimitrakopoulos](https://github.com/gdimitrak)

Past: Nikolaos Altanis

# License
Fast-Float4HLS is licensed with the MIT License. You are completely free to re-distribute your work derived from Fast-Float4HLS

