# Fast-Float4HLS
Fast-Float4HLS is a C++ header only library for floating point arithmetic operations developed to be used for High Level Synthesis implemetations. 

Inside the header file <fast-float.h> is described a templated datatype 'fast_float<M,E>' together with a set of operations. The datatype supports different floating point representations depending on the definitions of the mantissa M and exponent E widths through the template.
For example:
* fast_float<23,8>  is used for the single precision representation.
* fast_float<52,11> is used for the double precision representation.
* fast_float<7,8>  is used for the representation of Brain Floating Point, bfloat16, developed by Google Brain.

Fast-Float4HLS is based on ac_fixed library that is available in [HLSLibs](https://github.com/hlslibs/ac_types).

The execution of the given examples require the inclusion of mc_scverify library that is available in [ac_simutils](https://github.com/hlslibs/ac_simutils/tree/master/include).

The installation of the required libraries can be done by downloading the code from the links above or by running the ```set_libs.sh``` script, that downloads the libraries on the working directory, or by using the following commands:

```console
git clone http://github.com/hlslibs/ac_types.git
git clone http://github.com/hlslibs/ac_simutils.git
```

# Supported Operators

* Addition
* Multiplication
* Multiply-Add
* Dot Product

Addition and multiplication are overloaded also to + and * operators. The same holds for comparisons and assignment operators

## Denormals
The operations of addition, multiplication and MAC support denormalized values through a template parameter ''DENORMALS''. By default, the operators do not support denormalized values, while converting them into zero values.

Currently the overloaded operators +, -, *, +=, -=, *= compute the corresponding operation without support for denormalized values.

# Pending Features

* Allow for possible increase the output precision of dot product.
* support rounding on dot product
* Implement Division

# Contributors

Currently active: [Dionysios Filippas](https://github.com/dionisisfil), [Giorgos Dimitrakopoulos](https://github.com/gdimitrak)

Past: Nikolaos Altanis

# License
Fast-Float4HLS is licensed with the MIT License. You are completely free to re-distribute your work derived from Fast-Float4HLS

