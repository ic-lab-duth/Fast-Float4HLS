# Fast-Float
Fast-Float is a C++ library for floating point arithmetic operations meant to be used by High Level Synthesis tools.

# Operators

Addition, Multiplication, MAC, Dot Product, comparisong and assignment operators

## Denormals
The operations of addition, multiplication and MAC support denormalized values through a template parameter ''DENORMALS''. By default, the operators do not support denormalized values, while converting them into zero values.

The operators +, -, *, +=, -=, *= compute the corresponding operation without support for denormalized values.