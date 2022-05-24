# Matrix multiplication implementation

This example provides three different implementations matrix-matrix multiplications using the Fast-Float4HLS library. 

> matrix_mul_basic.cpp : Contains an implementation using the basic operators + and *.

> matrix_mul_fma.cpp   : Contains an implementation using the fma implementation.

> matrix_mul_dot.cpp    : Contains an implementation using the dot product implementation. 

All implementations are based on the Fast-Float4HLS library. The templatized C++ functions are implemented for HLS using the Catapult HLS tool. The sizes of the two input matrices to be multiplied can be selected through the template parameters ''N'', ''M'' and ''K''. 

To synthesize the design on Catapult HLS use the *go_hls.tcl* script. The given examples synthesize the multiplication of two matrices of size 4x4 and 4x2 accordingly. When changing the size of the matrices through the template parameter 

```c++
const int N=4, M=4, K=2;
matXmat<N,M,K>(A,B,OM);
```

you should also change the value of the corresponding variables at the top of the *go_hls.tcl* TCL script.

```tcl
set N 4
set M 4
set K 2
```

Furthermore, through the SET_IMPL variable you can select the example to be synthesized.

```tcl
set IMPLEMENTATIONS_LIST {basic fma dot}
set SEL_IMPL 2
```

To synthesize the design launch catapult in the directory of the example.

```bash
catapult -f go_hls.tcl
```

