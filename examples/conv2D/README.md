# 2D Convolution Engine implementation

This is an example implementation of a 2D convolution engine using the Fast-Float4HLS library. The outputs are being calculated using the dot product function of the Fast-Float4HLS library. The implementation can be found in the conv2D.cpp file together with the testbench main function. The conv2D is a templatized C++ function implemented for HLS using the Catapult HLS tool. The template parameters define the size of the input features map (''H'' and ''W'') as well as the size of the kernel (''K'').

To synthesize the design on Catapult HLS use the *go_hls.tcl* script. The given example synthesizes a 2D convolution engine that convolves a 5x5 kernel with a 14x14 input. When changing the template parameters in the source code, that define the design

```c++
const int H = 14;
const int W = 14;
const int K = 5;
conv2D<H,W,K>(img,flt,out);
```

you should also change the corresponding variables at the top of the *go_hls.tcl* TCL script.

```tcl
set IF_HEIGHT 14
set IF_WIDTH  14
set KERNEL_S  5
```

Apart from the input features and kernel sizes, the tcl scripts uses the sizes of the floating point representation fields to help with the architecture definition of the line buffers.

```tcl
set MANW 7
set EXPW 8
```

To synthesize the design launch catapult in the directory of the example.

```bash
catapult -f go_hls.tcl
```

