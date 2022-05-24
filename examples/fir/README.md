# FIR Filter implementation

This is an example implementation of a Finite Impulse Response filter using the Fast-Float4HLS library. The filter is implemented using the fma function of the Fast-Float4HLS library and can be found in the fir.cpp file together with the testbench main function. The fir function is a templatized C++ function implemented for HLS using the Catapult HLS tool. The order of the filter can be selected through the template parameter ''N''. 

To synthesize the design on Catapult HLS use the *go_hls.tcl* script. The given example synthesizes a 5-order FIR filter. When changing the order of the filter through the template parameter 

```c++
const int TAPS = 5;
fir<TAPS>(In,Coeff,Out);
```

you should also change the value of the variable ''TAPS'' at the top of the *go_hls.tcl* TCL script.

```tcl
set TAPS 5
```

To synthesize the design launch catapult in the directory of the example.

```bash
catapult -f go_hls.tcl
```

