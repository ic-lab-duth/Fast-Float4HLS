# Bubble Sort implementation

This is an example implementation of a Bubble Sorting algorithm using the Fast-Float4HLS library. The implementation can be found in the bubble_sort.cpp file together with the testbench main function. The bubbleSort function is a templatized C++ function implemented for HLS using the Catapult HLS tool. The size of the input array to be sorted can be selected through the template parameter ''N''. 

To synthesize the design on Catapult HLS use the *go_hls.tcl* script. The given example synthesizes a bubble sort algorithm that sorts the 10 elements of an array. When changing the size of the array through the template parameter 

```c++
const int NSIZE = 10;
bubbleSort<NSIZE>(inA);
```

you should also change the value of the variable ''ARRAY_SIZE'' at the top of the *go_hls.tcl* TCL script.

```tcl
set ARRAY_SIZE 10
```

To synthesize the design launch catapult in the directory of the example.

```bash
catapult -f go_hls.tcl
```

