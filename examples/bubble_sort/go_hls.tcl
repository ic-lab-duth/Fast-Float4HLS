# Define template parameters for top level
#   ARRAY_SIZE       : The size of the array to be sorted

set ARRAY_SIZE 10

# Definition of top level module
set TOP "bubbleSort<${ARRAY_SIZE}>"

# Definition path to top (tb) cpp file
set FPATH "./bubble_sort.cpp"

# Definition project name
set PROJECT "bubble_sort_example"

# GO_HLS

project new -name $PROJECT
solution new -state initial
solution options defaults
solution options set /ComponentLibs/SearchPath . -append
solution options set /Interface/DefaultClockPeriod 2
solution options set /Input/CppStandard c++11
solution options set /Input/CompilerFlags { -DHLS_CATAPULT=1 }
solution options set /Input/SearchPath { ../../ }
solution options set /Flows/QuestaSIM/SCCOM_OPTS {-O3 -x c++ -Wall -Wno-unused-label -Wno-unknown-pragmas}
solution options set /Flows/SCVerify/USE_CCS_BLOCK true
flow package require /SCVerify
flow package require /QuestaSIM
solution file add $FPATH -type C++
directive set -DESIGN_GOAL latency
go new
solution design set $TOP -top
go analyze
directive set -CLOCKS {clk {-CLOCK_PERIOD 2.0 -CLOCK_UNCERTAINTY 0.0 -CLOCK_HIGH_TIME 1.0}}
solution library add nangate-45nm_beh -- -rtlsyntool OasysRTL -vendor Nangate -technology 045nm
solution library add ccs_sample_mem
go libraries
go assembly
go architect
go allocate
go extract

#COSIM
flow run /SCVerify/launch_make ./scverify/Verify_concat_sim_rtl_v_msim.mk {} SIMTOOL=msim sim

