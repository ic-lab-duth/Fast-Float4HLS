# Define template parameters for top level
#   IF_HEIGHT  : Height of the input features map
#   IF_WIDTH   : Width  of the input features map
#   KERNEL_S   : Size of the Kernel (K x K) 

set IF_HEIGHT 14
set IF_WIDTH  14
set KERNEL_S  5

# Define mantissa and exponent width of the datatype
# used for the synthesis
set MANW 7
set EXPW 8

# Define top level module
set TOP "conv2D<${IF_HEIGHT}, ${IF_HEIGHT}, ${KERNEL_S}>"
set TOP_HIER "conv2D<${IF_HEIGHT},${IF_HEIGHT},${KERNEL_S}>"

# Define path to top (tb) cpp file
set FPATH "./conv2D.cpp"

# Define project name
set PROJECT "conv2D_example"

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
directive set -CLOCKS {clk {-CLOCK_PERIOD 8.0 -CLOCK_UNCERTAINTY 0.0 -CLOCK_HIGH_TIME 4.0}}
solution library add nangate-45nm_beh -- -rtlsyntool OasysRTL -vendor Nangate -technology 045nm
solution library add ccs_sample_mem
go libraries
go assembly
directive set "/${TOP_HIER}/core/lb.mantissa:rsc" -PACKING_MODE compact
directive set "/${TOP_HIER}/core/lb.exponent"     -RESOURCE lb.mantissa:rsc
directive set "/${TOP_HIER}/core/lb.exponent"     -BASE_BIT $MANW
directive set "/${TOP_HIER}/core/lb.exponent:rsc" -MAP_TO_MODULE {}
directive set "/${TOP_HIER}/core/lb.sign"         -RESOURCE lb.mantissa:rsc
directive set "/${TOP_HIER}/core/lb.sign"         -BASE_BIT [expr $MANW+$EXPW]
directive set "/${TOP_HIER}/core/lb.sign:rsc"     -MAP_TO_MODULE {}
directive set "/${TOP_HIER}/core/lb.mantissa:rsc" -BLOCK_SIZE $IF_WIDTH
go architect
ignore_memory_precedences -from *write_mem(lb* -to *read_mem(lb*
go allocate
go extract

#COSIM
flow run /SCVerify/launch_make ./scverify/Verify_concat_sim_rtl_v_msim.mk {} SIMTOOL=msim sim

