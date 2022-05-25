# Define template parameters for top level
#   DAT_TYPE   : The name of the defined datatype in the source code
#   NODEs      : The number of the graph's nodes
#   INP_FEAT   : The size of the input featured vector
#   OUT_F_L1   : The size of the output features vector of the first layer
#   OUT_F_L2   : The size of the output features vector of the second layer
#   NON_ZERO   : The number of non-zero values in the adjacency matrix

set DAT_TYPE btype

set NODES    3327
set INP_FEAT 3703
set OUT_F_L1 21
set OUT_F_L2 6
set NON_ZERO 12431

# Define top level module
set TOP "gcn<${DAT_TYPE}, ${NODES}, ${INP_FEAT}, ${OUT_F_L1}, ${OUT_F_L2}, ${NON_ZERO}>"

# Define path to top (tb) cpp file
set FPATH "./gcn_tb.cpp"

# Define project name
set PROJECT "gcn_example"

# GO_HLS

project new -name $PROJECT
solution new -state initial
solution file add ./fast_float.h
solution file add ./gcn.h
solution file add ./gcn_tb.cpp
solution file add ./helper.h
solution file add ./defs.h
solution file add ./matrix.h
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
directive set -CLOCKS {clk {-CLOCK_PERIOD 10.0 -CLOCK_UNCERTAINTY 0.0 -CLOCK_HIGH_TIME 5.0}}
solution library add nangate-45nm_beh -- -rtlsyntool OasysRTL -vendor Nangate -technology 045nm
solution library add ccs_sample_mem
go libraries
go assembly
go architect
go allocate
go extract

#COSIM
#flow run /SCVerify/launch_make ./scverify/Verify_concat_sim_rtl_v_msim.mk {} SIMTOOL=msim sim


