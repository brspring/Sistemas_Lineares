#!/bin/bash

METRICA="FLOPS_DP"
CPU=3

LIKWID_HOME=/home/soft/likwid
CFLAGS="-I${LIKWID_HOME}/include -DLIKWID_PERFMON"
LFLAGS="-L${LIKWID_HOME}/lib -llikwid"

PROGRAMA="perfSL"

echo "performance" > /sys/devices/system/cpu/cpufreq/policy3/scaling_governor

likwid-perfctr -C ${CPU} -g ${METRICA} -m ${PROGRAMA}

echo "powersave" > /sys/devices/system/cpu/cpufreq/policy3/scaling_governor

"FLOPS_DP" | grep -v "AVX"