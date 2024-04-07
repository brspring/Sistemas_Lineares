#!/bin/bash

METRICA="FLOPS_DP"
CPU=15
LIKWID_HOME=/home/soft/likwid

PROGRAMA="./perfSL"

make
# Redireciona a saída do likwid-perfctr para um arquivo
likwid-perfctr -C ${CPU} -g ${METRICA} -m ${PROGRAMA} > saida_likwid.txt

cat saida_likwid.txt | grep -A 20  "EG clássico"
cat saida_likwid.txt | grep "DP" | grep -v "AVX"
