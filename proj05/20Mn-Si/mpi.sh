#!/bin/bash
#
#$ -q all.q
#$ -pe mpi_16 16
#$ -N si
#$ -M impurity@oz
#$ -m bea
#$ -o run.out
#$ -e run.err
#$ -V
#$ -cwd
#

mpirun -np 16 python step0.py
