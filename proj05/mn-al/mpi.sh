#!/bin/bash
#
#$ -q all.q
#$ -pe mpi_32 32
#$ -N mn-al
#$ -M impurity@oz
#$ -m bea
#$ -o run.out
#$ -e run.err
#$ -V
#$ -cwd
#

mpirun -np 32 python step0.py
