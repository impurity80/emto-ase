#!/bin/bash
#
#$ -q all.q
#$ -pe mpi_32 32
#$ -N co-ni
#$ -M impurity@oz
#$ -m bea
#$ -o run.out
#$ -e run.err
#$ -V
#$ -cwd
#

mpirun -np 32 python step0.py
