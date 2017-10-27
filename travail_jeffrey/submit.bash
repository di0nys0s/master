#!/bin/bash 

#PBS -N HelloWorld 
#PBS -A cjt-923-aa 
#PBS -l walltime=300 
#PBS -l nodes=2:ppn=8
module load blas-libs/mkl/11.0   
module load compilers/intel/2013
./cscanE
