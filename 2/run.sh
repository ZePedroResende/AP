#!/bin/bash
module load cmake/3.10.0
module load gcc/5.3.0
module load gnu/openmpi_eth/1.8.4
cd $HOME/AP2/2/build
rm -rf *
cmake ..
make

mpirun -np 18 -mca btl self,sm,tcp rb
