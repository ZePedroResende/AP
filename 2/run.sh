#!/bin/bash
module load gcc/5.3.0
module load gnu/openmpi_eth/1.8.4
cd $HOME/AP/2/build
rm -rf *
cmake ..
make

mpirun -np 2 -mca btl self,sm,tcp rb
