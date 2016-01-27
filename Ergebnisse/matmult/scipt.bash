#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --partition=sandy
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:0

module load gcc
module load intel

gcc -O0 -std=c99 ./main.c
./a.out > gccO0.csv
gcc -O1 -std=c99 ./main.c
./a.out > gccO1.csv
gcc -O2 -std=c99 ./main.c
./a.out > gccO2.csv
gcc -O3 -std=c99 ./main.c
./a.out > gccO3.csv
icc -O0 -std=c99 ./main.c
./a.out > iccO0.csv
icc -O1 -std=c99 ./main.c
./a.out > iccO1.csv
icc -O2 -std=c99 ./main.c
./a.out > iccO2.csv
icc -O3 -std=c99 ./main.c
./a.out > iccO3.csv

