#! /bin/bash
#SBATCH -N 2
#SBATCH --ntasks 32
mpirun ./a2 --m 2688 --n 4096 --epsilon 0.001 --max-iterations 1000
mpirun ./a2 --m 2688 --n 4096 --epsilon 0.001 --max-iterations 2000
mpirun ./a2 --m 1152 --n 1152 --epsilon 0.001 --max-iterations 1000