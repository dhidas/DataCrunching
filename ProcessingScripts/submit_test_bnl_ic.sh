module load gcc/6.4.0 python/3.6.5 /sdcc/covid19/dhidas/modulefiles/covid19
sbatch -t 1-00:00:00 -A covid-19 -p long -J covidt1 -N 2 -n 72 --wrap="srun python3 run_mpi.py ena+db-test.can"
