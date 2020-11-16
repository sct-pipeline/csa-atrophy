#SBATCH --account=def-jcohen
#SBATCH --time=0-12:00        # time (DD-HH:MM)
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32    # number of OpenMP processes
#SBATCH --mem=128G
#SBATCH --mail-user=paul.bautin@polymtl.ca
#SBATCH --mail-type=ALL
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd $SCRATCH
