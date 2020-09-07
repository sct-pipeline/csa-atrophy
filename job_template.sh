#SBATCH --account=def-jcohen
#SBATCH --time=0-08:00        # time (DD-HH:MM)
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32    # number of OpenMP processes
#SBATCH --mem=128G
#SBATCH --mail-user=paul.bautin@polymtl.ca
#SBATCH --mail-type=ALL
#SBATCH --mail-type=END

cd $SCRATCH
