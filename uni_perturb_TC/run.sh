#!/bin/sh

mkdir -p output_files
nThreads=32   # Adjusted to match ntasks

for l in $(seq 32 32 64)
do
  for clusterId in $(seq 1 12)
  do
    cat << EOS | sbatch --verbose
#!/bin/tcsh
#SBATCH --job-name=tc_noisy
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=32               # Adjusted from 64 to 32
#SBATCH --cpus-per-task=1
#SBATCH --mem=120G
#SBATCH --time=10:00:00

#SBATCH --output=output_files/job_%A_%a.out
#SBATCH --error=output_files/job_%A_%a.err

module load matlab
cd \$PWD

matlab -batch "main $l $clusterId $nThreads; exit"
EOS
  done
done
