#!/bin/sh

mkdir -p output_files
nThreads=16   # Adjusted to match ntasks

for l in $(seq 128 128)
do
  for clusterId in $(seq 21 25)
  do
    cat << EOS | sbatch --verbose
#!/bin/tcsh
#SBATCH --job-name=tc_noisy
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=12               # Adjusted from 64 to 32
#SBATCH --cpus-per-task=1
#SBATCH --mem=240G
#SBATCH --time=23:00:00

#SBATCH --output=output_files/job_%A_%a.out
#SBATCH --error=output_files/job_%A_%a.err

module load matlab
cd \$PWD

matlab -batch "main $l $clusterId $nThreads; exit"
EOS
  done
done
