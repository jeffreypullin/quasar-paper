#!/bin/bash

#SBATCH -J quasar-analysis-run
#SBATCH -A mrc-bsu-sl2-cpu
#SBATCH -p icelake-himem
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --mail-type=NONE
#SBATCH --ntasks=2
#SBATCH --no-requeue

. /etc/profile.d/modules.sh
module purge
module load rhel8/default-icl

export JAVA_HOME=/home/jp2045/software/jdk-23.0.1
export PATH=/home/jp2045/software/jdk-23.0.1/bin:$PATH

workdir="$SLURM_SUBMIT_DIR"

./nextflow run main.nf -resume
