#!/bin/bash

export JAVA_HOME=/home/jp2045/software/jdk-23.0.1
export PATH=/home/jp2045/software/jdk-23.0.1/bin:$PATH

submit_flag=''

print_usage() {
  printf "Usage: -s to run in a slurm script, otherwise runs nextflow locally\n"
}

while getopts 's' flag; do
  case "${flag}" in
    s) submit_flag='true' ;;
    *) print_usage
       exit 1 ;;
  esac
done

if [[ "$submit_flag" = "true" ]]
then
  sbatch slurm-run.sh
else
  ./nextflow run main.nf -resume
fi

