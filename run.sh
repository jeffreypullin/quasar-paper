#!/bin/bash

export JAVA_HOME=/home/jp2045/software/jdk-23.0.1
export PATH=/home/jp2045/software/jdk-23.0.1/bin:$PATH

./nextflow run main.nf -resume
