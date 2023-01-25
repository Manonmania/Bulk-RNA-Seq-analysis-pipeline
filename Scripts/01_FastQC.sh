#!/bin/bash -l

spack load fastqc@0.11.9

mkdir FastQC

for my_file in rawREADS/*.fastq.gz
        do
                fastqc -o FastQC $my_file
        done
