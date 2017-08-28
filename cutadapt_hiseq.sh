#!/bin/bash

# cutadapt finally supports paired-end trimming!
#
# provide arguments as follows: 
#   cutadapt_hiseq.sh <6-char tag> <read1> <read2>
# 
# make sure read 1 contains "R1" and read 2 "R2".
#
# optional flags:
#   -q FRONT,BACK   quality (from front, from back). trim at q=20
#   -m LENGTH       retain reads of minimum length
#   --trim-n        removes flanking Ns (NNNATCGNNN --> ATCG)

outfile=`echo ${2} | sed 's/.fastq$//' | sed 's/.fq$//'`.trim.fastq
logfile=`echo ${2} | sed 's/.fastq$//' | sed 's/.fq$//'`.trim.log

cutadapt -O 10 -q 20,20 -m 25 --trim-n -b GATCGGAAGAGCACACGTCTGAACTCCAGTCAC${1}ATCTCGTATGCCGTCTTCTGCTTG -b AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -B GATCGGAAGAGCACACGTCTGAACTCCAGTCAC${1}ATCTCGTATGCCGTCTTCTGCTTG -B AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT ${2} ${3} -o ${outfile} -p ${outfile/R1/R2} > ${logfile}
