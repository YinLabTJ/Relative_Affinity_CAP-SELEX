#!/bin/bash

# An example of computing relative affinity for our published CAP-SELEX matrices
# To run this pipeline, you need a batch file with TF_pairs information (File Batch_test is an example)and a folder (e.g. ./SELEX) that stores sequence file
# The sequence files should contain DNA sequence only, not fastq format. And the name of the sequence files should be consistent with the batch file.
# Take "Batch_test.txt" and folder “SELEX" as an example

#export PATH="/data/software/Relative_Affinity_CAP-SELEX:$PATH"  #spacek40 path

#1. run spacek40 for local_max information
perl script/step1_spacek40.pl Batch_test.txt SELEX

#2. computing relative affinity
#filter_list was created with MI pipeline using the same batch_file
for i in `grep -v HGNC Batch_test.txt | awk '{print $1"_"$3}'`
do
	perl script/step2.relative_affinity.pl Batch_test.txt $i 10 SELEX MI_output/filter_list
done

#the output tables and figures should be in ./output
