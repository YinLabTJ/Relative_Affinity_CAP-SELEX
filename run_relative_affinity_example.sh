#!/bin/bash

# An example of computing relative affinity for our published CAP-SELEX matrices

export PATH="/data/software/relative_affinity/upload_pipeline_v2:$PATH"  #spacek40 path

#1. run spacek40 for local_max information
perl script/step1_spacek40.pl Batch_test.txt SELEX data/Curated_Prey_Final.txt

#2. computing relative affinity
#filter_list was created with MI pipeline using the same batch_file
for i in `grep -v HGNC Batch_test.txt | awk '{print $1"_"$3}'`
do
	perl script/step2.relative_affinity.pl Batch_test.txt $i 10 ./SELEX data/Curated_Prey_Final.txt MI_output/filter_list
done

#the result will be in ./output
