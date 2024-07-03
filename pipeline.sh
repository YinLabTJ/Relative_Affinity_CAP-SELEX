export PATH="/data/software/spacek40:$PATH" #spacek40 path
perl step0_get_list.pl Batch_test.txt
perl step1_spacek40_generate.pl /data/workdata/SELEX # /data/workdata/SELEX is the directory of seq files
sh step1_spacek40_sample.sh
sh step1_spacek40_pair.sh
perl step2.generate.pl Batch_test.txt filter.8mer.xls /data/workdata/SELEX #filter.8mer.xls come from the early step of MI pipeline
sh step2.sh
perl step3.create_R_input.pl Batch_test.txt
perl step4.run.pl Batch_test.txt
