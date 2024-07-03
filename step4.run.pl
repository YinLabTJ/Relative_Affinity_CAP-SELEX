use Cwd;
use strict;
my $batch_file=shift;
my $dir=cwd;
open BATCH,"$batch_file";
while(<BATCH>){
	chomp;
	my @b=split;
	next unless($b[0]=~/_/);
	my $tf_pair=$b[0]."_".$b[2];
	if(-e "output/$tf_pair/Cor_input2.xls"){
		`perl step4.sort_and_merge.pl $tf_pair`;
		`cd $dir/output/$tf_pair && /usr/bin/Rscript $dir/CAP-SELEX_vs_HT-SELEX.R`;
	}
	`cd $dir`;
}
close IN; close OUT; close OUT2; close OUT3; close BATCH;
