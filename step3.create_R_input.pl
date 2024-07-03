use strict;
my $batch_file=shift;
open BATCH,"$batch_file";
while(<BATCH>){
	chomp;
	my @b=split;
	next unless($b[0]=~/_/);
	my $tf_pair=$b[0]."_".$b[2];
	if(-e "output/$tf_pair/relative_affinity.xls"){
		open IN,"output/$tf_pair/relative_affinity.xls";
		open OUT,">output/$tf_pair/Cor_input.xls";
		open OUT2,">output/$tf_pair/Cor_input2.xls";
		open OUT3,">output/$tf_pair/Cor_input3.xls";
		print OUT "TF1_and_TF2\tTF1_or_TF2\n";
		print OUT2 "TF1_and_TF2\tTF1_or_TF2\n";
		print OUT3 "TF1_and_TF2\tTF1\tTF2\n";
		while(<IN>){
			chomp;
			my @t=split /\t/,$_;
			if($t[0] eq $b[0]){
				my $tf1tf2=(($t[3]/$t[4])/($t[5]/$t[6]))**(1/3);
				my $tf1=(($t[7]/$t[8])/($t[9]/$t[10]))**(1/3);
				my $tf2=(($t[12]/$t[13])/($t[14]/$t[15]))**(1/3);
				if($tf1>=$tf2){
					print OUT "$tf1tf2\t$tf1\n";
					print OUT2 "$tf1tf2\t$tf1\t$t[1]\tTF1\n";
					print OUT3 "$tf1tf2\t$tf1\t$tf2\t$t[1]\tTF1\n";
				}else{
					print OUT "$tf1tf2\t$tf2\n";
					print OUT2 "$tf1tf2\t$tf2\t$t[1]\tTF2\n";
					print OUT3 "$tf1tf2\t$tf1\t$tf2\t$t[1]\tTF2\n";
				}	
			}
		}
	}
}
close IN; close OUT; close OUT2; close OUT3; close BATCH;
