use strict;
my $batchfile=shift; my $filter_pair=shift; my $seq_file_dir=shift;
open IN,"$filter_pair";
open OUT,">step2.sh";
while(<IN>){
	chomp;
	my @t=split;
	next if(@t>7 && $t[7]=~/Type/);
	print OUT "perl step2.work_3b0.pl $batchfile $t[0] $t[1] 10 $seq_file_dir\n";
}
close IN; close OUT;
