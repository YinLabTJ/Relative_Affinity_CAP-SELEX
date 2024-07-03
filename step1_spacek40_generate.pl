use strict;
my $seq_file_dir=shift;
`mkdir -p spacek40_out`;
open OUT,">step1_spacek40_sample.sh";
open IN,"single.list";
while(<IN>){
	chomp;
	my @t=split;
	print OUT "spacek40 -40N -nogaps $seq_file_dir/$t[1]0u_sig.seq /data/workdata/SELEX/$t[1]$t[2]$t[3]_sig.seq 8 8 1 > spacek40_out/$t[0].out\n";
}
close IN;
close OUT;
open OUT,">step1_spacek40_pair.sh";
open IN,"pair.list";
while(<IN>){
	chomp;
	my @t=split;
	print OUT "spacek40 -40N -nogaps $seq_file_dir/$t[1]0u_sig.seq /data/workdata/SELEX/$t[1]$t[2]3u_sig.seq 10 10 1 > spacek40_out/$t[0]\_$t[1].out\n";
}
close IN;
close OUT;


