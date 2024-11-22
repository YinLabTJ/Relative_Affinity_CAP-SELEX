#!/usr/bin/perl
#Running spacek40 for local_max information
use strict;
use File::Spec;
my $cwd=File::Spec->rel2abs(__FILE__);
my ($vol,$shdir,$file)=File::Spec->splitpath($cwd);

die "Usage: perl step1_spacek40.pl batch_file seq_files_dir\n\tnames of sequence files should be consistent with Batch_file\n" unless(@ARGV==2);

my $batch_file=shift; my $seq_file_dir=shift; 
my $monomer_batch_file="$shdir/../data/Curated_Prey_Final.txt"; #information of individual TFs

my %barcode; my %barcodes; my %tf;

#Loading HT-SELEX information
open IN,"$monomer_batch_file";
<IN>;
while(<IN>){
	chomp;
	my @t=split /\t/,$_;
	if($t[0]=~/_/){
		my @t2=split /_/,$t[0];
		$barcode{$t2[0]}=$t[1]."\t".$t[2]."\t".$t[3];
		$barcode{$t2[1]}=$t[4]."\t".$t[5]."\t".$t[6];
	}else{
		$barcode{$t[0]}=$t[1]."\t".$t[2]."\t".$t[3];
	}
}
close IN;

#Loading CAP-SELEX batch file
open IN,"$batch_file";
<IN>;
while(<IN>){
	chomp;
	my @t=split;
	next unless($t[0]=~/_/);
	my $pair=$t[0];
	my @t2=split /_/,$t[0];
	$barcodes{$pair}.="," if(defined($barcodes{$pair}));
	$barcodes{$pair}.=$t[2]."\t".$t[3]."\t".$t[6];
	$tf{$t2[0]}=1;
	$tf{$t2[1]}=1;
}
close IN;

`mkdir -p spacek40_out`;
open OUT,">spacek40_out/running_spacek40.sh";

foreach my $key(sort keys %barcodes){
	my @t=split /_/,$key;
	if(defined($tf{$t[0]}) && defined($tf{$t[1]})){
		my @pair_info=split /,/,$barcodes{$key};
		foreach my $pair_info(@pair_info){
			my @col=split /\t/,$pair_info;
			my $treat_seq; my $control_seq;
			if($col[2]=~/b/){
				my @c=split /b/,$col[2];
				$treat_seq="$seq_file_dir/$col[0]$col[1]$c[0]u_sig.seq";
				if($c[1]==0){
					$control_seq="$seq_file_dir/$col[0]0u_sig.seq";
				}else{
					$control_seq="$seq_file_dir/$col[0]$col[1]$c[1]u_sig.seq";
				}
			}else{
				$treat_seq="$seq_file_dir/$col[0]$col[1]$col[2]u_sig.seq";
				$control_seq="$seq_file_dir/$col[0]0u_sig.seq";
			}
			die "Treatment file of $key doesn't exist\n" unless(-e $treat_seq);
			die "Control file of $key doesn't exist\n" unless(-e $control_seq);
			print OUT "spacek40 -40N -nogaps $control_seq $treat_seq 10 10 1 > spacek40_out/$key\_$col[0]$col[1].out\n";
		}
	}
}
close OUT;

`bash spacek40_out/running_spacek40.sh`;
`rm spacek40_out/running_spacek40.sh`;


