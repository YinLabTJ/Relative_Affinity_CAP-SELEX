use strict;
my $batch_file=shift;
open IN,"Curated_Prey_Final.txt";
<IN>;
my %barcode; my %barcodes; my %single;
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
open IN,"$batch_file";
<IN>;
while(<IN>){
	chomp;
	my @t=split;
	next unless($t[0]=~/_/);
	my $pair=$t[0];
	my @t2=split /_/,$t[0];
	$barcodes{$pair}.="," if(defined($barcodes{$pair}));
	$barcodes{$pair}.=$t[2]."\t".$t[3];
	$single{$t2[0]}=1;
	$single{$t2[1]}=1;
}
close IN;
open OUT,">single.list";
foreach my $key(sort keys %barcode){
	print OUT "$key\t$barcode{$key}\n" if(defined$single{$key});
}
close OUT;
open OUT,">pair.list";
foreach my $key(sort keys %barcodes){
	my @t=split /,/,$barcodes{$key};
	foreach my $a(@t){
        	print OUT "$key\t$a\n";
	}
}
close OUT;
