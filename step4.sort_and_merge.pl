use Cwd;
use strict;
my $pair=shift;
open IN,"output/$pair/Cor_input2.xls";
my %distance_to_x_equal_y;
my %distance_to_x_equal_y2;
<IN>;
while(<IN>){
	chomp;
	my @t=split;
	$distance_to_x_equal_y2{$t[2]}=($t[0]-$t[1])/(2**0.5);
	next if($t[0]<0.5);
	$distance_to_x_equal_y{$t[2]}=($t[0]-$t[1])/(2**0.5);
}
close IN;

my $count=0; my %furthest; my %pop;
open OUT,">output/$pair/cluster1.out";
foreach my $key(sort{$distance_to_x_equal_y{$b}<=>$distance_to_x_equal_y{$a}} keys %distance_to_x_equal_y){
	my $revcom=reverse($key);
	$revcom=~tr/ATCG/TAGC/;
	unless(defined($furthest{$revcom})){
		$count++;
		$furthest{$key}=$count;
		print OUT "$key\t$furthest{$key}\n";
	}
	last if($count>=30);
	#print "$key\t$distance_to_x_equal_y{$key}\n";
}
close OUT;

my %top100;
$count=0; 
foreach my $key(sort{$distance_to_x_equal_y2{$b}<=>$distance_to_x_equal_y2{$a}} keys %distance_to_x_equal_y2){
	my $revcom=reverse($key);
	$revcom=~tr/ATCG/TAGC/;
	unless(defined($top100{$revcom})){
                $count++;
	}
	$top100{$key}=1;
	$top100{$revcom}=1;
	last if($count>=100);
}

`grep local_max spacek40_out/$pair.out | cut -f6 > output/$pair/tmp`;
my %local_max;
open IN,"output/$pair/tmp";
while(<IN>){
	chomp;
	$local_max{$_}=1;
}
`rm output/$pair/tmp`;
close IN;


open IN,"output/$pair/Cor_input2.xls";
open OUT,">output/$pair/Cor_input_color.xls";
<IN>;
print OUT "CAP_SELEX\tHT_SELEX\tClass\n";
while(<IN>){
	chomp;
	my @t=split;
	if(defined($top100{$t[2]})){
		print OUT "$t[0]\t$t[1]\tTop100\n";
	}else{
		print OUT "$t[0]\t$t[1]\t$t[3]\n";
	}
	if(defined($local_max{$t[2]})){
		print OUT "$t[0]\t$t[1]\tlocal_max\n";
	}
}
close IN; close OUT;
