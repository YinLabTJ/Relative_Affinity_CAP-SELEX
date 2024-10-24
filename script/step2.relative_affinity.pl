#!/usr/bin/perl
#computing relative affinity
use strict;
use File::Spec;

die "Usage: perl step2.relative_affinity.pl batch_file TF_pair_and_barcode k_length seq_file_dir monomer_batch_file filter_list_file\n" unless(@ARGV==6);

my $cwd=File::Spec->rel2abs(__FILE__);
my ($vol,$dir,$file)=File::Spec->splitpath($cwd);

my $batch_file=shift; my $pair_barcode=shift; my $k=shift; my $seq_file_dir=shift; my $monomer_batch_file=shift; my $filter_list=shift;
my %single_c; my %single_t; my %pair_c; my %pair_t;

my @p=split /_/,$pair_barcode;
die "TF_pair $pair_barcode  should follow such format: TF1_TF2_barcode\n" unless(@p==3);

open IN,"$filter_list";
while(<IN>){
	chomp;
	my @t=split;
	if("$t[0]_$t[1]" eq $pair_barcode && $t[7]=~/Type/){
		die "Individual motifs of TFs: $t[0] are too similiar.\n";
	}
}
close IN; close OUT;


`mkdir -p output/$pair_barcode`;

#create all k-mer combination for fifth-order markov
my %allkmer;
my $max=4**$k;
for(my $i=0;$i<$max;$i++){
	my $seq="";
	my $dividend=$i;
        for(my $j=1;$j<=$k;$j++){
        	my $divisor=4**$j;
                my $remainder=$dividend % $divisor;
                if($remainder/(4**($j-1))==0){
                	$seq.="A";
                }elsif($remainder/(4**($j-1))==1){
                        $seq.="C";
                }elsif($remainder/(4**($j-1))==2){
                        $seq.="G";
                }elsif($remainder/(4**($j-1))==3){
                        $seq.="T";
                }else{
                        print "$remainder=$dividend % $divisor\n";
                }
                $dividend=$dividend-$remainder;
	}
        $allkmer{$seq}=1;
}

open OUT,">output/$pair_barcode/relative_affinity_xyplot.input";
print OUT "CAP_SELEX\tHT_SELEX\tClass\n";

open IN,"$monomer_batch_file";
while(<IN>){
	chomp;
	my @t=split;
	my @t2=split /_/,$t[0];
	$single_c{$t2[0]}="$seq_file_dir/$t[1]0_sig.seq";
	$single_t{$t2[0]}="$seq_file_dir/$t[1]$t[2]$t[3]_sig.seq";
	if(@t2>1){
		$single_c{$t2[1]}="$seq_file_dir/$t[4]0_sig.seq";
		$single_t{$t2[1]}="$seq_file_dir/$t[4]$t[5]$t[6]_sig.seq";
	}
}
close IN;
my $top_tf1tf2; my $top_tf1; my $top_tf2;
my %distance_to_x_equal_y2; my %distance_to_x_equal_y; my %higher_monomer;
my %pair_rel_aff; my %tf1_rel_aff; my %tf2_rel_aff;
my $spacek_out;
open IN,"$batch_file";
<IN>; #headline
while(<IN>){
	chomp;
	my @t=split;
	my $combine=$t[0]."_".$t[2];
	my $cycle_num;
	if($t[6]=~/b/){
		my @c=split /b/,$t[6];
		$c[1]=~s/u//;
		$cycle_num=$c[0]-$c[1];
	}else{
		$cycle_num=1;
	}

	next unless($combine eq $pair_barcode);
	$pair_c{$combine}="$seq_file_dir/$t[2]0_sig.seq";
	$pair_t{$combine}="$seq_file_dir/$t[2]$t[3]3u_sig.seq";
	$spacek_out=$pair_barcode.$t[3].".out";

	my @tf=split /_/,$t[0];
	print "Not enough TFs in $t[0]\n" if(@tf<2);
	
	#compute the frequency of the most abundant sequence in cycle r 
	$top_tf1tf2=&get_top_fc($pair_c{$combine},$pair_t{$combine},$k,1,%allkmer);
	$top_tf1=&get_top_fc($single_c{$tf[0]},$single_t{$tf[0]},$k,1,%allkmer);
	$top_tf2=&get_top_fc($single_c{$tf[1]},$single_t{$tf[1]},$k,1,%allkmer);

	
	#kmer frequncey of cycle 0 is simulted using 5th-order markov because of high diversity
	my ($tf1_top_c,%tf1_c)=&fifth_order_markov($single_c{$tf[0]},$k,$top_tf1,%allkmer);
	my ($tf1_top_t,%tf1_t)=&kmer_in_seq($single_t{$tf[0]},$k,$top_tf1,%allkmer);
	my ($tf2_top_c,%tf2_c)=&fifth_order_markov($single_c{$tf[1]},$k,$top_tf2,%allkmer);
        my ($tf2_top_t,%tf2_t)=&kmer_in_seq($single_t{$tf[1]},$k,$top_tf2,%allkmer);
	
	my ($tf1tf2_top_c,%tf1tf2_c)=&fifth_order_markov($pair_c{$combine},$k,$top_tf1tf2,%allkmer);
        my ($tf1tf2_top_t,%tf1tf2_t)=&kmer_in_seq($pair_t{$combine},$k,$top_tf1tf2,%allkmer);
	foreach my $key(keys %allkmer){
		$tf1tf2_t{$key}=0 if(!defined($tf1tf2_t{$key}));
		$tf1_t{$key}=0 if(!defined($tf1_t{$key}));
		$tf2_t{$key}=0 if(!defined($tf2_t{$key}));
		$tf1tf2_top_c=1	if($tf1tf2_top_c==0);
		$tf1_top_c=1 if($tf1_top_c==0);
		$tf2_top_c=1 if($tf2_top_c==0);
		$tf1tf2_c{$key}=1 if($tf1tf2_c{$key}==0);
		$tf1_c{$key}=1 if($tf1_c{$key}==0);
		$tf2_c{$key}=1 if($tf2_c{$key}==0);

		
		#relative affinity
		$pair_rel_aff{$key}=(($tf1tf2_t{$key}/$tf1tf2_top_t)/($tf1tf2_c{$key}/$tf1tf2_top_c))**(1/$cycle_num);
		$tf1_rel_aff{$key}=(($tf1_t{$key}/$tf1_top_t)/($tf1_c{$key}/$tf1_top_c))**(1/$cycle_num);
		$tf2_rel_aff{$key}=(($tf2_t{$key}/$tf2_top_t)/($tf2_c{$key}/$tf2_top_c))**(1/$cycle_num);


		if($tf1_rel_aff{$key}>$tf2_rel_aff{$key}){
			$distance_to_x_equal_y2{$key}=($pair_rel_aff{$key}-$tf1_rel_aff{$key})/(2**0.5);
			$higher_monomer{$key}="TF1";
			$distance_to_x_equal_y{$key}=$distance_to_x_equal_y2{$key} if($pair_rel_aff{$key}>=0.5);
		}else{
			$distance_to_x_equal_y2{$key}=($pair_rel_aff{$key}-$tf2_rel_aff{$key})/(2**0.5);
			$higher_monomer{$key}="TF2";
			$distance_to_x_equal_y{$key}=$distance_to_x_equal_y2{$key} if($pair_rel_aff{$key}>=0.5);
		}
	}
}
close IN;

my %top100;
my $top_count=0;
foreach my $key(sort{$distance_to_x_equal_y2{$b}<=>$distance_to_x_equal_y2{$a}} keys %distance_to_x_equal_y2){
	my $revcom=reverse($key);
	$revcom=~tr/ATCG/TAGC/;
	unless(defined($top100{$revcom})){
		$top_count++;
	}
	$top100{$key}=1;
	$top100{$revcom}=1;
	last if($top_count>=100);
}


open IN,"spacek40_out/$spacek_out";
my %local_max;
while(<IN>){
	chomp;
	if(/local_max/){
		my @t=split;
		$local_max{$t[5]}=1;
	}
}
close IN;

my $higher_monomer_ra;
foreach my $key(sort keys %pair_rel_aff){
	if($higher_monomer{$key} eq "TF1"){
		$higher_monomer_ra=$tf1_rel_aff{$key};
	}else{
		$higher_monomer_ra=$tf2_rel_aff{$key};
	}
	if(defined($top100{$key})){
		print OUT "$pair_rel_aff{$key}\t$higher_monomer_ra\tTop100\n";
	}else{
		print OUT "$pair_rel_aff{$key}\t$higher_monomer_ra\t$higher_monomer{$key}\n";
	}
	if(defined($local_max{$key})){
		print OUT "$pair_rel_aff{$key}\t$higher_monomer_ra\tlocal_max\n";
	}	
}

close OUT;

`Rscript $dir/CAP-SELEX_vs_HT-SELEX.R $pair_barcode`;
`rm Rplots.pdf` if(-e "Rplots.pdf");


############################ subroutines ###################

sub kmer_in_seq(){
	my ($seqfile,$k,$top_in_tf,%fc_top100)=@_;
	open SEQ,"$seqfile";
	my %tf_top100; my $tf_top=0;
	while(<SEQ>){
		chomp;
		my $len=length($_);
                my $i;
		for($i=0;$i<=$len-$k;$i++){
                        my $kmer=substr($_,$i,$k);
			my $kmer2=reverse($kmer); $kmer2=~tr/ATCG/TAGC/;
                        if($kmer eq $top_in_tf){
                                $tf_top++;
                        }
			if($kmer2 eq $top_in_tf){
				$tf_top++;
			}
                        if(defined($fc_top100{$kmer})){
                                $tf_top100{$kmer}++;
                        }
			if(defined($fc_top100{$kmer2})){
				$tf_top100{$kmer2}++;
			}
                }
        }
        close SEQ;
	return ($tf_top,%tf_top100);
}


sub fifth_order_markov(){
        my ($seqfile,$k,$top_in_tf,%allkmer)=@_;
        open SEQ,"$seqfile";
        my %sixmer; my %double; my %first; my %markov; my %exp; my %predict;
        while(<SEQ>){
                chomp;
                my $len=length($_);
                for(my $i=0;$i<=$len-6;$i++){
                        my $subseq=substr($_,$i,6);
                        $sixmer{$subseq}++;
                }
        }
        close SEQ;
        foreach my $key(sort keys %sixmer){
                my $seq1=substr($key,0,5);
                my $seq2=substr($key,1,5);
                my $pair=$seq1."_".$seq2;
                $double{$pair}+=$sixmer{$key};
                $first{$seq1}+=$sixmer{$key};
        }
        foreach my $key(sort keys %double){
                my @t=split /_/,$key;
                $markov{$key}=$double{$key}/$first{$t[0]};
        }
        open SEQ,"$seqfile";
        while(<SEQ>){
                chomp;
                my @t=split;
                my $len=length($_);
                for(my $i=0;$i<=40-$k;$i++){
                        my $subseq=substr($_,$i,5);
                        my $subseq=substr($_,$i,$k);
                        $exp{$subseq}++;
                }

        }
        
	my $tf_top=0; my %tf_top100;
	foreach my $key(keys %allkmer){
                my $start=substr($key,0,5);
                $predict{$key}=$first{$start};
                for(my $i=0;$i<=$k-6;$i++){
                        my $seq1=substr($key,$i,5);
                        my $seq2=substr($key,$i+1,5);
                        my $markov_key=$seq1."_".$seq2;
                        $predict{$key}=$predict{$key}*$markov{$markov_key};
                }
		my $kmer2=reverse($key); $kmer2=~tr/ATCG/TAGC/;
                if($key eq $top_in_tf){
			$tf_top+=$predict{$key};
                }
                if($kmer2 eq $top_in_tf){
                        $tf_top+=$predict{$key};
                }
                $tf_top100{$key}+=$predict{$key};
		$tf_top100{$kmer2}+=$predict{$key};
        }
	close SEQ;
	return ($tf_top,%tf_top100);
}

sub get_top_fc(){
        my ($cycle0_file,$cycle3_file,$k,$RevCom,%allkmer)=@_;
        my %cycle0_count; my %cycle3_count; my %fc; my $top_fc_kmer;
	open SEQ,"$cycle0_file";
	my %sixmer; my %double; my %first; my %markov; my %exp; my %predict;
        while(<SEQ>){
                chomp;
                my $len=length($_);
                for(my $i=0;$i<=$len-6;$i++){
                        my $subseq=substr($_,$i,6);
                        $sixmer{$subseq}++;
                }
        }
        close SEQ;
        foreach my $key(sort keys %sixmer){
                my $seq1=substr($key,0,5);
                my $seq2=substr($key,1,5);
                my $pair=$seq1."_".$seq2;
                $double{$pair}+=$sixmer{$key};
                $first{$seq1}+=$sixmer{$key};
        }
        foreach my $key(sort keys %double){
                my @t=split /_/,$key;
                $markov{$key}=$double{$key}/$first{$t[0]};
        }
        open SEQ,"$cycle0_file";
        while(<SEQ>){
                chomp;
                my @t=split;
                my $len=length($_);
                for(my $i=0;$i<=40-$k;$i++){
                        my $subseq=substr($_,$i,5);
                        my $subseq=substr($_,$i,$k);
                        $exp{$subseq}++;
                }
        }
	my $tf_top=0;
        foreach my $key(keys %allkmer){
                my $start=substr($key,0,5);
                $predict{$key}=$first{$start};
                for(my $i=0;$i<=$k-6;$i++){
                        my $seq1=substr($key,$i,5);
                        my $seq2=substr($key,$i+1,5);
                        my $markov_key=$seq1."_".$seq2;
                        $predict{$key}=$predict{$key}*$markov{$markov_key};
                }
                my $kmer2=reverse($key); $kmer2=~tr/ATCG/TAGC/;
                $cycle0_count{$key}+=$predict{$key};
                $cycle0_count{$kmer2}+=$predict{$key};
        }
        close SEQ;


	
	#signal file
	open SEQ,"$cycle3_file";
	while(<SEQ>){
		chomp;
		my @t=split //,$_;
		for(my $i-0;$i<=40-$k;$i++){
			my $kmer=substr($_,$i,$k);
			$cycle3_count{$kmer}++;
			if($RevCom==1){
				my $rc_seq=$kmer;
				$rc_seq=reverse($rc_seq);
				$rc_seq=~tr/ATCG/TAGC/;
				$cycle3_count{$rc_seq}++;
			}
		}
	}
	close SEQ;
	foreach my $key(keys %cycle3_count){
		next if($key=~/N/);
		if(!defined($cycle0_count{$key}) || $cycle0_count{$key}==0){
			$fc{$key}=$cycle3_count{$key};
		}else{
			$fc{$key}=$cycle3_count{$key}/$cycle0_count{$key};
		}
	}
	foreach my $key(sort{$fc{$b}<=>$fc{$a}} keys %fc){
		$top_fc_kmer=$key;
		last;
	}
	return ($top_fc_kmer);
}
