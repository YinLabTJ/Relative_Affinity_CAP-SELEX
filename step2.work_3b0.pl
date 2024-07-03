use strict;
my $batch_file=shift; my $pair=shift; my $barcode=shift; my $k=shift; my $seq_file_dir=shift;
my %single0; my %single3; my %pair0; my %pair3;

`mkdir -p output`;
my $combine2=$pair."_".$barcode;
`mkdir -p output/$combine2`;

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

open IN,"Curated_Prey_Final.txt";
open OUT,">output/$combine2/relative_affinity.xls";
open OUT2,">output/$combine2/relative_affinity2.xls";
while(<IN>){
	chomp;
	my @t=split;
	my @t2=split /_/,$t[0];
	$single0{$t2[0]}="$seq_file_dir/$t[1]0_sig.seq";
	$single3{$t2[0]}="$seq_file_dir/$t[1]$t[2]3u_sig.seq";
	if(@t2>1){
		$single0{$t2[1]}="$seq_file_dir/$t[4]0_sig.seq";
		$single3{$t2[1]}="$seq_file_dir/$t[4]$t[5]3u_sig.seq";
	}
}
close IN;
my $top_tf1tf2; my $top_tf1; my $top_tf2;
open IN,"$batch_file";
<IN>; #tytle line
while(<IN>){
	chomp;
	my @t=split;
	my $combine1=$t[0]."_".$t[2];
	$pair0{$combine1}="$seq_file_dir/$t[2]0_sig.seq";
	$pair3{$combine1}="$seq_file_dir/$t[2]$t[3]3u_sig.seq";
	next unless($combine1 eq $combine2);
	$top_tf1tf2=&get_top_fc($pair0{$combine1},$pair3{$combine1},$k,1,%allkmer);

	my @te=split /_/,$t[0];
	print "Not enough TFs in $t[0]\n" if(@te<2);
	$top_tf1=&get_top_fc($single0{$te[0]},$single3{$te[0]},$k,1,%allkmer);
	
	$top_tf2=&get_top_fc($single0{$te[1]},$single3{$te[1]},$k,1,%allkmer);


	#TF1 cycle0
	print "$t[0] start\n";
	print "$pair0{$combine1}\n$pair3{$combine1}\n$single0{$te[0]}\n$single3{$te[0]}\n$single0{$te[1]}\n$single3{$te[1]}\n";
	my ($tf1_top0,%tf1_0)=&fifth_order_markov($single0{$te[0]},$k,$top_tf1,%allkmer);
	my ($tf1_top3,%tf1_3)=&kmer_in_seq($single3{$te[0]},$k,$top_tf1,%allkmer);
	my ($tf2_top0,%tf2_0)=&fifth_order_markov($single0{$te[1]},$k,$top_tf2,%allkmer);
        my ($tf2_top3,%tf2_3)=&kmer_in_seq($single3{$te[1]},$k,$top_tf2,%allkmer);
	
	my ($tf1tf2_top_0,%tf1tf2_0)=&fifth_order_markov($pair0{$combine1},$k,$top_tf1tf2,%allkmer);
        my ($tf1tf2_top_3,%tf1tf2_3)=&kmer_in_seq($pair3{$combine1},$k,$top_tf1tf2,%allkmer);
	foreach my $key(keys %allkmer){
		print OUT "$t[0]\t";
		print OUT2 "$t[0]\t";
		print OUT "$key\t$top_tf1\t";
		print OUT2 "$key\t$top_tf1tf2\t";
		if(!defined($tf1tf2_3{$key})){
			$tf1tf2_3{$key}=0;
		}
		if(!defined($tf1_3{$key})){
			$tf1_3{$key}=0;
		}
		if(!defined($tf2_3{$key})){
                        $tf2_3{$key}=0;
                }
		if($tf1tf2_top_0==0){
			$tf1tf2_top_0=1;
		}
		if($tf1_top0==0){
                        $tf1_top0=1;
                }
		if($tf2_top0==0){
                        $tf2_top0=1;
                }
		if($tf1tf2_0{$key}==0){
			$tf1tf2_0{$key}=1;
		}
		if($tf1_0{$key}==0){
                        $tf1_0{$key}=1;
                }
		if($tf2_0{$key}==0){
                        $tf2_0{$key}=1;
                }
		print OUT "$tf1tf2_3{$key}\t$tf1tf2_top_3\t$tf1tf2_0{$key}\t$tf1tf2_top_0\t";
		print OUT2 "$tf1tf2_3{$key}\t$tf1tf2_top_3\t$tf1tf2_0{$key}\t$tf1tf2_top_0\t";
		print OUT "$tf1_3{$key}\t$tf1_top3\t$tf1_0{$key}\t$tf1_top0\t";
		print OUT2 "$top_tf1\t$tf1_3{$key}\t$tf1_top3\t$tf1_0{$key}\t$tf1_top0\t";
		print OUT "$top_tf2\t";
		print OUT2 "$top_tf2\t";
		print OUT "$tf2_3{$key}\t$tf2_top3\t$tf2_0{$key}\t$tf2_top0\n";
		print OUT2 "$tf2_3{$key}\t$tf2_top3\t$tf2_0{$key}\t$tf2_top0\n";
	}
	print "$t[0] end\n";

}
close IN;
close OUT; close OUT2;

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
	 print "$top_fc_kmer\t$fc{$top_fc_kmer}\n";
	return ($top_fc_kmer);
}
