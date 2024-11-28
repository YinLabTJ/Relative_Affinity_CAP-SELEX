#include <time.h> 
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#define MAX(a,b) (((a)>(b))?(a):(b))

char *VERSION = "spacek40 v0.163 NOV 31 2014";
/* multinomial 2 bug fixed */
/* from 0.144 lower memory use for kmer counting */

char *COMMAND;

/* GLOBAL VARIABLES */
float pseudocount = 0.0001;
short int Nlength = 20;
short int max_Nlength = 40;
short int max_width_of_pwm = 83; /* Nlength * 2 + 1 */
long int max_number_of_sequences = 10000000;
short int number_of_files = 2;
double pvalue_cache[1001][1001];
short int head_to_tail = 0;
short int head_to_head = 1;
short int tail_to_tail = 2;
short int print_local_max = 0;
short int align_matches = 0;
short int pwm_align = 0;
short int contacts = 0;
short int rna = 0;
double local_max_min_percent = 0.1;
short int nocall = 0;
short int methylCGcompare = 0;

__uint128_t mask_ULL[42][42]; /* PRIMARY mask_ULL FOR EACH SEPARATE NUCLEOTIDE, VALUES GIVEN IN MAIN PROGRAM */
__uint128_t lowmask_ULL[42];  /* LOW MASK FOR EACH KMER, VALUES GIVEN IN MAIN PROGRAM */
__uint128_t highmask_ULL[42]; /* HIGH MASK FOR EACH KMER, VALUES GIVEN IN MAIN PROGRAM */

char *tf_kmers[] = {"AATCAA","ACATGT","ACCGCA","ACCGGA","AGATAA","AGGTCA","ATCGAT","CAACAG","CACCTG","CACGCA","CACGTC","CACGTG","CACTTA","CAGCTG","CATAAA","CATATG","CATTCC","CCATAT","CCATTA","CCCGCC","CCCGGA","CCGGAT","CCGTTA","CGAAAC","CGTAAA","CTAGTG","CTGTCA","GAAACC","GAACAA","GACCAC","GACGTC","GAGGAA","GCCACG","GGCAAC","GGCGCC","GGGGAA","GGTACA","GGTGTG","GTCACG","TAAACA","TAATTA","TACGTA","TATGCA","TGACAG","TGCATA","TGCCAA","TGCGCA","TGCGGG","TGCTGA","TGGAAA","TTCTAG","ATGCCC","GTCCGC","GTGAAA","CCGCCA","TCGCGA","CGACCA","CGCTGT","ACCCAC","ACCGGT","CCATGG","ATCAAA","AACGAT","TTCGAA","AAATAG","TCTAGA","TGCCCT","CACGCC","GATGCA","TGACTC","TGAGTC","TGACAC","TCCCCA","TAAACG","TAATTG","CAATAA","ACATGA","CAAGGT","GTCCAA","AAGTCA","AGTTCA", "END"};
char *tf_names[] = {"PBX","P53","RUNX","ETSI","GATA","NucRes1","CUT","SCRT","TCF4","EGR","CREB3","EboxHLH","NKX","NHLH","HOX9to12","bHLHcat","TEAD","SRF","PITX","HighE2F","ETSII","SPDEF","MYB","IRF","HOX13","CTCF","MEIS","IRF2","SOX","GLI","cre_bZIP","ETSIII","SP_KLF","RFX","LowE2F","NFKB1","AR","Tbox","PAX","FOX","Homeo","TEF_GMEB","POU","MEIS","CUT2","NFI","CEBP","GCM","MAF","NFAT","HSF","HIC","HINFP","PRD_IRF","YY","ZBED","ZBTB7","ZIC","ZNF143","GRHL","EBF","TCF7_LEF","HOMEZ","HSFY","MEF","SMAD3","TFAP","SREBF","CRE_CEBP","prebZIP1","prebZIP2","oddFox_a","FOXO","BARHLa","BARHLb","HOX9to12b","IRX","ESRR","HNF4","NRE","VDR", "END"};
char *pwm_row_ids[] = {"Top_Strand_Base_Count_A","Top_Strand_Base_Count_C","Top_Strand_Base_Count_G","Top_Strand_Base_Count_T","Top_Strand_Maj_Groove_Phosphate", "Top_Strand_Min_Groove_Phosphate", "Top_Strand_Maj_Groove_Sugar", "Top_Strand_Min_Groove_Sugar", "Top_Strand_Maj_Groove_Base", "Top_Strand_Min_Groove_Base", "Bottom_Strand_Min_Groove_Base", "Bottom_Strand_Maj_Groove_Base", "Bottom_Strand_Min_Groove_Sugar", "Bottom_Strand_Maj_Groove_Sugar", "Bottom_Strand_Min_Groove_Phosphate", "Bottom_Strand_Maj_Groove_Phosphate"};
long int *tf_kmer_values;
short int number_of_tf_kmers;

char *yesnocall[] = {"No","Yes", "ND"};
char *dna_iupac =    "ACGTRYMKWSBDHVN";
char *dna_rc_iupac = "TGCAYRKMWSVHDBN";
char *dna_bitiupac= "?ACMGRSVTWYHKDBN";
char *rna_iupac =    "ACGURYMKWSBDHVN";
char *rna_rc_iupac = "UGCAYRKMWSVHDBN";
char *rna_bitiupac= "?ACMGRSVUWYHKDBN";

char *dnaforward = "ACGTN";
char *rnaforward = "ACGUN";
char *dnareverse = "TGCAN";
char *rnareverse = "UGCAN";
char *forward_lc = "acgtn";
short int iupac_length = 15;
long int **flank_kmer_count;

/* GLOBAL FLAGS */ 
short int count_both_instances_of_palindromic_hit = 0;
short int count_unequal_hits_only = 0;
short int count_only_forward_instance_of_palindromic_hit = 0;
short int count_only_reverse_instance_of_palindromic_hit = 0;
short int prefer_forward_strand = 0;
short int prefer_reverse_strand = 0;
short int last_count_was_forward = 0;
short int noname = 0;
short int paths = 0;
short int barcodelogo = 0;
short int gray_bars = 0;
short int barcodelogolabels = 0;
double scale_bars = 0;
short int max_scale_bar = 0;
short int circles = 1;
signed short int flank_kmer_pos = -100;


/* STRUCTURE DECLARATIONS AND INITIALIZATION SUBROUTINES */
struct flags {short int extendedoutput; short int remove_non_unique; short int print_counts; short int print_p_values; short int print_input_sequences; 
short int output_all_gap_lengths_in_one_line; short int print_nucleotides; short int print_frequencies; short int count_also_spaced_kmers; 
short int only_palindromes; short int dinucleotide_properties; short int kmer_table; short int information_content_output;
short int even_background; short int complex_background; short int flank_with_pwm; short int flank; short int kmer_count;};

/* RGB COLOR */
struct rgb_color {short int red; short int green; short int blue;};

/* ALIGNMENT SCORE */
struct alignscore {long int **score; long int **count; long int **direct_repeat; long int **inverted_repeat;};
short int alignscore_structure_init (struct alignscore *a)
{
short int counter;
short int strand;
short int gap;
(*a).score = malloc(sizeof(long int *) * (Nlength * 2) + 5);
(*a).count = malloc(sizeof(long int *) * (Nlength * 2) + 5);
(*a).direct_repeat = malloc(sizeof(long int *) * Nlength + 5);
(*a).inverted_repeat = malloc(sizeof(long int *) * Nlength + 5);
for(counter = 0; counter < Nlength * 2; counter++)
{
(*a).score[counter] = malloc(sizeof(long int) * 3 + 5);
(*a).count[counter] = malloc(sizeof(long int) * 3 + 5);
for(strand = 0; strand < 2; strand++)
{
(*a).score[counter][strand] = 0;
(*a).count[counter][strand] = 0;
}
}
for(counter = 0; counter < Nlength; counter++)
{
(*a).direct_repeat[counter] = malloc(sizeof(long int) * Nlength + 5);
(*a).inverted_repeat[counter] = malloc(sizeof(long int) * Nlength + 5);
for(gap = 0; gap < Nlength; gap++)
{
(*a).direct_repeat[counter][gap] = 0;
(*a).inverted_repeat[counter][gap] = 0;
}
}
}

/* ORIENTED MATCH */
struct oriented_match {short int position; signed short int strand; double score; short int id;};
short int oriented_match_init (struct oriented_match *i)
{  
    (*i).position = 0;
    (*i).strand = 0;
    (*i).score = 0;
    (*i).id = 0;
    return(0);
}

/* BIT REPRESENTATION OF IUPAC */
struct bitiupac_structure {short int length; __uint128_t *base; __uint128_t sequence_value_ULL;};
short int bitiupac_structure_init (struct bitiupac_structure *s)
{
(*s).base = malloc(sizeof(__uint128_t)*4+5);
(*s).base[0] = 0;
(*s).base[1] = 0;
(*s).base[2] = 0;
(*s).base[3] = 0;
(*s).sequence_value_ULL = 0;
(*s).length = 0;
}

/* STATISTICS FOR KMER-MONO CORRELATIONS */
struct similarity_stats {double correlation; long int sum_total_deviation; long int sum_relative_deviation;};

struct sumtable {__uint128_t *sums; long int *counts; long int *max_min_counts; long int *max_max_counts; __uint128_t *max_max_kmer; __uint128_t *max_min_kmer; long int ***cloud_counts;};
short int sumtable_init(struct sumtable *s, short int kmer_length)
{
    (*s).sums = malloc(sizeof(__uint128_t)*20+5);
    (*s).counts = malloc(sizeof(long int)*20+5);    
    (*s).max_min_counts = malloc(sizeof(long int)*20+5); 
    (*s).max_min_kmer = malloc(sizeof(__uint128_t)*20+5); 
    (*s).max_max_counts = malloc(sizeof(long int)*20+5); 
    (*s).max_max_kmer = malloc(sizeof(__uint128_t)*20+5);
    (*s).cloud_counts = malloc(sizeof(long int *)*kmer_length+5);      
    short int counter;
    long int counter2;
    for (counter = 0; counter < 20; counter++)
    {
    (*s).sums[counter] = 0;   
    (*s).counts[counter] = 0;
    (*s).max_max_counts[counter] = 0;
    (*s).max_max_kmer[counter] = 0;
    (*s).max_min_counts[counter] = 0;
    (*s).max_min_kmer[counter] = 0;
    }
    for (counter = 0; counter < kmer_length; counter++)
    {
    (*s).cloud_counts[counter] = malloc(sizeof(long int *)*pow(4,kmer_length)+5);  
    for (counter2 = 0; counter2 <  pow(4,kmer_length); counter2++) 
    {
    (*s).cloud_counts[counter][counter2] = malloc(sizeof(long int) * 6 + 5);
    (*s).cloud_counts[counter][counter2][0] = 0;
    (*s).cloud_counts[counter][counter2][1] = 0;
    (*s).cloud_counts[counter][counter2][2] = 0;
    (*s).cloud_counts[counter][counter2][3] = 0;
    (*s).cloud_counts[counter][counter2][4] = 0;
    (*s).cloud_counts[counter][counter2][5] = 0;
    }
    }
}

/* ADJACENT DINUC MARKOV MODEL */
struct adjacent_dinucleotide_model {char *name; short int width; double **fraction; double **mononuc_fraction;};
short int adjacent_dinucleotide_model_init(struct adjacent_dinucleotide_model *d, char *name, short int width)
{
short int dinucleotide;
short int first;
short int second;
(*d).width = width+1;
(*d).name = malloc(1000);
strcpy((*d).name, name);
(*d).fraction = malloc(sizeof(double *) * 16 + 5);
(*d).mononuc_fraction = malloc(sizeof(double *) * 4 + 5);
for (first = 0; first < 16; first++)
{
(*d).fraction[first] = malloc(sizeof(double) * (width+1) + 5);
for (second = 0; second <= width; second++) (*d).fraction[first][second] = 0;
}
for (first = 0; first < 4; first++)
{
(*d).mononuc_fraction[first] = malloc(sizeof(double) * (width+1) + 5);
for (second = 0; second <= width; second++) (*d).mononuc_fraction[first][second] = 0;
}
return(0);    
}

/* BASE DEPENDENCY MATRIX */
struct base_dependency_matrix {char *name; short int width; long int ***incidence; double **total_relative_deviation; double **count_statistic_expected_total_relative_deviation; double **information_content; double **eo_correlation; double **permutated_correlation; double max_expected_relative_deviation; double max_relative_deviation;};
short int base_dependency_matrix_init(struct base_dependency_matrix *m, char *name, short int width)
{
short int dinucleotide;
short int first;
short int second;
(*m).width = width;
(*m).max_relative_deviation = 0;
(*m).max_expected_relative_deviation = 0;
(*m).name = malloc(1000);
strcpy((*m).name, name);
(*m).incidence = malloc(sizeof(long int *) * width + 5);
(*m).information_content = malloc(sizeof(double *) * width + 5);
(*m).total_relative_deviation = malloc(sizeof(double *) * width + 5);
(*m).count_statistic_expected_total_relative_deviation = malloc(sizeof(double *) * width + 5);
(*m).eo_correlation = malloc(sizeof(double *) * width + 5);
(*m).permutated_correlation = malloc(sizeof(double *) * width + 5);
for (first = 0; first < width; first++)
{
(*m).information_content[first] = malloc(sizeof(double) * width + 5);
(*m).total_relative_deviation[first] = malloc(sizeof(double) * width + 5);
(*m).count_statistic_expected_total_relative_deviation[first] = malloc(sizeof(double) * width + 5);
(*m).eo_correlation[first] = malloc(sizeof(double) * width + 5);
(*m).permutated_correlation[first] = malloc(sizeof(double) * width + 5);
(*m).incidence[first] = malloc(sizeof(long int *) * width + 5);
for (second = 0; second < width; second++)
{
(*m).incidence[first][second] = malloc(sizeof(long int) * 16 + 5);
(*m).information_content[first][second] = 0;
(*m).total_relative_deviation[first][second] = 0;
(*m).count_statistic_expected_total_relative_deviation[first][second] = 0;
(*m).eo_correlation[first][second] = 0;
(*m).permutated_correlation[first][second] = 0;
for (dinucleotide = 0; dinucleotide < 16; dinucleotide++) (*m).incidence[first][second][dinucleotide] = 0;
}
}
return(0);
}

/* DINUCLEOTIDE PROPERTIES */
struct dinucleotide_properties {short int number_of_dinucleotide_properties; char **dinucleotide_property_string; double **dinucleotide_property;};
short int dinucleotide_properties_init(struct dinucleotide_properties *d)
{
(*d).number_of_dinucleotide_properties = 16;
short int AA = 0;
short int AC = 1;
short int AG = 2;
short int AT = 3;
short int CA = 4;
short int CC = 5;
short int CG = 6;
short int CT = 7;
short int GA = 8;
short int GC = 9;
short int GG = 10;
short int GT = 11;
short int TA = 12;
short int TC = 13;
short int TG = 14;
short int TT = 15;
short int PMID = 16;

(*d).dinucleotide_property_string = malloc(((*d).number_of_dinucleotide_properties + 1) * sizeof(char *) + 5);
(*d).dinucleotide_property = malloc(((*d).number_of_dinucleotide_properties + 1) * sizeof(double *) + 5);
short int counter;
short int counter2;
    
for (counter = 0; counter <= (*d).number_of_dinucleotide_properties; counter++) 
{
(*d).dinucleotide_property_string[counter] = malloc(200);
(*d).dinucleotide_property[counter] = malloc(20 * sizeof(double) + 5);
for(counter2 = 0; counter2 < 17; counter2++) (*d).dinucleotide_property[counter][counter2] = -1.1111;
}

short int Twist = 0;
strcpy((*d).dinucleotide_property_string[Twist], "Twist");
(*d).dinucleotide_property[Twist][AA]=38.9;
(*d).dinucleotide_property[Twist][AC]=31.12;
(*d).dinucleotide_property[Twist][AG]=32.15;
(*d).dinucleotide_property[Twist][AT]=33.81;
(*d).dinucleotide_property[Twist][CA]=41.41;
(*d).dinucleotide_property[Twist][CC]=34.96;
(*d).dinucleotide_property[Twist][CG]=32.91;
(*d).dinucleotide_property[Twist][CT]=32.15;
(*d).dinucleotide_property[Twist][GA]=41.31;
(*d).dinucleotide_property[Twist][GC]=38.5;
(*d).dinucleotide_property[Twist][GG]=34.96;
(*d).dinucleotide_property[Twist][GT]=31.12;
(*d).dinucleotide_property[Twist][TA]=33.28;
(*d).dinucleotide_property[Twist][TC]=41.31;
(*d).dinucleotide_property[Twist][TG]=41.41;
(*d).dinucleotide_property[Twist][TT]=38.9;
(*d).dinucleotide_property[Twist][PMID]=8996793;

short int Stacking_energy = 1;
strcpy((*d).dinucleotide_property_string[Stacking_energy], "Stacking energy");
(*d).dinucleotide_property[Stacking_energy][AA]=-12;
(*d).dinucleotide_property[Stacking_energy][AC]=-11.8;
(*d).dinucleotide_property[Stacking_energy][AG]=-11.5;
(*d).dinucleotide_property[Stacking_energy][AT]=-10.6;
(*d).dinucleotide_property[Stacking_energy][CA]=-12.3;
(*d).dinucleotide_property[Stacking_energy][CC]=-9.5;
(*d).dinucleotide_property[Stacking_energy][CG]=-13.1;
(*d).dinucleotide_property[Stacking_energy][CT]=-11.5;
(*d).dinucleotide_property[Stacking_energy][GA]=-11.4;
(*d).dinucleotide_property[Stacking_energy][GC]=-13.2;
(*d).dinucleotide_property[Stacking_energy][GG]=-9.5;
(*d).dinucleotide_property[Stacking_energy][GT]=-11.8;
(*d).dinucleotide_property[Stacking_energy][TA]=-11.2;
(*d).dinucleotide_property[Stacking_energy][TC]=-11.4;
(*d).dinucleotide_property[Stacking_energy][TG]=-12.3;
(*d).dinucleotide_property[Stacking_energy][TT]=-12;
(*d).dinucleotide_property[Stacking_energy][PMID]=9199773;

short int Rise = 2;
strcpy((*d).dinucleotide_property_string[Rise], "Rise");
(*d).dinucleotide_property[Rise][AA]=3.16;
(*d).dinucleotide_property[Rise][AC]=3.41;
(*d).dinucleotide_property[Rise][AG]=3.63;
(*d).dinucleotide_property[Rise][AT]=3.89;
(*d).dinucleotide_property[Rise][CA]=3.23;
(*d).dinucleotide_property[Rise][CC]=4.08;
(*d).dinucleotide_property[Rise][CG]=3.6;
(*d).dinucleotide_property[Rise][CT]=3.63;
(*d).dinucleotide_property[Rise][GA]=3.47;
(*d).dinucleotide_property[Rise][GC]=3.81;
(*d).dinucleotide_property[Rise][GG]=4.08;
(*d).dinucleotide_property[Rise][GT]=3.41;
(*d).dinucleotide_property[Rise][TA]=3.21;
(*d).dinucleotide_property[Rise][TC]=3.47;
(*d).dinucleotide_property[Rise][TG]=3.23;
(*d).dinucleotide_property[Rise][TT]=3.16;
(*d).dinucleotide_property[Rise][PMID]=8996793;

short int Bend = 3;
strcpy((*d).dinucleotide_property_string[Bend], "Bend");
(*d).dinucleotide_property[Bend][AA]=3.07;
(*d).dinucleotide_property[Bend][AC]=2.97;
(*d).dinucleotide_property[Bend][AG]=2.31;
(*d).dinucleotide_property[Bend][AT]=2.6;
(*d).dinucleotide_property[Bend][CA]=3.58;
(*d).dinucleotide_property[Bend][CC]=2.16;
(*d).dinucleotide_property[Bend][CG]=2.81;
(*d).dinucleotide_property[Bend][CT]=2.31;
(*d).dinucleotide_property[Bend][GA]=2.51;
(*d).dinucleotide_property[Bend][GC]=3.06;
(*d).dinucleotide_property[Bend][GG]=2.16;
(*d).dinucleotide_property[Bend][GT]=2.97;
(*d).dinucleotide_property[Bend][TA]=6.74;
(*d).dinucleotide_property[Bend][TC]=2.51;
(*d).dinucleotide_property[Bend][TG]=3.58;
(*d).dinucleotide_property[Bend][TT]=3.07;
(*d).dinucleotide_property[Bend][PMID]=8996793;

short int Tip = 4;
strcpy((*d).dinucleotide_property_string[Tip], "Tip");
(*d).dinucleotide_property[Tip][AA]=1.76;
(*d).dinucleotide_property[Tip][AC]=2;
(*d).dinucleotide_property[Tip][AG]=0.9;
(*d).dinucleotide_property[Tip][AT]=1.87;
(*d).dinucleotide_property[Tip][CA]=-1.64;
(*d).dinucleotide_property[Tip][CC]=0.71;
(*d).dinucleotide_property[Tip][CG]=0.22;
(*d).dinucleotide_property[Tip][CT]=0.9;
(*d).dinucleotide_property[Tip][GA]=1.35;
(*d).dinucleotide_property[Tip][GC]=2.5;
(*d).dinucleotide_property[Tip][GG]=0.71;
(*d).dinucleotide_property[Tip][GT]=2;
(*d).dinucleotide_property[Tip][TA]=6.7;
(*d).dinucleotide_property[Tip][TC]=1.35;
(*d).dinucleotide_property[Tip][TG]=-1.64;
(*d).dinucleotide_property[Tip][TT]=1.76;
(*d).dinucleotide_property[Tip][PMID]=8996793;

short int Inclination = 5;
strcpy((*d).dinucleotide_property_string[Inclination], "Inclination");
(*d).dinucleotide_property[Inclination][AA]=-1.43;
(*d).dinucleotide_property[Inclination][AC]=-0.11;
(*d).dinucleotide_property[Inclination][AG]=-0.92;
(*d).dinucleotide_property[Inclination][AT]=0;
(*d).dinucleotide_property[Inclination][CA]=1.31;
(*d).dinucleotide_property[Inclination][CC]=-1.11;
(*d).dinucleotide_property[Inclination][CG]=0;
(*d).dinucleotide_property[Inclination][CT]=0.92;
(*d).dinucleotide_property[Inclination][GA]=-0.33;
(*d).dinucleotide_property[Inclination][GC]=0;
(*d).dinucleotide_property[Inclination][GG]=1.11;
(*d).dinucleotide_property[Inclination][GT]=0.11;
(*d).dinucleotide_property[Inclination][TA]=0;
(*d).dinucleotide_property[Inclination][TC]=0.33;
(*d).dinucleotide_property[Inclination][TG]=-1.31;
(*d).dinucleotide_property[Inclination][TT]=1.43;
(*d).dinucleotide_property[Inclination][PMID]=8996793;

short int Major_Groove_Width = 6;
strcpy((*d).dinucleotide_property_string[Major_Groove_Width], "Major_Groove_Width");
(*d).dinucleotide_property[Major_Groove_Width][AA]=12.15;
(*d).dinucleotide_property[Major_Groove_Width][AC]=12.37;
(*d).dinucleotide_property[Major_Groove_Width][AG]=13.51;
(*d).dinucleotide_property[Major_Groove_Width][AT]=12.87;
(*d).dinucleotide_property[Major_Groove_Width][CA]=13.58;
(*d).dinucleotide_property[Major_Groove_Width][CC]=15.49;
(*d).dinucleotide_property[Major_Groove_Width][CG]=14.42;
(*d).dinucleotide_property[Major_Groove_Width][CT]=13.51;
(*d).dinucleotide_property[Major_Groove_Width][GA]=13.93;
(*d).dinucleotide_property[Major_Groove_Width][GC]=14.55;
(*d).dinucleotide_property[Major_Groove_Width][GG]=15.49;
(*d).dinucleotide_property[Major_Groove_Width][GT]=12.37;
(*d).dinucleotide_property[Major_Groove_Width][TA]=12.32;
(*d).dinucleotide_property[Major_Groove_Width][TC]=13.93;
(*d).dinucleotide_property[Major_Groove_Width][TG]=13.58;
(*d).dinucleotide_property[Major_Groove_Width][TT]=12.15;
(*d).dinucleotide_property[Major_Groove_Width][PMID]=8996793;

short int Major_Groove_Depth = 7;
short int Major_Groove_Size = 8;
short int Major_Groove_Distance = 9;
short int Minor_Groove_Width = 10;
short int Minor_Groove_Depth = 11;
short int Minor_Groove_Size = 12;
short int Minor_Groove_Distance = 13;
short int Persistance_Length = 14;
short int Melting_Temperature = 15;
short int None = 16;
strcpy((*d).dinucleotide_property_string[Major_Groove_Depth], "Major_Groove_Depth");
strcpy((*d).dinucleotide_property_string[Major_Groove_Size], "Major_Groove_Size");
strcpy((*d).dinucleotide_property_string[Major_Groove_Distance], "Major_Groove_Distance");
strcpy((*d).dinucleotide_property_string[Minor_Groove_Width], "Minor_Groove_Width");
strcpy((*d).dinucleotide_property_string[Minor_Groove_Depth], "Minor_Groove_Depth");
strcpy((*d).dinucleotide_property_string[Minor_Groove_Size], "Minor_Groove_Size");
strcpy((*d).dinucleotide_property_string[Minor_Groove_Distance], "Minor_Groove_Distance");
strcpy((*d).dinucleotide_property_string[Persistance_Length], "Persistance_Length");
strcpy((*d).dinucleotide_property_string[Melting_Temperature], "Melting_Temperature");
strcpy((*d).dinucleotide_property_string[None], "None");
    
(*d).dinucleotide_property[Major_Groove_Depth][AA]=9.12;
(*d).dinucleotide_property[Major_Groove_Depth][AC]=9.41;
(*d).dinucleotide_property[Major_Groove_Depth][AG]=8.96;
(*d).dinucleotide_property[Major_Groove_Depth][AT]=8.96;
(*d).dinucleotide_property[Major_Groove_Depth][CA]=8.67;
(*d).dinucleotide_property[Major_Groove_Depth][CC]=8.45;
(*d).dinucleotide_property[Major_Groove_Depth][CG]=8.81;
(*d).dinucleotide_property[Major_Groove_Depth][CT]=8.96;
(*d).dinucleotide_property[Major_Groove_Depth][GA]=8.76;
(*d).dinucleotide_property[Major_Groove_Depth][GC]=8.67;
(*d).dinucleotide_property[Major_Groove_Depth][GG]=8.45;
(*d).dinucleotide_property[Major_Groove_Depth][GT]=9.41;
(*d).dinucleotide_property[Major_Groove_Depth][TA]=9.6;
(*d).dinucleotide_property[Major_Groove_Depth][TC]=8.76;
(*d).dinucleotide_property[Major_Groove_Depth][TG]=8.67;
(*d).dinucleotide_property[Major_Groove_Depth][TT]=9.12;
(*d).dinucleotide_property[Major_Groove_Depth][PMID]=8996793;
(*d).dinucleotide_property[Major_Groove_Size][AA]=3.98;
(*d).dinucleotide_property[Major_Groove_Size][AC]=3.98;
(*d).dinucleotide_property[Major_Groove_Size][AG]=4.7;
(*d).dinucleotide_property[Major_Groove_Size][AT]=4.7;
(*d).dinucleotide_property[Major_Groove_Size][CA]=3.98;
(*d).dinucleotide_property[Major_Groove_Size][CC]=3.98;
(*d).dinucleotide_property[Major_Groove_Size][CG]=4.7;
(*d).dinucleotide_property[Major_Groove_Size][CT]=4.7;
(*d).dinucleotide_property[Major_Groove_Size][GA]=3.26;
(*d).dinucleotide_property[Major_Groove_Size][GC]=3.26;
(*d).dinucleotide_property[Major_Groove_Size][GG]=3.98;
(*d).dinucleotide_property[Major_Groove_Size][GT]=3.98;
(*d).dinucleotide_property[Major_Groove_Size][TA]=3.26;
(*d).dinucleotide_property[Major_Groove_Size][TC]=3.26;
(*d).dinucleotide_property[Major_Groove_Size][TG]=3.98;
(*d).dinucleotide_property[Major_Groove_Size][TT]=3.98;
(*d).dinucleotide_property[Major_Groove_Size][PMID]=7897660;
(*d).dinucleotide_property[Major_Groove_Distance][AA]=3.38;
(*d).dinucleotide_property[Major_Groove_Distance][AC]=3.03;
(*d).dinucleotide_property[Major_Groove_Distance][AG]=3.36;
(*d).dinucleotide_property[Major_Groove_Distance][AT]=3.02;
(*d).dinucleotide_property[Major_Groove_Distance][CA]=3.79;
(*d).dinucleotide_property[Major_Groove_Distance][CC]=3.38;
(*d).dinucleotide_property[Major_Groove_Distance][CG]=3.77;
(*d).dinucleotide_property[Major_Groove_Distance][CT]=3.36;
(*d).dinucleotide_property[Major_Groove_Distance][GA]=3.4;
(*d).dinucleotide_property[Major_Groove_Distance][GC]=3.04;
(*d).dinucleotide_property[Major_Groove_Distance][GG]=3.38;
(*d).dinucleotide_property[Major_Groove_Distance][GT]=3.03;
(*d).dinucleotide_property[Major_Groove_Distance][TA]=3.81;
(*d).dinucleotide_property[Major_Groove_Distance][TC]=3.4;
(*d).dinucleotide_property[Major_Groove_Distance][TG]=3.79;
(*d).dinucleotide_property[Major_Groove_Distance][TT]=3.38;
(*d).dinucleotide_property[Major_Groove_Distance][PMID]=7897660;
(*d).dinucleotide_property[Minor_Groove_Width][AA]=5.3;
(*d).dinucleotide_property[Minor_Groove_Width][AC]=6.04;
(*d).dinucleotide_property[Minor_Groove_Width][AG]=5.19;
(*d).dinucleotide_property[Minor_Groove_Width][AT]=5.31;
(*d).dinucleotide_property[Minor_Groove_Width][CA]=4.79;
(*d).dinucleotide_property[Minor_Groove_Width][CC]=4.62;
(*d).dinucleotide_property[Minor_Groove_Width][CG]=5.16;
(*d).dinucleotide_property[Minor_Groove_Width][CT]=5.19;
(*d).dinucleotide_property[Minor_Groove_Width][GA]=4.71;
(*d).dinucleotide_property[Minor_Groove_Width][GC]=4.74;
(*d).dinucleotide_property[Minor_Groove_Width][GG]=4.62;
(*d).dinucleotide_property[Minor_Groove_Width][GT]=6.04;
(*d).dinucleotide_property[Minor_Groove_Width][TA]=6.4;
(*d).dinucleotide_property[Minor_Groove_Width][TC]=4.71;
(*d).dinucleotide_property[Minor_Groove_Width][TG]=4.79;
(*d).dinucleotide_property[Minor_Groove_Width][TT]=5.3;
(*d).dinucleotide_property[Minor_Groove_Width][PMID]=8996793;
(*d).dinucleotide_property[Minor_Groove_Depth][AA]=9.03;
(*d).dinucleotide_property[Minor_Groove_Depth][AC]=8.79;
(*d).dinucleotide_property[Minor_Groove_Depth][AG]=8.98;
(*d).dinucleotide_property[Minor_Groove_Depth][AT]=8.91;
(*d).dinucleotide_property[Minor_Groove_Depth][CA]=9.09;
(*d).dinucleotide_property[Minor_Groove_Depth][CC]=8.99;
(*d).dinucleotide_property[Minor_Groove_Depth][CG]=9.06;
(*d).dinucleotide_property[Minor_Groove_Depth][CT]=8.98;
(*d).dinucleotide_property[Minor_Groove_Depth][GA]=9.11;
(*d).dinucleotide_property[Minor_Groove_Depth][GC]=8.98;
(*d).dinucleotide_property[Minor_Groove_Depth][GG]=8.99;
(*d).dinucleotide_property[Minor_Groove_Depth][GT]=8.79;
(*d).dinucleotide_property[Minor_Groove_Depth][TA]=9;
(*d).dinucleotide_property[Minor_Groove_Depth][TC]=9.11;
(*d).dinucleotide_property[Minor_Groove_Depth][TG]=9.09;
(*d).dinucleotide_property[Minor_Groove_Depth][TT]=9.03;
(*d).dinucleotide_property[Minor_Groove_Depth][PMID]=8996793;
(*d).dinucleotide_property[Minor_Groove_Size][AA]=2.98;
(*d).dinucleotide_property[Minor_Groove_Size][AC]=3.26;
(*d).dinucleotide_property[Minor_Groove_Size][AG]=3.98;
(*d).dinucleotide_property[Minor_Groove_Size][AT]=3.26;
(*d).dinucleotide_property[Minor_Groove_Size][CA]=3.7;
(*d).dinucleotide_property[Minor_Groove_Size][CC]=3.98;
(*d).dinucleotide_property[Minor_Groove_Size][CG]=4.7;
(*d).dinucleotide_property[Minor_Groove_Size][CT]=3.98;
(*d).dinucleotide_property[Minor_Groove_Size][GA]=2.98;
(*d).dinucleotide_property[Minor_Groove_Size][GC]=3.26;
(*d).dinucleotide_property[Minor_Groove_Size][GG]=3.98;
(*d).dinucleotide_property[Minor_Groove_Size][GT]=3.26;
(*d).dinucleotide_property[Minor_Groove_Size][TA]=2.7;
(*d).dinucleotide_property[Minor_Groove_Size][TC]=2.98;
(*d).dinucleotide_property[Minor_Groove_Size][TG]=3.7;
(*d).dinucleotide_property[Minor_Groove_Size][TT]=2.98;
(*d).dinucleotide_property[Minor_Groove_Size][PMID]=7897660;
(*d).dinucleotide_property[Minor_Groove_Distance][AA]=2.94;
(*d).dinucleotide_property[Minor_Groove_Distance][AC]=4.22;
(*d).dinucleotide_property[Minor_Groove_Distance][AG]=2.79;
(*d).dinucleotide_property[Minor_Groove_Distance][AT]=4.2;
(*d).dinucleotide_property[Minor_Groove_Distance][CA]=3.09;
(*d).dinucleotide_property[Minor_Groove_Distance][CC]=2.8;
(*d).dinucleotide_property[Minor_Groove_Distance][CG]=3.21;
(*d).dinucleotide_property[Minor_Groove_Distance][CT]=2.79;
(*d).dinucleotide_property[Minor_Groove_Distance][GA]=2.95;
(*d).dinucleotide_property[Minor_Groove_Distance][GC]=4.24;
(*d).dinucleotide_property[Minor_Groove_Distance][GG]=2.8;
(*d).dinucleotide_property[Minor_Groove_Distance][GT]=4.22;
(*d).dinucleotide_property[Minor_Groove_Distance][TA]=2.97;
(*d).dinucleotide_property[Minor_Groove_Distance][TC]=2.95;
(*d).dinucleotide_property[Minor_Groove_Distance][TG]=3.09;
(*d).dinucleotide_property[Minor_Groove_Distance][TT]=2.94;
(*d).dinucleotide_property[Minor_Groove_Distance][PMID]=7897660;
(*d).dinucleotide_property[Persistance_Length][AA]=35;
(*d).dinucleotide_property[Persistance_Length][AC]=60;
(*d).dinucleotide_property[Persistance_Length][AG]=60;
(*d).dinucleotide_property[Persistance_Length][AT]=20;
(*d).dinucleotide_property[Persistance_Length][CA]=60;
(*d).dinucleotide_property[Persistance_Length][CC]=130;
(*d).dinucleotide_property[Persistance_Length][CG]=85;
(*d).dinucleotide_property[Persistance_Length][CT]=60;
(*d).dinucleotide_property[Persistance_Length][GA]=60;
(*d).dinucleotide_property[Persistance_Length][GC]=85;
(*d).dinucleotide_property[Persistance_Length][GG]=130;
(*d).dinucleotide_property[Persistance_Length][GT]=60;
(*d).dinucleotide_property[Persistance_Length][TA]=20;
(*d).dinucleotide_property[Persistance_Length][TC]=60;
(*d).dinucleotide_property[Persistance_Length][TG]=60;
(*d).dinucleotide_property[Persistance_Length][TT]=35;
(*d).dinucleotide_property[Persistance_Length][PMID]=3627268;
(*d).dinucleotide_property[Melting_Temperature][AA]=54.5;
(*d).dinucleotide_property[Melting_Temperature][AC]=97.73;
(*d).dinucleotide_property[Melting_Temperature][AG]=58.42;
(*d).dinucleotide_property[Melting_Temperature][AT]=57.02;
(*d).dinucleotide_property[Melting_Temperature][CA]=54.71;
(*d).dinucleotide_property[Melting_Temperature][CC]=85.97;
(*d).dinucleotide_property[Melting_Temperature][CG]=72.55;
(*d).dinucleotide_property[Melting_Temperature][CT]=58.42;
(*d).dinucleotide_property[Melting_Temperature][GA]=86.44;
(*d).dinucleotide_property[Melting_Temperature][GC]=136.12;
(*d).dinucleotide_property[Melting_Temperature][GG]=85.97;
(*d).dinucleotide_property[Melting_Temperature][GT]=97.73;
(*d).dinucleotide_property[Melting_Temperature][TA]=36.73;
(*d).dinucleotide_property[Melting_Temperature][TC]=86.44;
(*d).dinucleotide_property[Melting_Temperature][TG]=54.71;
(*d).dinucleotide_property[Melting_Temperature][TT]=54.5;
(*d).dinucleotide_property[Melting_Temperature][PMID]=-1;
(*d).dinucleotide_property[None][AA]=1;
(*d).dinucleotide_property[None][AC]=1;
(*d).dinucleotide_property[None][AG]=1;
(*d).dinucleotide_property[None][AT]=1;
(*d).dinucleotide_property[None][CA]=1;
(*d).dinucleotide_property[None][CC]=1;
(*d).dinucleotide_property[None][CG]=1;
(*d).dinucleotide_property[None][CT]=1;
(*d).dinucleotide_property[None][GA]=1;
(*d).dinucleotide_property[None][GC]=1;
(*d).dinucleotide_property[None][GG]=1;
(*d).dinucleotide_property[None][GT]=1;
(*d).dinucleotide_property[None][TA]=1;
(*d).dinucleotide_property[None][TC]=1;
(*d).dinucleotide_property[None][TG]=1;
(*d).dinucleotide_property[None][TT]=1;
  
for(counter = 0; counter < (*d).number_of_dinucleotide_properties; counter++) for(counter2 = 0; counter2 < 17; counter2++) if((*d).dinucleotide_property[counter][counter2] == -1.1111) printf("\nUNINITIALIZED DINUCLEOTIDE PROPERTY %i,%i",counter,counter2);
    
}

/* COUNT PAIR */
struct countpair {long int sequence_value; short int 
    contains_CpG; short int 
    add_logo; short int gap_position; short int gap_width; long int x_count; long int y_count; short int x_local_max; short int y_local_max; long int sum;};

/* SEQUENCE INCIDENCE TABLE */
struct sequence_incidence_table {__uint128_t sequence_value_ULL; long int incidence;};

/* KMER INCIDENCE TABLE */
struct kmer_incidence_table {short int kmer_length; long int kmer; long int incidence; float max_enrichment; float local_max_enrichment; long int local_max_max_incidence; short int max_gap_length; short int local_max_max_gap_length; short int local_max; short int any_local_max; short int preferred;};

/* KMER LIST */
struct kmer_list {char **kmer; __uint128_t *kmer_value; short int kmer_length; short int number_of_kmers;};
short int kmer_list_init (struct kmer_list *k, char **kmers, short int kmer_length, __uint128_t *kmer_value, short int query_sequence_length, short int number_of_kmers)
{
unsigned long int query_sequence_value;
unsigned long int position;
unsigned long int position_value;
unsigned long int nucleotide_value;
short int counter;
short int counter2;
char *forward;
if (rna == 0) forward = dnaforward;
else forward = rnaforward;
    
(*k).kmer = kmers;
(*k).kmer_length = kmer_length;
(*k).number_of_kmers = number_of_kmers;
/* COUNTS KMER VALUES */
for (counter = 0; counter < number_of_kmers; counter++)
{
	for(query_sequence_value = 0, position = 0, position_value = pow(4,query_sequence_length-1); position < query_sequence_length; position++, position_value /= 4)
	{
		for (nucleotide_value = 0; nucleotide_value < 4 && kmers[counter][position] != forward[nucleotide_value]; nucleotide_value++);
		if(nucleotide_value == 4) {printf("\nERROR IN QUERY SEQUENCE\n"); exit (1);}
		query_sequence_value += position_value * nucleotide_value;
	}
(*k).kmer_value[counter] = query_sequence_value;
}
return(0);
}

/* KMER MATCHES */
struct kmer_matches {char **kmer; __uint128_t *kmer_value; short int kmer_length; short int number_of_kmers; long int number_of_sequences; short int **number_of_matches;};
short int kmer_matches_init (struct kmer_list *k, struct kmer_matches *i, short int number_of_kmers, long int number_of_sequences)
{
short int counter;
short int counter2;
short int initial_value = 0;
(*i).number_of_kmers = number_of_kmers;
(*i).number_of_sequences = number_of_sequences;
(*i).kmer_length = (*k).kmer_length;
(*i).kmer = malloc(sizeof(double *) * 5 + 5);
(*i).number_of_matches = malloc(sizeof(double *) * 5 + 5);
for (counter = 0; counter < 5; counter++) 
{
(*i).number_of_matches[counter] = malloc(sizeof(double) * (*i).number_of_sequences + 5);
for (counter2 = 0; counter2 < number_of_sequences; counter2++) (*i).number_of_matches[counter][counter2] = initial_value;
}
return(0);
}


/* COUNT PWM */
struct count_pwm {char *name; short int width; long int max_counts; double **incidence;};
short int count_pwm_init (struct count_pwm *i, char *name, short int width, double initial_value)
{
short int maximum_width = max_width_of_pwm;
short int counter;
short int counter2;
(*i).name = malloc(1000);
strcpy ((*i).name, name);
(*i).width = width;
(*i).max_counts = initial_value;
(*i).incidence = malloc(sizeof(double *) * (5 + contacts * 12) + 5);
for (counter = 0; counter < 5 + contacts * 12; counter++) 
{
(*i).incidence[counter] = malloc(sizeof(double) * maximum_width + 5);
for (counter2 = 0; counter2 < maximum_width; counter2++) (*i).incidence[counter][counter2] = initial_value;
}
return(0);
}
short int count_pwm_free (struct count_pwm *i)
{
short int counter;
free((*i).name);
for (counter = 0; counter < 5; counter++) free((*i).incidence[counter]);
free((*i).incidence);
}
                          
struct pairwise_correlation {short int first_base; short int second_base; double delta_ic; short int max_dinucleotide; short int min_dinucleotide; double min_fold_change; double max_fold_change;};

/* NORMALIZED PWM */
struct normalized_pwm {char *name; char *seed; short int width; long int max_counts; double *information_content; short int *original_position; double *position_score; long int *total_counts_for_column; double **fraction; struct pairwise_correlation *pairwise_correlation; short int negative_values_allowed; struct oriented_match match;};
short int normalized_pwm_init (struct normalized_pwm *i, char *name, short int width, double initial_value)
{
short int maximum_width = max_width_of_pwm;
short int counter;
short int counter2;
(*i).negative_values_allowed = 0;
(*i).pairwise_correlation = malloc(sizeof(struct pairwise_correlation) * 10 + 5);
for (counter = 0; counter < 10; counter++)
{
(*i).pairwise_correlation[counter].first_base = 0;
(*i).pairwise_correlation[counter].delta_ic = 0;
(*i).pairwise_correlation[counter].second_base = 0;
(*i).pairwise_correlation[counter].max_dinucleotide = 0;
(*i).pairwise_correlation[counter].min_dinucleotide = 0;
(*i).pairwise_correlation[counter].max_fold_change = 0;
(*i).pairwise_correlation[counter].min_fold_change = 0;
}
(*i).name = malloc(100);
strcpy ((*i).name, name);
(*i).seed = malloc(1000);
strcpy ((*i).seed, "UNKNOWN");
(*i).width = width;
(*i).max_counts = initial_value;
(*i).fraction = malloc(sizeof(double *) * (5 + contacts * 12) + 5);
(*i).information_content = malloc(sizeof(double) * maximum_width + 5);
(*i).position_score = malloc(sizeof(double) * maximum_width + 5);
(*i).original_position = malloc(sizeof(short int) * maximum_width + 5);
(*i).total_counts_for_column = malloc(sizeof(long int) * maximum_width + 5);

for (counter = 0; counter < 5 + contacts * 12; counter++) 
{
(*i).fraction[counter] = malloc(sizeof(double) * maximum_width + 5);
for (counter2 = 0; counter2 < maximum_width; counter2++) (*i).fraction[counter][counter2] = initial_value;
}
for (counter2 = 0; counter2 < maximum_width; counter2++)
{
(*i).information_content[counter2] = 0;
(*i).position_score[counter2] = 0;
(*i).original_position[counter2] = counter2;
(*i).total_counts_for_column[counter2] = 0;
}
(*i).match.position = 0;
(*i).match.strand = 0;
(*i).match.score = 0;
return(0);
}
short int normalized_pwm_free (struct normalized_pwm *i)
{
short int counter;
free((*i).name);
free((*i).pairwise_correlation);
free((*i).information_content);
free((*i).position_score);
free((*i).total_counts_for_column);    
for (counter = 0; counter < 5; counter++) free((*i).fraction[counter]);
free((*i).fraction);
}

/* DINUCLEOTIDE PROPERTIES MATRIX */
struct dinucleotide_properties_matrix {char *name; short int width; short int number_of_dinucleotide_properties; long int max_counts; long int **count; double **score; char **dinucleotide_property_string; short int query_sequence_length;};
short int dinucleotide_properties_matrix_init (struct dinucleotide_properties_matrix *i, struct dinucleotide_properties *d, char *name, short int width, double initial_value, short int query_sequence_length)
{
(*i).query_sequence_length = query_sequence_length;
(*i).number_of_dinucleotide_properties = (*d).number_of_dinucleotide_properties;
short int maximum_width = width;
short int counter;
short int counter2;
(*i).name = malloc(1000);
(*i).max_counts = 0;
strcpy ((*i).name, name);
(*i).width = width;
(*i).score = malloc(sizeof(double *) * ((*i).number_of_dinucleotide_properties + 1) + 5);
(*i).count = malloc(sizeof(long int *) * ((*i).number_of_dinucleotide_properties + 1) + 5);
(*i).dinucleotide_property_string = malloc(sizeof(char *) * ((*i).number_of_dinucleotide_properties + 1) + 5);
for (counter = 0; counter <= (*i).number_of_dinucleotide_properties; counter++) 
{
(*i).score[counter] = malloc(sizeof(double) * maximum_width + 5);
(*i).count[counter] = malloc(sizeof(long int) * maximum_width + 5);
(*i).dinucleotide_property_string[counter] = malloc(1000);
strcpy((*i).dinucleotide_property_string[counter], (*d).dinucleotide_property_string[counter]);
for (counter2 = 0; counter2 < width; counter2++) 
{
(*i).count[counter][counter2] = 0;
(*i).score[counter][counter2] = initial_value;
}
}
return(0);
}

/* COUNT CONNECTING MATRIX */
struct count_connecting_matrix {char *name; short int width; short int height; long int number_of_total_matches; long int **incidence;};
short int count_connecting_matrix_init (struct count_connecting_matrix *i, char *name, short int width, double initial_value)
{
short int maximum_width = max_width_of_pwm;
short int counter;
short int counter2;
(*i).number_of_total_matches = 0;
(*i).name = malloc(100);
strcpy ((*i).name, name);
(*i).width = width;
(*i).incidence = malloc(sizeof(double *) * 5 + 5);
for (counter = 0; counter < 5; counter++) 
{
(*i).incidence[counter] = malloc(sizeof(double) * maximum_width + 5);
for (counter2 = 0; counter2 < maximum_width; counter2++) (*i).incidence[counter][counter2] = initial_value;
}
return(0);
}

/* NORMALIZED CONNECTING MATRIX */
struct normalized_connecting_matrix {char *name; short int width; short int height; long int *orientation_count; double *orientation_fraction; double **fraction;};
short int normalized_connecting_matrix_init (struct normalized_connecting_matrix *i, char *name, short int width, double initial_value)
{
short int maximum_width = max_width_of_pwm;
short int counter;
short int counter2;
(*i).name = malloc(100);
strcpy ((*i).name, name);
(*i).width = width;
(*i).fraction = malloc(sizeof(double *) * 5 + 5);
(*i).orientation_fraction = malloc(sizeof(double) * 5 + 5);
(*i).orientation_count = malloc(sizeof(long int) * 5 + 5);
for (counter = 0; counter < 4; counter++) 
{
(*i).fraction[counter] = malloc(sizeof(double) * maximum_width + 5);
for (counter2 = 0; counter2 < maximum_width; counter2++) (*i).fraction[counter][counter2] = initial_value;
}
for (counter2 = 0; counter2 < 4; counter2++)
{
(*i).orientation_fraction[counter2] = 0;
(*i).orientation_count[counter2] = 0;
}
return(0);
}


struct hit_position_matrix {char *name; short int width; long int **incidence; double **fraction;};
short int hit_position_matrix_init (struct hit_position_matrix *i, char *name, short int width, double initial_value)
{
short int maximum_width = max_Nlength;
short int counter;
short int counter2;
(*i).name = malloc(100);
strcpy ((*i).name, name);
(*i).width = width;
(*i).fraction = malloc(sizeof(double *) * 3 + 5);
(*i).incidence = malloc(sizeof(long int *) * 3 + 5);
for (counter = 0; counter < 2; counter++) 
{
(*i).fraction[counter] = malloc(sizeof(double) * maximum_width + 5);
(*i).incidence[counter] = malloc(sizeof(long int) * maximum_width + 5);
for (counter2 = 0; counter2 < maximum_width; counter2++) 
{
(*i).fraction[counter][counter2] = initial_value;
(*i).incidence[counter][counter2] = initial_value;
}
}
return(0);
}


struct match {short int *position; double *score;};
short int match_init (struct match *i, short int width)
{
short int maximum_width = Nlength+10;
short int counter;
short int counter2;
(*i).position = malloc(sizeof(short int) * maximum_width + 5);
(*i).score = malloc(sizeof(double) * maximum_width + 5);
for (counter = 0; counter < maximum_width; counter++)
{
(*i).position[counter] = 0;
(*i).score[counter] = 0;
}
return(0);
}




/* SVG_TILE */
struct svg_tile {short int x; short int y; short int height; short int width; short int red; short int green; short int blue; short int stroke;};

/* SUBROUTINES */

/* SUBROUTINE THAT ESCAPES OFFENDING CHARACTERS FROM STRINGS FOR INCLUSION IN SVG COMMENTS */
char *svgsafe(char *string)
{
    signed long int counter;
    signed long int counter2;
    
    for(counter = strlen(string); counter >= 0; counter--)
    {
        if(string[counter] == '-' && string[counter+1] == '-') string[counter] = '_';
    }
    return(string);
}

/* SUBROUTINE THAT PRINTS A LINE-TILE TO SVG HEATMAP */
short int Linetile (FILE *outfile, struct base_dependency_matrix *dep, struct svg_tile *t, short int add_line)
{
fprintf(outfile, " <rect x=\"%i\" y=\"%i\" width=\"%i", (*t).x * (*t).width, (*t).y * (*t).height,  (*t).width);
fprintf(outfile, "\" height=\"");
fprintf(outfile, "%i", (*t).height);
fprintf(outfile, "\" style=\"");
fprintf(outfile, "fill:rgb(");
fprintf(outfile, "%i,%i,%i);stroke-width:%i;stroke:rgb(0,0,0)\"/>\n", (*t).red, (*t).green, (*t).blue, (*t).stroke);
return(0);
}

/* SUBROUTINE THAT REVERSES BASE CHARACTER */
char Reverse_base (char base)
{
short int counter;
char reverse_nucleotide;
char *nucleotide_iupac = dna_iupac;
char *nucleotide_rc_iupac = dna_rc_iupac;
if (rna == 1)
{
nucleotide_iupac = rna_iupac;
nucleotide_rc_iupac = rna_rc_iupac;
}
    
for (counter = 0; counter < strlen(nucleotide_iupac); counter++)
{
if(base == nucleotide_iupac[counter])
{
reverse_nucleotide = nucleotide_rc_iupac[counter];
break;
}
}
if (base == nucleotide_iupac[counter]) return(reverse_nucleotide);
else return ('*');
}

/* CHECKS IF STRING IS AN IUPAC PALINDROME, ALLOWS EXCLUSION OF TWO BASES FOR DINUCLEOTIDE ANALYSES */
Is_this_string_iupac_palindrome (char *test_string, signed short int exclude1, signed short int exclude2)
{
    char *nucleotide_iupac = dna_iupac;
    char *nucleotide_rc_iupac = dna_rc_iupac;
    if (rna == 1)
    {
        nucleotide_iupac = rna_iupac;
        nucleotide_rc_iupac = rna_rc_iupac;
    }
   
    /* printf("\nCalls IupacP %s", test_string); */
    signed short int string_position = strlen(test_string)-1;
    signed short int string_length = string_position;
    signed short int iupac_position = 0;
    signed short int rc_iupac_position = 0;
    for (;string_position >= 0 && iupac_position == rc_iupac_position; string_position--)
    {
        if (string_position == exclude1-1 || string_position == exclude2-1 || string_position == string_length + 1 - exclude1 || string_position == string_length + 1 - exclude2) continue;
        for(iupac_position = strlen(nucleotide_iupac); iupac_position >= 0 && test_string[string_position] != nucleotide_iupac[iupac_position]; iupac_position--);
        for(rc_iupac_position = strlen(nucleotide_rc_iupac); rc_iupac_position >= 0, test_string[string_length-string_position] != nucleotide_rc_iupac[rc_iupac_position]; rc_iupac_position--);
        /* printf("\nIupacP %i,%c,%i,%c", iupac_position, iupac[iupac_position], rc_iupac_position, rc_iupac[rc_iupac_position]); */
    }
    if (string_position == -1) return (1);
    else return (0);
}

/* SUBROUTINE THAT PLOTS SVG GREEN_BLACK_RED HEATMAP TO outfile AT x,y USING DATA IN n SCALED TO GREEN min_score BLACK black_score RED max_score. IF SCALE VARIABLES SET TO -1 DETERMINES MIN CENTER AND/OR MAX */
short int Svg_heatmap (FILE *outfile, short int x, short int y, long int min_score, float black_score, long int max_score, struct count_connecting_matrix *n)
{

short int current_pwm = 0;
short int top_position;
short int counter;
short int nucleotide_value;
short int tile_width = 20;
short int tile_height = 30;
short int stroke = 1;

double scale_end_minimum = 1000;
short int number_of_scale_tiles = 6;
double current_scale_tile_value;
double scale_tile_increment;
short int tile = 0;
    
short int x_position;
short int row;

char **orientation;
orientation = malloc(sizeof(char *) * 5 + 5);
orientation[0] = " HT ";
orientation[1] = " HH ";
orientation[2] = " TT ";
short int red;
short int green;
short int blue = 0;

char *font = "Courier";

fprintf(outfile, "<g id=\"group_%s\" transform=\"translate(%i, %i)\" >", (*n).name, x, y);

/* DEFINES MIN MAX AND CENTER IF NOT GIVEN */
if (max_score == -1) for(row = 0; row < 3; row++) for(x_position = 0; x_position < (*n).width; x_position++) if ((*n).incidence[row][x_position] > max_score && (*n).incidence[row][x_position] != 0) max_score = (*n).incidence[row][x_position];
if (min_score == -1) {min_score = max_score; for(row = 0; row < 3; row++) for(x_position = 0; x_position < (*n).width; x_position++) if ((*n).incidence[row][x_position] < min_score && (*n).incidence[row][x_position] != 0) min_score = (*n).incidence[row][x_position];}
if (black_score == -1) black_score = (max_score+min_score+0.0001) / 2;

current_scale_tile_value = min_score;
scale_tile_increment = (double) (max_score-min_score) / (number_of_scale_tiles-1);
    
/* PRINTS COORDINATES AND ORIENTATION LABELS */
for(row = 0; row < 3; row++) {fprintf(outfile, " <text  x=\"15\" y=\"%i\" fill = \"black\" stroke = \"black\" font-size=\"30\" font-family = \"%s\">%s</text>\n", (row + 1) * tile_height + 25, font, orientation[row]);}
for(x_position = 0; x_position < (*n).width; x_position++) 
{
fprintf(outfile, " <text  x=\"%i\" y=\"22\" fill = \"black\" stroke = \"black\" font-size=\"14\" font-family = \"%s\">%i</text>\n", x_position * tile_width + 95 - (x_position > 9) * 5, font, x_position);
}

    
/* GENERATES HEATMAP */
for(row = 0; row < 3; row++)
{
for(x_position = 0; x_position < (*n).width; x_position++)
{
/* DETERMINES COLOR */
red = 0;
green = 0;
if ((*n).incidence[row][x_position] > black_score) red = ((*n).incidence[row][x_position] - black_score) * 256 / (max_score-black_score);
else green = (black_score - (*n).incidence[row][x_position]) * 256 / (black_score-min_score);

/* PRINTS OUT TILE */
fprintf(outfile, " <rect x=\"%i\" y=\"%i\" width=\"%i", x_position * tile_width + 90, (row + 1) * tile_height, tile_width);
fprintf(outfile, "\" height=\"");
fprintf(outfile, "%i", tile_height);
fprintf(outfile, "\" style=\"");
fprintf(outfile, "fill:rgb(");
fprintf(outfile, "%i,%i,%i);stroke-width:%i;stroke:rgb(0,0,0)\"/>\n", red, green, blue, stroke);
}
    
/* PRINTS SCALE TILES */
    for(tile = 0; tile < number_of_scale_tiles / 3; tile++)
{
    red = 0;
    green = 0;
    if (current_scale_tile_value > black_score) red = (current_scale_tile_value - black_score) * 256 / (max_score-black_score);
    else green = (black_score - current_scale_tile_value) * 256 / (black_score-min_score);
    /* PRINTS OUT TILE */
    fprintf(outfile, " <rect x=\"%i\" y=\"%i\" width=\"%i", x_position * tile_width + 100, (row + 1) * tile_height + tile * tile_height/2, tile_width/2);
    fprintf(outfile, "\" height=\"");
    fprintf(outfile, "%i", tile_height/2);
    fprintf(outfile, "\" style=\"");
    fprintf(outfile, "fill:rgb(");
    fprintf(outfile, "%i,%i,%i);stroke-width:%i;stroke:rgb(0,0,0)\"/>\n", red, green, blue, stroke); 
    fprintf(outfile, " <text  x=\"%i\" y=\"%i\" fill = \"black\" stroke = \"gray\" font-size=\"14\" font-family = \"%s\">%.0f (%.1f%%)</text>\n",  x_position * tile_width + 115 + ((current_scale_tile_value<9.5)+(current_scale_tile_value<99.5)+(current_scale_tile_value<999.5)+(current_scale_tile_value<9999.5))*8, (row + 1) * tile_height + tile * tile_height/2 + tile_height/2 - 4, font,current_scale_tile_value, 200 * (double) current_scale_tile_value / (*n).number_of_total_matches);
    current_scale_tile_value += scale_tile_increment;
}
}

/* PRINTS OUT NAME OF CONNECTING MATRIX AND OTHER DATA */
fprintf(outfile, "<text  x=\"");
fprintf(outfile, "%i", x_position * tile_width + 250);
fprintf(outfile, "\" y=\"");
fprintf(outfile, "%i", row * tile_height);
fprintf(outfile, "\" fill = \"black\" stroke = \"black\" font-size=\"30\" font-family = \"");
fprintf(outfile, "%s", font);
fprintf(outfile, "\" >");
fprintf(outfile, "%s", (*n).name);
fprintf(outfile, "</text>\n");
fprintf(outfile, "</g>\n"); 
return(0);
}

/* SUBROUTINE THAT CALCULATES MEANS FOR SCORES IN DINUCLEOTIDE PROPERTIES MATRIX n */
short int Average_dinuceotide_properties_matrix (struct dinucleotide_properties_matrix *n)
{
short int row;
short int x_position;
for(row = 0; row < (*n).number_of_dinucleotide_properties; row++) for(x_position = 0; x_position < (*n).width; x_position++)
if ((*n).count[row][x_position] != 0) 
{
(*n).score[row][x_position] = (double) ((*n).score[row][x_position] / (double) (*n).count[row][x_position]);
(*n).count[row][x_position] = 1;
}
else (*n).count[row][x_position] = 0;
}

/* SUBROUTINE THAT PLOTS SVG GREEN_BLACK_RED BASE DEPENDENCY (DEVIATION) HEATMAP */
short int Svg_base_dependency_heatmap (FILE *outfile, char *seed, short int x, short int y, struct base_dependency_matrix *n, struct base_dependency_matrix *ne, struct normalized_pwm **p)
{
    
    char *tempstring;
    tempstring = malloc(1000);
    tempstring[0] = '\0';
    
    short int current_pwm = 0;
    short int top_position;
    short int counter;
    short int counter2;
    short int nucleotide_value;
    short int tile_width = 20;
    short int tile_height = 15;
    double inner_tile_size = 0.5;
    double circle_radius;
    short int stroke = 1;
    long int total_counts;
    
    short int nucleotide1;
    short int nucleotide2;
    
    short int x_position;
    short int row;
    
    short int red = 0;
    short int green = 0;
    short int blue = 0;
    double max_score;
    double min_score;
    double black_score;
    double scale_end_minimum = 0.7;
    short int number_of_scale_tiles = 10;
    double current_scale_tile_value;
    double scale_tile_increment;
    
    short int x_displacement = 8;
    short int right_crop = Nlength/2 - 6;
    short int left_crop = Nlength/2 - 6;
    short int top_crop = Nlength-5;
    short int horizontal_seed_offset = Nlength/2-2;
    short int vertical_seed_offset = 0;
    char *font = "Courier";
    char *less_than_expected_color = "Blue";
    char *more_than_expected_color = "Yellow";
    char *circle_color;
    double current_ic;
 
    /* STARTS PRINTING */
    fprintf(outfile, "<g id=\"group_%s\" transform=\"translate(%i, %i)\" >", (*n).name, x, y);
        
    /* DEFINES MIN MAX AND CENTER */
    for (row=top_crop, min_score = 1E20, max_score = -1E20; row <= (*n).width - (right_crop+left_crop); row++) for(counter = left_crop; counter < (*n).width - right_crop; counter++) 
    {
        current_ic = (*n).total_relative_deviation[row][counter];
        if (current_ic != 0 && current_ic <= 20 && current_ic >= -20)
        {
            if (current_ic > max_score) max_score = current_ic;
            if (current_ic < min_score) min_score = current_ic;
        }
    }

    /* SCALES MIN TO 0 */
    min_score = 0;
    /* SCALES HEATMAP TO ABSOLUTE SCALE IF MAX IS NOT HIGHER THAN LIMIT */
    if(max_score < scale_end_minimum) 
    {
    max_score = scale_end_minimum;
    min_score = 0;
    }
    black_score = (max_score+min_score) / 2;
    scale_tile_increment = (max_score-min_score) / (number_of_scale_tiles-1);
    current_scale_tile_value = min_score;
    
    /* GENERATES HEATMAP */
    for(counter = 0, counter2 = 0, row = top_crop; row <= (*n).width - (right_crop+left_crop); row++)
    {
        for(x_position = left_crop + x_displacement; x_position < (*n).width - right_crop; x_position++)
        {
            /* DETERMINES COLOR */
            red = 0;
            green = 0;
            current_ic = (*n).total_relative_deviation[row][x_position] - (*n).count_statistic_expected_total_relative_deviation[row][x_position];
            if (current_ic > black_score) red = (int) ((current_ic - black_score) * 256 / (max_score-black_score));
            else green = (int) ((black_score - current_ic) * 256 / (black_score-min_score));
            /* printf("\nRow, Col, Red, Green %i, %i, %i,%i", row, x_position, red, green); */
            /* PRINTS OUT TILE */
            fprintf(outfile, " <rect x=\"%i\" y=\"%i\" width=\"%i", (x_position-x_displacement) * tile_width + 160, (row + 1 - top_crop) * tile_height, tile_width);
            fprintf(outfile, "\" height=\"");
            fprintf(outfile, "%i", tile_height);
            fprintf(outfile, "\" style=\"");
            fprintf(outfile, "fill:rgb(");
            fprintf(outfile, "%i,%i,%i);stroke-width:%i;stroke:rgb(0,0,0)\"/>\n", red, green, blue, stroke);

            current_ic = (*n).total_relative_deviation[row][x_position];
            if (current_ic > black_score) red = (int) ((current_ic - black_score) * 256 / (max_score-black_score));
            else green = (int) ((black_score - current_ic) * 256 / (black_score-min_score));
            /* printf("\nRow, Col, Red, Green %i, %i, %i,%i", row, x_position, red, green); */
            /* PRINTS OUT INNER TILE */
            fprintf(outfile, " <rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f", ((x_position-x_displacement+(1-inner_tile_size)/2) * tile_width + 160), ((row + 1 - top_crop+(1-inner_tile_size)/2) * tile_height), (tile_width*inner_tile_size));
            fprintf(outfile, "\" height=\"");
            fprintf(outfile, "%.1f", (tile_height*inner_tile_size));
            fprintf(outfile, "\" style=\"");
            fprintf(outfile, "fill:rgb(");
            fprintf(outfile, "%i,%i,%i);stroke-width:%i;stroke:rgb(0,0,0)\"/>\n", red, green, blue, 0);            
           
            /* PRINTS OUT CIRCLES */            
            for(total_counts = 0, nucleotide1 = 0; nucleotide1 < 4; nucleotide1++) for(nucleotide2 = 0; nucleotide2 < 4; nucleotide2++) total_counts += (*n).incidence[row][x_position][nucleotide1*4+nucleotide2];
            for(nucleotide1 = 0; nucleotide1 < 4; nucleotide1++) for(nucleotide2 = 0; nucleotide2 < 4; nucleotide2++)
            {
                if(total_counts>0)
                {
                circle_radius = ((float) (abs((*n).incidence[row][x_position][nucleotide1*4+nucleotide2]-(*ne).incidence[row][x_position][nucleotide1*4+nucleotide2])*tile_width))/(2* (float) total_counts);
                circle_color = more_than_expected_color;
                if ((*n).incidence[row][x_position][nucleotide1*4+nucleotide2] < (*ne).incidence[row][x_position][nucleotide1*4+nucleotide2]) circle_color = less_than_expected_color;
                if(circle_radius > 0.1)
                {
            fprintf(outfile, " <circle cx=\"%.2f\" cy=\"%.2f\" r=\"%.2f\" fill=\"%s\" stroke-width=\"%.2f\"/>\n", ((x_position-x_displacement+(float)(nucleotide2+1)/5) * tile_width + 160), ((row + 1 - top_crop + (float)(nucleotide1+1)/5) * tile_height), circle_radius,circle_color, (float) 0);
                }
                }
            }
            
            /* PRINTS HORIZONTAL SEED */
            if(row == top_crop + 4 && counter2 < strlen(seed))
            {
                fprintf(outfile, " <text  x=\"%i\" y=\"%i\" fill = \"white\" stroke = \"white\" font-size=\"14\" font-family = \"%s\">%c</text>\n", (horizontal_seed_offset + x_position-x_displacement) * tile_width + 165, (row - top_crop) * tile_height + 10, font,seed[counter2]);
                counter2++;
            }
        }
        /* PRINTS VERTICAL SEED */
        if (row-top_crop > 3 && counter < strlen(seed)) 
        {
        fprintf(outfile, " <text  x=\"%i\" y=\"%i\" fill = \"gray\" stroke = \"gray\" font-size=\"14\" font-family = \"%s\">%c</text>\n", (x_position-x_displacement-8-Nlength/2+(int) strlen(seed)) * tile_width + 165, (vertical_seed_offset + row - top_crop + 1) * tile_height + 10, font,seed[counter]);
        counter++;
        }
        /* PRINTS SCALE TILE */
        if(row-top_crop > 0 && number_of_scale_tiles >= 0)
        {
            /* DETERMINES COLOR */
            red = 0;
            green = 0;
            if (number_of_scale_tiles == 0) 
            {
            current_scale_tile_value = (*n).max_expected_relative_deviation;
            strcpy(tempstring, "max error: ");
            }
            if (current_scale_tile_value > black_score) red = (int) ((current_scale_tile_value - black_score) * 256 / (max_score-black_score));
            else green = (int) ((black_score - current_scale_tile_value) * 256 / (black_score-min_score));
            /* PRINTS OUT TILE */
            fprintf(outfile, " <rect x=\"%i\" y=\"%i\" width=\"%i", (x_position-x_displacement) * tile_width + 180, (row + 1 - top_crop) * tile_height + (number_of_scale_tiles == 0) * 20, tile_width);
            fprintf(outfile, "\" height=\"");
            fprintf(outfile, "%i", tile_height);
            fprintf(outfile, "\" style=\"");
            fprintf(outfile, "fill:rgb(");
            fprintf(outfile, "%i,%i,%i);stroke-width:%i;stroke:rgb(0,0,0)\"/>\n", red, green, blue, stroke);
            x_position++;
            /* PRINTS OUT VALUE */
            fprintf(outfile, " <text  x=\"%i\" y=\"%i\" fill = \"black\" stroke = \"gray\" font-size=\"14\" font-family = \"%s\">%s%.1f%%</text>\n", (x_position-x_displacement) * tile_width + 190 + (current_scale_tile_value<0.1)*8, (vertical_seed_offset + row - top_crop + 1) * tile_height + 10 + (number_of_scale_tiles == 0) * 20, font, tempstring, current_scale_tile_value * 100);
            number_of_scale_tiles--;
            current_scale_tile_value += scale_tile_increment;
        }
        }
    
    /* PRINTS OUT NAME OF MATRIX AND OTHER DATA */
    fprintf(outfile, "<text  x=\"");
    fprintf(outfile, "%i", x_position * tile_width + 350);
    fprintf(outfile, "\" y=\"");
    fprintf(outfile, "%i", (row - 2 - top_crop) * tile_height);
    fprintf(outfile, "\" fill = \"black\" stroke = \"black\" font-size=\"30\" font-family = \"");
    fprintf(outfile, "%s", font);
    fprintf(outfile, "\" >");
    fprintf(outfile, "%s", (*n).name);
    fprintf(outfile, "</text>\n");
    fprintf(outfile, "</g>\n"); 
    free(tempstring);
    return(y+ (row + 20 - top_crop) * tile_height);
    
}

/* SUBROUTINE THAT PLOTS SVG GREEN_BLACK_RED DINUCLEOTIDE PROPERTIES HEATMAP TO outfile AT x,y USING DATA IN n SCALED TO GREEN min_score BLACK black_score RED max_score. IF SCALE VARIABLES SET TO -1 DETERMINES MIN CENTER AND/OR MAX */
short int Svg_dinucleotide_heatmap (FILE *outfile, short int x, short int y, struct dinucleotide_properties_matrix *n)
{

short int current_pwm = 0;
short int top_position;
short int counter;
short int nucleotide_value;
short int tile_width = 20;
short int tile_height = 15;
short int stroke = 1;

short int x_position;
short int row;

short int red;
short int green;
short int blue = 0;

double max_score;
double min_score;
double black_score;
short int x_displacement = 8;
short int right_crop = 7;
short int left_crop = (*n).query_sequence_length - 2;
char *font = "Courier";

/* CALCULATES AVERAGE SCORES */
Average_dinuceotide_properties_matrix (n);
 
/* STARTS PRINTING */
fprintf(outfile, "<g id=\"group_%s\" transform=\"translate(%i, %i)\" >", (*n).name, x, y);

/* PRINTS COORDINATES AND ORIENTATION LABELS
for(row = 0; row < (*n).number_of_dinucleotide_properties; row++) {fprintf(outfile, " <text  x=\"700\" y=\"%i\" fill = \"black\" stroke = \"black\" font-size=\"14\" font-family = \"%s\">%s</text>\n", (row + 1) * tile_height + 10, font, (*n).dinucleotide_property_string[row]);}
/* for(x_position = 2; x_position < (*n).width; x_position++) 
{
fprintf(outfile, " <text  x=\"%i\" y=\"15\" fill = \"black\" stroke = \"black\" font-size=\"14\" font-family = \"%s\">%i</text>\n", x_position * tile_width + 155 - (x_position > 9) * 5, font, x_position);
} */

/* GENERATES HEATMAP */
for(row = 0; row < (*n).number_of_dinucleotide_properties; row++)
{
for(x_position = left_crop + x_displacement; x_position < (*n).width - right_crop; x_position++)
{
/* DETERMINES COLOR */

/* DEFINES MIN MAX AND CENTER IF NOT GIVEN */
for(counter = left_crop + x_displacement, min_score = 1E20, max_score = -1E20; counter < (*n).width - right_crop; counter++) 
{
/*    printf("\n==%i",row);
    printf(",%i",counter);
    printf("\t%.1f", (*n).score[row][counter]);
    fflush(stdout); */
if ((*n).score[row][counter] != 0)
{
if ((*n).score[row][counter] > max_score) max_score = (*n).score[row][counter];
if ((*n).score[row][counter] < min_score) min_score = (*n).score[row][counter];
}
}
black_score = (max_score+min_score) / 2;
red = 0;
green = 0;
if ((*n).score[row][x_position] > black_score) red = (short int) (((*n).score[row][x_position] - black_score) * 256 / (max_score-black_score));
else green = (short int) ((black_score - (*n).score[row][x_position]) * 256 / (black_score-min_score));

/* PRINTS OUT TILE */
fprintf(outfile, " <rect x=\"%i\" y=\"%i\" width=\"%i", (x_position-x_displacement) * tile_width + 150, (row + 1) * tile_height, tile_width);
fprintf(outfile, "\" height=\"");
fprintf(outfile, "%i", tile_height);
fprintf(outfile, "\" style=\"");
fprintf(outfile, "fill:rgb(");
fprintf(outfile, "%i,%i,%i);stroke-width:%i;stroke:rgb(0,0,0)\"/>\n", red, green, blue, stroke);
}
fprintf(outfile, " <text  x=\"%i\" y=\"%i\" fill = \"black\" stroke = \"black\" font-size=\"14\" font-family = \"%s\">%s</text>\n", (x_position-x_displacement) * tile_width + 170, (row + 1) * tile_height + 10, font, (*n).dinucleotide_property_string[row]);
}

/* PRINTS OUT NAME OF MATRIX AND OTHER DATA */
fprintf(outfile, "<text  x=\"");
fprintf(outfile, "%i", x_position * tile_width + 350);
fprintf(outfile, "\" y=\"");
fprintf(outfile, "%i", row * tile_height);
fprintf(outfile, "\" fill = \"black\" stroke = \"black\" font-size=\"30\" font-family = \"");
fprintf(outfile, "%s", font);
fprintf(outfile, "\" >");
fprintf(outfile, "%s", (*n).name);
fprintf(outfile, "</text>\n");
fprintf(outfile, "</g>\n"); 
return(y+ (row + 2) * tile_height);
}


/* SUBROUTINE THAT DRAWS A DINUCLEOTIDE LINE */
short int Dinucleotide_line (FILE *outfile, struct normalized_pwm *n, short int x, short int y, short int line_number)
{
char *dinucleotide_string = "GG";
char *dinucleotide_string_lc = "tt";
char *font = "Courier";
char *forward;
if (rna == 0) forward = dnaforward;
else forward = rnaforward;

short int nucleotide_width = 20;
short int left;
left = x + nucleotide_width * ((*n).pairwise_correlation[line_number].second_base + 0.5);
short int right;
right = x + nucleotide_width * ((*n).pairwise_correlation[line_number].first_base + 0.5);

/* PRINTS OUT POLYLINE */
fprintf(outfile, "<polyline points =\"");
fprintf(outfile, "%i,%i %i,%i %i,%i %i,%i\" ",left,y,left,y-10*line_number,right,y-10*line_number,right,y);
fprintf(outfile, "fill = \"none\" stroke = \"black\" stroke-width = \"3\"/>");

/* PRINTS TEXT */
/* CHANGE IN INFORMATION CONTENT */
fprintf(outfile, " <text  x=\"%i\" y=\"%i\" fill = \"black\" stroke = \"blue\" font-size=\"12\" font-family = \"%s\" >", (right+left) / 2 - 30,y-10*line_number-6, font);
fprintf(outfile, "%.1f%%", (*n).pairwise_correlation[line_number].delta_ic * 100);
fprintf(outfile, "</text>\n");
/* MAX */
fprintf(outfile, " <text  x=\"%i\" y=\"%i\" fill = \"black\" stroke = \"black\" font-size=\"12\" font-family = \"%s\" >", right+5,y-10*line_number, font);
fprintf(outfile, "%c%c %.2f", forward[((*n).pairwise_correlation[line_number].max_dinucleotide & 12) >> 2], forward[((*n).pairwise_correlation[line_number].max_dinucleotide & 3)], (*n).pairwise_correlation[line_number].max_fold_change);
fprintf(outfile, "</text>\n");
/* MIN */
fprintf(outfile, " <text  x=\"%i\" y=\"%i\" fill = \"black\" stroke = \"black\" font-size=\"12\" font-family = \"%s\" >", right+5,y-10*line_number+10, font);
fprintf(outfile, "%c%c %.2f", forward_lc[((*n).pairwise_correlation[line_number].min_dinucleotide & 12) >> 2], forward_lc[((*n).pairwise_correlation[line_number].min_dinucleotide & 3)], (*n).pairwise_correlation[line_number].min_fold_change);
fprintf(outfile, "</text>\n");

return (0);
}

/* SUBROUTINE THAT ADDS NUCLEOTIDE PATHS */
short int Add_nucleotide_paths(FILE *outfile)
{
    fprintf(outfile, "<g display=\"none\"> <path style=\"stroke:none\" fill=\"green\" id=\"A\" d=\"m 6.8910692,-2.720526 -3.7670144,0 -0.3782159,1.030503 0.3782159,0 q 0.5597597,0 0.7942476,0.230832 0.2344959,0.230832 0.2420559,0.618304 0,0.370975 -0.2344879,0.601807 Q 3.6913745,-0.008248 3.1240548,0 L 1.0363116,0 Q 0.47655181,0 0.24205591,-0.230832 0.00756801,-0.461664 0,-0.849135 0,-1.228359 0.24962391,-1.467439 0.49923981,-1.706511 1.0741276,-1.690023 l 2.435703,-6.611698 -1.0136156,0 q -0.5597598,0 -0.7942477,-0.230832 -0.2344959,-0.23084 -0.2420559,-0.618304 0,-0.379223 0.2344879,-0.610055 0.2344959,-0.230832 0.8018157,-0.23908 l 3.3509827,0.0082 3.0862308,8.301721 q 0.544632,0 0.718608,0.131904 0.34796,0.272056 0.34796,0.717232 0,0.370975 -0.234496,0.601807 -0.234488,0.230832 -0.794248,0.23908 l -2.0877433,0 q -0.5597597,0 -0.7942556,-0.230832 -0.2344879,-0.230832 -0.2420559,-0.618303 0,-0.370976 0.2344959,-0.601816 0.2344879,-0.230832 0.8018156,-0.239072 l 0.3782159,0 -0.3706559,-1.030503 z M 6.2481095,-4.410549 4.999998,-7.757618 l -1.2556715,3.347069 2.503783,0 z\" />\n");
    fprintf(outfile, "<path style=\"stroke:none\" id=\"C\" fill=\"blue\" d=\"m 7.9562048,-9.347881 q 0.1733519,-0.204239 0.3740799,-0.306359 0.2007359,-0.10212 0.4470878,-0.10212 0.410584,0 0.666056,0.243512 0.255472,0.243519 0.2646,0.824815 l 0,1.366839 q 0,0.581295 -0.25548,0.824815 -0.255472,0.243512 -0.675176,0.251368 -0.3740879,0 -0.6021918,-0.180672 Q 7.9470768,-6.606355 7.8375889,-7.10125 7.7737489,-7.43117 7.582117,-7.611842 7.2080291,-7.965338 6.5419654,-8.177434 5.8759096,-8.389529 5.1916059,-8.389529 q -0.8394156,0 -1.5419754,0.314215 -0.7025517,0.314216 -1.2408715,1.0212 -0.5383197,0.706983 -0.5383197,1.681046 l 0,1.044767 q 0,1.162599 0.9762716,1.940279 0.9762796,0.777679 2.7281029,0.777679 1.0401436,0 1.7609513,-0.243512 0.4197038,-0.141399 0.8941596,-0.557735 0.2919679,-0.251368 0.4561998,-0.322072 0.16424,-0.07072 0.374088,-0.07856 0.374088,0 0.656936,0.24352 0.282848,0.243512 0.282848,0.57344 0,0.329927 -0.383208,0.706983 -0.556576,0.54988 -1.4324878,0.864095 Q 7.0072932,0 5.5839418,0 3.9233584,0 2.591239,-0.589152 1.5145994,-1.060471 0.75729571,-2.073814 0,-3.087158 0,-4.296885 L 0,-5.38878 Q 0,-6.645643 0.67518371,-7.729682 1.3503675,-8.813729 2.563871,-9.410736 q 1.2135035,-0.597008 2.554743,-0.589152 0.8120396,0 1.5145913,0.157112 0.7025598,0.157104 1.3229995,0.494887 z\" />\n");
    fprintf(outfile, "<path style=\"stroke:none\" id=\"G\" fill=\"orange\" d=\"m 9.2454565,-3.181462 0,2.160255 Q 8.0399288,-0.432048 7.1986567,-0.219952 6.3573851,-0.007856 5.3686815,0 4.007034,0 2.8448665,-0.369208 1.9342188,-0.652 1.413843,-1.084055 0.89346721,-1.516103 0.44248341,-2.356638 -0.00850839,-3.197173 1.6360563e-4,-4.249797 l 0,-1.162607 q 0,-1.610374 1.08411159437,-2.875101 1.4570394,-1.712487 3.9635024,-1.712487 0.7285277,0 1.3789914,0.133544 0.6504638,0.133544 1.2488955,0.400624 0.3555843,-0.290648 0.6938322,-0.290648 0.3902798,0 0.6331198,0.24352 0.24284,0.243519 0.24284,0.824823 l 0,1.044775 q 0,0.581304 -0.24284,0.824823 -0.24284,0.24352 -0.6417918,0.251368 -0.3122239,0 -0.5377202,-0.172816 Q 7.6496485,-6.857811 7.5368966,-7.258442 7.4241526,-7.659066 7.2680407,-7.808322 7.0338728,-8.03613 6.444113,-8.208945 5.8543613,-8.381769 5.1084896,-8.389625 q -1.0580876,0 -1.8386473,0.392775 -0.5550638,0.290648 -1.0233996,1.013352 -0.4683358,0.722703 -0.4683358,1.571094 l 0,1.162607 q 0,1.296143 0.8239197,1.963863 0.8239276,0.667711 2.6365589,0.675567 1.2228715,0 2.2289271,-0.416343 l 0,-1.154752 -1.8126313,0 q -0.6417917,0 -0.9106476,-0.219951 -0.2688639,-0.219952 -0.2775359,-0.58916 0,-0.361352 0.2688559,-0.581303 0.2688639,-0.219952 0.9193276,-0.227808 l 3.1569272,0.008 q 0.6418,0 0.910656,0.219952 0.268856,0.219951 0.2775365,0.581303 0,0.2828 -0.1821365,0.494896 -0.182128,0.212096 -0.572408,0.314215 z\" />\n");
    if (rna == 0) fprintf(outfile, "<path style=\"stroke:none\" id=\"T\" fill=\"red\" d=\"m 5.9510296,-8.301721 0,6.611698 1.2994395,0 q 0.6967997,0 0.9886956,0.230832 0.2919039,0.230832 0.3013199,0.618304 0,0.370975 -0.2919039,0.601807 Q 7.9566848,-0.008248 7.2504691,0 L 2.7306949,0 Q 2.0338952,0 1.7419993,-0.230832 1.4500954,-0.461664 1.4406794,-0.849135 q 0,-0.370976 0.2919039,-0.601816 0.2918959,-0.230832 0.9981116,-0.239072 l 1.2900155,0 0,-6.611698 -2.0903912,0 0,1.599334 q 0,0.610056 -0.2636559,0.865624 -0.2636479,0.255567 -0.70621569,0.263807 -0.4237278,0 -0.6873757,-0.255567 Q 0.00941601,-6.084083 0,-6.702387 l 0,-3.297605 9.9999965,0.0082 0,3.289357 q 0,0.610056 -0.263656,0.865624 -0.263648,0.255567 -0.706216,0.263807 -0.4237279,0 -0.6873758,-0.255567 -0.2636559,-0.25556 -0.2730719,-0.873864 l 0,-1.599334 -2.1186472,0 z\" />\n</g>\n");
    else fprintf(outfile, "<path style=\"stroke:none\" id=\"U\" fill=\"red\" d=\"m 1.0966549,-8.8225952 0,5.1428531 c 0,2.3025166 1.5799292,3.68066910008 3.9219361,3.68066910008 2.6394046,0 3.9405131,-1.31092600008 3.9405131,-3.98318990008 l 0,-4.8403323 c 0.037178,0 0.092936,0 0.130111,0 0.6691451,0 0.9107811,-0.1344542 0.9107811,-0.5882348 0,-0.386554 -0.204461,-0.588235 -0.650557,-0.588235 -0.1487001,0 -0.278812,0 -0.42751,0 l -2.0074259,0 c -0.6691449,0 -0.985131,0.06722 -0.985131,0.588235 0,0.4873935 0.353161,0.5882348 1.1895921,0.5882348 0.1486989,0 0.3159855,0 0.4832717,0 l 0,4.6554585 c 0,1.9999961 -0.6691472,2.9915941 -2.4907082,2.9915941 -0.9107805,0 -1.710039,-0.3865551 -2.1375483,-1.0588229 C 2.5650578,-2.8898264 2.5650574,-3.5452893 2.5650574,-4.318397 l 0,-4.5041982 c 0.1486989,0 0.2788108,0 0.3903349,0 0.7992565,0 1.1338299,-0.1008413 1.1338299,-0.5882348 0,-0.571428 -0.3717481,-0.588235 -1.1338299,-0.588235 l -1.8587374,0 C 0.37174808,-9.999065 0,-9.965445 0,-9.41083 c 0,0.4537806 0.2602238,0.5882348 0.87360628,0.5882348 0.07435,0 0.14869912,0 0.22304862,0 z\" />\n</g>\n");
            
}

struct plot_point {double x; double y;};

short int Add_kmercount_logo_to_Svg (FILE *outfile, struct plot_point *plot_pos, long int current_kmer, long int *****results, short int file_number, short int current_kmer_length, short int current_gap_position, short int current_gap_length, short int local_max, char *legend, short int logo_id, short int count_also_spaced_kmers)
{
short int position;
short int nucleotide_value;
long int total_count;
long int current_count;
long int mutant_kmer;
short int counter;
long int x_position = (*plot_pos).x;
long int y_position = (*plot_pos).y;
char *font = "Courier";
short int tile_height = 30;
short int font_size = 16;
char *forward;
if (rna == 0) forward = dnaforward;
else forward = rnaforward;
long int font_position = 0;
double pwm_position_width = 10;
double pwm_position_height = 25;
double y_margin = 0;
double x_margin = 0;
double swap;
double **order;
    
order = malloc(sizeof(double *) * 4 + 5);
for (counter = 0; counter < 3; counter++) order[counter] = malloc(sizeof(double) * 6 + 5);
char *font_style[] = {"regular", "bold"};
char *font_color[] = {"none", "gray", "black"};
    
fprintf(outfile, "<g> <title>SEED:");
for(position = current_kmer_length-1; position > -1 ; position--) 
{
if(current_kmer_length - position - 1 == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) fprintf(outfile, "n");
fprintf(outfile, "%c", forward[(current_kmer & mask_ULL[1][position]) >> (position * 2)]);
}
fprintf(outfile, "</title>");
    
for (position = current_kmer_length-1; position >= 0; position--)
{
    
/* IDENTIFIES THE FOUR KMER COUNTS FOR THIS POSITION */
for(total_count = 0, nucleotide_value = 0; nucleotide_value < 4; nucleotide_value++)
{
mutant_kmer = ((current_kmer & (~(mask_ULL[1][position]))) | (nucleotide_value << (position * 2))) & lowmask_ULL[current_kmer_length-1];
current_count = results[file_number][current_kmer_length][current_gap_position][current_gap_length][mutant_kmer];
total_count += current_count;
}
        for(nucleotide_value = 0; nucleotide_value < 4; nucleotide_value++)
        {
            mutant_kmer = ((current_kmer & (~(mask_ULL[1][position]))) | (nucleotide_value << (position * 2))) & lowmask_ULL[current_kmer_length-1];
            order[0][nucleotide_value] = nucleotide_value;
            order[1][nucleotide_value] = ((double) results[file_number][current_kmer_length][current_gap_position][current_gap_length][mutant_kmer]) / total_count;
        }
        
        /* BUBBLE SORT */
        for (counter = 0; counter < 3; counter++) 
        {
            for(nucleotide_value = counter; nucleotide_value < 4; nucleotide_value++)
            {
                if (order[1][counter] < order[1][nucleotide_value]) 
                {
                    swap = order[0][counter];
                    order[0][counter] = order[0][nucleotide_value];
                    order[0][nucleotide_value] = swap;
                    swap = order[1][counter];
                    order[1][counter] = order[1][nucleotide_value];
                    order[1][nucleotide_value] = swap;
                }
            }
        }
    
        for(font_position = y_position + 3.5, nucleotide_value = 0; nucleotide_value < 4; nucleotide_value++)
        {
            /* DRAWS CG SEED HIGHLIGHT RECTANGLE */
            if (nucleotide_value == 0 && position > 0 && position != current_gap_position-1) if (((current_kmer & mask_ULL[2][position-1]) >> ((position-1) * 2)) == 6) 
            {
                fprintf(outfile, "<rect x=\"%li\" y=\"%.2f\" width=\"%.2f\" height=\"%.2f\" style=\"fill-opacity:0.1;fill:lime;stroke-width:1;stroke:black\"/>\n", x_position, font_position-1.5, pwm_position_width * 2, pwm_position_height+2);
            }
            /* PRINTS OUT SCALED PATH NUCLEOTIDES */
            fprintf(outfile, "<use xlink:href=\"#%c\" ", forward[(int) order[0][nucleotide_value]]);
            fprintf(outfile, " transform=\"translate(%li,%f) scale(%f,%f)\" visibility=\"visible\" opacity=\"%.2f\"/>\n", x_position, font_position + (order[1][nucleotide_value] * pwm_position_height), ((double) pwm_position_width) / 10, order[1][nucleotide_value] * ((double) pwm_position_height) / 10, 1.0);
            font_position += (order[1][nucleotide_value] * pwm_position_height);
        }
        x_position += pwm_position_width;

    if(current_kmer_length - position == current_gap_position && count_also_spaced_kmers != 0) 
        {
            /* PRINTS OUT Ns */
            fprintf(outfile, "<text  transform=\"scale(1, 2)\" x=\"%li\" y=\"%.2f\" fill = \"%s\" stroke = \"%s\" font-size=\"%i\" font-family = \"%s\">", x_position, (((double) y_position + tile_height - 10) / 2), font_color[(current_gap_length != 0) * 2], "none", font_size/2, font);
            if(current_gap_length < 10) fprintf(outfile, "-");
            fprintf(outfile,"-%iN-</text>\n", current_gap_length);
            x_position += (font_size/2-3) * 5;
        }
    }
    if (legend[0] != '\0') fprintf(outfile, "<text  x=\"%li\" y=\"%li\" text-anchor = \"end\" fill = \"%s\" font-weight=\"%s\" font-size=\"%i\" font-family = \"%s\">%li  %s</text>\n", x_position + 80, y_position + tile_height - 10, font_color[local_max+1], font_style[local_max], font_size, font, results[file_number][current_kmer_length][current_gap_position][current_gap_length][current_kmer], legend);
    (*plot_pos).x = x_position + x_margin;
    (*plot_pos).y = y_position + pwm_position_height + y_margin; 

    for (counter = 0; counter < 3; counter++) free(order[counter]);
    free(order);
    fprintf(outfile, "</g>"); 
    return(0);
}


/* SUBROUTINE THAT GENERATES A LOGO FILE FOR ADM THAT ADDS EDGES TO SHOW PREFERENTIAL PATHS */
short int Svg_logo_ADM (char *filename, long int offset, long int yoffset, struct adjacent_dinucleotide_model *a, double mononucleotide_frequency_cutoff, double log_fold_change_cutoff, double absolute_deviation_cutoff, double gray_dinucleotide_cutoff)
{
    FILE *outfile;
    outfile = fopen (filename, "w");

    short int top_position;
    short int counter;
    short int first;
    short int second;
    double starty;
    double endy;
    short int nucleotide_value;
    short int pwm_position = 0;
    short int warning = 0;
    short int shift_right = 0;
    double **order;
    order = malloc(sizeof(double *) * 4 + 5);
    for (counter = 0; counter < 3; counter++) order[counter] = malloc(sizeof(double) * 6 + 5);
    double **position_memory;
    position_memory = malloc(sizeof(double *) * 4 + 5);
    for (counter = 0; counter < 3; counter++) {position_memory[counter] = malloc(sizeof(double) * 6 + 5);}
    for (counter = 0; counter < 4; counter++) position_memory[1][counter] = 0;
    
    double swap;
    double font_position;
    char *forward;
    if (rna == 0) forward = dnaforward;
    else forward = rnaforward;
    
    char **colors;
    colors = malloc(sizeof(char *) * 7 + 5);
    for (counter = 0; counter < 7; counter++) colors[counter] = malloc(sizeof(char) * 20 + 5);
    strcpy(colors[0], "green");
    strcpy(colors[1], "blue");
    strcpy(colors[2], "orange");
    strcpy(colors[3], "red");
    strcpy(colors[4], "black");
    strcpy(colors[5], "black");
    
    char *font = "Courier";
    fprintf(outfile, "<?xml version=\"1.0\" standalone=\"no\"?> <!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">");
    fprintf(outfile, "<!--%s : command %s -->\n", svgsafe(VERSION), svgsafe(COMMAND));
    
        if(noname == 1)
        {
            fprintf(outfile, "<svg width=\"");
            fprintf(outfile, "%i", (*a).width * 20);
            fprintf(outfile, "\" height=\"100\" x=\"0\" y=\"%i\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">", 300);
        }
        
		else fprintf(outfile, "<svg width=\"2000\" height=\"100\" x=\"0\" y=\"%i\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">", 300);
 
    Add_nucleotide_paths(outfile); /* Adds nucleotide paths */
    
    /* GENERATES LOGO */
top_position = 20;
    
        
        for(pwm_position = 0; pwm_position < (*a).width; pwm_position++)
        {
            
            /* USES ALPHABETIC ORDER (NO SORT), AND SHIFTS BASE GRAPHICAL POSITION MEMORY */
            for (counter = 0; counter < 4; counter++) {order[0][counter] = counter; order[1][counter] = (*a).mononuc_fraction[counter][pwm_position];}
            
            /* PRINTS NUCLEOTIDES */
            for(font_position = 0, nucleotide_value = 0; nucleotide_value < 4; nucleotide_value++)
            {
                if (order[1][nucleotide_value] > 0)
                {
                    /* PRINTS OUT SCALED NUCLEOTIDES */
                    fprintf(outfile, "<use xlink:href=\"#%c\" ", forward[(int) order[0][nucleotide_value]]);
                    fprintf(outfile, " transform=\"translate(%li,%f) scale(2,%f)\" visibility=\"visible\" />\n", pwm_position * 20 + offset, font_position + (order[1][nucleotide_value] * 100) + yoffset, order[1][nucleotide_value] * 10);                    
                    font_position += (order[1][nucleotide_value] * 100);
                }
            }
            
            
            /* DETERMINES IF ADJACENT BASES AFFECT EACH OTHER MORE THAN CUTOFF */
            for(shift_right = 0, first = 0; first < 4; first++) if((*a).mononuc_fraction[first][pwm_position] > mononucleotide_frequency_cutoff) for(second = 0; second < 4; second++) if((*a).mononuc_fraction[second][pwm_position+1] > mononucleotide_frequency_cutoff && (((log10((*a).fraction[first*4+second][pwm_position] / ((*a).mononuc_fraction[second][pwm_position+1] + pseudocount)))) > log_fold_change_cutoff || (*a).fraction[first*4+second][pwm_position] - (*a).mononuc_fraction[second][pwm_position+1] > absolute_deviation_cutoff)) shift_right = 1;
            
            if (shift_right == 1)
            {
            for(starty = 0, first = 0; first < 4; first++) 
            {
            starty += (*a).mononuc_fraction[first][pwm_position] * 50;
            if((*a).mononuc_fraction[first][pwm_position] > mononucleotide_frequency_cutoff) 
            {
                printf("\nBase %c at position %i over cutoff", forward[first], pwm_position);
            for(endy = 0, second = 0; second < 4; second++)
            {
            endy += (*a).mononuc_fraction[second][pwm_position+1] * 50;
                /* DRAWS LINE TO CONNECT PREFERENTIAL PAIR AND SETS SHIFT RIGHT FLAG IF DEVIATION FROM PWM IS DETECTED */
                if((*a).mononuc_fraction[second][pwm_position+1] > mononucleotide_frequency_cutoff) 
                {          
                    printf("\nDinucleotide %c%c at position %i over both cutoffs; cond %.3f vs mono %.3f: log fold %.3f", forward[first], forward[second], pwm_position, (*a).fraction[first*4+second][pwm_position], ((*a).mononuc_fraction[second][pwm_position+1] + pseudocount), (log10((*a).fraction[first*4+second][pwm_position] / ((*a).mononuc_fraction[second][pwm_position+1] + pseudocount))));
                    if(((log10((*a).fraction[first*4+second][pwm_position] / ((*a).mononuc_fraction[second][pwm_position+1] + pseudocount)))) > log_fold_change_cutoff || (*a).fraction[first*4+second][pwm_position] - (*a).mononuc_fraction[second][pwm_position+1] > absolute_deviation_cutoff) {printf("\t ***over LOGFOLD or ABSDEV cutoff"); shift_right = 1;
                        fprintf(outfile, "<polyline points =\"");
                        fprintf(outfile, "%li,%.0f %li,%.0f\" ", pwm_position * 20 + offset +22, starty, pwm_position * 20 + offset +38, endy); 
                        fprintf(outfile, "fill = \"none\" stroke = \"black\" stroke-width = \"%f\"/>\n", (*a).fraction[first*4+second][pwm_position] * (*a).mononuc_fraction[first][pwm_position] * 15);}
                        else 
                            if ((*a).fraction[first*4+second][pwm_position] * (*a).mononuc_fraction[first][pwm_position] > gray_dinucleotide_cutoff) 
                            {
                            fprintf(outfile, "<polyline points =\"");
                            fprintf(outfile, "%li,%.0f %li,%.0f\" ", pwm_position * 20 + offset +22, starty, pwm_position * 20 + offset +38, endy); 
                            fprintf(outfile, "fill = \"none\" stroke = \"lightgrey\" stroke-width = \"%f\"/>\n", (*a).fraction[first*4+second][pwm_position] * (*a).mononuc_fraction[first][pwm_position] * 15);   
                            }
                    
                        
                }
                endy += (*a).mononuc_fraction[second][pwm_position+1] * 50;
            }
            }
            starty += (*a).mononuc_fraction[first][pwm_position] * 50;
            }
            offset += 20;
            }
            
        }
        
        /* PRINTS OUT NAME OF LOGO AND OTHER DATA */
        fprintf(outfile, "<text  x=\"");
        fprintf(outfile, "%li", pwm_position * 20 + 20 + offset);
        fprintf(outfile, "\" y=\"");
        fprintf(outfile, "%i", top_position + 30);
        fprintf(outfile, "\" fill = \"black\" stroke = \"%s\" font-size=\"30\" font-family = \"", colors[4 - 2 * warning + 3 * (warning == 2)]);
        fprintf(outfile, "%s", font);
        fprintf(outfile, "\" >");
        fprintf(outfile, "%s", (*a).name);
        fprintf(outfile, "</text>\n");
    
    fprintf(outfile, "</svg>");
    fclose(outfile);
    return(0);
}

/* SUBROUTINE THAT FINDS BEST AND SECOND BEST SCORE AND POSITION FOR KMER IN A WIDER PWM  RETURN VALUES PLACED IN STRUCTURES o AND o2*/
short int Kmerposandscore (double **pwm, signed short int pwm_width, long long int kmer_sequence_value, signed short int kmer_length, struct oriented_match *o, struct oriented_match *o2)
{
    signed short int nucleotide;
    signed short int rev_nucleotide;
    signed short int kmer_position;
    signed short int current_strand;
    double score = 0;
    double revscore = 0;
    double best_score = -100;
    double second_best_score = -100;
    short int counter;
    signed short int pwm_position;
    signed short int current_position;
    long long int reverse_complement_kmer_sequence_value = 0;
    long long int current_reverse_complement_kmer_sequence_value;
    long long int current_kmer_sequence_value = kmer_sequence_value;
    
    for (counter = 0; counter < kmer_length; counter++) 
    {
        reverse_complement_kmer_sequence_value <<= 2;
        reverse_complement_kmer_sequence_value |= (3 ^ (current_kmer_sequence_value & 3));
        current_kmer_sequence_value >>= 2;
    }
    
    for (pwm_position = 0; pwm_position < pwm_width - kmer_length; pwm_position++)
    {
        
        current_position = pwm_width-1-pwm_position;
        for(kmer_position = 0, score= 0, revscore = 0, current_kmer_sequence_value = kmer_sequence_value, current_reverse_complement_kmer_sequence_value = reverse_complement_kmer_sequence_value; kmer_position < kmer_length; kmer_position++, current_kmer_sequence_value >>= 2, current_reverse_complement_kmer_sequence_value >>= 2) 
        {
            nucleotide = (current_kmer_sequence_value & 3);
            rev_nucleotide = (current_reverse_complement_kmer_sequence_value & 3);
            score += pwm[nucleotide][current_position-kmer_position];
            revscore += pwm[rev_nucleotide][current_position-kmer_position];
        }
        current_strand = 1;
        if (revscore == score) current_strand = 0;
        if (revscore > score) {score = revscore; current_strand = -1;}
        if (score >= best_score) /* SCORES PWM AT THIS POSITION */
        {
            second_best_score = best_score;
            best_score = score;
            (*o2).position = (*o).position;
            (*o2).strand = (*o2).position;
            (*o).position = current_position;
            (*o).strand = current_strand;
        }
        else if (score >= second_best_score)
        {
            second_best_score = score;
            (*o2).position = current_position;
            (*o2).strand = current_strand;            
        }
    }
    (*o).score = best_score;
    (*o2).score = second_best_score;
    return (0);
}

/* FUNCTION FOR QSORTING ACCORDING TO SCORE */
int Sort_according_to_score (const void *a, const void *b)
{
	struct oriented_match *leftseq;
	struct oriented_match *rightseq;
	leftseq = (struct oriented_match *) a;
	rightseq = (struct oriented_match *) b;
	{
		if ((*leftseq).score == (*rightseq).score) return(0);
		if ((*leftseq).score > (*rightseq).score) return(-1); 
		else return (1);
	}
}

/* SUBROUTINE THAT GENERATES FORWARD SEQUENCE VALUE as long int */
long int Generate_sequence_value (char *searchstring)
{
short int query_sequence_length = strlen(searchstring);
long int query_sequence_value ;
short int position;
short int nucleotide_value;
long int position_value;
    
for(query_sequence_value = 0, position = 0, position_value = pow(4,query_sequence_length-1); position < query_sequence_length; position++, position_value /= 4)
{
    for (nucleotide_value = 0; nucleotide_value < 4 && searchstring[position] != dnaforward[nucleotide_value]; nucleotide_value++);
    if(nucleotide_value == 4) {printf("\nERROR IN QUERY SEQUENCE\n"); exit (1);}
    query_sequence_value += position_value * nucleotide_value;
}
return(query_sequence_value);
}

/* SUBROUTINE THAT GENERATES FORWARD SEQUENCE VALUE as ULL*/
__uint128_t sequence_value_ULL (char *current_sequence)
{
    __uint128_t current_sequence_value_ULL;
    short int position;
    short int nucleotide_value;
    __uint128_t position_value;
    __uint128_t left_position_value = pow(4,Nlength-2);
	/* FORWARD STRAND */
	for(current_sequence_value_ULL = 0, position = 0, position_value = left_position_value; position < Nlength-1; position++, position_value /= 4)
	{
		for (nucleotide_value = 0; nucleotide_value < 4 && current_sequence[position] != dnaforward[nucleotide_value]; nucleotide_value++);
		if(nucleotide_value == 4) return ('\0');
		current_sequence_value_ULL += position_value * nucleotide_value;
	}
    return (current_sequence_value_ULL);
}

/* SUBROUTINE THAT DRAWS KMER LINES */
short int Kmerlines (FILE *outfile, char *current_pwm_name, double **pwm, short int pwm_width, short int x, short int y, short int kmer_length, short int number_of_lines_to_draw, double cutoff)
{
short int line_number = 1;
short int current_line;
short int counter;
double score;
short int dash_score;
short int current_kmer_number;
short int nucleotide_width = 20;
short int left = 0;
short int oldleft = 0;
short int right;
short int matches_rejected = 0;
char *font = "Helvetica";
char *style = "italic";
    
char **colors;
colors = malloc(100);
colors[0] = "black";
colors[1] = "gray";
colors[2] = "gray";
colors[3] = "gray";
colors[4] = "gray";
colors[5] = "gray";
colors[6] = "gray";
colors[7] = "gray";
colors[8] = "gray";
colors[9] = "gray";
colors[10] = "gray";
    
char **dash;
dash = malloc(100);
dash[0] = "1, 5";
dash[1] = "2, 5";
dash[2] = "5, 5";
dash[3] = "10, 5";
dash[4] = "20, 5";
dash[5] = "200, 1";
    
signed short int offset=-1;
short int max_allowed_overlap = 3;
    
struct oriented_match *match;
match = malloc (sizeof(struct oriented_match) * 2 * (number_of_tf_kmers + 1) + 5);
for(counter = 0; counter < number_of_tf_kmers * 2 + 1; counter++) oriented_match_init(&match[counter]);

/* finds two best scores for all representative kmers and sorts them */
for (current_kmer_number = 0; current_kmer_number < number_of_tf_kmers; current_kmer_number++) 
{
match[current_kmer_number*2].id = current_kmer_number; 
match[current_kmer_number*2+1].id = current_kmer_number;
score = Kmerposandscore (pwm, pwm_width, tf_kmer_values[current_kmer_number],kmer_length, &match[current_kmer_number*2], &match[current_kmer_number*2+1]); 
/* printf("\nTF %i %s score is %.2f at position %i strand %i ", current_kmer_number, tf_names[current_kmer_number], match[current_kmer_number].score, match[current_kmer_number].position, match[current_kmer_number].strand); */
}
qsort (match, number_of_tf_kmers * 2, sizeof(struct oriented_match), Sort_according_to_score);
printf("\nBest characteristic kmer match for %s is number %i: %s with score %.2f at position %i strand %i\n", current_pwm_name, match[0].id, tf_names[match[0].id], match[0].score, match[0].position, match[0].strand) ; 

left = 0;
/* draws lines */
for (current_line = 0; current_line < number_of_lines_to_draw + matches_rejected; current_line++)
{
oldleft = left;
if (match[current_line].score < cutoff) break;
    
if (match[current_line].strand == 1) left = x + nucleotide_width * (match[current_line].position-5) + 0.5;
else left = x + nucleotide_width * (match[current_line].position-kmer_length+1) + 0.5;
right = left + nucleotide_width * kmer_length;

/* CHECKS IF OVERLAP TOO HIGH, IF YES, TAKES NEXT MODEL UP TO 10 */
if (abs(oldleft - left) < ((kmer_length - max_allowed_overlap) * nucleotide_width)) {matches_rejected++; if (matches_rejected > 8) continue; else break;}
    
/* PRINTS OUT LINE */
if(match[current_line].score > 5) dash_score = 5;
else dash_score = (short int) ((match[current_line].score-3) * 2.5);
if (dash_score < 0) dash_score = 0;
    
if (nocall == 0) 
{
fprintf(outfile, "<polyline points =\"%i,%i %i,%i\" fill = \"none\" stroke = \"%s\" stroke-width = \"1\" stroke-dasharray = \"%s\"/>",left,y+current_line*offset-2*line_number,right,y+current_line*offset-2*line_number, colors[current_line], dash[dash_score]);

/* PRINTS ARROWHEAD */
if (match[current_line].strand == 1) fprintf(outfile, "<polyline points =\"%i,%i %i,%i,  %i,%i\" fill = \"none\" stroke = \"%s\" stroke-width = \"1\"/>",right-4, y+current_line*offset-2*line_number-3, right,y+current_line*offset-2*line_number, right-4, y+current_line*offset-2*line_number+3, colors[current_line]);
if (match[current_line].strand == -1) fprintf(outfile, "<polyline points =\"%i,%i %i,%i,  %i,%i\" fill = \"none\" stroke = \"%s\" stroke-width = \"1\"/>",left+4, y+current_line*offset-2*line_number-3, left,y+current_line*offset-2*line_number, left+4, y+current_line*offset-2*line_number+3, colors[current_line]);
    
/* PRINTS TEXT */
fprintf(outfile, " <text  x=\"%li\" y=\"%i\" fill = \"%s\" font-size=\"11\" font-family = \"%s\" font-style = \"%s\">", (right+left) / 2 - strlen(tf_names[match[current_line].id])*5, y-2*line_number-4+current_line*offset, colors[current_line], font, style);
fprintf(outfile, "%s", tf_names[match[current_line].id]);
fprintf(outfile, "</text>\n");
}
    
}

free(match);
free(colors);
}

/* SUBROUTINE THAT CHECKS IF CpG DINUCLEOTIDE FREQUENCY IS DIFFERENT */
short int CG_difference_call (struct base_dependency_matrix *bd_background, struct base_dependency_matrix *bd_signal, short int position, short int **colorcode)
{
if (position > Nlength * 2 - 2) return (0);
double difference_cutoff = 0.2;

double difference;
short int difference_found = 0;
    
short int dinucleotide;
short int current_position;

colorcode[0][1] = colorcode[1][1];
colorcode[1][1] = 0;
colorcode[0][2] = colorcode[1][2];
colorcode[1][2] = 0;

double background_total_count = 0;
double signal_total_count = 0;
    
for (dinucleotide = 0; dinucleotide < 16; dinucleotide++)
{
background_total_count += (*bd_background).incidence[position+1][position][dinucleotide];
signal_total_count += (*bd_signal).incidence[position+1][position][dinucleotide];
}
/* printf("\nTotal counts: %.1f, %.1f", background_total_count, signal_total_count); */
    
difference = (double) (*bd_signal).incidence[position+1][position][6] / background_total_count - (double) (*bd_background).incidence[position+1][position][6] / signal_total_count;
/* printf("\nCG difference at position %i is %.2f", position, difference); */
if (difference > difference_cutoff) {colorcode[0][1] = 1; colorcode[1][2] = 1; difference_found = 1;} // CpG higher
if (difference + difference_cutoff < 0) {colorcode[0][1] = -1; colorcode[1][2] = -1; difference_found = 1;} // CpG lower     
difference = (double) (*bd_signal).incidence[position+1][position][9] / background_total_count - (double) (*bd_background).incidence[position+1][position][9] / signal_total_count;
/* if (difference > difference_cutoff) {colorcode[0][2] = 1; colorcode[1][1] = 1; difference_found = 1;} // GpC higher, control (G)
if (difference + difference_cutoff < 0) {colorcode[0][2] = -1; colorcode[1][1] = -1; difference_found = 1;} // GpC lower, control (G) */
if (difference_found == 1) return (1);
return(0);  
}

/* SUBROUTINE THAT GENERATES AN SVG LOGO FILE */
short int Svg_logo (char *filename, short int number_of_pwms, struct normalized_pwm **n, struct count_connecting_matrix **cm, struct dinucleotide_properties_matrix *d, struct base_dependency_matrix *bd_background, struct base_dependency_matrix *bd_signal, struct base_dependency_matrix *bd_expected, struct alignscore *all_hits_align_scores, double lambda, short int warning)
{
FILE *outfile;
outfile = fopen (filename, "w");
if (warning > 2 || warning < 0) warning = 2;

short int circles = 0;
short int top_backbone_contacts = 0;
short int top_base_contacts = 0;
short int bottom_backbone_contacts = 0;
short int bottom_base_contacts = 0;
short int current_contact_position;
short int contacts_from_pwm;
short int contacts_defined;
short int color_start;
signed short int dot_offset;
short int CG_difference;    
short int current_pwm = 0;
short int top_position;
short int counter;
short int nucleotide_value;
short int height_scale;
signed short int start_rectangle;
double rectangle_cutoff = 0.25;
short int pwm_position = 0;
short int offset = 0;

short int **colorcode;
colorcode = malloc(sizeof(short int *) * 2 + 5);
colorcode[0] = malloc(sizeof(short int) * 4 + 5);
colorcode[1] = malloc(sizeof(short int) * 4 + 5);
for (counter = 0; counter < 4; counter++) {colorcode[0][counter] = 0; colorcode[1][counter] = 0;}
        
double **order;
order = malloc(sizeof(double *) * 4 + 5);
for (counter = 0; counter < 3; counter++) order[counter] = malloc(sizeof(double) * 6 + 5);
double swap;
double font_position;
char *nucleotide_char;
char *nucleotide_iupac; 
char *nucleotide_bitiupac;
    
if (rna == 1) {nucleotide_char = rnaforward; nucleotide_iupac = dna_iupac; nucleotide_bitiupac = dna_bitiupac;}
else {nucleotide_char = dnaforward; nucleotide_iupac = rna_iupac; nucleotide_bitiupac = rna_bitiupac;}

short int iupac_bits = 0;
short int rgbcolors[5][4];
rgbcolors[0][0] = 0;
rgbcolors[0][1] = 128;
rgbcolors[0][2] = 0;
rgbcolors[1][0] = 0;
rgbcolors[1][1] = 0;
rgbcolors[1][2] = 255;
rgbcolors[2][0] = 255;
rgbcolors[2][1] = 165;
rgbcolors[2][2] = 0;
rgbcolors[3][0] = 255;
rgbcolors[3][1] = 0;
rgbcolors[3][2] = 0;
rgbcolors[4][0] = 211;
rgbcolors[4][1] = 211;
rgbcolors[4][2] = 211;

char **colors;
colors = malloc(sizeof(char *) * 16 + 5);    
for (counter = 0; counter < 15; counter++) colors[counter] = malloc(sizeof(char) * 20 + 5);
strcpy(colors[0], "green");
strcpy(colors[1], "blue");
strcpy(colors[2], "orange");
strcpy(colors[3], "red");
strcpy(colors[4], "black");
strcpy(colors[5], "lightgray");
strcpy(colors[6], "white");
strcpy(colors[7], "cornflowerblue");
strcpy(colors[8], "none");
strcpy(colors[9], "orchid");
strcpy(colors[10], "midnightblue");
strcpy(colors[11], "aliceblue");
strcpy(colors[12], "saddlebrown");
strcpy(colors[13], "moccasin");
    
/* INITIALIZES TF KMER SEQUENCE VALUES */
for (number_of_tf_kmers = 0; strcmp(tf_kmers[number_of_tf_kmers], "END") != 0 ; number_of_tf_kmers++);
tf_kmer_values = malloc(sizeof(char *) * number_of_tf_kmers + 5);
for (counter = 0; counter < number_of_tf_kmers; counter++)
{
tf_kmer_values[counter] = Generate_sequence_value(tf_kmers[counter]);
/*  printf("\n%s %li", tf_kmers[counter], tf_kmer_values[counter]); */
}
    
    
char *font = "Courier";
fprintf(outfile, "<?xml version=\"1.0\" standalone=\"no\"?>\n<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n");
fprintf(outfile, "<!--%s : command %s -->\n", svgsafe(VERSION), svgsafe(COMMAND));

if (cm == '\0') 
{
	if(noname == 1 && number_of_pwms == 1)
	{
		fprintf(outfile, "<svg width=\"");
		fprintf(outfile, "%i", (*n)[0].width * 20);
		fprintf(outfile, "\" height=\"%i\" x=\"0\" y=\"%i\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n", 100 + 130 * contacts, 300 * current_pwm);
	}
	
		else fprintf(outfile, "<svg width=\"2000\" height=\"%i\" x=\"0\" y=\"%i\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n", number_of_pwms * 100 + 130 * contacts, 300 * current_pwm);
}
	else fprintf(outfile, "<svg width=\"2000\" height=\"5000\" x=\"0\" y=\"%i\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n", 300 * current_pwm);
    
if (paths == 1) Add_nucleotide_paths(outfile); /* Adds nucleotide paths */
    
/* GENERATES LOGO */
for(top_position = 20; current_pwm < number_of_pwms; current_pwm++)
{
    
if (current_pwm == 3) if (d != '\0')
{
 offset = Svg_dinucleotide_heatmap(outfile, 0, current_pwm * 120 - 40, d) - 350;
 Dinucleotide_line(outfile, n[current_pwm], 0, current_pwm * 120 - 120,1);
 offset = Svg_base_dependency_heatmap(outfile, (*n[current_pwm]).seed, 0, current_pwm * 120 + offset, bd_signal, bd_expected, n) - 500;
}

if (cm != '\0') fprintf(outfile, "<g id=\"group%i\" transform=\"translate(%i, %i)\" >", current_pwm, 0 + contacts * 20, current_pwm * 120 + offset);
else
fprintf(outfile, "<g id=\"group%i\" transform=\"translate(%i, %i)\" >", current_pwm, ((*(n[current_pwm])).match.position + contacts) * 20, current_pwm * 95 + offset);
    
/* DRAWS KMER LINE */ 
if (current_pwm > 0) Kmerlines (outfile, (*(n[current_pwm])).name, (*n[current_pwm]).fraction, (*n[current_pwm]).width, 0, 0, 6, 2, 2.0);
else printf("\n");
    
for(pwm_position = 0; pwm_position < (*(n[current_pwm])).width; pwm_position++)
{

/* DETERMINES ORDER BY BUBBLE SORT*/
double positive_sum = 0;
double negative_sum = 0;
for (counter = 0; counter < 4; counter++) 
{
order[0][counter] = counter; 
order[1][counter] = (*(n[current_pwm])).fraction[counter][pwm_position];
if(order[1][counter] > 0) positive_sum = positive_sum + order[1][counter];
else negative_sum = negative_sum + order[1][counter];
}
for (counter = 0; counter < 3; counter++) 
{
for(nucleotide_value = counter; nucleotide_value < 4; nucleotide_value++)
{
if (order[1][counter] == (0.0 / 0.0) || order[1][counter] == (1.0 / 0.0)) order[1][counter] = 0;
if (order[1][counter] < order[1][nucleotide_value]) 
{
swap = order[0][counter];
order[0][counter] = order[0][nucleotide_value];
order[0][nucleotide_value] = swap;
swap = order[1][counter];
order[1][counter] = order[1][nucleotide_value];
order[1][nucleotide_value] = swap;
}
}
}

/* CHECKS IF THERE IS DIFFERENCE IN CG OR GC FREQUENCY BETWEEN CONTROL AND SAMPLE */
if (methylCGcompare == 1) CG_difference = CG_difference_call (bd_background, bd_signal, pwm_position, colorcode);
    
for(font_position = 0, nucleotide_value = 0, iupac_bits = 0; nucleotide_value < 4; nucleotide_value++)
{
if (order[1][nucleotide_value] > 0 || (order[1][nucleotide_value] < 0 && (*(n[current_pwm])).negative_values_allowed == 1))
{
if (paths == 0 && barcodelogo == 0)
{
// if (current_pwm == 3) CG_difference = 1;

/* PRINTS OUT SCALED FONT NUCLEOTIDES */
fprintf(outfile, " <text  x=\"%i\" y=\"%f", pwm_position * 20, font_position / (order[1][nucleotide_value] * 4.5) + top_position  - order[1][0] * 0.9);
/* else fprintf(outfile, "%f", font_position); */
fprintf(outfile, "\" fill = \"%s\" stroke=\"", colors[(int) order[0][nucleotide_value]]);
if (colorcode[0][(int) order[0][nucleotide_value]] == 0 || methylCGcompare == 0 || current_pwm > 5) fprintf(outfile, "%s", colors[(int) order[0][nucleotide_value]]);
else
{
if (colorcode[0][(int) order[0][nucleotide_value]] == 1) fprintf(outfile, "%s", "black\" stroke-width = \"1.5");
if (colorcode[0][(int) order[0][nucleotide_value]] == -1) fprintf(outfile, "%s", "black\" opacity = \"0.25");
}
fprintf(outfile, "\" font-size=\"30\" font-family = \"");
fprintf(outfile, "%s\" ", font);
if(contacts == 1 && pwm_align == 1 && current_pwm < 2 && (*(n[current_pwm])).position_score[pwm_position] < rectangle_cutoff) fprintf(outfile, "opacity = \"0.1\" ");
fprintf(outfile, "transform = \"scale(1, ");
fprintf(outfile, "%f", order[1][nucleotide_value] * 4.5);
fprintf(outfile, ")\" >");
fprintf(outfile, "%c", nucleotide_char[(int) order[0][nucleotide_value]]);
fprintf(outfile, "</text>\n");
font_position += (order[1][nucleotide_value] * 90);
}
else if (barcodelogo == 1)
{
    /* PRINTS OUT BARCODE LOGO RECTANGLE */
    if (max_scale_bar == 1) height_scale = 0; else height_scale = nucleotide_value;
    fprintf(outfile, " <rect x=\"%.3f\" y=\"%.2f\" width=\"%.3f", pwm_position * 20 + font_position, (double) top_position + 25 * (1-scale_bars*order[1][height_scale]), order[1][nucleotide_value] * 15);
    fprintf(outfile, "\" height=\"");
    fprintf(outfile, "%.2f", 50*(scale_bars == 0) + 50 * (scale_bars*order[1][height_scale]));
    fprintf(outfile, "\" style=\"");
    if (gray_bars == 0) fprintf(outfile, "fill:%s;stroke-width:%i;stroke:%s\"/>\n", colors[(int) order[0][nucleotide_value]], 1, colors[(int) order[0][nucleotide_value]]);
    else
    {
        fprintf(outfile, "fill:rgb(%.0f,%.0f,%.0f);stroke-width:%i;stroke:rgb(0,0,0)\"/>\n", order[1][nucleotide_value]*rgbcolors[(int) order[0][nucleotide_value]][0]+(1-order[1][nucleotide_value])*rgbcolors[4][0], order[1][nucleotide_value]*rgbcolors[(int) order[0][nucleotide_value]][1]+(1-order[1][nucleotide_value])*rgbcolors[4][1], order[1][nucleotide_value]*rgbcolors[(int) order[0][nucleotide_value]][2]+(1-order[1][nucleotide_value])*rgbcolors[4][2], 0);
    }
    if(order[1][nucleotide_value] > order[1][0] / 2) iupac_bits |= 1 << (int) order[0][nucleotide_value];
    if ((nucleotide_value == 3 || (nucleotide_value < 4 && order[1][nucleotide_value+1] == 0)) && barcodelogo == 1 && barcodelogolabels == 1) fprintf(outfile, " <text  x=\"%i\" y=\"%i\" fill = \"white\" stroke =\"white\" font-size=\"15\" >%c</text>\n", pwm_position * 20 + 2, top_position + 30, nucleotide_bitiupac[iupac_bits]);
    font_position += (order[1][nucleotide_value] * 15);

}
else if (paths == 1)
{
    /* PRINTS OUT SCALED PATH NUCLEOTIDES */
    fprintf(outfile, "<use xlink:href=\"#%c\" ", nucleotide_char[(int) order[0][nucleotide_value]]);
    if ((*(n[current_pwm])).negative_values_allowed == 0) 
    {
    fprintf(outfile, " transform=\"translate(%i,%f) scale(2,%f)\" visibility=\"visible\" />\n", pwm_position * 20 + offset, font_position + (order[1][nucleotide_value] * 100), order[1][nucleotide_value] * 10);
    font_position += (order[1][nucleotide_value] * 100);
    }
    else 
    {
    /* printf("\n%f",order[1][nucleotide_value]); */
    if (order[1][nucleotide_value] < 0) order[1][nucleotide_value] = -order[1][nucleotide_value];
    fprintf(outfile, " transform=\"translate(%i,%f) scale(2,%f)\" visibility=\"visible\" />\n", pwm_position * 20 + offset, (1 - positive_sum) * 50 + font_position + (order[1][nucleotide_value] * 50), order[1][nucleotide_value] * 5);
    font_position += (order[1][nucleotide_value] * 50);
    }

}

}
}


    
/* PRINTS BASE CONTACTS */
if (contacts == 1 && ((pwm_align == 0) || (pwm_align == 1 && current_pwm == 2)))
{
current_contact_position = pwm_position;
dot_offset = 0;
    if (pwm_align == 0) contacts_from_pwm = current_pwm;
    else {contacts_from_pwm = 0; dot_offset--;}
    
    for (color_start = 7; contacts_from_pwm <= current_pwm - pwm_align; contacts_from_pwm++, color_start += 2, dot_offset += 2)
{
    if (pwm_align == 1) current_contact_position = pwm_position + MAX((*(n[0])).width,(*(n[1])).width) - (*(n[contacts_from_pwm])).match.position;
    /* printf("\nPosition : %i", current_contact_position); */
    contacts_defined = 1;
    top_backbone_contacts = 0;
    top_base_contacts = 0;
    bottom_backbone_contacts = 0;
    bottom_base_contacts = 0;
    if (current_contact_position < 0 || current_contact_position > (*(n[contacts_from_pwm])).width) contacts_defined = 0;
    if (contacts_defined == 1)
    {
    top_backbone_contacts = (*(n[contacts_from_pwm])).fraction[4][current_contact_position]+(*(n[contacts_from_pwm])).fraction[5][current_contact_position]+(*(n[contacts_from_pwm])).fraction[6][current_contact_position]+(*(n[contacts_from_pwm])).fraction[7][current_contact_position];
    top_base_contacts = (*(n[contacts_from_pwm])).fraction[8][current_contact_position]+(*(n[contacts_from_pwm])).fraction[9][current_contact_position];
    }
    /* TOP BACKBONE ELLIPSE */
    fprintf(outfile, " <ellipse cx=\"%i\" cy=\"%i\" rx=\"%i\" ry=\"%i\" style=\"opacity:0.5;fill:%s;stroke-width:%i;stroke:%s\"/>\n", pwm_position * 20, top_position + 90, 8, 12, colors[color_start-(top_backbone_contacts==0)], 1, "saddlebrown");
    /* CONNECTOR LINES */
    fprintf(outfile, " <polyline points =\"%i,%i %i,%i %i,%i\" style=\"opacity:1;fill:%s;stroke-width:%i;stroke:%s\"/>\n", pwm_position * 20 - 5 + 10 * (pwm_position == 0), top_position + 110, pwm_position * 20, top_position + 103, pwm_position * 20 + 5, top_position + 110, "none", 1, "gray");
    /* TOP BACKBONE DOTS */
    for (circles = 0; circles < top_backbone_contacts; circles++) {
        fprintf(outfile, " <circle cx=\"%i\" cy=\"%.2f\" r=\"%i\" style=\"fill:%s;stroke-width:%.2f;stroke:%s\"/>\n", pwm_position * 20 + dot_offset, top_position + 84 + (double)circles/top_backbone_contacts * 20 + dot_offset, 2, colors[color_start+4-(circles<(*(n[contacts_from_pwm])).fraction[4][current_contact_position]+(*(n[contacts_from_pwm])).fraction[6][current_contact_position])], 0.5, "black"); }
    /* TOP BASE RECTANGLE */
    fprintf(outfile, " <rect x=\"%i\" y=\"%i\" width=\"%i\" height=\"%i\" style=\"opacity:0.5;fill:%s;stroke-width:%i;stroke:%s\"/>\n", pwm_position * 20 + 2, top_position + 110, 16, 30, colors[color_start-(top_base_contacts==0)], 1, "black");
    /* TOP BASE DOTS */
    for (circles = 0; circles < top_base_contacts; circles++) {
        fprintf(outfile, " <circle cx=\"%i\" cy=\"%.2f\" r=\"%i\" style=\"fill:%s;stroke-width:%.2f;stroke:%s\"/>\n", pwm_position * 20 + 10 + dot_offset, top_position + 118 + (double)circles/top_base_contacts * 20 + dot_offset, 2, colors[color_start+4-(circles<(*(n[contacts_from_pwm])).fraction[8][current_contact_position])], 0.5, "black"); }
    
    if (contacts_defined == 1)
    {
    bottom_backbone_contacts = (*(n[contacts_from_pwm])).fraction[12][current_contact_position]+(*(n[contacts_from_pwm])).fraction[13][current_contact_position]+(*(n[contacts_from_pwm])).fraction[14][current_contact_position]+(*(n[contacts_from_pwm])).fraction[15][current_contact_position];
    bottom_base_contacts = (*(n[contacts_from_pwm])).fraction[10][current_contact_position]+(*(n[contacts_from_pwm])).fraction[11][current_contact_position];
    }
    /* BOTTOM BACKBONE ELLIPSE */
    fprintf(outfile, " <ellipse cx=\"%i\" cy=\"%i\" rx=\"%i\" ry=\"%i\" style=\"opacity:0.5;fill:%s;stroke-width:%i;stroke:%s\"/>\n", pwm_position * 20 + 20, top_position + 195, 8, 12, colors[color_start-(bottom_backbone_contacts==0)], 1, "saddlebrown");
    /* CONNECTOR LINES */
    fprintf(outfile, " <polyline points =\"%i,%i %i,%i %i,%i\" style=\"opacity:1;fill:%s;stroke-width:%i;stroke:%s\"/>\n", pwm_position * 20 + 15, top_position + 175, pwm_position * 20 + 20, top_position + 182, pwm_position * 20 + 25 - 10 * (pwm_position == (*(n[current_pwm])).width - 1), top_position + 175, "none", 1, "gray");
    /* BOTTOM BACKBONE DOTS */
    for (circles = 0; circles < bottom_backbone_contacts; circles++) {
        fprintf(outfile, " <circle cx=\"%i\" cy=\"%.2f\" r=\"%i\" style=\"fill:%s;stroke-width:%.2f;stroke:%s\"/>\n", pwm_position * 20 + 20 + dot_offset, top_position + 201 - (double)circles/bottom_backbone_contacts * 20 - dot_offset, 2, colors[color_start+4-(circles<(*(n[contacts_from_pwm])).fraction[15][current_contact_position]+(*(n[contacts_from_pwm])).fraction[13][current_contact_position])], 0.5, "black"); }
    /* BOTTOM BASE RECTANGLE */
    fprintf(outfile, " <rect x=\"%i\" y=\"%i\" width=\"%i\" height=\"%i\" style=\"opacity:0.5;fill:%s;stroke-width:%i;stroke:%s\"/>\n", pwm_position * 20 + 2, top_position + 145, 16, 30, colors[color_start-(bottom_base_contacts==0)], 1, "black");
    /* BOTTOM BASE DOTS */
    for (circles = 0; circles < bottom_base_contacts; circles++) {
        fprintf(outfile, " <circle cx=\"%i\" cy=\"%.2f\" r=\"%i\" style=\"fill:%s;stroke-width:%.2f;stroke:%s\"/>\n", pwm_position * 20 + 10 + dot_offset, top_position + 167 - (double)circles/bottom_base_contacts * 20 - dot_offset, 2, colors[color_start+4-(circles<(*(n[contacts_from_pwm])).fraction[11][current_contact_position])], 0.5, "black"); }
}
}
    
    
/* PRINTS REPEAT TILES */
if (align_matches == 1 && current_pwm == 2)
{
/* REPEAT */
fprintf(outfile, "<rect x=\"%i\" y=\"%.2f\" width=\"%i\" height=\"%.2f\" style=\"fill:cornflowerblue;stroke:none;stroke-width:2;opacity:0.5\" />", (pwm_position + 5) * 20, top_position + 70 - 0.5 * ((double) all_hits_align_scores[1].score[pwm_position][0])/((double) all_hits_align_scores[1].count[pwm_position][0]+0.01), 20, ((double) all_hits_align_scores[1].score[pwm_position][0])/((double) all_hits_align_scores[1].count[pwm_position][0]+0.01));
/* INVERTED REPEAT */
fprintf(outfile, "<rect x=\"%i\" y=\"%.2f\" width=\"%i\" height=\"%.2f\" style=\"fill:midnightblue;stroke:none;stroke-width:2;opacity:0.5\" />", (pwm_position + 5) * 20, top_position + 80 - 0.5 * ((double) all_hits_align_scores[1].score[pwm_position][1])/((double) all_hits_align_scores[1].count[pwm_position][1]+0.01), 20, ((double) all_hits_align_scores[1].score[pwm_position][1])/((double) all_hits_align_scores[1].count[pwm_position][1]+0.01));
    
}
}

    
    
if((*(n[current_pwm])).negative_values_allowed == 1)
{   /* PRINTS OUT POLYLINE */
    fprintf(outfile, "<polyline points =\"");
    fprintf(outfile, "%i,%i %i,%i \" ", offset,50,pwm_position * 20 + offset,50);
    fprintf(outfile, "fill = \"none\" stroke = \"black\" stroke-width = \"0.5\"/>");
}
    
/* PRINTS OUT NAME OF LOGO AND OTHER DATA */
if (noname == 0)
{
fprintf(outfile, "<text  x=\"");
fprintf(outfile, "%i", pwm_position * 20 + 20);
fprintf(outfile, "\" y=\"");
fprintf(outfile, "%i", top_position + 30);
fprintf(outfile, "\" fill = \"black\" stroke = \"%s\" font-size=\"30\" font-family = \"", colors[4 - 2 * warning + 3 * (warning == 2)]);
fprintf(outfile, "%s", font);
fprintf(outfile, "\" >");
fprintf(outfile, "%s", (*(n[current_pwm])).name);
if (current_pwm == 0 || current_pwm == 4) fprintf (outfile, " lambda=%.3f", lambda);
fprintf(outfile, "</text>\n");
}
    
/* PRINTS ALIGNMENT BOX IF TWO LOGOS ARE ALIGNED */
if (cm == '\0' && contacts == 0)
{
for(start_rectangle = -1, pwm_position = 0; pwm_position < (*(n[current_pwm])).width; pwm_position++)
{
/* printf("\nBase %i score = %.2f", pwm_position, (*(n[current_pwm])).position_score[pwm_position]); */
if ((*(n[current_pwm])).position_score[pwm_position] > rectangle_cutoff && start_rectangle == -1) {start_rectangle = pwm_position; /* printf("\nStarts rectangle at %i", pwm_position); */ continue;}
if (start_rectangle != -1 && ((*(n[current_pwm])).position_score[pwm_position] < rectangle_cutoff || pwm_position == (*(n[current_pwm])).width)) 
{
/* printf("\nEnds rectangle at %i", pwm_position); */
fprintf(outfile, "<rect x=\"%i\" y=\"%i\" width=\"%i\" height=\"%i\" style=\"fill:cornflowerblue;stroke:none;stroke-width:2;opacity:0.1\" />", start_rectangle * 20 - 1, top_position - 18 - 95 * (current_pwm == 2), (pwm_position-start_rectangle) * 20, 187);
fprintf(outfile, "<rect x=\"%i\" y=\"%i\" width=\"%i\" height=\"%i\" style=\"fill:none;stroke:midnightblue;stroke-width:2;opacity:1\" />", start_rectangle * 20 - 1, top_position - 18 - 95 * (current_pwm == 2), (pwm_position-start_rectangle) * 20, 187);
start_rectangle = -1;
}
}
if (start_rectangle != -1) 
{
fprintf(outfile, "<rect x=\"%i\" y=\"%i\" width=\"%i\" height=\"%i\" style=\"fill:cornflowerblue;stroke:none;stroke-width:2;opacity:0.1\" />", start_rectangle * 20 - 1, top_position - 18 - 95 * (current_pwm == 2), (pwm_position-start_rectangle) * 20, 187);    
fprintf(outfile, "<rect x=\"%i\" y=\"%i\" width=\"%i\" height=\"%i\" style=\"fill:none;stroke:midnightblue;stroke-width:2;opacity:1\" />", start_rectangle * 20 - 1, top_position - 18 - 95 * (current_pwm == 2), (pwm_position-start_rectangle) * 20, 187); 
}
}
    
fprintf(outfile, "</g>\n");
}

if (cm != '\0')
{
Svg_heatmap(outfile, 180, current_pwm * 120 + 50 + offset, -1, -1, -1, cm[0]);
Svg_heatmap(outfile, 180, current_pwm * 120 + 200 + offset, -1, -1, -1, cm[1]);
Svg_heatmap(outfile, 180, current_pwm * 120 + 350 + offset, 0, 1, -1, cm[2]);
}

free(tf_kmer_values);
free(colors);
free(order);
    
fprintf(outfile, "</svg>\n");
fclose(outfile);
return(0);
}

/* SUBROUTINE THAT GENERATES PWM FROM IUPAC */
short int Iupac_to_pwm(struct normalized_pwm *n, char *searchstring)
{
short int counter;
short int pwm_position;
short int current_match_position;
short int nucleotide_value;
(*n).width =  strlen(searchstring);

char **canbe;
canbe = malloc(sizeof(char *) * 4 + 5);
for (counter = 0; counter < 4; counter++) canbe[counter] = malloc(200);
strcpy (canbe[0], "100010101001111");
strcpy (canbe[1], "010001100110111");
strcpy (canbe[2], "001010010111011");
strcpy (canbe[3], "000101011011101");

char *conversion_string;
conversion_string = malloc(10);
conversion_string[1] = '\0';

char *nucleotide_iupac = dna_iupac;
if (rna == 1) nucleotide_iupac = rna_iupac;

/* BUILDS PWM */
for(pwm_position = 0; pwm_position < (*n).width; pwm_position++)
{
printf("\n");
for(current_match_position = 0; (current_match_position < iupac_length) & (searchstring[pwm_position] != nucleotide_iupac[current_match_position]); current_match_position++);
if(current_match_position == iupac_length) {printf("\nDefective IUPAC");exit(1);}
for (nucleotide_value = 0; nucleotide_value < 4; nucleotide_value++) 
{
conversion_string[0] = canbe[nucleotide_value][current_match_position];
(*n).fraction[nucleotide_value][pwm_position] = atof(conversion_string) - 1;
printf("%s\t", conversion_string);
}
}

for (counter = 0; counter < 4; counter++) free(canbe[counter]);
free(conversion_string);
return(0);
}

/* SUBROUTINE THAT GENERATES PWM FOR MULTINOMIAL */
short int Multinomial_pwm(struct normalized_pwm *n, char *searchstring)
{
short int counter;
short int position;
short int nucleotide_value;
(*n).width =  strlen(searchstring);
char *forward;
if (rna == 0) forward = dnaforward;
else forward = rnaforward;

for (position = 0; position < (*n).width; position++)
{
for (nucleotide_value = 0; nucleotide_value < 4; nucleotide_value++) (*n).fraction[nucleotide_value][position] = - 1;
for (nucleotide_value = 0; nucleotide_value < 5 && searchstring[position] != forward[nucleotide_value]; nucleotide_value++);
if (nucleotide_value == 5) return ('\0');
if (nucleotide_value == 4) for (nucleotide_value = 0; nucleotide_value < 4; nucleotide_value++) (*n).fraction[nucleotide_value][position] = 0; /* N */
else (*n).fraction[nucleotide_value][position] = 0;
}
return(0);
}

/* NOT READY */
int Invert_pwm(struct normalized_pwm *bqp, struct normalized_pwm *qp, short int multinomial)
{
return(0);
}

/* SUBROUTINE THAT CONVERTS COUNT PWM TO NORMALIZED PWM (ZEROES NEGATIVE VALUES, ROWS IN EACH COLUMN ADD TO 1) */
short int Count_to_normalized_pwm (struct normalized_pwm *n, struct count_pwm *c)
{
short int counter;
short int position;
double total_nucleotides = 0;
double count_value;
(*n).width = (*c).width;
for (position = 0; position < (*c).width; position++)
{
for (counter = 0, total_nucleotides = 0; counter < 4; counter++) if ((*c).incidence[counter][position] > 0) total_nucleotides += (*c).incidence[counter][position];
for ((*n).information_content[position] = 0, (*n).total_counts_for_column[position] = 0, counter = 0; counter < 4; counter++) 
{
count_value = (*c).incidence[counter][position];
if (count_value < 0) count_value = 0;
(*n).fraction[counter][position] = count_value / total_nucleotides;
(*n).information_content[position] -= -log((*n).fraction[counter][position])*(*n).fraction[counter][position];
(*n).total_counts_for_column[position] += count_value;
}
}
return (0);
}

/* SUBROUTINE THAT COPIES NORMALIZED PWM */
short int Copy_normalized_pwm (struct normalized_pwm *n, struct normalized_pwm *c)
{
	short int counter;
	short int position;
	double total_nucleotides = 0;
	(*n).width = (*c).width;
	(*n).max_counts = (*c).max_counts;
	for (position = 0; position < (*c).width; position++)
	{
		for (counter = 0; counter < 4 + contacts * 12; counter++) (*n).fraction[counter][position] = (*c).fraction[counter][position];
		(*n).information_content[position] = (*c).information_content[position];
		(*n).total_counts_for_column[position] = (*c).total_counts_for_column[position];
	}
	return (0);
}

/* SUBROUTINE THAT SUBTRACTS NORMALIZED c FROM NORMALIZED PWM n */
short int Subtract_normalized_pwm (struct normalized_pwm *n, struct normalized_pwm *c, short int offset)
{
	short int counter;
	short int position;
	double total_nucleotides = 0;
	if((*c).width - offset < (*n).width) (*n).width = (*c).width - offset;
	for (position = 0; position < (*n).width; position++)
	{
		for (counter = 0; counter < 4; counter++) (*n).fraction[counter][position] -= (*c).fraction[counter][position+offset];
	}
	return (0);
}

/* SUBROUTINE THAT REVERSE COMPLEMENTS NORMALIZED PWM (NOT TESTED) */
short int Reverse_complement_normalized_pwm (struct normalized_pwm *rc, struct normalized_pwm *f)
{
	short int counter;
	short int position;
	short int width = (*f).width-1;
	double total_nucleotides = 0;
	(*rc).width = (*f).width;
	(*rc).max_counts = (*f).max_counts;
    strcpy((*rc).name, "rc_");
    strcat((*rc).name, (*f).name);
    
	for (position = 0; position < (*f).width; position++)
	{
		/* REVERSE COMPLEMENTS BASE FREQUENCIES */
        for (counter = 0; counter < 4; counter++) (*rc).fraction[3-counter][width-position] = (*f).fraction[counter][position];
		(*rc).information_content[width-position] = (*f).information_content[position];
		(*rc).original_position[width-position] = (*f).original_position[position];
        (*rc).position_score[width-position] = (*f).position_score[position];
        (*rc).total_counts_for_column[width-position] = (*f).total_counts_for_column[position];
        
        /* REVERSE COMPLEMENTS BASE CONTACTS */
        if (contacts == 1) for (counter = 4; counter < 16; counter++) (*rc).fraction[19-counter][width-position] = (*f).fraction[counter][position]; 
            
	}
	return (0);
}

/* SUBROUTINE THAT RENORMALIZES NORMALIZED PWM (ROWS IN EACH COLUMN ADD TO 1) */
short int Normalize_pwm (struct normalized_pwm *n)
{
	short int counter;
	short int position;
	double total_nucleotides = 0;
	double normalized_value = 0;
	for (position = 0; position < (*n).width; position++)
	{
		for (counter = 0, total_nucleotides = 0; counter < 4; counter++) 
        {
        if ((*n).fraction[counter][position] > 0) total_nucleotides += (*n).fraction[counter][position]; 
        else if ((*n).negative_values_allowed == 1) total_nucleotides += -(*n).fraction[counter][position];
        }
        
		for (counter = 0; counter < 4; counter++) 
		{
		normalized_value = ((double) (*n).fraction[counter][position]) / total_nucleotides;
		if ((normalized_value < 0) && ((*n).negative_values_allowed == 0)) normalized_value = 0;
		(*n).fraction[counter][position] = normalized_value;
		}
	}
	return (0);
}

/* GENERATE NORMALIZED CONNECTING MATRIX FROM COUNT CONNECTING MATRIX */
short int Count_to_normalized_connecting_matrix (struct normalized_connecting_matrix *n, struct count_connecting_matrix *c)
{
short int counter;
short int position;
long int total_occurrence = 0;
double normalized_value;
(*n).width = (*c).width;
for (counter = 0; counter < 3; counter++)
{
for (position = 1, total_occurrence = 0; position < (*c).width; position++) if ((*c).incidence[counter][position] > 0) total_occurrence += (*c).incidence[counter][position];
(*n).orientation_count[counter] = total_occurrence;
for (position = 1; position < (*c).width; position++) 
{
normalized_value = (double) (*c).incidence[counter][position] / (double) total_occurrence;
if(normalized_value < 0) normalized_value = 0;
(*n).fraction[counter][position] = normalized_value;
}
}
for (counter = 0, total_occurrence = 0; counter < 3; counter++) total_occurrence += (*n).orientation_count[counter];
for (counter = 0; counter < 3; counter++) {(*n).orientation_fraction[counter] = (double) (*n).orientation_count[counter] / (double) total_occurrence; (*n).fraction[counter][0] = (*n).orientation_fraction[counter];}
return (0);
}

/* FUNCTION FOR QSORTING ACCORDING TO SEQUENCE VALUE (ULL) */
int Sort_according_to_sequence_value_ULL (const void *a, const void *b)
{
	struct sequence_incidence_table *leftseq;
	struct sequence_incidence_table *rightseq;
	leftseq = (struct sequence_incidence_table *) a;
	rightseq = (struct sequence_incidence_table *) b;
	
	/* if ((*leftseq).incidence < 100 && (*rightseq).incidence < 100) */
	{
		if ((*leftseq).sequence_value_ULL == (*rightseq).sequence_value_ULL) return(0);
		if ((*leftseq).sequence_value_ULL > (*rightseq).sequence_value_ULL) return(1); 
		else return (-1);
	}
	/* if ((*leftseq).incidence > (*rightseq).incidence) return(1); 
	else return(-1); */
}

/* FUNCTION FOR QSORTING ACCORDING TO DOUBLE */
int Sort_according_to_double (const void *a, const void *b)
{
	double *leftseq;
	double *rightseq;
	leftseq = (double *) a;
	rightseq = (double *) b;
	
	/* if ((*leftseq).incidence < 100 && (*rightseq).incidence < 100) */
	{
		if ((*leftseq) == (*rightseq)) return(0);
		if ((*leftseq) > (*rightseq)) return(1); 
		else return (-1);
	}
	/* if ((*leftseq).incidence > (*rightseq).incidence) return(1); 
     else return(-1); */
}

/* FUNCTION FOR QSORTING ACCORDING TO PREFERRED KMER RANK */
int Sort_according_to_preferred (const void *a, const void *b)
{
	struct kmer_incidence_table *leftseq;
	struct kmer_incidence_table *rightseq;
	leftseq = (struct kmer_incidence_table *) a;
	rightseq = (struct kmer_incidence_table *) b;
	
	/* if ((*leftseq).incidence < 100 && (*rightseq).incidence < 100) */
	{
		if ((*leftseq).preferred == (*rightseq).preferred) return(0);
		if ((*leftseq).preferred < (*rightseq).preferred) return(1); 
		else return (-1);
	}
	/* if ((*leftseq).incidence > (*rightseq).incidence) return(1); 
     else return(-1); */
}

/* FUNCTION FOR QSORTING ACCORDING TO INCIDENCE */
int Sort_according_to_incidence (const void *a, const void *b)
{
	struct kmer_incidence_table *leftseq;
	struct kmer_incidence_table *rightseq;
	leftseq = (struct kmer_incidence_table *) a;
	rightseq = (struct kmer_incidence_table *) b;
	
	/* if ((*leftseq).incidence < 100 && (*rightseq).incidence < 100) */
	{
		if ((*leftseq).incidence == (*rightseq).incidence) return(0);
		if ((*leftseq).incidence > (*rightseq).incidence) return(-1); 
		else return (1);
	}
	/* if ((*leftseq).incidence > (*rightseq).incidence) return(1); 
     else return(-1); */
}

/* FUNCTION FOR QSORTING ACCORDING TO HIGH INCIDENCE LOCAL MAX */
int Sort_according_to_local_max (const void *a, const void *b)
{
	struct kmer_incidence_table *leftseq;
	struct kmer_incidence_table *rightseq;
	leftseq = (struct kmer_incidence_table *) a;
	rightseq = (struct kmer_incidence_table *) b;
	{
		if ((*leftseq).local_max == (*rightseq).local_max) return(0);
		if ((*leftseq).local_max > (*rightseq).local_max) return(-1); 
		else return (1);
	}
}

/* FUNCTION FOR QSORTING ACCORDING TO ENRICHMENT */
int Sort_according_to_enrichment (const void *a, const void *b)
{
	struct kmer_incidence_table *leftseq;
	struct kmer_incidence_table *rightseq;
	leftseq = (struct kmer_incidence_table *) a;
	rightseq = (struct kmer_incidence_table *) b;
	{
		if ((*leftseq).max_enrichment == (*rightseq).max_enrichment) return(0);
		if ((*leftseq).max_enrichment > (*rightseq).max_enrichment) return(-1); 
		else return (1);
	}
}

/* FUNCTION FOR QSORTING ACCORDING TO MAX LOCAL MAX ENRICHMENT */
int Sort_according_to_local_max_enrichment (const void *a, const void *b)
{
	struct kmer_incidence_table *leftseq;
	struct kmer_incidence_table *rightseq;
	leftseq = (struct kmer_incidence_table *) a;
	rightseq = (struct kmer_incidence_table *) b;
	{
		if ((*leftseq).local_max_enrichment == (*rightseq).local_max_enrichment) return(0);
		if ((*leftseq).local_max_enrichment > (*rightseq).local_max_enrichment) return(-1); 
		else return (1);
	}
}

/* FUNCTION FOR QSORTING ACCORDING TO MAX LOCAL MAX */
int Sort_according_to_local_max_max_incidence (const void *a, const void *b)
{
	struct kmer_incidence_table *leftseq;
	struct kmer_incidence_table *rightseq;
	leftseq = (struct kmer_incidence_table *) a;
	rightseq = (struct kmer_incidence_table *) b;
	{
		if ((*leftseq).local_max_max_incidence == (*rightseq).local_max_max_incidence) return(0);
		if ((*leftseq).local_max_max_incidence > (*rightseq).local_max_max_incidence) return(-1); 
		else return (1);
	}
}

/* FUNCTION FOR QSORTING ACCORDING TO ANY LOCAL MAX */
int Sort_according_to_any_local_max (const void *a, const void *b)
{
	struct kmer_incidence_table *leftseq;
	struct kmer_incidence_table *rightseq;
	leftseq = (struct kmer_incidence_table *) a;
	rightseq = (struct kmer_incidence_table *) b;
	{
		if ((*leftseq).any_local_max == (*rightseq).any_local_max) return(0);
		if ((*leftseq).any_local_max > (*rightseq).any_local_max) return(-1); 
		else return (1);
	}
}

/* FUNCTION FOR QSORTING ACCORDING TO COUNT */
int Sort_according_to_count (const void *a, const void *b)
{
	struct countpair *leftseq;
	struct countpair *rightseq;
	leftseq = (struct countpair *) a;
	rightseq = (struct countpair *) b;
	
	/* if ((*leftseq).incidence < 100 && (*rightseq).incidence < 100) */
	{
		if ((*leftseq).sum == (*rightseq).sum) return(0);
		if ((*leftseq).sum > (*rightseq).sum) return(-1); 
		else return (1);
	}
	/* if ((*leftseq).incidence > (*rightseq).incidence) return(1); 
     else return(-1); */
}


/* FUNCTION FOR QSORT */
int Numeric_sort_long_ints (const void *a, const void *b)
{
	long int *leftseq;
	long int *rightseq;
	leftseq = (long int *) a;
	rightseq = (long int *) b;
	
	/* if ((*leftseq).incidence < 100 && (*rightseq).incidence < 100) */
		if (*leftseq == *rightseq) return(0);
		if (*leftseq > *rightseq) return(1); 
		else return (-1);
}

/* SUBROUTINE THAT ADDS A SEQUENCE TO PWM */
short int Add_to_pwm (struct count_pwm *p, short int match_position, __uint128_t sequence_value_ULL, struct normalized_pwm *background_pwm)
{
/* printf("\ncalls pwm"); */
short int pwm_position;
short int background_position;
(*p).max_counts++;
    
if (background_pwm == '\0')
{
for (pwm_position = (*p).width / 2 - 2 + Nlength - match_position; pwm_position > (*p).width / 2 - 2 - match_position + 1; pwm_position--, sequence_value_ULL >>= 2) (*p).incidence[sequence_value_ULL & 3][pwm_position]++;
}
else
{
for (background_position = Nlength - 2, pwm_position = (*p).width / 2 - 2 + Nlength - match_position; pwm_position > (*p).width / 2 - 2 - match_position + 1; pwm_position--, background_position--, sequence_value_ULL >>= 2) 
{
(*p).incidence[sequence_value_ULL & 3][pwm_position] += 0.25 / (*background_pwm).fraction[sequence_value_ULL & 3][background_position];
/* if (pwm_position == 10 && sequence_value_ULL & 3 == 3)
{
printf("\n10T");
} */

/* printf ("\nMatch position %i, Background_position %i, Nucleotide value %i, Background value %f", match_position, background_position, sequence_value_ULL & 3, (*background_pwm).fraction[sequence_value_ULL & 3][background_position]); */
}
}
return(0);
}

/* SUBROUTINE THAT ADDS A SEQUENCE TO PWM BUT EXCLUDES CONSENSUS NUCLEOTIDE IF MUTATION AT OTHER SITE (MULTINOMIAL N DISTRIBUTION) */
short int Multinomial_add_to_pwm (struct count_pwm *p, struct normalized_pwm *qp, short int match_position, double score, double cut_off, __uint128_t sequence_value_ULL, struct normalized_pwm *background_pwm)
{
/* printf("\ncalls pwm"); */
short int pwm_position;
short int background_position;
short int query_position = 0;    
    
(*p).max_counts++;
if (background_pwm == '\0')
{
for (pwm_position = (*p).width / 2 - 2 + Nlength - match_position; pwm_position > (*p).width / 2 - 2 - match_position + 1; pwm_position--, sequence_value_ULL >>= 2) (*p).incidence[sequence_value_ULL & 3][pwm_position]++;
}
else
{
for (background_position = Nlength - 2, pwm_position = (*p).width / 2 - 2 + Nlength - match_position; pwm_position > (*p).width / 2 - 2 - match_position + 1; pwm_position--, background_position--, sequence_value_ULL >>= 2) 
{
/* CHECK IF NUCLEOTIDE WITHIN QUERY */
query_position = background_position - match_position + 1;

if (query_position >= 0 && query_position < (*qp).width)
{
/* EXCLUDE IF SCORE IS JUST AT CUT-OFF AND NUCLEOTIDE IS THE CONSENSUS BASE */
if ((*qp).fraction[sequence_value_ULL & 3][query_position] != 0 || score > cut_off + 1)
{
(*p).incidence[sequence_value_ULL & 3][pwm_position] += 0.25 / (*background_pwm).fraction[sequence_value_ULL & 3][background_position];
/* printf ("\nBase %c counted due to PWM mismatch at position %i, query_position %i", forward[sequence_value_ULL & 3], background_position, query_position); */
}
/* else printf ("\nBase %c EXCLUDED at background position %i, query_position %i, PWM score %.1f, score %.1f, cut-off %f", forward[sequence_value_ULL & 3], background_position, query_position, (*qp).fraction[sequence_value_ULL & 3][query_position], score, cut_off); */
}
else 
{
/* printf ("\nBase %c counted: outside query at position %i, query_position %i", forward[sequence_value_ULL & 3], background_position, query_position); */
(*p).incidence[sequence_value_ULL & 3][pwm_position] += 0.25 / (*background_pwm).fraction[sequence_value_ULL & 3][background_position];
}
/* printf ("\nMatch position %i, Background_position %i, Nucleotide value %i, Background value %f", match_position, background_position, sequence_value_ULL & 3, (*background_pwm).fraction[sequence_value_ULL & 3][background_position]); */
}
}
return(0);
}

/* SUBROUTINE THAT ADDS A DINUCLEOTIDE PROPERTIES TO DINUCLEOTIDE PROPERTIES MATRIX (DOES NOT NORMALIZE, ADDS TO SUM AND INCREASES COUNTS BY ONE */
short int Add_to_dinucleotide_matrix (struct dinucleotide_properties_matrix *p, struct dinucleotide_properties *d, struct base_dependency_matrix *dep, short int match_position, __uint128_t sequence_value_ULL)
{
/* SCORES DINUCLEOTIDE PROPERTIES */

/* printf("\ncalls add_to_di"); */
short int pwm_position;
short int gap_position;
short int property;
short int background_position;
__uint128_t current_sequence_value_ULL;
__uint128_t deleted_sequence_value_ULL;
(*p).max_counts++;
for (property = 0; property < (*p).number_of_dinucleotide_properties; property++)
{
for (pwm_position = (*p).width / 2 - 2 + Nlength - match_position, current_sequence_value_ULL = sequence_value_ULL; pwm_position > (*p).width / 2 - 2 - match_position + 2; pwm_position--, current_sequence_value_ULL >>= 2) 
{
(*p).score[property][pwm_position] += (*d).dinucleotide_property[property][current_sequence_value_ULL & 15];
(*p).count[property][pwm_position]++;
/* printf("\n%i, %i, %i", property, pwm_position, (*p).count[property][pwm_position]); */
}
}

/* COUNTS DINUCLEOTIDES WITH ALL POSSIBLE GAPS */
for (pwm_position = (*p).width / 2 - 2 + Nlength - match_position, current_sequence_value_ULL = sequence_value_ULL; pwm_position > (*p).width / 2 - 2 - match_position + 2; pwm_position--, current_sequence_value_ULL >>= 2) 
{
for (gap_position = pwm_position - 1, deleted_sequence_value_ULL = current_sequence_value_ULL; gap_position > (*p).width / 2 - 2 - match_position + 1; gap_position--, deleted_sequence_value_ULL >>= 2) 
{
((*dep).incidence)[pwm_position][gap_position][(current_sequence_value_ULL & 3) + (deleted_sequence_value_ULL & 12)]++;
/* printf("\n%i, %i, %i", pwm_position, gap_position, ((current_sequence_value_ULL & 3) + (deleted_sequence_value_ULL & 12))); */
}
}

return(0);
}

/* SUBROUTINE THAT ADDS A DINUCLEOTIDE PROPERTIES TO DINUCLEOTIDE PROPERTIES MATRIX (DOES NOT NORMALIZE, ADDS TO SUM AND INCREASES COUNTS BY ONE */
short int Multinomial_add_to_dinucleotide_matrix (struct dinucleotide_properties_matrix *p, struct dinucleotide_properties *d, struct base_dependency_matrix *dep, struct normalized_pwm *qp, short int match_position, double score, double cut_off, __uint128_t sequence_value_ULL)
{
    /* SCORES DINUCLEOTIDE PROPERTIES */
    
    /* printf("\ncalls add_to_di, score =%f", score); */
    short int pwm_position;
    short int gap_position;
    short int property;
    short int counter;
    short int *total_score;
    total_score = malloc(sizeof(short int) * Nlength + 5);
    short int background_position;
    short int background_position2;
    short int current_matches_to_consensus;
    __uint128_t current_sequence_value_ULL;
    __uint128_t deleted_sequence_value_ULL;
    short int query_position;
    
/* CALCULATES TOTAL SCORE FOR EACH QUERY POSITION (to identify Ns)*/
    for (query_position = 0; query_position < (*qp).width; query_position++) for (counter = 0, total_score[query_position] = 0; counter < 4; counter++) total_score[query_position] += (*qp).fraction[counter][query_position];
   
    query_position = 0;
    (*p).max_counts++;
    
    for (property = 0; property < (*p).number_of_dinucleotide_properties; property++)
    {
        for (background_position = Nlength - 2, pwm_position = (*p).width / 2 - 2 + Nlength - match_position, current_sequence_value_ULL = sequence_value_ULL; pwm_position > (*p).width / 2 - 2 - match_position + 2; background_position--, pwm_position--, current_sequence_value_ULL >>= 2) 
        {
            /* FIGURES OUT IF THE DINUCLEOTIDE COMPONENTS MATCH CONSENSUS, IF THEY DO, DOES NOT COUNT IF SCORE NEAR CUTOFF; DOES NOT CONSIDER A MATCH TO N AS MATCH TO 'CONSENSUS' */
            current_matches_to_consensus = 0;
            query_position = background_position - match_position + 1;
            if (query_position >= 0 && query_position < (*qp).width) current_matches_to_consensus = ((*qp).fraction[current_sequence_value_ULL & 3][query_position] == 0 && total_score[query_position] != 0);
            query_position = background_position - match_position + 2;
            if (query_position >= 0 && query_position < (*qp).width) current_matches_to_consensus += ((*qp).fraction[(current_sequence_value_ULL & 12) >> 2][query_position] == 0 && total_score[query_position] != 0);
            
            if(score > cut_off + current_matches_to_consensus)
            {
            (*p).score[property][pwm_position] += (*d).dinucleotide_property[property][current_sequence_value_ULL & 15];
            (*p).count[property][pwm_position]++;
            }
            /* printf("\n%i, %i, %i", property, pwm_position, (*p).count[property][pwm_position]); */
        }
    }
    
    /* COUNTS DINUCLEOTIDES WITH ALL POSSIBLE GAPS */
    for (background_position = Nlength - 2, pwm_position = (*p).width / 2 - 2 + Nlength - match_position, current_sequence_value_ULL = sequence_value_ULL; pwm_position > (*p).width / 2 - 2 - match_position + 2; background_position--, pwm_position--, current_sequence_value_ULL >>= 2) 
    {
        for (background_position2 = background_position - 1, gap_position = pwm_position - 1, deleted_sequence_value_ULL = current_sequence_value_ULL; gap_position > (*p).width / 2 - 2 - match_position + 1; background_position2--, gap_position--, deleted_sequence_value_ULL >>= 2) 
        {
            /* EXCLUDES DINUCLEOTIDE IF SCORE IS AT BORDER OF CUTOFF AND MONONUCLEOTIDES MATCH CONSENSUS, DOES NOT CONSIDER A MATCH TO N */
            current_matches_to_consensus = 0;
            query_position = background_position - match_position + 1;
            if (query_position >= 0 && query_position < (*qp).width) current_matches_to_consensus = ((*qp).fraction[current_sequence_value_ULL & 3][query_position] == 0 && total_score[query_position] != 0);
            query_position = background_position2 - match_position + 1;
            if (query_position >= 0 && query_position < (*qp).width) current_matches_to_consensus += ((*qp).fraction[(deleted_sequence_value_ULL & 12) >> 2][query_position] == 0 && total_score[query_position] != 0);
            /* printf("\n%c%c %i %i  \tm%i s%f %f", forward[(deleted_sequence_value_ULL & 12) >> 2], forward[current_sequence_value_ULL & 3], background_position, background_position2, current_matches_to_consensus, score, cut_off);
            printf(" %i, %i, :%i-%i", pwm_position, gap_position, background_position - match_position + 1, query_position); */ 
            if (score > cut_off + current_matches_to_consensus) 
            {
            /* printf(" counted"); */
            ((*dep).incidence)[pwm_position][gap_position][(current_sequence_value_ULL & 3) + (deleted_sequence_value_ULL & 12)]++;
            }
            
        }
    }
    free(total_score);
    return(0);
}


/* SUBROUTINE THAT CONVERTS A PWM TO A LOG FOLD AFFINITY PWM (RELATIVE AFFINITY IS 10 TO THE POWER OF THE SUM OF THE POSITION WEIGHTS) */
short int Log_fold_affinity_pwm (struct normalized_pwm *p)
{
short int line = 0;
short int pwm_position = 0;
double column_maximum;
for (pwm_position = 0; pwm_position < (*p).width; pwm_position++)
{
for(column_maximum = 0, line = 0; line <= 3; line++) if ((*p).fraction[line][pwm_position] > column_maximum) column_maximum = (*p).fraction[line][pwm_position];
for(line = 0; line <= 3; line++) 
{
if ((*p).fraction[line][pwm_position] != 0) (*p).fraction[line][pwm_position] = log10((*p).fraction[line][pwm_position] / column_maximum);
else (*p).fraction[line][pwm_position] = log10(pseudocount);
}
}
for (line = 0; line < 4; line++) {printf("\n"); for (pwm_position = 0; pwm_position < (*p).width; pwm_position++) printf("\t%f", (*p).fraction[line][pwm_position]);}
return(0);
}


/* SUBROUTINE THAT CONVERTS A PWM TO A LOG2 RATIO PWM (DEFINES SPECIFICITY OF EACH NUCLEOTIDE AT EACH POSITION) */
short int Log2_ratio_pwm (struct normalized_pwm *p)
{
    short int line = 0;
    short int pwm_position = 0;
    for (pwm_position = 0; pwm_position < (*p).width; pwm_position++)
    {
        for(line = 0; line <= 3; line++) 
        {
            (*p).fraction[line][pwm_position] = log2(((*p).fraction[line][pwm_position] + pseudocount) / (pseudocount + 1 - (*p).fraction[line][pwm_position]));
        }
    }
    for (line = 0; line < 4; line++) {printf("\n"); for (pwm_position = 0; pwm_position < (*p).width; pwm_position++) printf("\t%f", (*p).fraction[line][pwm_position]);}
    return(0);
}

/* SUBROUTINE THAT CONVERTS A PWM TO A LOGe RATIO PWM (DEFINES SPECIFICITY OF EACH NUCLEOTIDE AT EACH POSITION) */
short int Logeminus10_ratio_pwm (struct normalized_pwm *p)
{
    short int line = 0;
    short int pwm_position = 0;
    for (pwm_position = 0; pwm_position < (*p).width; pwm_position++)
    {
        for(line = 0; line <= 3; line++) 
        {
            (*p).fraction[line][pwm_position] = log(((*p).fraction[line][pwm_position] + pseudocount) / (pseudocount + 1 - (*p).fraction[line][pwm_position])) - 10;
        }
    }
    for (line = 0; line < 4; line++) {printf("\n"); for (pwm_position = 0; pwm_position < (*p).width; pwm_position++) printf("\t%f", (*p).fraction[line][pwm_position]);}
    return(0);
}


/* SUBROUTINE THAT CONVERTS A PWM TO A LOG RATIO PWM (DEFINES SPECIFICITY OF EACH NUCLEOTIDE AT EACH POSITION) */
short int Log_ratio_pwm (struct normalized_pwm *p)
{
short int line = 0;
short int pwm_position = 0;
for (pwm_position = 0; pwm_position < (*p).width; pwm_position++)
{
for(line = 0; line <= 3; line++) 
{
(*p).fraction[line][pwm_position] = log10(((*p).fraction[line][pwm_position] + pseudocount) / (pseudocount + 1 - (*p).fraction[line][pwm_position]));
}
}
for (line = 0; line < 4; line++) {printf("\n"); for (pwm_position = 0; pwm_position < (*p).width; pwm_position++) printf("\t%f", (*p).fraction[line][pwm_position]);}
return(0);
}


/* SUBROUTINE THAT LOADS A PWM AND NORMALIZES IT */
short int Load_pwm (struct normalized_pwm *p, char *filename, short int normalize)
{
long int counter;
char text1;
short int line = 0;
short int pwm_position = 0;
char *current_string;
current_string = malloc(200);
FILE *pwmfile;
if ((pwmfile = fopen(filename, "r")) == '\0') {printf("\nNo File: %s", filename); exit (2);}
for(line = 0; line <= 3 + contacts * 12;)
{
	for(counter = 0; counter < 30; counter++)
	{
		text1 = getc(pwmfile);
		if (text1 == EOF || text1 == '\n' || text1 == '\t')
		{
		current_string[counter] = '\0'; 
		if (counter > 0 && (current_string[0] == '0' || current_string[0] == '1' || current_string[0] == '2' || current_string[0] == '3' || current_string[0] == '4' || current_string[0] == '5' || current_string[0] == '6' || current_string[0] == '7' || current_string[0] == '8' || current_string[0] == '9' || current_string[0] == ' ' || current_string[0] == '-')) 
		{
		(*p).fraction[line][pwm_position] = atof(current_string);
		/* printf("\n%f", (*p).fraction[line][pwm_position]); */
		pwm_position++; 
		}
		if (text1 == '\n' || text1 == EOF) {(*p).width = pwm_position; line++; pwm_position = 0;}
		break;
		}
		current_string[counter]=text1;
		/* printf ("%c", text1); */
	}
}
free (current_string);
if (normalize == 1) Normalize_pwm(p);
// for (line = 0; line < 4 + contacts * 12; line++) {printf("\n"); for (pwm_position = 0; pwm_position < (*p).width; pwm_position++) printf("\t%f", (*p).fraction[line][pwm_position]);}
if (text1 == EOF && line != 3) return(1);
else return (0);
}

/* SUBROUTINE THAT FINDS EXACT MATCH OF QUERY SEQUENCE IN TEST SEQUENCE */
short int Findexactmatch (__uint128_t query_sequence_ULL, short int query_sequence_length, __uint128_t test_sequence_ULL, struct match *match)
{
short int number_of_matches = 0;
short int current_position;
__uint128_t mask2_ULL;

/* printf ("\nExact match for %i using %i as query", test_sequence_ULL, query_sequence_ULL); */
for (current_position = 0, mask2_ULL = pow(4,query_sequence_length)-1; current_position < Nlength - query_sequence_length; query_sequence_ULL <<= 2, mask2_ULL <<= 2, current_position++)
{
if ((test_sequence_ULL & mask2_ULL) == query_sequence_ULL)
{
number_of_matches++;
(*match).position[number_of_matches] = Nlength - current_position - query_sequence_length;
/* printf ("\nMATCH at %i: both values %i", match.position[number_of_matches], query_sequence_ULL); */
}
/* else printf ("\n%llu - No match at %i: %llu is not %llu", test_sequence_ULL, Nlength - current_position - query_sequence_length, test_sequence_ULL & mask2_ULL, query_sequence_ULL); */
}
(*match).position[0] = number_of_matches; 
return (number_of_matches);
}

/* SUBROUTINE THAT FINDS MATCH OF QUERY PWM (score above cut_off) IN TEST SEQUENCE */
short int Findpwmmatch (struct normalized_pwm *p, double cut_off, __uint128_t test_sequence_ULL, struct match *match)
{
short int number_of_matches = 0;
short int nucleotide;
short int current_position;
double score;
signed short int pwm_position;

for (current_position = 0; current_position < Nlength - (*p).width; current_position++)
{
/* GENERATES PWM SCORE */
for(score = 0, pwm_position = (*p).width-1; pwm_position >= 0; pwm_position--) 
{
nucleotide = (test_sequence_ULL & mask_ULL[1][current_position+pwm_position]) >> (2 * (pwm_position+current_position));
/* printf("\nNucleotide %c at sequence_position %i and pwm_position %i", forward[nucleotide], current_position, pwm_position); */
score += (*p).fraction[nucleotide][(*p).width-1-pwm_position];
}

/* printf ("\nScore %f at position %i", score, Nlength-current_position-(*p).width); */
if (score >= cut_off) /* SCORES PWM AT THIS POSITION */
{
number_of_matches++;
(*match).score[number_of_matches] = score;
(*match).position[number_of_matches] = Nlength - current_position - (*p).width;
/* printf ("\nMATCH at %i: score %f higher than cut-off %f", match.position[number_of_matches], score, cut_off); */
}
/* else printf ("\nNo match at %i: %i is not %i", Nlength - current_position - query_sequence_length, query_sequence_ULL, test_sequence_ULL & query_sequence_ULL); */
}
(*match).position[0] = number_of_matches; 
return (number_of_matches);
}

/* REVERSE COMPLEMENTS SEQUENCE VALUE */
__uint128_t Reverse_complement_sequence_value (__uint128_t kmer_sequence_value, short int kmer_length)
{
__uint128_t rc_seq_value = 0;
short int counter;
for (counter = 0; counter < kmer_length; counter++, rc_seq_value <<= 2, kmer_sequence_value >>= 2) rc_seq_value |= (3 ^ (kmer_sequence_value & 3));
return(rc_seq_value >> 2);
}

/* REVERSE COMPLEMENTS SEQUENCE VALUE (long int) */
long int Reverse_complement_sequence_value_li (long int kmer_sequence_value, short int kmer_length)
{
    long int rc_seq_value = 0;
    short int counter;
    for (counter = 0; counter < kmer_length; counter++, rc_seq_value <<= 2, kmer_sequence_value >>= 2) rc_seq_value |= (3 ^ (kmer_sequence_value & 3));
    return(rc_seq_value >> 2);
}

/* PRINTS KMER */
void Kmerprint (__uint128_t print_sequence_value, short int kmer_length)
{
char *forward;
if (rna == 0) forward = dnaforward;
else forward = rnaforward;
    
for (; kmer_length > 0; kmer_length--, print_sequence_value >>= 2) printf("%c",forward[print_sequence_value & 3]);
}

/* SUBROUTINE THAT FINDS BEST SCORE AND POSITION FOR KMER IN A WIDER PWM  */
double Kmerscore (struct normalized_pwm *p, __uint128_t kmer_sequence_value, short int kmer_length, short int *kmermatch_position)
{
short int nucleotide;
short int rev_nucleotide;
short int kmer_position;
double score = 0;
double revscore = 0;
double best_score = -100;
signed short int pwm_position;
__uint128_t reverse_complement_kmer_sequence_value;

reverse_complement_kmer_sequence_value = Reverse_complement_sequence_value(kmer_sequence_value, kmer_length);

/* printf("\nFORWARD_SEQ: "); Kmerprint(kmer_sequence_value, kmer_length); */
    
for (pwm_position = 0; pwm_position <= (*p).width - kmer_length; pwm_position++)
{
/* GENERATES KMER SCORE */
for(score = 0, revscore = 0, kmer_position = kmer_length-1; kmer_position >= 0; kmer_position--) 
{
nucleotide = (kmer_sequence_value & mask_ULL[1][kmer_position]) >> (2 * kmer_position);
rev_nucleotide = (reverse_complement_kmer_sequence_value & mask_ULL[1][kmer_position]) >> (2 * kmer_position);
/* printf("\nNucleotide %c at kmer_position %i and pwm_position %i, score addition %.2f", forward[nucleotide], kmer_position, pwm_position, (*p).fraction[nucleotide][(*p).width-1-pwm_position-kmer_position]); */
score += (*p).fraction[nucleotide][(*p).width-1-pwm_position-kmer_position];
revscore += (*p).fraction[rev_nucleotide][(*p).width-1-pwm_position-kmer_position];
}
if (revscore > score) score = revscore;
/* printf ("\nScore %f at position %i", score, (*p).width-1-pwm_position-kmer_position); */
if (score >= best_score) /* SCORES PWM AT THIS POSITION */
{
best_score = score;
*kmermatch_position = kmer_position;
/* printf ("\nMATCH at %i: score %f higher than cut-off %f", match.position[number_of_matches], score, cut_off); */
}
/* else printf ("\nNo match at %i: %i is not %i", Nlength - current_position - query_sequence_length, query_sequence_ULL, test_sequence_ULL & query_sequence_ULL); */
}

return (best_score);
}


/* faster SUBROUTINE THAT FINDS BEST SCORE AND POSITION FOR KMER IN A WIDER PWM  */
double fastKmerscore (double **pwm, signed short int pwm_width, long long int kmer_sequence_value, signed short int kmer_length)
{
    signed short int nucleotide;
    signed short int rev_nucleotide;
    signed short int kmer_position;
    double score = 0;
    double revscore = 0;
    double best_score = -100;
    short int counter;
    signed short int pwm_position;
    signed short int current_position;
    long long int reverse_complement_kmer_sequence_value = 0;
    long long int current_reverse_complement_kmer_sequence_value;
    long long int current_kmer_sequence_value = kmer_sequence_value;
    
    for (counter = 0; counter < kmer_length; counter++) 
    {
        reverse_complement_kmer_sequence_value <<= 2;
        reverse_complement_kmer_sequence_value |= (3 ^ (current_kmer_sequence_value & 3));
        current_kmer_sequence_value >>= 2;
    }
    /* printf("\nFORWARD_SEQ: "); Kmerprint(kmer_sequence_value, kmer_length);
     printf("\nREVERSE_SEQ: "); Kmerprint(reverse_complement_kmer_sequence_value, kmer_length); */
    
    for (pwm_position = 0; pwm_position < pwm_width - kmer_length; pwm_position++)
    {
        
        current_position = pwm_width-1-pwm_position;
        for(kmer_position = 0, score= 0, revscore = 0, current_kmer_sequence_value = kmer_sequence_value, current_reverse_complement_kmer_sequence_value = reverse_complement_kmer_sequence_value; kmer_position < kmer_length; kmer_position++, current_kmer_sequence_value >>= 2, current_reverse_complement_kmer_sequence_value >>= 2) 
        {
            nucleotide = (current_kmer_sequence_value & 3);
            //           nucleotide = ((nucleotide & 3) << 2) + ((nucleotide & 12) >> 2);
            rev_nucleotide = (current_reverse_complement_kmer_sequence_value & 3);
            //         rev_nucleotide = ((rev_nucleotide & 3) << 2) + ((rev_nucleotide & 12) >> 2);
            // printf("\nDinucleotide %c%c at kmer_position %i and pwm_position %i, score addition %.2f", forward[(nucleotide & 12) >> 2], forward[nucleotide & 3], kmer_position, pwm_position, (*a).fraction[nucleotide][(*a).width-1-pwm_position-kmer_position]);
            score += pwm[nucleotide][current_position-kmer_position];
            revscore += pwm[rev_nucleotide][current_position-kmer_position];
        }
        if (revscore > score) score = revscore;
        /* printf ("\nScore %f at position %i", score, (*p).width-1-pwm_position-kmer_position); */
        if (score >= best_score) /* SCORES PWM AT THIS POSITION */
        {
            best_score = score;
            /* printf ("\nMATCH at %i: score %f higher than cut-off %f", match.position[number_of_matches], score, cut_off); */
        }
        /* else printf ("\nNo match at %i: %i is not %i", Nlength - current_position - query_sequence_length, query_sequence_ULL, test_sequence_ULL & query_sequence_ULL); */
    }
    
    return (best_score);
}

/* SUBROUTINE THAT FINDS BEST SCORE AND POSITION FOR A GAPPED KMER IN A WIDER PWM  */
double fastgappedKmerscore (double **pwm, signed short int pwm_width, long long int kmer_sequence_value, long long int kmer2_sequence_value, signed short int kmer_length, signed short int gap)
{
    signed short int nucleotide;
    signed short int rev_nucleotide;
    signed short int kmer_position;
    double score = 0;
    double revscore = 0;
    double best_score = -100;
    short int counter;
    signed short int pwm_position;
    signed short int current_position;
    signed short int current_gap;
    long long int reverse_complement_kmer_sequence_value = 0;
    long long int current_reverse_complement_kmer_sequence_value;
    long long int current_kmer_sequence_value = kmer_sequence_value;
    long long int reverse_complement_kmer2_sequence_value = 0;
    long long int current_reverse_complement_kmer2_sequence_value;
    long long int current_kmer2_sequence_value = kmer2_sequence_value;
    
    /* GENERATES REVERSE COMPLEMENTS */
    for (counter = 0; counter < kmer_length; counter++)
    {
        reverse_complement_kmer_sequence_value <<= 2;
        reverse_complement_kmer2_sequence_value <<= 2;
        reverse_complement_kmer_sequence_value |= (3 ^ (current_kmer2_sequence_value & 3));
        reverse_complement_kmer2_sequence_value |= (3 ^ (current_kmer_sequence_value & 3));
        current_kmer_sequence_value >>= 2;
        current_kmer2_sequence_value >>= 2;
    }
    
    for (pwm_position = 0; pwm_position < pwm_width - kmer_length * 2 - gap; pwm_position++)
    {
        /* printf("\nFORWARD_SEQ: "); Kmerprint(kmer_sequence_value, kmer_length); */
        
        current_position = pwm_width-1-pwm_position;
        for(kmer_position = 0, score= 0, revscore = 0, current_kmer_sequence_value = kmer_sequence_value, current_reverse_complement_kmer_sequence_value = reverse_complement_kmer_sequence_value, current_kmer2_sequence_value = kmer2_sequence_value, current_reverse_complement_kmer2_sequence_value = reverse_complement_kmer2_sequence_value; kmer_position < kmer_length; kmer_position++, current_kmer_sequence_value >>= 2, current_reverse_complement_kmer_sequence_value >>= 2, current_kmer2_sequence_value >>= 2, current_reverse_complement_kmer2_sequence_value >>= 2) 
        {
            /* KMER1 */
            nucleotide = (current_kmer_sequence_value & 3);
            rev_nucleotide = (current_reverse_complement_kmer_sequence_value & 3);
            score += pwm[nucleotide][current_position-kmer_position];
            revscore += pwm[rev_nucleotide][current_position-kmer_position];
            /* KMER2 */
            nucleotide = (current_kmer2_sequence_value & 3);
            rev_nucleotide = (current_reverse_complement_kmer2_sequence_value & 3);
            score += pwm[nucleotide][current_position-kmer_position-gap-kmer_length];
            revscore += pwm[rev_nucleotide][current_position-kmer_position-gap-kmer_length];
        }
        if (revscore > score) score = revscore;
        /* printf ("\nScore %f at position %i", score, (*p).width-1-pwm_position-kmer_position); */
        if (score >= best_score) /* SCORES PWM AT THIS POSITION */
        {
            best_score = score;
            /* printf ("\nMATCH at %i: score %f higher than cut-off %f", match.position[number_of_matches], score, cut_off); */
        }
        /* else printf ("\nNo match at %i: %i is not %i", Nlength - current_position - query_sequence_length, query_sequence_ULL, test_sequence_ULL & query_sequence_ULL); */
    }
    
    return (best_score);
}

/* SUBROUTINE THAT FINDS BEST SCORE AND POSITION FOR A GAPPED KMER IN A WIDER PWM  */
double vectorfastgappedKmerscore (double **a, signed short int pwm_width, long long int kmer_sequence_value, long long int kmer2_sequence_value, signed short int kmer_length, signed short int gap)
{
    signed short int nucleotide;
    signed short int rev_nucleotide;
    signed short int kmer_position;
    double score = 0;
    double revscore = 0;
    double best_score = -100;
    short int counter;
    double *pwm = a[0];
    double *pwm_position;
    double *pwm_position_after_gap;
    signed short int current_position;
    signed short int current_gap;
    long long int reverse_complement_kmer_sequence_value = 0;
    long long int current_reverse_complement_kmer_sequence_value;
    long long int current_kmer_sequence_value = kmer_sequence_value;
    long long int reverse_complement_kmer2_sequence_value = 0;
    long long int current_reverse_complement_kmer2_sequence_value;
    long long int current_kmer2_sequence_value = kmer2_sequence_value;
    
    /* GENERATES REVERSE COMPLEMENTS */
    for (counter = 0; counter < kmer_length; counter++)
    {
        reverse_complement_kmer_sequence_value <<= 2;
        reverse_complement_kmer2_sequence_value <<= 2;
        reverse_complement_kmer_sequence_value |= (3 ^ (current_kmer2_sequence_value & 3));
        reverse_complement_kmer2_sequence_value |= (3 ^ (current_kmer_sequence_value & 3));
        current_kmer_sequence_value >>= 2;
        current_kmer2_sequence_value >>= 2;
    }
    
    for (counter = 0; counter < pwm_width - kmer_length * 2 - gap; pwm += 4, pwm_position++)
    {
        /* printf("\nFORWARD_SEQ: "); Kmerprint(kmer_sequence_value, kmer_length); */
        
        pwm_position = pwm;
        pwm_position_after_gap = pwm+gap*4;
        for(kmer_position = 0, score= 0, revscore = 0, current_kmer_sequence_value = kmer_sequence_value, current_reverse_complement_kmer_sequence_value = reverse_complement_kmer_sequence_value, current_kmer2_sequence_value = kmer2_sequence_value, current_reverse_complement_kmer2_sequence_value = reverse_complement_kmer2_sequence_value; kmer_position < kmer_length; pwm_position += 4,pwm_position_after_gap += 4, kmer_position++, current_kmer_sequence_value >>= 2, current_reverse_complement_kmer_sequence_value >>= 2, current_kmer2_sequence_value >>= 2, current_reverse_complement_kmer2_sequence_value >>= 2) 
        {
            /* KMER1 */
            nucleotide = (current_kmer_sequence_value & 3);
            rev_nucleotide = (current_reverse_complement_kmer_sequence_value & 3);
            score += pwm_position[nucleotide];
            revscore += pwm_position[rev_nucleotide];
            /* KMER2 */
            nucleotide = (current_kmer2_sequence_value & 3);
            rev_nucleotide = (current_reverse_complement_kmer2_sequence_value & 3);
            score += pwm_position_after_gap[nucleotide];
            revscore += pwm_position_after_gap[rev_nucleotide];
        }
        if (revscore > score) score = revscore;
        /* printf ("\nScore %f at position %i", score, (*p).width-1-pwm_position-kmer_position); */
        if (score >= best_score) /* SCORES PWM AT THIS POSITION */
        {
            best_score = score;
            /* printf ("\nMATCH at %i: score %f higher than cut-off %f", match.position[number_of_matches], score, cut_off); */
        }
        /* else printf ("\nNo match at %i: %i is not %i", Nlength - current_position - query_sequence_length, query_sequence_ULL, test_sequence_ULL & query_sequence_ULL); */
    }
    
    return (best_score);
}



/* SUBROUTINE THAT ADDS ZERO INFORMATION CONTENT (EQUAL) FLANKS TO A PWM */
short int Add_flanks(struct normalized_pwm *p, short int length, double fill_value)
{
short int nucleotide;
short int pwm_position = (*p).width + 2 * length;
for ( ; pwm_position >= (*p).width + length; pwm_position--) for (nucleotide = 0; nucleotide < 4; nucleotide++) (*p).fraction[nucleotide][pwm_position] = fill_value;
for ( ; pwm_position >= length; pwm_position--) for (nucleotide = 0; nucleotide < 4; nucleotide++) (*p).fraction[nucleotide][pwm_position] = (*p).fraction[nucleotide][pwm_position-length];
for ( ; pwm_position >= 0; pwm_position--) for (nucleotide = 0; nucleotide < 4; nucleotide++) (*p).fraction[nucleotide][pwm_position] = fill_value;
(*p).width = (*p).width + 2 * length;
return(0);
}

/* SUBROUTINE THAT PRINTS SEQUENCES */
short int Seqprint (__uint128_t sequence_int, short int length)
{
	signed short int position;
    char *forward;
    if (rna == 0) forward = dnaforward;
    else forward = rnaforward;
    
	if (length > Nlength) return (1);
    for(position = length-1; position > -1 ; position--) printf("%c", forward[(sequence_int & mask_ULL[1][position]) >> (position * 2)]);
	return(0);
}

/* SUBROUTINE THAT SORTS AND INDEXES UNIQUE SEQUENCES */
long int Sort_and_index (struct sequence_incidence_table **sorted_list, struct sequence_incidence_table *unordered_list, short int current_sorted_list, short int number_of_unordered_sequences, struct sequence_incidence_table **sorted_index)
{
	/* printf("%i\t%i\ti=%i\t%i\n", current_sorted_list, number_of_unordered_sequences, sorted_list[current_sorted_list][0].sequence_value_ULL, sorted_list[current_sorted_list][1].sequence_value_ULL); */
	short int shortcounter;
	short int new_sorted_list = 1 - current_sorted_list;
	short int old_index_value = 5000;
	short int new_index_value = 5001;
	long int counter;
	long int number_of_sequences_to_merge;
	struct sequence_incidence_table *cursor;
	struct sequence_incidence_table *insertion_point;
	struct sequence_incidence_table *low_incidence;
	cursor = sorted_list[current_sorted_list];
	insertion_point = sorted_list[new_sorted_list];
	
	/* printf("\nFirst values are: current list %i and future list %i", cursor[0].sequence_value_ULL, insertion_point[0].sequence_value_ULL); */

	/* printf("\nUNORDERED LIST:");
	for (shortcounter = 0; shortcounter < 100; shortcounter++) if (unordered_list[shortcounter].incidence > 1) {printf("\nSequence %i occurs %i times", unordered_list[shortcounter].sequence_value_ULL, unordered_list[shortcounter].incidence);}  */

	/* SORTS UNORDERED LIST */
	qsort (unordered_list, number_of_unordered_sequences, sizeof(struct sequence_incidence_table), Sort_according_to_sequence_value_ULL);
	
	/* printf("\nSORTED UNORDERED LIST:");
	for (shortcounter = 0; shortcounter < 100; shortcounter++) printf("\t%i", unordered_list[shortcounter].incidence);  */

		
	/* FINDS BORDER BETWEEN LOW INCIDENCE AND HIGH INCIDENCE SEQUENCES */
	for(low_incidence = unordered_list; low_incidence < unordered_list + 50 && (*low_incidence).incidence > 100; low_incidence++);
	
	number_of_sequences_to_merge = (sorted_index[4097] - cursor) + (unordered_list - low_incidence) + number_of_unordered_sequences;
	
	/* printf("\nMerges %i sequences", number_of_sequences_to_merge); */
	
	/* MERGES THE TWO SORTED LISTS */
	for (counter = 0, shortcounter = 0; counter + shortcounter < number_of_sequences_to_merge; )
	{
		if (shortcounter < number_of_unordered_sequences && (cursor[counter].sequence_value_ULL > low_incidence[shortcounter].sequence_value_ULL | cursor + counter >= sorted_index[4097])) 
		{
		insertion_point[counter+shortcounter].sequence_value_ULL = low_incidence[shortcounter].sequence_value_ULL;
		insertion_point[counter+shortcounter].incidence = low_incidence[shortcounter].incidence;
		shortcounter++;
		continue;
		}
		insertion_point[counter+shortcounter].sequence_value_ULL = cursor[counter].sequence_value_ULL;
		insertion_point[counter+shortcounter].incidence = cursor[counter].incidence;
		counter++;
	}
	
	/* printf("\nSORTED LIST:");
	for (counter = 0; counter < number_of_sequences_to_merge; counter++) printf("\t%i", insertion_point[counter].sequence_value_ULL); */
	
	/* INDEXES LIST */
	for (counter = 0; counter < number_of_sequences_to_merge; counter++)
	{
		new_index_value = sorted_list[new_sorted_list][counter].sequence_value_ULL / 65536;
		if (old_index_value == new_index_value) continue;
		old_index_value = sorted_list[new_sorted_list][counter].sequence_value_ULL / 65536;
		sorted_index[old_index_value] = sorted_list[new_sorted_list]+counter;
		/* printf("\nINDEX %i VALUE: %i \t sequence %i, %i", old_index_value, counter, (*sorted_index[old_index_value]).sequence_value_ULL, (*sorted_index[old_index_value]).sequence_value_ULL / 65536); */
	}
	sorted_index[4097] = sorted_list[new_sorted_list] + number_of_sequences_to_merge;
	return (low_incidence - unordered_list);
}


/* SUBROUTINE THAT CALLS WINFLAT */
double Winflat (long int sample, long int background, long int total_sample, long int total_background)
{
	
	if (sample < 1000 && background < 1000)
	{
		if (pvalue_cache[sample][background] != 2) return (pvalue_cache[sample][background]);
	}
	
	FILE *winflat;
	short int counter;
	short int line;
	char text1;
	char *addstring;
	addstring = malloc (255);
	char *system_command;
	system_command = malloc (2000);
	strcpy (system_command, "winflat -xvalue ");
	sprintf(addstring, "%li", sample);
	strcat (system_command, addstring);
	strcat (system_command, " -yvalue ");
	sprintf(addstring, "%li", background);
	strcat (system_command, addstring);
	strcat (system_command, " -diff ");
	sprintf(addstring, "%li", total_sample);
	strcat (system_command, addstring);
	strcat (system_command, " ");
	sprintf(addstring, "%li", total_background);
	strcat (system_command, addstring);
	strcat (system_command, " | sed s/.*=/\\ /");
	/* printf("%s", system_command); */
	winflat = popen(system_command, "r");
	for(line = 1; line <= 2; line++) /* captures output from line 2 */
	{
		for(counter = 0; counter < 100; counter++)
		{
			text1 = getc(winflat);
			if (text1 == EOF) return(2);
			if (text1 == '\n') break;
			addstring[counter]=text1;
		}
	}
	pclose(winflat);
	addstring[counter] = '\0';
	if (sample < 1000 && background < 1000) pvalue_cache[sample][background] = atof(addstring);
	return (atof(addstring));
}

/* TESTS IF THE CURRENT KMER IS A PALINDROME OR DIRECT REPEAT, RETURNS NON-ZERO VALUE IF IT IS, 0 IF IT IS NOT                                               */
/* ODD NUMBERED PALINDROMES OR REPEATS WITH ANY BASE IN THE MIDDLE RETURN 3 OR 1, AND PERFECT PALINDROMES OR REPEATS RETURN 4 OR 2, RESPECTIVELY             */
short int Is_this_sequence_dimer (long int current_kmer, short int current_kmer_length)
{
short int counter = current_kmer_length * 2 - 2;
short int pos = 0;
short int mask = 3;
short int shift = current_kmer_length + current_kmer_length % 2;  /* ADDS REMAINDER TO ALLOW ODD PALINDROMES AND REPEATS */
/* TESTS IF SEQUENCE IS A PALINDROME */
for (; (counter >= current_kmer_length) && (((current_kmer >> pos) & mask) + ((current_kmer >> counter) & mask) == 3); counter--, counter--, pos++, pos++);
if (counter < current_kmer_length) 
{
if (shift == current_kmer_length) return (4);                     /* PERFECT PALINDROME */
else return (3);                                                  /* PALINDROME WITH CENTER NUCLEOTIDE */
}
/* TESTS IF SEQUENCE IS A DIRECT REPEAT */
mask = pow(2, (current_kmer_length - current_kmer_length % 2)) - 1;
if (((current_kmer ^ (current_kmer >> shift)) & mask) == 0) 
{
if (shift == current_kmer_length) return (2);                     /* PERFECT REPEAT */
else return (1);                                                  /* REPEAT WITH CENTER NUCLEOTIDE */
}
return (0);
}
	
Exit_with_error (char *error_code, double *error_values)
{
short int counter;
printf("\n\n**** ERROR %s", error_code);
for (counter = 1; counter < error_values[0]; counter++) printf("\t%f", error_values[counter]);
exit(1);
}

/* CALCULATES INFORMATION CONTENT FOR INTEGER DISTRIBUTION */
double Information_content (long int *counts, long int start, long int end)
{
long int length;
length = end - start + 1;
long int counter;
long int total = 0;
double fraction;
double ic = log2(length);
    
for (counter = start; counter < length; counter++) total += counts[counter];
for (counter = start; counter < length; counter++)
{        
fraction = (double) counts[counter] / (double) total;
if (counts[counter] != 0) ic += fraction * log2 (fraction);
/* printf("\n%i\t%i\t%f", counter, counts[counter], ic); */
}
return(ic);
}

/* SUBROUTINE THAT DETERMINES MAXIMUM LENGTH OF ALIGNMENT TO SELF FOR EACH POSITION IN A QUERY SEQUENCE */
short int *Alignscore (__uint128_t query_sequence_ULL, struct alignscore *a, short int *max_alignment_length)
{
    short int current_position2;
    short int current_position;
    short int gap;
    __uint128_t mask2_ULL;
    __uint128_t test_sequence_ULL = query_sequence_ULL;
    __uint128_t original_sequence_ULL = query_sequence_ULL;
    short int query_sequence_length;
    short int strand;
    short int counter; 
    for (strand = 0; strand < 2; strand++) for(query_sequence_length = Nlength * 2 - 1; query_sequence_length >= 0; query_sequence_length--) 
    {
    (*a).score[query_sequence_length][strand] = 0;
    (*a).count[query_sequence_length][strand] = 0;
    }
    for(query_sequence_length = Nlength-1; query_sequence_length >= 0; query_sequence_length--) for (gap = 0; gap < Nlength; gap++)
    {
    (*a).direct_repeat[query_sequence_length][gap] = 0;
    (*a).inverted_repeat[query_sequence_length][gap] = 0;
    }
    
    for (strand = 0; strand < 2; strand++)
    {
        for(query_sequence_length = Nlength - 1; query_sequence_length > 0; query_sequence_length--)
        {
            for (current_position2 = 0; current_position2 < Nlength - query_sequence_length; current_position2++)
            {
                query_sequence_ULL = (test_sequence_ULL & mask_ULL[query_sequence_length][current_position2]) >> (2 * current_position2);
                for (current_position = 0, mask2_ULL = pow(4,query_sequence_length)-1; current_position < Nlength - query_sequence_length; 
                query_sequence_ULL <<= 2, mask2_ULL <<= 2, current_position++)
                {
                    /* Seqprint(test_sequence_ULL, Nlength-1);
                    printf("\tcompared against\t");
                    Seqprint(query_sequence_ULL, Nlength-1);
                    printf("\tlength,pos2,pos,strand=%i,%i,%i,%i\t",query_sequence_length, current_position2, current_position, strand);
                    Seqprint(mask_ULL[query_sequence_length][current_position2], Nlength-1);
                    printf("\n"); */
                    if(current_position2 == current_position && strand == 0)
                    {
                        /* printf("Original position excluded\n"); */
                        continue;
                    }
                    if ((original_sequence_ULL & mask2_ULL) == query_sequence_ULL)
                    {
                        if (strand == 0) (*a).direct_repeat[query_sequence_length][abs(current_position-current_position2)]++;
                        else (*a).inverted_repeat[query_sequence_length][abs(Nlength-1-current_position-current_position2)]++;
                        
                        for (counter = current_position; counter < current_position + query_sequence_length; counter++)
                        {
                            if (query_sequence_length > (*a).score[counter][strand]) 
                            {
                                /* printf("\npos,strand %i,%i highest match = %i", counter, strand, query_sequence_length); */
                                (*a).score[counter][strand] = query_sequence_length;
                                (*a).count[counter][strand]++;
                            }
                        }
                        if(query_sequence_length > max_alignment_length[strand]) max_alignment_length[strand] = query_sequence_length;
                    }
                }
            }
        }
        test_sequence_ULL = Reverse_complement_sequence_value(test_sequence_ULL, Nlength-1);
    }
    return (max_alignment_length);
}

/* SUBROUTINE THAT ADDS ALIGNMENT SCORES RELATIVE TO MATCH POSITIONS */
short int *Add_to_alignscore (__uint128_t query_sequence_ULL, struct alignscore *t, struct alignscore *a, short int match_position, short int *max_align_score)
{
    Alignscore (query_sequence_ULL, a, max_align_score);
    short int position;
    short int strand;
    short int gap;
    for(strand = 0; strand < 2; strand++)
    {
        for(position = 0; position < Nlength; position++)
        {
            (*t).score[position+Nlength-1-match_position][strand] += (*a).score[position][strand];
            (*t).count[position+Nlength-1-match_position][strand] += (*a).count[position][strand];
            for(gap = 0; gap < Nlength; gap++)
            {
                (*t).direct_repeat[position][gap] += (*a).direct_repeat[position][gap];
                (*t).inverted_repeat[position][gap] += (*a).inverted_repeat[position][gap];
            }
        }
    }
    return (max_align_score);
}


/* CALCULATES HAMMING DISTANCE */
short int bitHamming (__uint128_t sequence1_ULL, __uint128_t sequence2_ULL, short int sequence2_length)
{
    __uint128_t difference_ULL;
    short int hamming_distance = 0;
    short int reverse_complement_hamming_distance = 0;
    difference_ULL = sequence1_ULL ^ sequence2_ULL;
    for (hamming_distance = 0; difference_ULL; hamming_distance += ((difference_ULL & 3) != 0), difference_ULL >>= 2);
    difference_ULL = sequence1_ULL ^ Reverse_complement_sequence_value(sequence2_ULL, sequence2_length);
    for (reverse_complement_hamming_distance = 0; difference_ULL; reverse_complement_hamming_distance += ((difference_ULL & 3) != 0), difference_ULL >>= 2);
    if (hamming_distance < reverse_complement_hamming_distance) return (hamming_distance);
    else return (reverse_complement_hamming_distance);
}

/* CALCULATES EDIT DISTANCE WITHOUT INDELS USING FORWARD AND REVCOMP */
short int bitEditdist (__uint128_t sequence1_ULL, __uint128_t sequence2_ULL, short int sequence2_length)
{
    __uint128_t difference_ULL;
    __uint128_t current_sequence1_ULL = sequence1_ULL;
    __uint128_t current_sequence2_ULL = sequence2_ULL;
    short int minimum_edit_distance = sequence2_length+2;
    short int current_edit_distance = 0;
    short int current_reverse_complement_edit_distance = 0;
    short int counter;
    
    /* DELETIONS OF SEQUENCE 1, FORWARD */
    for(counter = 0; counter < sequence2_length; counter++, current_sequence1_ULL >>= 2)
    {
    difference_ULL = (current_sequence1_ULL ^ current_sequence2_ULL) & lowmask_ULL[sequence2_length-counter-1];
    for (current_edit_distance = counter; difference_ULL; current_edit_distance += ((difference_ULL & 3) != 0), difference_ULL >>= 2);
    if (current_edit_distance < minimum_edit_distance) minimum_edit_distance = current_edit_distance;    
    }

    /* DELETIONS OF SEQUENCE 2, FORWARD */
    current_sequence1_ULL = sequence1_ULL;
    for(counter = 0; counter < sequence2_length; counter++, current_sequence2_ULL >>= 2)
    {
    difference_ULL = (current_sequence1_ULL ^ current_sequence2_ULL) & lowmask_ULL[sequence2_length-counter-1];
    for (current_edit_distance = counter; difference_ULL; current_edit_distance += ((difference_ULL & 3) != 0), difference_ULL >>= 2);
    if (current_edit_distance < minimum_edit_distance) minimum_edit_distance = current_edit_distance;    
    }

    /* DELETIONS OF SEQUENCE 2, REVERSE */
    current_sequence2_ULL = sequence2_ULL;
    current_sequence1_ULL = Reverse_complement_sequence_value(current_sequence1_ULL, sequence2_length);
    for(counter = 0; counter < sequence2_length; counter++, current_sequence2_ULL >>= 2)
    {
    difference_ULL = (current_sequence1_ULL ^ current_sequence2_ULL) & lowmask_ULL[sequence2_length-counter-1];
    for (current_edit_distance = counter; difference_ULL; current_edit_distance += ((difference_ULL & 3) != 0), difference_ULL >>= 2);
    if (current_edit_distance < minimum_edit_distance) minimum_edit_distance = current_edit_distance;    
    }

    /* DELETIONS OF SEQUENCE 1, REVERSE */
    current_sequence2_ULL = sequence2_ULL;
    for(counter = 0; counter < sequence2_length; counter++, current_sequence1_ULL >>= 2)
    {
    difference_ULL = (current_sequence1_ULL ^ current_sequence2_ULL) & lowmask_ULL[sequence2_length-counter-1];
    for (current_edit_distance = counter; difference_ULL; current_edit_distance += ((difference_ULL & 3) != 0), difference_ULL >>= 2);
    if (current_edit_distance < minimum_edit_distance) minimum_edit_distance = current_edit_distance;    
    }

    return(minimum_edit_distance);
}

/* CALCULATES EDIT DISTANCE WITHOUT INDELS, RETURNS MINIMUM EDIT DIST, DOES NOT CONSIDER REVERSE COMPLEMENTS */
short int bitEditds (__uint128_t sequence1_ULL, __uint128_t sequence2_ULL, short int gap1, short int gap2, short int sequence2_length)
{
    __uint128_t difference_l_ULL;
    __uint128_t difference_r_ULL;    
    short int halfsite_length = sequence2_length / 2;    

    short int minimum_edit_distance = sequence2_length+2;
    short int current_edit_distance = 0;
    short int current_reverse_complement_edit_distance = 0;
    short int counter;
    signed short int gap_difference = gap2 - gap1;
    
    __uint128_t current_sequence1l_ULL = (sequence1_ULL >> (halfsite_length * 2)) & lowmask_ULL[halfsite_length];
    __uint128_t current_sequence1r_ULL = sequence1_ULL & lowmask_ULL[halfsite_length];
    __uint128_t current_sequence2l_ULL = (sequence2_ULL >> (halfsite_length * 2)) & lowmask_ULL[halfsite_length];
    __uint128_t current_sequence2r_ULL = sequence2_ULL & lowmask_ULL[halfsite_length];
    
    
    /* DELETIONS OF SEQUENCES 1, FORWARD */
    /* GAP DIFFERENCE SHIFT */
    if(gap_difference)
    {
    current_sequence1r_ULL >>= 2;
    }
    /* LOOP FOR TOTAL SHIFT */
    for(counter = 0; counter < 2 - gap_difference; counter++, current_sequence1l_ULL >>= 2, current_sequence1r_ULL >>= 2)
    {
        difference_l_ULL = (current_sequence1l_ULL ^ current_sequence2l_ULL) & lowmask_ULL[sequence2_length-counter-1];
        difference_r_ULL = (current_sequence1r_ULL ^ current_sequence2r_ULL) & lowmask_ULL[sequence2_length-counter-1-gap_difference];
        for (current_edit_distance = counter; difference_l_ULL | difference_r_ULL; current_edit_distance += ((difference_l_ULL & 3) != 0) + ((difference_r_ULL & 3) != 0), difference_l_ULL >>= 2, difference_r_ULL >>= 2);
        if (current_edit_distance < minimum_edit_distance) minimum_edit_distance = current_edit_distance;
    }
    
    current_sequence1l_ULL = (sequence1_ULL >> halfsite_length) & lowmask_ULL[halfsite_length];
    current_sequence1r_ULL = sequence1_ULL & lowmask_ULL[halfsite_length];
    /* DELETIONS OF SEQUENCES 2, FORWARD */
    /* GAP DIFFERENCE SHIFT */
    if(gap_difference)
    {
    current_sequence1r_ULL >>= 2;
    }
    /* LOOP FOR TOTAL SHIFT */
    for(counter = 0; counter < 2 - gap_difference; counter++, current_sequence2l_ULL >>= 2, current_sequence2r_ULL >>= 2)
    {
        difference_l_ULL = (current_sequence1l_ULL ^ current_sequence2l_ULL) & lowmask_ULL[sequence2_length-counter-1];
        difference_r_ULL = (current_sequence1r_ULL ^ current_sequence2r_ULL) & lowmask_ULL[sequence2_length-counter-1-gap_difference];
        for (current_edit_distance = counter; difference_l_ULL | difference_r_ULL; current_edit_distance += ((difference_l_ULL & 3) != 0) + ((difference_r_ULL & 3) != 0), difference_l_ULL >>= 2, difference_r_ULL >>= 2);
        if (current_edit_distance < minimum_edit_distance) minimum_edit_distance = current_edit_distance; 
    }
    
    return(current_edit_distance);
}


/* DETERMINES IF GAPPED KMERS ARE WITHIN 1 EDIT DISTANCE (only for even gapped kmers) DOES NOT WORK */
short int GappedkmerEditdist (__uint128_t sequence1_ULL, __uint128_t sequence2_ULL, short int gap1, short int gap2, short int sequence2_length)
{
    short int forward_dist;
    short int reverse_dist;
    
    /* DISTANCES BETWEEN HALF-SITES (forward) */
    forward_dist = bitEditds (sequence1_ULL, sequence2_ULL, gap1, gap2, sequence2_length);
    /* DISTANCES BETWEEN HALF-SITES (reverse) */
    sequence1_ULL = Reverse_complement_sequence_value(sequence1_ULL, sequence2_length);
    reverse_dist = bitEditds (sequence1_ULL, sequence2_ULL, gap1, gap2, sequence2_length);
    
    if (forward_dist < reverse_dist) return (forward_dist);
    else return(reverse_dist);
}

/* CALCULATES HAMMING DISTANCE USING BITCOMB (does not work) */
short int combHamming8 (long int sequence1, long int sequence2)
{
short int hamming_distance;
unsigned int v; // count the number of bits set in v
unsigned int vr; // count the number of bits set in vr
    
/* EORS and COMBS to get bits */
v = sequence1 ^ sequence2;
v = (v & 21845) | ((v & 43690) >> 1);
vr = sequence1 ^ Reverse_complement_sequence_value(sequence2, 8);
vr = (vr & 21845) | ((vr & 43690) >> 1);

/* COUNTS BITS */

unsigned int c; // c accumulates the total bits set in v
unsigned int cr; // cr accumulates the total bits set in vr
for (c = 0; v; c++) v &= v - 1; // clear the least significant bit set
for (cr = 0; vr; cr++) vr &= vr - 1; // clear the least significant bit set
if(cr > c) return(c); else return(cr);
}

/* CALCULATES SUM OF LOW PAIR COUNTS FOR EACH EDIT DISTANCE */
short int Lowpaircounts (struct sumtable *s, long int *signal_counts, long int cutoff, short int sequence_length)
{
signed long long int current_sequence1_ULL;
__uint128_t current_sequence2_ULL; 
__uint128_t last_kmer_value_ULL;
__uint128_t current_lowcount_kmer_ULL;
__uint128_t current_highcount_kmer_ULL;
long int *max_min_counts = (*s).max_min_counts;     
long int *max_max_counts = (*s).max_max_counts; 
long int ***cloud_counts = (*s).cloud_counts;
__uint128_t *max_min_kmer = (*s).max_min_kmer;  
__uint128_t *max_max_kmer = (*s).max_max_kmer;  
__uint128_t *sums = (*s).sums;
long int *counts = (*s).counts;
long int *current_count;
long int *current_sum;    
long int *local_max; 
local_max = malloc(sizeof(long int)*sequence_length+5);
current_count = malloc(sizeof(long int)*sequence_length+5);
current_sum = malloc(sizeof(long int)*sequence_length+5);
    
char *forward;
char *nucleotide_bitiupac;
if (rna == 0) {forward = dnaforward; nucleotide_bitiupac = dna_bitiupac;}
else {forward = rnaforward; nucleotide_bitiupac = rna_bitiupac;}
    
    
long int max_counts;    
short int counter;
long int counter2;
signed short int position;
short int current_edit_distance;
long int current_counts;
last_kmer_value_ULL = pow(4,sequence_length);
long int current_kmer_number = 0; 
long int max_kmer_number = 0; 
for(max_counts = 0, current_sequence1_ULL = 0; current_sequence1_ULL < last_kmer_value_ULL; current_sequence1_ULL++) if (signal_counts[current_sequence1_ULL] > max_counts) max_counts = signal_counts[current_sequence1_ULL];
cutoff = (long int) max_counts * local_max_min_percent;
/* #pragma omp parallel for schedule(static) default(shared) private(current_counts, current_edit_distance, max_min_counts, max_max_counts, max_min_kmer, max_max_kmer, current_sequence2_ULL) */
for(current_sequence1_ULL = 0; current_sequence1_ULL < last_kmer_value_ULL; current_sequence1_ULL++) 
{
if (signal_counts[current_sequence1_ULL] > cutoff) 
{
for(counter = 0 ; counter < sequence_length; counter++)
{
current_count[counter] = 0;
current_sum[counter] = 0;
local_max[counter] = 1;
}
for(current_sequence2_ULL = 0; current_sequence2_ULL < last_kmer_value_ULL; current_sequence2_ULL++)
{
    if(signal_counts[current_sequence1_ULL] > signal_counts[current_sequence2_ULL]) {current_lowcount_kmer_ULL = current_sequence2_ULL;current_highcount_kmer_ULL = current_sequence1_ULL;}
    else 
    {
        current_lowcount_kmer_ULL = current_sequence1_ULL; 
        current_highcount_kmer_ULL = current_sequence2_ULL;
    }
    current_counts = signal_counts[current_lowcount_kmer_ULL]; 
    if(current_counts > cutoff) 
    {
    current_edit_distance = bitEditdist (current_sequence1_ULL, current_sequence2_ULL, sequence_length);
    /* printf("%i-",current_hamming_distance); */
    if (current_counts > max_min_counts[current_edit_distance]) {max_min_counts[current_edit_distance] = current_counts; max_max_counts[current_edit_distance] = signal_counts[current_highcount_kmer_ULL]; max_min_kmer[current_edit_distance] = current_lowcount_kmer_ULL;max_max_kmer[current_edit_distance] = current_highcount_kmer_ULL;}
    current_sum[current_edit_distance] += current_counts;
    current_count[current_edit_distance]++;
    if(signal_counts[current_sequence1_ULL] < signal_counts[current_sequence2_ULL]) local_max[current_edit_distance] = 0;
    }
}
for(counter = 0 ; counter < sequence_length; counter++)
{
cloud_counts[counter][current_kmer_number][1] = (long int) current_sequence1_ULL;
if (counter > 0) 
{
cloud_counts[counter][current_kmer_number][2] = cloud_counts[counter-1][current_kmer_number][2] + current_count[counter];
cloud_counts[counter][current_kmer_number][3] = cloud_counts[counter-1][current_kmer_number][3] + current_sum[counter];
cloud_counts[counter][current_kmer_number][4] = local_max[counter] * cloud_counts[counter-1][current_kmer_number][4];
}
else 
{
cloud_counts[counter][current_kmer_number][2] = current_count[counter];
cloud_counts[counter][current_kmer_number][3] = current_sum[counter];
cloud_counts[counter][current_kmer_number][4] = local_max[counter];
}
counts[counter] += current_count[counter];
sums[counter] += current_sum[counter];
}
current_kmer_number++;
}
}

/* FINDS CLOUD MAX COUNTS */
for(current_highcount_kmer_ULL = 0, counter = 0 ; counter < sequence_length; counter++)
{
    for(counter2 = 0; counter2 < current_kmer_number; counter2++) if (cloud_counts[counter][counter2][3] > current_highcount_kmer_ULL) current_highcount_kmer_ULL = cloud_counts[counter][counter2][3];
    for(counter2 = 0; counter2 < current_kmer_number; counter2++) if (cloud_counts[counter][counter2][3] == current_highcount_kmer_ULL) cloud_counts[counter][counter2][5] = 1;
} 
    
/* PRINTS PAIRS */
    printf("\nEdit_distance\tNumber\tAverage_lowcount\tMax_lowcount\tMax_highcount\tMax_low_kmer\tMax_high_kmer");
    for(counter = 0; counter <= sequence_length; counter++) 
    {
    printf("\n%i\t%li\t%.1f\t%li\t%li\t", counter, (*s).counts[counter], ((long int) (*s).sums[counter])/((*s).counts[counter]+0.1), (*s).max_min_counts[counter], (*s).max_max_counts[counter]);
    if((*s).counts[counter] == 0) printf("-"); 
    else for(position = sequence_length-1; position > -1 ; position--) printf("%c", forward[((*s).max_min_kmer[counter] & mask_ULL[1][position]) >> (position * 2)]);
    printf("\t");
    if((*s).counts[counter] == 0) printf("-"); 
    else for(position = sequence_length-1; position > -1 ; position--) printf("%c", forward[((*s).max_max_kmer[counter] & mask_ULL[1][position]) >> (position * 2)]);
    printf("\t\teditdist");
    }
    printf("\n");

    printf("\n");
    
    char *max_string[3];
    max_string[0] = "";
    max_string[1] = "local_max";
    max_string[2] = "cloud_max";
    short int max_nucleotide;
    short int nucleotide;
    short int iupac_bits;

    
/* PRINTS CLOUD COUNTS (annotates also local and cloud max) */
 printf("\nEdit_distance\tkmer\tIUPAC\tSelf count\tNeighbors\tTotal count\tlocal_max\tcloud_max\tcc");
 for(counter = 0 ; counter < sequence_length; counter++)
 {
 for(counter2 = 0 ; counter2 < current_kmer_number; counter2++)
 {
/* PRINTS DISTANCE AND KMER */
     printf("\n%li\t", cloud_counts[counter][counter2][0]);
    for(position = sequence_length-1; position > -1 ; position--) printf("%c", forward[(cloud_counts[counter][counter2][1] & mask_ULL[1][position]) >> (position * 2)]);
    printf("\t");
/* PRINTS IUPAC */
     for(position = sequence_length-1; position > -1 ; position--) 
    {
        max_nucleotide = (cloud_counts[counter][counter2][1] & mask_ULL[1][position]) >> (position * 2);
        max_counts = signal_counts[cloud_counts[counter][counter2][1]];
        for(iupac_bits = 1 << max_nucleotide, nucleotide = 0; nucleotide < 4;nucleotide++)
        {
        if(nucleotide == max_nucleotide) continue;
        if(
            (signal_counts[
                           (cloud_counts[counter][counter2][1] & (~(mask_ULL[1][position]))) | (nucleotide << (position * 2))
                                                                                          ]
            ) * 2 > max_counts) iupac_bits |= (1 << nucleotide);    
        }
        printf("%c", nucleotide_bitiupac[iupac_bits]);
    }
/* PRINTS COUNT STATS */
     printf("\t%li\t%li\t%li\t\t%s\t%s\tcc%li",  cloud_counts[0][counter2][3], cloud_counts[counter][counter2][2],cloud_counts[counter][counter2][3], max_string[cloud_counts[counter][counter2][4] * (counter != 0)], max_string[cloud_counts[counter][counter2][5] * 2], cloud_counts[counter][counter2][0]);
 }
 }
    free(current_count);
    free(current_sum);
    free(local_max);
}

/* REMOVES PALINDROMIC MATCHES */
short int Remove_palindromic_matches(struct match **match, short int query_sequence_length)
{
    /* RETURNS IF GLOBAL FLAG SET TO COUNT ALL HITS */
    if (count_both_instances_of_palindromic_hit == 1) return(1);
    short int fw_match = 1;
    short int rev_match = 1;
    short int counter;
    float fw_score;
    float rev_score;
    /* printf("\nCALLS remove palindromic matches");
    printf("\nFW pos-scores: ");
    for (fw_match = 1; fw_match <= (*match)[0].position[0]; fw_match++) printf("%i-%f ", (*match)[0].position[fw_match],(*match)[0].score[fw_match]);
    printf("\nREV pos-scores: ");
    for (fw_match = 1; fw_match <= (*match)[1].position[0]; fw_match++) printf("%i-%f ", (*match)[1].position[fw_match],(*match)[1].score[fw_match]);
    */
    
    for (fw_match = 1; fw_match <= (*match)[0].position[0]; fw_match++) for (rev_match = 1; rev_match <= (*match)[1].position[0]; rev_match++)
    {
    /* CHECKS IF HITS ARE COMPLETELY OVERLAPPING ON OPPOSITE STRANDS */
    if((*match)[0].position[fw_match] == Nlength + 1 - query_sequence_length - (*match)[1].position[rev_match])
    {
        /* DOES NOT REMOVE MATCHES (CONTINUES LOOP) IF COUNT ONLY UNEQUAL FLAG IS SET AND SCORES ARE EQUAL */
        if (count_unequal_hits_only == 1 && (*match)[0].score[fw_match] == (*match)[1].score[rev_match]) continue;
        /* printf("\nPALINDROMIC MATCH AT FW %i REV %i FWPOS %i REVPOS %i QUERYLENGTH %i fwSCORE %f revSCORE %f", fw_match, rev_match, (*match)[0].position[fw_match], (*match)[1].position[rev_match], query_sequence_length, (*match)[0].score[fw_match], (*match)[1].score[rev_match]); */
        last_count_was_forward ^= 1;
        if (prefer_forward_strand == 1) last_count_was_forward = 0;
        if (prefer_reverse_strand == 1) last_count_was_forward = 1;
        if((((*match)[0].score[fw_match] < (*match)[1].score[rev_match] || ((*match)[0].score[fw_match] == (*match)[1].score[rev_match]) && last_count_was_forward == 1) && count_only_forward_instance_of_palindromic_hit == 0) || count_only_reverse_instance_of_palindromic_hit == 1)
        {
            /* printf("--WINNER %f", (*match)[1].score[rev_match]); */
            /* REMOVES HIT FROM FORWARD STRAND */
            for(counter = fw_match; counter < (*match)[0].position[0]; counter++) 
            {
                (*match)[0].position[counter] = (*match)[0].position[counter+1];
                (*match)[0].score[counter] = (*match)[0].score[counter+1];
            }
            (*match)[0].position[0]--;
            fw_match--;
            break;
        }
        else
        {
             /* REMOVES HIT FROM REVERSE STRAND */
            for(counter = rev_match; counter < (*match)[1].position[0]; counter++) 
            {
                (*match)[1].position[counter] = (*match)[1].position[counter+1];
                (*match)[1].score[counter] = (*match)[1].score[counter+1];
            }
            (*match)[1].position[0]--;
            rev_match--;
            break;
        }
    }
        
    }
    /* printf("\nFW pos-scores: ");
    for (fw_match = 1; fw_match <= (*match)[0].position[0]; fw_match++) printf("%i-%f ", (*match)[0].position[fw_match],(*match)[0].score[fw_match]);
    printf("\nREV pos-scores: ");
    for (fw_match = 1; fw_match <= (*match)[1].position[0]; fw_match++) printf("%i-%f ", (*match)[1].position[fw_match],(*match)[1].score[fw_match]);
    */
    
    return(0);
    
}

/* CALCULATES CORRELATION BETWEEN KMER COUNT AND SAME COUNT DERIVED FROM MONONUCLEOTIDE DISTRIBUTION */
double Kmer_mono_correlation(long int *observed_counts, double *expected_counts, short int kmer_length)
{
    double correlation;
    double score;
    double sizefactor;
    long int current_kmer;
    long int kmer_position;
    long int current_kmer_number;
    long int current_kmer_count;
    long int number_of_different_kmers = pow(4, kmer_length);
    long int total_kmer_count = 0;
    long int total_expected_kmer_count = 0;
    struct count_pwm p;
    count_pwm_init(&p, "c", kmer_length+1, 0);
    struct normalized_pwm n;
    normalized_pwm_init(&n, "n", kmer_length+1, 0);

    double sum_x_squared = 0;
    double sum_y_squared = 0;
    double sum_xy = 0;
    
    /* CALCULATES MONONUCLEOTIDE DISTRIBUTION */
    for(current_kmer_number = 0; current_kmer_number < number_of_different_kmers; current_kmer_number++)
    {
    total_kmer_count += observed_counts[current_kmer_number];
    for(current_kmer = current_kmer_number, kmer_position = 0; kmer_position < kmer_length; kmer_position++, current_kmer >>= 2)
    {
    p.incidence[current_kmer & 3][kmer_position] += observed_counts[current_kmer_number];
    }
    }
    Count_to_normalized_pwm(&n,&p);
    
    /* CALCULATES EXPECTED KMER COUNTS FROM MONONUCLEOTIDE DISTRIBUTION */
    for(current_kmer_number = 0; current_kmer_number < number_of_different_kmers; current_kmer_number++)
    {
    /* printf("\n%li: ", current_kmer_number); */
    for(current_kmer = current_kmer_number, score = 1, kmer_position = 0; kmer_position < kmer_length; kmer_position++, current_kmer >>= 2)
    {
    /* printf("%.2f * ", n.fraction[current_kmer & 3][kmer_position]); */
    score *= n.fraction[current_kmer & 3][kmer_position];
    }
    /* printf("1 = %.6f", score); */
    expected_counts[current_kmer_number] = 1000000 * score;
    total_expected_kmer_count += expected_counts[current_kmer_number];
    }    
    /* NORMALIZES */
    sizefactor = (double) total_kmer_count / (double) total_expected_kmer_count;
    for(current_kmer_number = 0; current_kmer_number < number_of_different_kmers; current_kmer_number++) expected_counts[current_kmer_number] = expected_counts[current_kmer_number] * sizefactor;
    
    /* CALCULATES UNCENTERED CORRELATION */
    for(current_kmer_number = 0; current_kmer_number < number_of_different_kmers; current_kmer_number++)
    {
        /* printf("\n%li\t%li", observed_counts[current_kmer_number], expected_counts[current_kmer_number]); */
        sum_x_squared += (double) pow(observed_counts[current_kmer_number], 2);
        sum_y_squared += (double) pow(expected_counts[current_kmer_number], 2);
        sum_xy += (double) observed_counts[current_kmer_number] * expected_counts[current_kmer_number]; 
    }
    correlation = sum_xy/(sqrt(sum_y_squared)*sqrt(sum_x_squared));
    
    count_pwm_free(&p);
    normalized_pwm_free(&n);
    return(correlation);
}

/* SUBROUTINE THAT CALCULATES INFORMATION CONTENT OF DINUCLEOTIDE DISTRIBUTIONS FOR A DINUCLEOTIDE BASE DEPENDENCY MATRIX */
short int Add_information_content (struct base_dependency_matrix *dep)
{
    short int first_base;
    short int second_base;
    short int dinucleotide;
    double information_content;
    short int ic_flag;
    double total;
    for (first_base = 0; first_base < (*dep).width; first_base++)
    {
        for (second_base = 0; second_base < (*dep).width; second_base++) 
        {
            if (first_base == second_base) continue;
            for (total = 0, dinucleotide = 0; dinucleotide < 16; dinucleotide++) total += (double) (*dep).incidence[first_base][second_base][dinucleotide];
            for (information_content = 0, ic_flag = 0, dinucleotide = 0; dinucleotide < 16; dinucleotide++) if ((*dep).incidence[first_base][second_base][dinucleotide] > 0) {information_content += ((*dep).incidence[first_base][second_base][dinucleotide] / total) * log2((*dep).incidence[first_base][second_base][dinucleotide] / total); ic_flag = 1;}
            if (ic_flag == 1) (*dep).information_content[first_base][second_base] = 4 + information_content;
        }
    }
    return(0);
}

/* SUBROUTINE THAT CALCULATES TOTAL RELATIVE DEVIATION OF DINUCLEOTIDE DISTRIBUTIONS FOR A DINUCLEOTIDE BASE DEPENDENCY MATRIX */
short int Add_total_relative_deviation (struct base_dependency_matrix *o, struct base_dependency_matrix *e)
{
    short int first_base;
    short int second_base;
    short int dinucleotide;
    short int end_trim = 5;
    double total_relative_deviation;
    double count_statistic_expected_total_relative_deviation;
    double total;
    long int *observed_error;
    observed_error = malloc(16 * sizeof(long int) + 5);
    double *propagated_error;
    propagated_error = malloc(16 * sizeof(double) + 5);
    
    for (first_base = 0; first_base < (*o).width; first_base++)
    {
        for (second_base = 0; second_base < (*o).width; second_base++) 
        {
            if (first_base == second_base) continue;
            for (total = 0, dinucleotide = 0; dinucleotide < 16; dinucleotide++) 
            {
            total += (double) (*o).incidence[first_base][second_base][dinucleotide];
            observed_error[dinucleotide] = (int) sqrt((*o).incidence[first_base][second_base][dinucleotide]);
            }
            for (total_relative_deviation = 0, count_statistic_expected_total_relative_deviation = 0, dinucleotide = 0; dinucleotide < 16; dinucleotide++) if ((*o).incidence[first_base][second_base][dinucleotide] >= 0) 
            {
                total_relative_deviation += abs((*o).incidence[first_base][second_base][dinucleotide] - (*e).incidence[first_base][second_base][dinucleotide])/ total; 
                Kmer_mono_correlation(observed_error,propagated_error, 2);
                count_statistic_expected_total_relative_deviation += propagated_error[dinucleotide] / total;
            }
            (*o).total_relative_deviation[first_base][second_base] = total_relative_deviation;
            (*o).count_statistic_expected_total_relative_deviation[first_base][second_base] = count_statistic_expected_total_relative_deviation;
            if(first_base > end_trim && second_base > end_trim && first_base < (*o).width - end_trim && second_base < (*o).width - end_trim)
            {
            if(total_relative_deviation > (*o).max_relative_deviation) (*o).max_relative_deviation = total_relative_deviation;
            if(count_statistic_expected_total_relative_deviation > (*o).max_expected_relative_deviation) (*o).max_expected_relative_deviation = count_statistic_expected_total_relative_deviation;
            }
        }
    }
    free(observed_error);
    free(propagated_error);
    return(0);
}

/* SUBROUTINE THAT CALCULATES EXPECTED SPACED DINUCLEOTIDE COUNTS BASED ON A PWM */
short int Expected_dinucleotides_pwm (struct base_dependency_matrix *dep, struct normalized_pwm *p, char *seed, long int counts)
{
    short int first_base;
    short int second_base;
    short int first_base2;
    short int second_base2;
    short int dinucleotide;
    signed short int first_pos_relative_to_seed;
    signed short int second_pos_relative_to_seed;
    short int seed_length = strlen(seed);
    double current_counts;
    
    for (first_base = (*p).width - 1; first_base >= 0; first_base--)
    {
        for (second_base = first_base - 1; second_base >= 0; second_base--) 
        {
            for (dinucleotide = 0; dinucleotide < 16; dinucleotide++) 
            {
                current_counts = (double) (counts * (*p).fraction[dinucleotide & 3][first_base] * (*p).fraction[(dinucleotide & 12) >> 2][second_base]);
                if (current_counts == (0.0 / 0.0) || current_counts == (1.0 / 0.0)) current_counts = 0;
                (*dep).incidence[first_base][second_base][dinucleotide] = current_counts;
            }
        }
    }
    
    
    for (first_base = (*p).width - 1; first_base >= 0; first_base--)
    {
        for (second_base = first_base - 1; second_base >= 0; second_base--) 
        {
            /* EVENS DINUCLEOTIDES IF SEED IS PALINDROMIC WITHOUT THE DINUCLEOTIDE */
            second_base2 = (*p).width - first_base + seed_length/2 + 3;
            first_base2 = (*p).width - second_base + seed_length/2 + 3;
            first_pos_relative_to_seed = first_base-Nlength+2;
            second_pos_relative_to_seed = second_base-Nlength+2;
            
            /* (second_pos_relative_to_seed == seed_length + 1 - first_pos_relative_to_seed) && */
            
            if (((second_pos_relative_to_seed == seed_length + 1 - first_pos_relative_to_seed) && first_pos_relative_to_seed > 0 && first_pos_relative_to_seed <= seed_length) && (second_pos_relative_to_seed > 0 && second_pos_relative_to_seed <= seed_length))
            {
                
                printf("\n%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%s\t_6:7", first_base2, second_base2, first_base, second_base, first_base2-Nlength+2, second_base2-Nlength+2, first_pos_relative_to_seed, second_pos_relative_to_seed, seed); 
                if(Is_this_string_iupac_palindrome(seed, first_pos_relative_to_seed, second_pos_relative_to_seed))
                {
                    printf("\tPALINDROME %li", (*dep).incidence[first_base][second_base][0]);
                    (*dep).incidence[first_base][second_base][0] = ((*dep).incidence[first_base][second_base][0] + (*dep).incidence[first_base2][second_base2][15]) / 2; /* AA & TT */
                    /* if (first_pos_relative_to_seed == 6 && second_pos_relative_to_seed == 1) (*dep).incidence[first_base][second_base][0] = -1000000 * first_base; */
                    (*dep).incidence[first_base2][second_base2][15] = (*dep).incidence[first_base][second_base][0]; /* AA & TT */
                    (*dep).incidence[first_base][second_base][1] = ((*dep).incidence[first_base][second_base][1] + (*dep).incidence[first_base2][second_base2][11]) / 2; /* AC & GT */
                    (*dep).incidence[first_base2][second_base2][11] = (*dep).incidence[first_base][second_base][1]; /* AC & GT */    
                    (*dep).incidence[first_base][second_base][2] = ((*dep).incidence[first_base][second_base][2] + (*dep).incidence[first_base2][second_base2][7]) / 2; /* AG & CT */
                    (*dep).incidence[first_base2][second_base2][7] = (*dep).incidence[first_base][second_base][2] ; /* AG & CT */
                    
                    (*dep).incidence[first_base][second_base][3] = ((*dep).incidence[first_base][second_base][3] + (*dep).incidence[first_base2][second_base2][3]) / 2; /* AT & self */
                    (*dep).incidence[first_base2][second_base2][3] = (*dep).incidence[first_base][second_base][3] ; /* AT & self */
                    
                    (*dep).incidence[first_base][second_base][4] = ((*dep).incidence[first_base][second_base][4] + (*dep).incidence[first_base2][second_base2][14]) / 2; /* CA & TG */
                    (*dep).incidence[first_base2][second_base2][14] = (*dep).incidence[first_base][second_base][4]; /* CA & TG */
                    (*dep).incidence[first_base][second_base][5] = ((*dep).incidence[first_base][second_base][5] + (*dep).incidence[first_base2][second_base2][10]) / 2; /* CC & GG */
                    (*dep).incidence[first_base2][second_base2][10] = (*dep).incidence[first_base][second_base][5]; /* CC & GG */ 
                    
                    (*dep).incidence[first_base][second_base][6] = ((*dep).incidence[first_base][second_base][6] + (*dep).incidence[first_base2][second_base2][6]) / 2; /* CG & self */
                    (*dep).incidence[first_base2][second_base2][6] = (*dep).incidence[first_base][second_base][6] ; /* CG & self */
                    
                    (*dep).incidence[first_base][second_base][8] = ((*dep).incidence[first_base][second_base][8] + (*dep).incidence[first_base2][second_base2][13]) / 2; /* GA & TC */
                    (*dep).incidence[first_base2][second_base2][13] = (*dep).incidence[first_base][second_base][8]; /* GA & TC */ 
                    
                    (*dep).incidence[first_base][second_base][9] = ((*dep).incidence[first_base][second_base][9] + (*dep).incidence[first_base2][second_base2][9]) / 2; /* GC & self */
                    (*dep).incidence[first_base2][second_base2][9] = (*dep).incidence[first_base][second_base][9] ; /* GC & self */
                    
                    (*dep).incidence[first_base][second_base][12] = ((*dep).incidence[first_base][second_base][12] + (*dep).incidence[first_base2][second_base2][12]) / 2; /* TA & self */
                    (*dep).incidence[first_base2][second_base2][12] = (*dep).incidence[first_base][second_base][12] ; /* TA & self */
                }
            }
        }
    }
    
    Add_information_content(dep);
    return(0);
}

/* SUBROUTINE THAT GENERATES A PERMUTATION OF 16 VALUES INTO dinucs FOR PERMUTATION OF DINUCLEOTIDES, RETURNS POINTER TO dinucs */
long int *Dinucleotide_permutation (long int *dinuc_counts)
{
	/* PERMUTATION OF DINUCLEOTIDES */
	unsigned int rseed = (unsigned int)time(NULL);
	srand (rseed);
	
	short int i;
	short int j;
	short int already_taken;
	short int *permutation;
	permutation = malloc(16*sizeof(short int)+5);
	long int *permutated_values;
	permutated_values = malloc(16*sizeof(long int)+5);
	short int *all;
	all = malloc(16*sizeof(short int)+5);
	long int misses = 0;
	long int test;
    
	/* LOOP TO GET ALL 16 RANDOMIZED POSITIONS */
	/*     for(test = 0; test < 1000000; test++)
	 { */
	for (i=0; i<16; i++) all[i] = -1;
	
	int current_random_number;
	
	for (i=0; i<16; i++)
	{
		for(already_taken = 1, misses = 0; already_taken == 1;)
		{
			already_taken = 0;
			current_random_number = rand () & 15;
			for(j=0; j < i; j++)
			{
				if (all[j] == current_random_number) 
				{
					misses++;    
					already_taken = 1;    
				}
			}
			
		}
		all[j] = current_random_number;
		permutation[i] = current_random_number;
		permutated_values[i] = dinuc_counts[current_random_number];
		/* printf ("%li\t%u\t%f\n", misses, permutation[i], permutated_values[i]); */
	}
	/* } */
	free (all);
	return (permutated_values);
}

/* SUBROUTINE THAT CALCULATES EXPECTED SPACED DINUCLEOTIDE COUNTS BY MAKING A MONONUCLEOTIDE MODEL */
short int Expected_dinucleotides (struct base_dependency_matrix *o, struct base_dependency_matrix *e, char *seed, long int counts)
{
    short int first_base;
    short int second_base;
    short int first_base2;
    short int second_base2;
    short int dinucleotide;
    signed short int first_pos_relative_to_seed;
    signed short int second_pos_relative_to_seed;
    short int seed_length = strlen(seed);
    double *current_counts;
	long int *permutated_counts;
    current_counts = malloc(16 * sizeof(double) + 5);
    
    for (first_base = (*o).width - 1; first_base >= 0; first_base--)
    {
        for (second_base = first_base - 1; second_base >= 0; second_base--) 
        {
			(*o).eo_correlation[first_base][second_base] = Kmer_mono_correlation((*o).incidence[first_base][second_base],current_counts, 2);
            for (dinucleotide = 0; dinucleotide < 16; dinucleotide++) (*e).incidence[first_base][second_base][dinucleotide] = (long int) current_counts[dinucleotide];
			permutated_counts = Dinucleotide_permutation((*o).incidence[first_base][second_base]);
			(*o).permutated_correlation[first_base][second_base] = Kmer_mono_correlation(permutated_counts,current_counts, 2);
			free(permutated_counts);
		}
    }
    Add_information_content(o);
    Add_total_relative_deviation(o,e);
	free(current_counts);
    return(0);
}  

/* SUBROUTINE THAT CALCULATES EXPECTED DINUCLEOTIDE COUNTS FOR TWO MISMATCHES BASED ON SINGLE MISMATCHES */
double Expected_mismatched_dinucleotides(short int dinucleotide, long int *counts) 
{
    short int counter;
    short int consensus = 0;
    double total = 0;
    double max = 0;
    double predicted_value = 0;
    for (counter = 0; counter < 16; counter++) 
    {
    total += counts[counter];
    if (counts[counter] > max) 
    {
    max = counts[counter];
    consensus = counter;
    }
    }
    if(consensus == dinucleotide) return(-0.1);
    if(((consensus & 3) == (dinucleotide & 3)) | ((consensus & 12) == (dinucleotide & 12))) return(-1.1);
    predicted_value = counts[(dinucleotide & 12) | (consensus & 3)] * counts[(dinucleotide & 3) | (consensus & 12)] / max;
    return(predicted_value);
}

/* SUBROUTINE THAT DETERMINES IF GAP IS AT EITHER OF THE CENTER POSITIONS, IF count_also is set to != 1 returns true */
short int Centergap (short int count_also_spaced_kmers, short int kmer_length, short int gap_position)
{
if (count_also_spaced_kmers != 1) return (1);
if (gap_position == kmer_length / 2) return (1); 
if (kmer_length % 2 == 1 && gap_position == kmer_length / 2 - 1) return (1); 
else return (0);
}

/* SUBROUTINE THAT DETERMINES IF A GAPPED KMER IS A LOCAL MAXIMUM */
short int Localmax(long int *****results, short int file_number, short int shortest_kmer, short int too_long_kmer, short int current_kmer_length,short int current_gap_position, short int current_gap_length, long int current_kmer, short int count_also_spaced_kmers, double kmer_length_difference_cutoff, char *kmerstring)
{
    long int kmer1_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][current_kmer];
    long int kmer2_incidence;
    long int compared_kmer = current_kmer;
    signed short int position;
    short int nucleotide;
    short int counter;
    short int higher_counts = 0;
    short int first_half;
    long int lowbit = 1;
    long int highbit = 2;
    short int shift;
    short int true_gap_position = current_gap_position;
    short int true_gap_length = current_gap_length;
    short int start;
    short int end;
    short int left = 0;
    short int right = 1;
    signed short int position2;
    char *LR = "LR";
    char *forward;
    if (rna == 0) forward = dnaforward;
    else forward = rnaforward;

    kmerstring[0] = '-';
    for(position2 = current_kmer_length-1, position = 0; position2 > -1 ; position2--, position++) 
    {
    if(current_kmer_length - position2 - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++, position++) kmerstring[position] = 'n'; 
    kmerstring[position] = forward[(compared_kmer & mask_ULL[1][position2]) >> (position2 * 2)];
    }
    kmerstring[position] = '\0';

    
    // if (print_local_max) printf("\n%s count: %li", kmerstring, kmer1_incidence); 
    
    /* HAMMING OF 1 */
    for(position=0; position < current_kmer_length; position++, lowbit <<= 2, highbit <<= 2) 
    {
        compared_kmer = lowbit ^ current_kmer;
        kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
        // if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for  Hamming1: ", kmer2_incidence); for(position2 = current_kmer_length-1; position2 > -1 ; position2--) {if(current_kmer_length - position2 - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position2]) >> (position2 * 2)]);}}
        if(kmer2_incidence > kmer1_incidence) return (0);
        compared_kmer = highbit ^ current_kmer;
        kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
        // if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for  Hamming1: ",kmer2_incidence); for(position2 = current_kmer_length-1; position2 > -1 ; position2--) {if(current_kmer_length - position2 - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position2]) >> (position2 * 2)]);}}
        if(kmer2_incidence > kmer1_incidence) return (0);
        compared_kmer = lowbit ^ compared_kmer;
        kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
        // if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for  Hamming1: ",kmer2_incidence); for(position2 = current_kmer_length-1; position2 > -1 ; position2--) {if(current_kmer_length - position2 - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position2]) >> (position2 * 2)]);}}
        if(kmer2_incidence > kmer1_incidence) return (0);
    }
    
    /* FULL SHIFT BY ONE */
    shift = (current_gap_length != 0);
    current_gap_position += shift;
    
    /* ONLY LOOK AT FULL SHIFT FOR UNGAPPED KMERS, OR FOR GAPPED KMERS IF GAP CAN SHIFT ALSO (ALL GAP POSITIONS HAVE BEEN COUNTED) */
    if (current_gap_length == 0 | (count_also_spaced_kmers == 2 && current_gap_position < current_kmer_length) | (current_kmer_length % 2 == 1 && count_also_spaced_kmers == 1 && current_gap_position == current_kmer_length / 2 + 1))
    {
    /* COMPARED FULL SHIFT RIGHT (same shift in gap if any) */
    compared_kmer = (current_kmer >> 2) & lowmask_ULL[current_kmer_length-1];
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
// if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for  right shift: ",kmer2_incidence); for(position = current_kmer_length-1; position > -1 ; position--) {if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position]) >> (position * 2)]);}}
    if(kmer2_incidence > kmer1_incidence) return (0);
    lowbit >>= 2; highbit >>= 2;
    compared_kmer = lowbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
// if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for  right shift: ",kmer2_incidence); for(position = current_kmer_length-1; position > -1 ; position--) {if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position]) >> (position * 2)]);}}
    if(kmer2_incidence > kmer1_incidence) return (0);
    compared_kmer = highbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
// if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for  right shift: ",kmer2_incidence); for(position = current_kmer_length-1; position > -1 ; position--) {if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position]) >> (position * 2)]);}}
    if(kmer2_incidence > kmer1_incidence) return (0);
    compared_kmer = lowbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
// if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for  right shift: ",kmer2_incidence); for(position = current_kmer_length-1; position > -1 ; position--) {if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position]) >> (position * 2)]);}}
    if(kmer2_incidence > kmer1_incidence) return (0);
    }
    
    current_gap_position -= shift;
    current_gap_position -= shift;
    if (current_gap_length == 0 | (count_also_spaced_kmers == 2 && current_gap_position > 0) | (current_kmer_length % 2 == 1 && count_also_spaced_kmers == 1 && current_gap_position == current_kmer_length / 2))
    {
    /* COMPARED FULL SHIFT LEFT (same shift in gap if any) */    
    compared_kmer = (current_kmer << 2) & lowmask_ULL[current_kmer_length-1];
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
// if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for  left shift: ",kmer2_incidence); for(position = current_kmer_length-1; position > -1 ; position--) {if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position]) >> (position * 2)]);}}
    if(kmer2_incidence > kmer1_incidence) return (0); 
    lowbit = 1; highbit = 2;
    compared_kmer = lowbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
// if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for  left shift: ",kmer2_incidence); for(position = current_kmer_length-1; position > -1 ; position--) {if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position]) >> (position * 2)]);}}
    if(kmer2_incidence > kmer1_incidence) return (0);
    compared_kmer = highbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
// if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for  left shift: ",kmer2_incidence); for(position = current_kmer_length-1; position > -1 ; position--) {if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position]) >> (position * 2)]);}}    
    if(kmer2_incidence > kmer1_incidence) return (0);
    compared_kmer = lowbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
// if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for  left shift: ",kmer2_incidence); for(position = current_kmer_length-1; position > -1 ; position--) {if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position]) >> (position * 2)]);}}
    if(kmer2_incidence > kmer1_incidence) return (0);
    }
    
    current_gap_position += shift;
    lowbit = 1; highbit = 2;
    if (count_also_spaced_kmers != 0)
    {
    if(current_gap_position == 0) current_gap_position = current_kmer_length / 2;
    
        
    /* COMPARE TO KMER WITH LONGER GAP */
    current_gap_length++;
    if (current_gap_length < Nlength - current_kmer_length && true_gap_length != 0)
    {
    compared_kmer = (current_kmer & highmask_ULL[current_kmer_length-current_gap_position]) | ((current_kmer << 2) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
    
    /* LOOP TO ANALYZE SHIFT OF EITHER HALF, STARTS WITH SECOND HALF */
    for(first_half = 0; first_half < 2; first_half++)
    {
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
// if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for  longer gap %c: ",kmer2_incidence, LR[first_half]); for(position = current_kmer_length-1; position > -1 ; position--) {if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position]) >> (position * 2)]);}}
    if(kmer2_incidence > kmer1_incidence) return (0);  
    compared_kmer = lowbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
// if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for  longer gap %c: ",kmer2_incidence, LR[first_half]); for(position = current_kmer_length-1; position > -1 ; position--) {if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position]) >> (position * 2)]);}}
    if(kmer2_incidence > kmer1_incidence) return (0);
    compared_kmer = highbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
// if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for  longer gap %c: ",kmer2_incidence, LR[first_half]); for(position = current_kmer_length-1; position > -1 ; position--) {if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position]) >> (position * 2)]);}}
    if(kmer2_incidence > kmer1_incidence) return (0);
    compared_kmer = lowbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
// if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for  longer gap %c: ",kmer2_incidence, LR[first_half]); for(position = current_kmer_length-1; position > -1 ; position--) {if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position]) >> (position * 2)]);}}
    if(kmer2_incidence > kmer1_incidence) return (0);
    
    /* SWITCHES TO FIRST HALF */
    compared_kmer = ((current_kmer >> 2) & highmask_ULL[current_kmer_length-current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
    lowbit <<= ((current_kmer_length - 1)*2); highbit <<= ((current_kmer_length - 1)*2);
    }
    }
    
        
    /* COMPARE TO KMER WITH SHORTER GAP */
    current_gap_length--;
    current_gap_length--;

    lowbit = 1; highbit = 2;
    lowbit <<= ((current_kmer_length - true_gap_position - 1)*2); highbit <<= ((current_kmer_length - true_gap_position - 1)*2);
    if (current_gap_length == 0) current_gap_position = 0;
    if (current_gap_length >= 0)
    {
        compared_kmer = (current_kmer & highmask_ULL[current_kmer_length - true_gap_position]) | ((current_kmer >> 2) & (lowmask_ULL[current_kmer_length-true_gap_position-1]));
        /* LOOP TO ANALYZE SHIFT OF EITHER HALF, STARTS WITH SECOND HALF */
        for(first_half = 0; first_half < 2; first_half++)
        {
            
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
// if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for  shorter gap: ",kmer2_incidence); for(position = current_kmer_length-1; position > -1 ; position--) {if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position]) >> (position * 2)]);}}
            if(kmer2_incidence > kmer1_incidence) return (0);  
            compared_kmer = lowbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
// if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for  shorter gap: ",kmer2_incidence); for(position = current_kmer_length-1; position > -1 ; position--) {if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position]) >> (position * 2)]);}}
            if(kmer2_incidence > kmer1_incidence) return (0);
            compared_kmer = highbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
// if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for  shorter gap: ",kmer2_incidence); for(position = current_kmer_length-1; position > -1 ; position--) {if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position]) >> (position * 2)]);}}
// printf("gpos, glen, kmer1, kmer2, higher_counts %i, %i, %i, %i, %i", current_gap_position, current_gap_length, kmer1_incidence, kmer2_incidence, higher_counts);
            if(kmer2_incidence > kmer1_incidence) return (0);
            compared_kmer = lowbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
// if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for  shorter gap: ",kmer2_incidence); for(position = current_kmer_length-1; position > -1 ; position--) {if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position]) >> (position * 2)]);}}
            if(kmer2_incidence > kmer1_incidence) return (0);
            
            /* SWITCHES TO FIRST HALF */
            compared_kmer = lowmask_ULL[current_kmer_length-1] & ((current_kmer << 2) & highmask_ULL[current_kmer_length - true_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-true_gap_position-1]));
            lowbit <<= 2; highbit <<= 2;
        }
    }
    
    
        /* COMPARES DIFFERENT GAP POSITIONS (SHIFTED BY ONE) */
        current_gap_length = true_gap_length;
        if ((count_also_spaced_kmers == 2 || (count_also_spaced_kmers == 1 && current_kmer_length % 2 == 1)) && true_gap_length > 0)
        {
        current_gap_position = true_gap_position + 1;
        lowbit = 1; highbit = 2;
        lowbit <<= ((current_kmer_length - current_gap_position)*2); highbit <<= ((current_kmer_length - current_gap_position)*2);
        compared_kmer = (current_kmer & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));

        /* LOOP TO ANALYZE SHIFT TO EITHER DIRECTION, STARTS WITH RIGHT BY ONE */
        for(first_half = 0; first_half < 2; first_half++)
        {
         //   printf("\nLOOP");
        if (current_gap_position < current_kmer_length && current_gap_position > 0 && (count_also_spaced_kmers == 2 || 
        (count_also_spaced_kmers == 1 && (current_gap_position == current_kmer_length / 2) || current_gap_position == current_kmer_length / 2 + 1)))
        {
        kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
// if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for  gap shift %c: ",kmer2_incidence, LR[1-first_half]); for(position = current_kmer_length-1; position > -1 ; position--) {if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position]) >> (position * 2)]);}}
        if(kmer2_incidence > kmer1_incidence) return (0);  
        compared_kmer = lowbit ^ compared_kmer;
        kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
// if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for  gap shift %c: ",kmer2_incidence, LR[1-first_half]); for(position = current_kmer_length-1; position > -1 ; position--) {if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position]) >> (position * 2)]);}}
        if(kmer2_incidence > kmer1_incidence) return (0);
        compared_kmer = highbit ^ compared_kmer;
        kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
// if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for  gap shift %c: ",kmer2_incidence, LR[1-first_half]); for(position = current_kmer_length-1; position > -1 ; position--) {if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position]) >> (position * 2)]);}}
        // printf("gpos, glen, kmer1, kmer2, higher_counts %i, %i, %i, %i, %i", current_gap_position, current_gap_length, kmer1_incidence, kmer2_incidence, higher_counts);
        if(kmer2_incidence > kmer1_incidence) return (0);
        compared_kmer = lowbit ^ compared_kmer;
        kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
// if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for  gap shift %c: ",kmer2_incidence, LR[1-first_half]); for(position = current_kmer_length-1; position > -1 ; position--) {if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position]) >> (position * 2)]);}}
        if(kmer2_incidence > kmer1_incidence) return (0);
        }
            
        /* SWITCHES LEFT BY ONE */
        current_gap_position--;
        current_gap_position--;
        compared_kmer = ((current_kmer) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
        lowbit <<= 2; highbit <<= 2;
        }
        }        
        
        start = 1;
        end = current_kmer_length;
        
        /* COMPARES UNGAPPED KMER TO ALL SINGLE GAPS IN ALL POSITIONS */
        if (count_also_spaced_kmers != 0 && true_gap_length == 0)
        {
            if(count_also_spaced_kmers == 1)
            {
                start = current_kmer_length / 2;
                end = start + 1 + (current_kmer_length % 2);
            }
        for(current_gap_position = start, current_gap_length = 1; current_gap_position < end; current_gap_position++) 
        {
        compared_kmer = (current_kmer & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer << 2) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
            /* LOOP TO ANALYZE SHIFT TO EITHER DIRECTION, STARTS WITH BEGINNING */
            for(lowbit = 1, highbit = 2, first_half = 0; first_half < 2; first_half++)
            {
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
// if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for  gapped: ",kmer2_incidence); for(position = current_kmer_length-1; position > -1 ; position--) {if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position]) >> (position * 2)]);}}
            if(kmer2_incidence > kmer1_incidence) return (0);  
            compared_kmer = lowbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
// if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for  gapped: ",kmer2_incidence); for(position = current_kmer_length-1; position > -1 ; position--) {if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position]) >> (position * 2)]);}}
            if(kmer2_incidence > kmer1_incidence) return (0);
            compared_kmer = highbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
// if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for  gapped: ",kmer2_incidence); for(position = current_kmer_length-1; position > -1 ; position--) {if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position]) >> (position * 2)]);}}
            // printf("gpos, glen, kmer1, kmer2, higher_counts %i, %i, %i, %i, %i", current_gap_position, current_gap_length, kmer1_incidence, kmer2_incidence, higher_counts);
            if(kmer2_incidence > kmer1_incidence) return (0);
            compared_kmer = lowbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
// if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for  gapped: ",kmer2_incidence); for(position = current_kmer_length-1; position > -1 ; position--) {if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position]) >> (position * 2)]);}}
            if(kmer2_incidence > kmer1_incidence) return (0);
           
        /* END */
        compared_kmer = ((current_kmer >> 2) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
            lowbit <<= ((current_kmer_length - 1)*2); highbit <<= ((current_kmer_length - 1)*2);
        }
        }

            
        }
        
    }
    
/* COMPARES KMER WITH ONE SHORTER */    
current_gap_position = true_gap_position;
current_gap_length = true_gap_length;
current_kmer_length--;
end = current_kmer_length;
if (current_kmer_length >= shortest_kmer)
{
/* IF NO GAP, INSERTS GAP AT ALL ALLOWED POSITIONS */
if(current_gap_length == 0)
{
    if (count_also_spaced_kmers != 0)
    {
        if(count_also_spaced_kmers == 1)
        {
            start = current_kmer_length / 2;
            end = start + 1 + (current_kmer_length % 2);
        }
        for(current_gap_position = start, current_gap_length = 1; current_gap_position < end; current_gap_position++) 
        {
            compared_kmer = ((current_kmer >> 2) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
                kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
// if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for  shorter with gap: ", ((long int) (kmer2_incidence * kmer_length_difference_cutoff))); for(position = current_kmer_length-1; position > -1 ; position--) {if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position]) >> (position * 2)]);}}
                if(kmer2_incidence  * kmer_length_difference_cutoff > kmer1_incidence) return (0);  
        }
        
        
    }
    
}

current_gap_position = true_gap_position;
current_gap_length = true_gap_length; 

if (count_also_spaced_kmers != 1) {left = 1; right = 1;}
if (current_gap_position == current_kmer_length / 2 && current_kmer_length % 2 == 0 && count_also_spaced_kmers == 1) {left = 1; right = 0;}
if (current_kmer_length % 2 == 1 && count_also_spaced_kmers == 1) left = 1;

    
/* LEFT PART */
if (current_gap_position < current_kmer_length) 
{
    if (left == 1 || true_gap_length == 0)
    {
compared_kmer = (current_kmer >> 2);
kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
// if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for  shorter left: ", ((long int) (kmer2_incidence * kmer_length_difference_cutoff))); for(position = current_kmer_length-1; position > -1 ; position--) {if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position]) >> (position * 2)]);}}
if(kmer2_incidence * kmer_length_difference_cutoff > kmer1_incidence) return (0);
    }
    }
/* RIGHT PART */

if (current_gap_position != 1 && right == 1) 
{
if (current_gap_position > 0) current_gap_position--;
    compared_kmer = (current_kmer & lowmask_ULL[current_kmer_length-1]);
kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
// if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for  shorter right: ", ((long int) (kmer2_incidence * kmer_length_difference_cutoff))); for(position = current_kmer_length-1; position > -1 ; position--) {if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position]) >> (position * 2)]);}}
if(kmer2_incidence * kmer_length_difference_cutoff > kmer1_incidence) return (0);
}
    
current_gap_position = true_gap_position;
/* LONGER GAP */
if(count_also_spaced_kmers != 0 && true_gap_position != 0)
{
current_gap_length++;
if (current_gap_position < current_kmer_length && left == 1) 
{
/* RIGHT BASE GAPPED */
compared_kmer = ((current_kmer >> 2) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
// if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for  shorter left with longer gap: ", ((long int) (kmer2_incidence * kmer_length_difference_cutoff))); for(position = current_kmer_length-1; position > -1 ; position--) {if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position]) >> (position * 2)]);}}
if(kmer2_incidence * kmer_length_difference_cutoff > kmer1_incidence) return (0);
}
    
/* LEFT BASE GAPPED */
if(current_gap_position > 1 && right == 1)
{
current_gap_position--;
compared_kmer = ((current_kmer >> 2) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
// if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for  shorter right with longer gap: ", ((long int) (kmer2_incidence * kmer_length_difference_cutoff))); for(position = current_kmer_length-1; position > -1 ; position--) {if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position]) >> (position * 2)]);}}
if(kmer2_incidence * kmer_length_difference_cutoff > kmer1_incidence) return (0);
}
current_gap_length--;
current_gap_position++;
}
 

    
    /* COMPARES HANGING SINGLE BASE TO UNGAPPED KMER */
    current_gap_position = true_gap_position;
    current_gap_length = true_gap_length;
    if(current_gap_position == 1)
    {
    current_gap_length = 0;
    current_gap_position = 0;
    compared_kmer = (current_kmer & lowmask_ULL[current_kmer_length-1]);  
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
// if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for : ", ((long int) (kmer2_incidence * kmer_length_difference_cutoff))); for(position = current_kmer_length-1; position > -1 ; position--) {if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position]) >> (position * 2)]);}}
    if(kmer2_incidence * kmer_length_difference_cutoff > kmer1_incidence) return (0);
    }
    else if(current_kmer_length-current_gap_position == 0)
    {
    current_gap_length = 0;
    current_gap_position = 0;
    compared_kmer = (current_kmer >> 2);   
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
// if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for : ", ((long int) (kmer2_incidence * kmer_length_difference_cutoff))); for(position = current_kmer_length-1; position > -1 ; position--) {if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position]) >> (position * 2)]);}}
    if(kmer2_incidence * kmer_length_difference_cutoff > kmer1_incidence) return (0);
     
    }
}
current_gap_position = true_gap_position;
current_gap_length = true_gap_length;
    
/* COMPARES KMER WITH ONE LONGER */
current_kmer_length++;
current_kmer_length++;
if (current_kmer_length < too_long_kmer)
{
compared_kmer = (current_kmer << 2);
/* LOOP TO ANALYZE SHIFT TO EITHER DIRECTION, STARTS WITH BEGINNING */
for(lowbit = 1, highbit = 2, first_half = 0; first_half < 2; first_half++)
{
if(count_also_spaced_kmers != 1 || true_gap_length == 0 || (first_half == 1 || (count_also_spaced_kmers == 1 && current_gap_position == current_kmer_length / 2)))
{
kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    // if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for longer : %i-", ((long int) (kmer2_incidence / kmer_length_difference_cutoff)), first_half); for(position = current_kmer_length-1; position > -1 ; position--) {if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position]) >> (position * 2)]);}}
if(kmer2_incidence > kmer1_incidence * kmer_length_difference_cutoff) return (0);  
compared_kmer = lowbit ^ compared_kmer;
kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
// if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for longer : ", ((long int) (kmer2_incidence / kmer_length_difference_cutoff))); for(position = current_kmer_length-1; position > -1 ; position--) {if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position]) >> (position * 2)]);}}
if(kmer2_incidence > kmer1_incidence * kmer_length_difference_cutoff) return (0);
compared_kmer = highbit ^ compared_kmer;
kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
// if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for longer : ", ((long int) (kmer2_incidence / kmer_length_difference_cutoff))); for(position = current_kmer_length-1; position > -1 ; position--) {if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position]) >> (position * 2)]);}}
// printf("gpos, glen, kmer1, kmer2, higher_counts %i, %i, %i, %i, %i", current_gap_position, current_gap_length, kmer1_incidence, kmer2_incidence, higher_counts);
if(kmer2_incidence > kmer1_incidence * kmer_length_difference_cutoff) return (0);
compared_kmer = lowbit ^ compared_kmer;
kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
// if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for longer : ", ((long int) (kmer2_incidence / kmer_length_difference_cutoff))); for(position = current_kmer_length-1; position > -1 ; position--) {if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position]) >> (position * 2)]);}}
if(kmer2_incidence > kmer1_incidence * kmer_length_difference_cutoff) return (0);
}
/* END */
compared_kmer = current_kmer;
lowbit <<= ((current_kmer_length - 1)*2); highbit <<= ((current_kmer_length - 1)*2);
if (true_gap_length == 0) continue;
current_gap_position++;
// printf("\n %i,%i,%i,%i,%i,%i, %i", true_gap_position, current_gap_position, current_kmer_length, count_also_spaced_kmers, current_gap_position,current_kmer_length / 2 + (current_kmer_length % 2), current_kmer_length % 2);
if (current_gap_position > current_kmer_length || (count_also_spaced_kmers == 1 && current_gap_position != (current_kmer_length / 2 + (current_kmer_length % 2)))) break; 
}
current_gap_position = true_gap_position;
current_gap_length = true_gap_length;
    
/* COMPARES TO ONE SHORTER GAP LENGTH */
if(count_also_spaced_kmers != 0 && true_gap_length >= 1)
{


current_gap_length--;

    
lowbit = 1; highbit = 2;
lowbit <<= ((current_kmer_length - current_gap_position-1)*2); highbit <<= ((current_kmer_length - current_gap_position-1)*2);
compared_kmer = ((current_kmer << 2) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
/* NO GAP LEFT */
if (current_gap_length == 0) current_gap_position = 0;
    
    /* LOOP TO ANALYZE ADDED BASE ON EITHER SIDE OF GAP, STARTS WITH RIGHT */
    for(first_half = 0; first_half < 2; first_half++)
    {
        if (current_gap_position < current_kmer_length && (current_gap_position == 0 || Centergap (count_also_spaced_kmers, current_kmer_length, current_gap_position)))
        {
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
// if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for longer with shorter gap: ", ((long int) (kmer2_incidence / kmer_length_difference_cutoff))); for(position = current_kmer_length-1; position > -1 ; position--) {if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position]) >> (position * 2)]);}}
            if(kmer2_incidence > kmer1_incidence * kmer_length_difference_cutoff) return (0);   
            compared_kmer = lowbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
// if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for longer with shorter gap : ", ((long int) (kmer2_incidence / kmer_length_difference_cutoff))); for(position = current_kmer_length-1; position > -1 ; position--) {if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position]) >> (position * 2)]);}}
            if(kmer2_incidence > kmer1_incidence * kmer_length_difference_cutoff) return (0); 
            compared_kmer = highbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
// if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for longer with shorter gap : ", ((long int) (kmer2_incidence / kmer_length_difference_cutoff))); for(position = current_kmer_length-1; position > -1 ; position--) {if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position]) >> (position * 2)]);}}
            // printf("gpos, glen, kmer1, kmer2, higher_counts %i, %i, %i, %i, %i", current_gap_position, current_gap_length, kmer1_incidence, kmer2_incidence, higher_counts);
            if(kmer2_incidence > kmer1_incidence * kmer_length_difference_cutoff) return (0); 
            compared_kmer = lowbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
// if (print_local_max) {printf("\n%s", kmerstring); printf("  compared to %li for longer with shorter gap : ", ((long int) (kmer2_incidence / kmer_length_difference_cutoff))); for(position = current_kmer_length-1; position > -1 ; position--) {if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(compared_kmer & mask_ULL[1][position]) >> (position * 2)]);}}
            if(kmer2_incidence > kmer1_incidence * kmer_length_difference_cutoff) return (0); 
        }
        
        /* SWITCHES ADDED BASE TO LEFT */
        current_gap_position++;
        if (current_gap_position > current_kmer_length || current_gap_position == 1) break;
        
        compared_kmer = ((current_kmer << 2) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
    }
current_gap_length++;
}        

    
}
current_gap_position = true_gap_position;
    
if(higher_counts == 0) return (1); else return (0);
}

struct alignment {short int strand; short int position;};

/* CONVERTS TWO-BIT SEQ TO BITIUPAC SEQ (INCLUDES 4bit wide SINGLE ULL REPRESENTATION IF LENGTH < 33) */
Bitseq_to_bitiupac(__uint128_t sequence_value_ULL, short int length, short int gap_position, short int gap_length, struct bitiupac_structure *s)
{
__uint128_t out_ULL = 0;
short int position;
short int current_base;

for(position = 0; position < length + gap_length; position++, sequence_value_ULL >>= 2)
{
if (position == gap_position) 
{
out_ULL <<= (gap_length << 2);
for (current_base = 0; current_base < 4; current_base++) (*s).base[current_base] <<= (gap_length << 2);
}
else 
{
out_ULL != (1 << ((sequence_value_ULL & 3) + (position << 2)));
(*s).base[sequence_value_ULL & 3] |= 1 << position;
}
}
if (position < 33) (*s).sequence_value_ULL = out_ULL;
(*s).length = length + gap_length;
}

/* REVERSE COMPLEMENTS BITIUPAC (IF LENGTH LESS THAN 33 ALSO REVCOMPS 4bit wide ULL VALUE */
void Reverse_complement_bitiupac(struct bitiupac_structure *s)
{
__uint128_t in_ULL = (*s).sequence_value_ULL;
__uint128_t out_ULL = 0;
short int position;
short int current_base;
    
__uint128_t *base_ULL;
base_ULL = malloc (sizeof(__uint128_t)*4+5);
for(current_base = 0; current_base < 4; current_base++) base_ULL[current_base] = 0;
    
for(position = 0; position < (*s).length; position++)
{
for(current_base = 0; current_base < 4; current_base++, base_ULL[current_base] <<= 1, (*s).base[current_base ^ 3] >>= 1)
{
base_ULL[current_base] != ((*s).base[current_base ^ 3] & 1);
}
for(current_base = 0; current_base < 4; current_base++) (*s).base[current_base] = base_ULL[current_base];
} 

if ((*s).length < 33)
{
for(position = 0; position < (*s).length; position++, in_ULL >>=4, out_ULL <<= 4) out_ULL |= (in_ULL ^ 15) & 15;
(*s).sequence_value_ULL = out_ULL >> 4;
}
free (base_ULL);
}

/* ADDS FLANKS TO BITIUPAC */
void Add_flanks_to_bitiupac(struct bitiupac_structure *s, short int flank_length, short int flank)
{
__uint128_t in_ULL = (*s).sequence_value_ULL;
__uint128_t out_ULL = 0;
short int position;
short int current_base;
 
if((*s).length + 2 * flank_length > 33)
{
for(position = 0; position < flank_length; position++) {out_ULL |= flank; out_ULL <<= 4;}
out_ULL <<= (((*s).length - 1) * 4);
out_ULL |= in_ULL;
for(position = 0; position < flank_length; position++) {out_ULL <<= 4; out_ULL |= flank;}
(*s).sequence_value_ULL = out_ULL;
}
else (*s).sequence_value_ULL = 0;
    
for(position = 0; position < flank_length; position++) for(current_base = 0; current_base < 4; current_base++,(*s).base[current_base] <<= 1) if ((flank & (1 << current_base)) != 0) (*s).base[current_base] |= 1;
    
for(current_base = 0; current_base < 4; current_base++) (*s).base[current_base] <<= (((*s).length - 1) * 4);
    
for(position = 0; position < flank_length; position++) for(current_base = 0; current_base < 4; current_base++, (*s).base[current_base] <<= 1) if ((flank & (1 << current_base)) != 0) (*s).base[current_base] |= 1;    
    
(*s).length += 2 * flank_length;
}

/* ALIGNS TWO BITIUPAC SEQUENCES, FINDS MINIMUM DISTANCE AND ALIGNMENT POSITION */
short int Align(struct alignment *a, struct bitiupac_structure *s1, struct bitiupac_structure *s2)
{    
    short int position_l;
    short int strand = 0;
    short int match_position = 0;
    short int min_score = 1000;
    short int score;
    short int current_base;
    __uint128_t  v;
    unsigned int c;
    __uint128_t l_ULL;
    __uint128_t s_ULL;
    
    struct bitiupac_structure *s;
    struct bitiupac_structure *l;
    if ((*s1).length > (*s2).length) {l = s1; s = s2;}
    else {l = s2; s = s1;}
    
    Add_flanks_to_bitiupac(l, (*s).length-1, 0); /* ADDS ? (ALWAYS MISMATCHING BASE) TO FLANKS */
    
    /* FORWARD */
    for (l_ULL = (*l).sequence_value_ULL, s_ULL = (*s).sequence_value_ULL, position_l = 0; position_l < (*l).length - (*s).length; position_l++, l_ULL >>= 1) 
    {
    for(c = 0, current_base = 0; current_base < 4; current_base++)
    {
    v = (*s).base[current_base] & (*l).base[current_base];
    for (; v; c++) v &= v - 1;  
    }
    score = (*s).length - v;
    if (score < min_score) {min_score = score; match_position = position_l; strand = 0;}    
    }
    /* REVERSE COMPLEMENT */
    Reverse_complement_bitiupac(l);
    for (l_ULL = (*l).sequence_value_ULL, s_ULL = (*s).sequence_value_ULL, position_l = 0; position_l < (*l).length - (*s).length; position_l++, l_ULL >>= 1) 
    {
        for(c = 0, current_base = 0; current_base < 4; current_base++)
        {
            v = (*s).base[current_base] & (*l).base[current_base];
            for (; v; c++) v &= v - 1;  
        }
        score = (*s).length - v;
        if (score < min_score) {min_score = score; match_position = position_l; strand = 0;}    
    }
    (*a).position = position_l;
    (*a).strand = strand;
    return(min_score);
}

/* SHAPE CALLER */
double Peak_shape(long int ***a, short int current_gap_position, long int current_kmer, short int length)
{
double shape;
double flanks = 1;
long int current_gap_length;
long int kmer_incidence;
long int max_value = 0;
long int max_value_position = -1;
/* FINDS GLOBAL MAXIMUM */
for (current_gap_length = 0; current_gap_length < length; current_gap_length++)
{
if (current_gap_length == 0) kmer_incidence = a[0][0][current_kmer];
else kmer_incidence = a[current_gap_position][current_gap_length][current_kmer];
    
if (kmer_incidence > max_value) {max_value = kmer_incidence; max_value_position = current_gap_length;}
}
    
double median;
double *normalized_values;
normalized_values = malloc(sizeof(double)*length+5);
for (current_gap_length = 0; current_gap_length < length; current_gap_length++)
{
if (current_gap_length == 0) kmer_incidence = a[0][0][current_kmer];
else kmer_incidence = a[current_gap_position][current_gap_length][current_kmer];
    
normalized_values[current_gap_length] = (double) kmer_incidence;   
}
    
qsort (normalized_values, length, sizeof(double), Sort_according_to_double);
median = normalized_values[length/2];

for (current_gap_length = 0; current_gap_length < length; current_gap_length++)
{
if (current_gap_length == 0) kmer_incidence = a[0][0][current_kmer];
else kmer_incidence = a[current_gap_position][current_gap_length][current_kmer];
    
normalized_values[current_gap_length] = (double) kmer_incidence - median;
}
    
shape = normalized_values[max_value_position];

if (max_value_position > 0 && max_value_position < length-2) flanks = normalized_values[max_value_position-1] + normalized_values[max_value_position+1];
if (max_value_position == 0) flanks = 2 * normalized_values[1];
if (max_value_position == length-1) flanks = 2 * normalized_values[length-2];
if (flanks < 0) flanks = 0;
shape /= flanks+0.1;
free(normalized_values);    
return shape;
}

/* RETURNS DISTANCE FROM TOP PEAK TO SECOND HIGHEST PEAK (LOCAL MAX MUST BE > cutoff * GLOBAL)*/
double *Peak_to_peak_distance(long int ***a, short int current_gap_position, long int current_kmer, short int length, double cutoff, double *returned_values)
{
long int current_gap_length;
long int kmer_incidence;
long int shorter_gap_incidence;
long int longer_gap_incidence;
    
long int global_max;
long int first_local_max_value = 0;
long int first_local_max_position = -1;
long int max_value = 0;
long int max_value_position = -1;
    
/* FINDS GLOBAL MAXIMUM */

for (current_gap_length = 0; current_gap_length < length; current_gap_length++)
{
if (current_gap_length == 0) kmer_incidence = a[0][0][current_kmer];
else kmer_incidence = a[current_gap_position][current_gap_length][current_kmer];
    
if (kmer_incidence > max_value) {max_value =  kmer_incidence; max_value_position = current_gap_length;}
}

/* FINDS FIRST LOCAL MAXIMUM BELOW GLOBAL MAX */
for (current_gap_length = 0; current_gap_length < length; current_gap_length++)
{
if (current_gap_length == 0) kmer_incidence = a[0][0][current_kmer];
else kmer_incidence = a[current_gap_position][current_gap_length][current_kmer];
    
if (current_gap_length == max_value_position) continue;
if (kmer_incidence > first_local_max_value) 
{
    
if (current_gap_length == 1) shorter_gap_incidence = a[0][0][current_kmer];
else if (current_gap_length > 1) shorter_gap_incidence = a[current_gap_position][current_gap_length-1][current_kmer];
    
if (current_gap_length > 1) if (kmer_incidence < shorter_gap_incidence) continue;
if (current_gap_length < length-2) if (kmer_incidence < a[current_gap_position][current_gap_length+1][current_kmer]) continue;
first_local_max_value = kmer_incidence;
first_local_max_position = current_gap_length;    
}
}
returned_values[0] = abs(max_value_position-first_local_max_position);
returned_values[1] = (double) first_local_max_value / ((double) max_value + 0.1);
if (returned_values[1] < cutoff) returned_values[0] = 0;
if (first_local_max_position == -1 || max_value_position == -1) returned_values[0] = 0;
return(returned_values);
}


short int First_digit_for_plot_scale (double max_value, double size_factor)
{
short int first_digit = ceil(max_value * size_factor);
if (first_digit == 7 || first_digit == 9) first_digit++;
return (first_digit);
}

short int Number_of_tickmarks_for_plot_scale (short int first_digit)
{
if (first_digit == 1 || first_digit == 6 || first_digit == 10) return (2);
if (first_digit == 4 || first_digit == 8) return (4);
return (first_digit);    
}

/* SSUBROUTINE THAT GENERATES AN XY-PLOT OF KMERS */
XYplot(char *filename, long int *****results, short int kmer_length, short int count_also_spaced_kmers, double kmer_length_difference_cutoff, char *xname, char *yname, char *kmerstring)
{

    /* OPENS SVG */
    FILE *outfile;
    outfile = fopen (filename, "w");
    
    fprintf(outfile, "<?xml version=\"1.0\" standalone=\"no\"?>\n<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n");
    fprintf(outfile, "<!--%s : command %s -->\n", svgsafe(VERSION), svgsafe(COMMAND));
    fprintf(outfile, "<svg width=\"2500\" height=\"5000\" x=\"0\" y=\"%i\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n", 0);
    
    Add_nucleotide_paths(outfile); /* ADDS NUCLEOTIDE PATHS */
    
    char *forward;
    if (rna == 0) forward = dnaforward;
    else forward = rnaforward;
    
    short int current_kmer_length = 0;
    long int current_kmer;
    short int current_gap_width;
    short int current_gap_position;
    long int current_pair;
    short int position;
    double max_x_value = 0;
    double max_y_value = 0;
    long int last_kmer = pow(4, kmer_length);
    long int *background_count;
    long int *signal_count;
    long int start_plotting_from_kmer_ranked;
    background_count = results[0][kmer_length][0][0];
    signal_count = results[1][kmer_length][0][0];  
    
    /* TRANSFERS COUNTS TO COUNTPAIR STRUCTURE */
    struct countpair *point;
    point = malloc (sizeof (struct countpair) * last_kmer * ((1 + (count_also_spaced_kmers != 0)) * (1 + (kmer_length % 2 == 1)) * (Nlength - kmer_length) + (count_also_spaced_kmers == 2) * (Nlength - kmer_length) * (kmer_length - 1)) + 5);
    
    for(current_pair = 0, current_gap_position = 0; current_gap_position <= (count_also_spaced_kmers == 1) * (kmer_length / 2 + kmer_length % 2) + (count_also_spaced_kmers == 2) * (kmer_length-1); current_gap_position++)
    {
    if (current_gap_position > 0 && (count_also_spaced_kmers == 1 && current_gap_position < (kmer_length / 2))) continue;
    for (current_gap_width = 0 + (current_gap_position != 0); current_gap_width <= (count_also_spaced_kmers != 0) * (Nlength - kmer_length - 2); current_gap_width++)
    {
    for (current_kmer = 0; current_kmer < last_kmer; current_kmer++, current_pair++) 
    {
    point[current_pair].add_logo = 0;
    point[current_pair].sequence_value = current_kmer;
    point[current_pair].gap_position = current_gap_position;
    point[current_pair].gap_width = current_gap_width;
    point[current_pair].x_count = results[0][kmer_length][current_gap_position][current_gap_width][current_kmer];
    point[current_pair].y_count = results[1][kmer_length][current_gap_position][current_gap_width][current_kmer];
    point[current_pair].sum = point[current_pair].x_count + point[current_pair].y_count;
    point[current_pair].x_local_max = Localmax(results, 0, kmer_length, kmer_length+1, kmer_length, point[current_pair].gap_position, point[current_pair].gap_width, point[current_pair].sequence_value, count_also_spaced_kmers, kmer_length_difference_cutoff, kmerstring);
    point[current_pair].y_local_max = Localmax(results, 1, kmer_length, kmer_length+1, kmer_length, point[current_pair].gap_position, point[current_pair].gap_width, point[current_pair].sequence_value, count_also_spaced_kmers, kmer_length_difference_cutoff, kmerstring);    
    }
    // printf("\nGap_position, Gap_length %i, %i", current_gap_position, current_gap_width); fflush(stdout);
    if(current_gap_position == 0) break;
    }
    }
    last_kmer = current_pair;
    qsort (point, last_kmer, sizeof(struct countpair), Sort_according_to_count); 
    
    /* STOPS WRITING A VERY LARGE FILE (MAX 65536 points) */   
    if (last_kmer > 65535) start_plotting_from_kmer_ranked = 65535;
    else start_plotting_from_kmer_ranked = last_kmer-1;
    
    double x_width = 500;
    double y_width = 500;
    double y_margin = 30;
    double x_margin = 100;
    double x_pos;
    double y_pos;
    double glyph_size = 5;
    double axis_stroke_width = 0.75;
    short int number_of_tickmarks;
    double last_x_tickmark_position;
    double last_y_tickmark_position;
    double current_tickmark_position;
    double scale;
    short int counter;
    short int first_x_digit;
    short int first_y_digit;
    short int x_local_max;
    short int y_local_max;
    short int current_tickmark;
    double tickmark_stroke_width = 0.75;
    double tickmark_length = 10;
    short int contains_CG;
    short int fill_glyph;
    char *glyphcolors[] = {"white", "black", "red"};
    short int localmax_added = 0;
    short int max_logos = 12;
    struct plot_point p;
    short int logos_added = 0;
    short int add_logo = 0;
    short int is_x = 0;
    short int localmax_point;
    short int start_new_pair = 2;
    char *x_or_null = "";
    char *y_or_null = "";
    double previous_x;
    double previous_y;
    double top_left_y;
    double top_left_x;
    short int boxes = 0;
    short int starting_logo = 0;
    
    fprintf(outfile, "<g id=\"xy_plot\">");
    
    /* FINDS MAX VALUES */
    for (current_kmer = 0; current_kmer < last_kmer; current_kmer++)
    {
    if (point[current_kmer].x_count > max_x_value) max_x_value = point[current_kmer].x_count;
    if (point[current_kmer].y_count > max_y_value) max_y_value = point[current_kmer].y_count;    
    }


    /* DRAWS X-AXIS */
    scale = pow(10.0, ceil(-log10(fabs(max_x_value))));
    first_x_digit = First_digit_for_plot_scale(max_x_value, scale);
    last_x_tickmark_position = first_x_digit / scale;
    // printf("\nX_maximum, end of scale = %.2f, %.2f", max_x_value, last_x_tickmark_position);
    fprintf(outfile, "<g> <title>%s</title> <line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" style=\"stroke:black;stroke-width:%.2f\"/> <text x=\"%.1f\" y=\"%.1f\" fill =\"black\" font-size = \"24\" text-anchor = \"middle\" font-style = \"italic\"> <tspan>%s (count x 10</tspan><tspan dy=\"-9\" dx=\"+2\" font-size = \"18\">%i</tspan><tspan dy=\"+9\" font-size = \"24\">)</tspan></text> </g>\n", xname, x_margin, y_width + y_margin, x_margin + x_width, y_width + y_margin, axis_stroke_width, x_margin + x_width / 2, y_width + y_margin + tickmark_length + 50, xname, 0 - (int) ceil(-log10(fabs(max_x_value))));

    /* DRAWS Y-AXIS */
    scale = pow(10.0, ceil(-log10(fabs(max_y_value))));
    first_y_digit = First_digit_for_plot_scale(max_y_value, scale);
    last_y_tickmark_position = first_y_digit / scale;
    // printf("\nY_maximum, end of scale = %.2f, %.2f", max_y_value, last_y_tickmark_position);
    fprintf(outfile, "<g> <title>%s</title> <line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" style=\"stroke:black;stroke-width:%.2f\"/> <text transform=\"rotate(-90, %.2f, %.2f)\" x=\"%.1f\" y=\"%.1f\" fill =\"black\" font-size = \"24\" text-anchor = \"middle\" font-style = \"italic\"> <tspan>%s (count x 10</tspan><tspan dy=\"-9\" dx=\"+2\" font-size = \"18\">%i</tspan><tspan dy=\"+9\" font-size = \"24\">)</tspan></text> </g>\n", yname, x_margin, y_width + y_margin, x_margin, y_margin, axis_stroke_width, x_margin - tickmark_length - 35, y_width / 2 + y_margin, x_margin - tickmark_length - 35, y_width / 2 + y_margin, yname, 0 - (int) ceil(-log10(fabs(max_y_value))));
    
    /* DRAWS Y-AXIS TICKMARKS */
    number_of_tickmarks = Number_of_tickmarks_for_plot_scale(first_y_digit);
    // printf("\nY_number of tickmarks = %i", number_of_tickmarks);
    for (current_tickmark = number_of_tickmarks; current_tickmark >= 0; current_tickmark--)
    {
    current_tickmark_position = y_margin + y_width * (last_y_tickmark_position - current_tickmark * last_y_tickmark_position / number_of_tickmarks) / last_y_tickmark_position;
    // printf("\nY_current tickmark %i position = %.0f", current_tickmark, current_tickmark_position);
    fprintf(outfile, "<g> <title>Y_tickmark_%i_value_%.0f</title> <line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" style=\"stroke:black;stroke-width:%.2f\"/> <text x=\"%.1f\" y=\"%.1f\" font-size = \"24\" fill =\"black\" >%i</text> </g>\n", current_tickmark, current_tickmark * last_y_tickmark_position / number_of_tickmarks, x_margin, current_tickmark_position, x_margin - tickmark_length, current_tickmark_position, tickmark_stroke_width, x_margin - tickmark_length - 15 - (current_tickmark == 2 && first_y_digit == 10) * 8, current_tickmark_position + 7, (int) current_tickmark * first_y_digit / number_of_tickmarks);   
    }
        
    /* DRAWS X-AXIS TICKMARKS */
    number_of_tickmarks = Number_of_tickmarks_for_plot_scale(first_x_digit);
    // printf("\nX_number of tickmarks = %i", number_of_tickmarks);
    for (current_tickmark = number_of_tickmarks; current_tickmark >= 0; current_tickmark--)
    {
    current_tickmark_position = x_margin + x_width * (last_x_tickmark_position - current_tickmark * last_x_tickmark_position / number_of_tickmarks) / last_x_tickmark_position;
    // printf("\nX_current tickmark %i position = %.0f", current_tickmark, current_tickmark_position);
    fprintf(outfile, "<g> <title>X_tickmark_%i_value_%.0f</title> <line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" style=\"stroke:black;stroke-width:%.2f\"/> <text x=\"%.1f\" y=\"%.1f\" font-size = \"24\" fill =\"black\" >%i</text> </g>\n", number_of_tickmarks - current_tickmark, (number_of_tickmarks - current_tickmark) * last_x_tickmark_position / number_of_tickmarks, current_tickmark_position, y_margin + y_width, current_tickmark_position, y_margin + y_width + tickmark_length, tickmark_stroke_width, current_tickmark_position - 6, y_margin + y_width + tickmark_length + 21, first_x_digit - (int) current_tickmark * first_x_digit / number_of_tickmarks);   
    }
    
    /* ADDS RESPONSIVE YELLOW MARKS TO BEHIND TOP LOCAL MAX POINTS */
    for (current_kmer = 0; current_kmer < last_kmer; current_kmer++)
    {
        add_logo = 0;
        if (point[current_kmer].y_local_max == 1) {add_logo = 1; is_x = 0;}
        if (point[current_kmer].x_local_max == 1) {add_logo = 1; is_x = 1;}
        if (add_logo == 1 && localmax_added < max_logos)
        {
            x_pos = x_margin + x_width * (double) point[current_kmer].x_count / last_x_tickmark_position;
            y_pos = y_margin + y_width - y_width * (double) point[current_kmer].y_count / last_y_tickmark_position;
            fprintf(outfile, "<circle id=\"point_outline_%li\" cx=\"%.2f\" cy=\"%.2f\" r=\"%.2f\" style=\"fill:yellow;opacity:0\"> <set attributeName=\"opacity\" to=\"0.5\" begin=\"logo_set_%li.mouseenter\" end=\"logo_set_%li.mouseleave\"/> </circle>", point[current_kmer].sum, x_pos, y_pos, glyph_size * 3, point[current_kmer].sequence_value, point[current_kmer].sequence_value);
            localmax_added++;
            point[current_kmer].add_logo = 1;
        }
    }
    
    /* DRAWS CIRCLES */
    for (localmax_point = localmax_added, current_kmer = start_plotting_from_kmer_ranked; current_kmer >= 0; current_kmer--)
    {        
        /* CHECKS FOR LOCAL MAX */
        y_local_max = 0;
        x_local_max = 0;
        if (point[current_kmer].y_local_max == 1) {y_local_max = 1; fill_glyph = 1;} else fill_glyph = 0;
        if (point[current_kmer].x_local_max == 1) {x_local_max = 1; fill_glyph = 1;}
        
        /* CHECKS FOR CG */
        for(contains_CG = 0, position = kmer_length-2; position > -1 ; position--) if (((point[current_kmer].sequence_value & mask_ULL[2][position]) >> (position * 2)) == 6 && kmer_length - position != point[current_kmer].gap_position) contains_CG = 1;
        x_pos = x_margin + x_width * (double) point[current_kmer].x_count / last_x_tickmark_position;
        y_pos = y_margin + y_width - y_width * (double) point[current_kmer].y_count / last_y_tickmark_position;
        if (point[current_kmer].add_logo == 1) 
        {
        fprintf(outfile, "<g id=\"maxpoint%i\">", localmax_point);
        localmax_point--;
        }
        else fprintf(outfile, "<g>");
        fprintf(outfile, "<title>");
        for(position = kmer_length-1; position > -1 ; position--) 
        {
        if(kmer_length - position - 1 == point[current_kmer].gap_position) for(counter = 0; counter < point[current_kmer].gap_width; counter++) fprintf(outfile, "n");
        fprintf(outfile, "%c", forward[(point[current_kmer].sequence_value & mask_ULL[1][position]) >> (position * 2)]);
        }
        fprintf(outfile, "</title> <circle cx=\"%.2f\" cy=\"%.2f\" r=\"%.2f\" style=\"fill:%s;stroke-width:%.2f;stroke:%s\"/> </g>\n", x_pos, y_pos, glyph_size, glyphcolors[(1+contains_CG) * fill_glyph], 0.75, glyphcolors[1+contains_CG]);
        /* DRAWS LINES FOR LOCAL MAX */
        if (x_local_max == 1) fprintf(outfile, "<g> <title>%s</title> <line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" style=\"stroke-width:%.2f;stroke:%s\"/> </g>\n", "x_local_max", x_pos, y_pos, x_pos + 2 * glyph_size, y_pos, 0.5, glyphcolors[1+contains_CG]); 
        if (y_local_max == 1) fprintf(outfile, "<g> <title>%s</title> <line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" style=\"stroke-width:%.2f;stroke:%s\"/> </g>\n", "y_local_max", x_pos, y_pos, x_pos, y_pos - 2 * glyph_size, 0.5, glyphcolors[1+contains_CG]); 
    }
    
    /* PRINTS BEST KMER LOGOS */
    fprintf(outfile, "</g><g id=\"kmerlogo_comparison\">");

    /* ADDS KMER LOGOS */
    p.y = y_margin;
    for (current_kmer = 0; current_kmer < last_kmer; current_kmer++)
    {
    add_logo = 0;
    if (point[current_kmer].y_local_max == 1) {add_logo = 1; is_x = 0;}
    if (point[current_kmer].x_local_max == 1) {add_logo = 1; is_x = 1;}

    /* if(add_logo == 1 || current_kmer < 20)
    {
        printf("\nLogo %li, sum %li\t", current_kmer, point[current_kmer].sum);
        if (add_logo == 1) printf ("Local Max\t"); else printf ("\t");
        for(position = kmer_length-1; position > -1 ; position--) 
        {
            if(kmer_length - position - 1 == point[current_kmer].gap_position) for(counter = 0; counter < point[current_kmer].gap_width; counter++) printf("n");
            printf("%c", forward[(point[current_kmer].sequence_value & mask_ULL[1][position]) >> (position * 2)]);
        }
    } */
        
    if(add_logo == 1)
    {
    
    if (start_new_pair == 2)
    {
    p.x = x_margin + x_width + 20;
    previous_y = p.y;
    top_left_y = p.y;
    top_left_x = p.x;
    starting_logo = logos_added;
    /* DRAWS HIGHLIGHT BOXES */
    for (boxes = 0; boxes <= 0 + (Is_this_sequence_dimer (point[current_kmer].sequence_value, kmer_length) != 4); boxes++)
    {
    fprintf(outfile, "<g id=\"highlight%i\"> <rect x=\"%.2f\" y=\"%.2f\" width=\"%.2f\" height=\"%.2f\" style=\"fill:yellow;opacity:0\"> <set attributeName=\"opacity\" to=\"0.3\" begin=\"maxpoint%i.mouseenter\" end=\"maxpoint%i.mouseleave\"/> </rect> </g> \n", starting_logo + boxes, top_left_x, top_left_y, 500.0, 56.0, starting_logo +1 + boxes, starting_logo +1 + boxes);
    }
    if (boxes == 1) {x_or_null = "x"; y_or_null = "y";}
    fprintf(outfile, "<g id=\"logo_set_%li\">", point[current_kmer].sequence_value);
    }
    previous_x = p.x;
    Add_kmercount_logo_to_Svg (outfile, &p, point[current_kmer].sequence_value, results, 0, kmer_length, point[current_kmer].gap_position, point[current_kmer].gap_width, point[current_kmer].x_local_max, x_or_null, logos_added, count_also_spaced_kmers);
        add_logo = 0;
    p.y += 2;
    p.x = previous_x;
    Add_kmercount_logo_to_Svg (outfile, &p, point[current_kmer].sequence_value, results, 1, kmer_length, point[current_kmer].gap_position, point[current_kmer].gap_width, point[current_kmer].y_local_max, y_or_null, logos_added, count_also_spaced_kmers);
        add_logo = 0;
    boxes++;    
    start_new_pair--;    
    if (Is_this_sequence_dimer (point[current_kmer].sequence_value, kmer_length) == 4) {p.y += 16; start_new_pair = 2;}
    if (start_new_pair == 0) {p.y += 16; start_new_pair = 2;}    
    logos_added++;
        
    if (start_new_pair == 2) 
    {
    fprintf(outfile, "</g>"); x_or_null = ""; y_or_null = "";
    }
    else 
    {
    p.x += 30; p.y = previous_y; x_or_null = "x"; y_or_null = "y";
    }
        
    if (logos_added >= max_logos) break;
    }
    }
    if (start_new_pair != 2) fprintf(outfile, "</g>");
    fprintf(outfile, "</g>");
    fprintf(outfile, "</svg>");
    fclose (outfile);
    free(point);
}

/* GENERATES AN SVG SUMMARY FILE FOR KMER COUNTS */
Kmer_svg(char *input_filename, long int *****results, short int file_number, short int shortest_kmer, short int too_long_kmer, short int count_also_spaced_kmers, double kmer_length_difference_cutoff, long int count_cutoff, long int background_counts, long int signal_counts, short int number_of_heatmap_rows, char *kmerstring, long int *preferred_kmers, short int preferred_length)
{
    char *font = "Courier";
    char *heatmap_name;
    heatmap_name = malloc(1000);
    char *tsv_filename;
    tsv_filename = malloc(1000);
    char *filename;
    filename = malloc(1000);
    strcpy(filename, input_filename);
    short int font_size = 28;
    double *peak_features;
    peak_features = malloc(sizeof(double)*3+5);
    double shape;
    char **arrow_color;
    arrow_color = malloc(sizeof(char *) * 3 + 5);
    arrow_color[0] = malloc(100);
    arrow_color[1] = malloc(100);
    strcpy(arrow_color[0], "black");
    strcpy(arrow_color[1], "red");
    char call = '-';
    
    float kmer_pseudocount;
    short int current_kmer_length = 0;
    short int current_gap_length = 0;
    short int current_gap_position = 0;
    long int current_kmer;
    long int previous_kmer;
    long int current_kmer_of_this_size; 
    if (filename[0] == '\0') sprintf(filename, "Kmer_summary%ito%i.svg", shortest_kmer, too_long_kmer-1);
    strcpy(tsv_filename, filename);
    strcat(tsv_filename, ".txt");
    
    long int number_of_kmers;
    if(shortest_kmer+1==too_long_kmer) number_of_kmers = pow(4, shortest_kmer);
    else number_of_kmers = pow(4, too_long_kmer);
    struct kmer_incidence_table *top_kmers;
    top_kmers = malloc(sizeof(struct kmer_incidence_table) * number_of_kmers +5);
    
    long int x_position = 0;
    short int current_heatmap_row;
    long int length_correction;
    if (too_long_kmer == 3) number_of_heatmap_rows = 16;
    // else number_of_heatmap_rows = 20;
    short int tile_height = 30;
    short int tile_width = 30;
    short int max_bar_width = 50;
    long int font_position = 0;
    double pwm_position_width = 10;
    double pwm_position_height = 25;
    short int local_maxes;
    short int any_local_maxes;
    long int local_max_max_score;
    float max_enrichment;
    float local_max_enrichment;
    double current_enrichment;
    double min_enrichment;
    long int y_position = tile_height;
    long int counter;
    signed short int position;
    short int red = 0;
    short int green = 255;
    short int blue = 0;
    short int blue_border = 0;
    short int yellow_border = 0;
    short int stroke = 1;
    short int nucleotide_value;
    short int heatmap;
    short int number_of_heatmaps = 7;
    long int mutant_kmer;
    double total_count;
    double swap;
    double **order;
    order = malloc(sizeof(double *) * 4 + 5);
    for (counter = 0; counter < 3; counter++) order[counter] = malloc(sizeof(double) * 6 + 5);
    
    double black_score;
    double max_score;
    double min_score;
    long int total_score;
    long int max_total_score = 0;
    long int current_count;
    long int current_background_count;
    float size_factor = (float) background_counts / (float) signal_counts;

    current_gap_length = 0;
    short int true_gap_position;
    short int first_shift = 0;
    short int second_shift = 0;
    
    short int mid1 = 3;
    short int mid2 = 3;
    short int ac1 = 0;
    short int ac2 = 1;
    short int yshift = 10;
    short int xshift = 20;
    
    char *forward;
    if (rna == 0) forward = dnaforward;
    else forward = rnaforward;
    
    /* OPENS SVG */
    FILE *outfile;
    outfile = fopen (filename, "w");
    FILE *tsv_outfile;
    tsv_outfile = fopen (tsv_filename, "w");
    
    fprintf(outfile, "<?xml version=\"1.0\" standalone=\"no\"?>\n<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n");
    fprintf(outfile, "<!--%s : command %s -->\n", svgsafe(VERSION), svgsafe(COMMAND));
    fprintf(outfile, "<svg width=\"2500\" height=\"5000\" x=\"0\" y=\"%i\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n", 0);
    
    Add_nucleotide_paths(outfile); /* Adds nucleotide paths */

    
    /* GENERATES KMER TABLE FOR IDENTIFICATION OF TOP KMERS (SUM OF ALL SPACINGS)*/
    for(current_kmer = 0, length_correction = 1, current_kmer_length = too_long_kmer-1; current_kmer_length >= shortest_kmer; current_kmer_length--, length_correction *= 4)
    {
    

        
    number_of_kmers = pow(4, current_kmer_length);
    true_gap_position = current_kmer_length / 2;
    kmer_pseudocount = 10;
    previous_kmer = current_kmer;
        
    for(current_kmer_of_this_size = 0; current_kmer_of_this_size < number_of_kmers; current_kmer_of_this_size++, current_kmer++)
    {
    (top_kmers[current_kmer]).preferred = 0;
    (top_kmers[current_kmer]).max_gap_length = 0;
    (top_kmers[current_kmer]).local_max_max_gap_length = -1;
    for(total_score = 0, any_local_maxes = 0, local_maxes = 0, max_score = 0, max_enrichment = 0, local_max_max_score = 0, local_max_enrichment = 0,current_gap_position = 0, current_gap_length = 0; current_gap_length < Nlength - current_kmer_length; current_gap_length++, current_gap_position = true_gap_position)
    {
        current_count = results[file_number][current_kmer_length][current_gap_position][current_gap_length][current_kmer_of_this_size];
        current_background_count = results[file_number-1][current_kmer_length][current_gap_position][current_gap_length][current_kmer_of_this_size];
        total_score += current_count;
        if (((float) current_count * size_factor + kmer_pseudocount) / ((float) current_background_count + kmer_pseudocount) > max_enrichment)
        {
            max_enrichment = ((float) current_count * size_factor + kmer_pseudocount) / ((float) current_background_count + kmer_pseudocount);
        }
        if (current_count / length_correction > max_score) 
        {
            max_score = current_count / length_correction;
            (top_kmers[current_kmer]).max_gap_length = current_gap_length;  
        }
        if(Localmax(results, file_number, shortest_kmer, too_long_kmer, current_kmer_length, current_gap_position, current_gap_length, current_kmer_of_this_size, count_also_spaced_kmers, kmer_length_difference_cutoff, kmerstring)) 
        {
        if (current_count > 0) 
        {
        any_local_maxes++;
        if (current_count > local_max_max_score) 
        {
        local_max_max_score = current_count;
        (top_kmers[current_kmer]).local_max_max_gap_length = current_gap_length;
        }
        if (((float) current_count * size_factor + kmer_pseudocount) / ((float) current_background_count + kmer_pseudocount) > local_max_enrichment)
        {
        local_max_enrichment = ((float) current_count * size_factor + kmer_pseudocount) / ((float) current_background_count + kmer_pseudocount);
        }
        }
        if (current_count > count_cutoff) 
        {
        local_maxes++; 
        }
        }
    }
        total_score /= length_correction;
        (top_kmers[current_kmer]).incidence = total_score; 
        if (total_score > max_total_score) max_total_score = total_score;
        (top_kmers[current_kmer]).max_enrichment = max_enrichment; 
        (top_kmers[current_kmer]).local_max_enrichment = local_max_enrichment; 
        (top_kmers[current_kmer]).local_max_max_incidence = local_max_max_score / length_correction;
        (top_kmers[current_kmer]).kmer = current_kmer_of_this_size;
        (top_kmers[current_kmer]).kmer_length = current_kmer_length;
        (top_kmers[current_kmer]).local_max = local_maxes;
        (top_kmers[current_kmer]).any_local_max = any_local_maxes;
    }
    
    /* CHECKS IF KMER IS ON THE PREFERRED LIST */
    if (current_kmer_length == preferred_length) for (counter = 0; counter < 10; counter++) (top_kmers[previous_kmer+preferred_kmers[counter]]).preferred = 10-counter;
        
    }
    
    number_of_kmers = current_kmer - 1;

    /* FINDS POSITIONS ON HEXBIN GRAPH 
    struct hex_tile {signed short int red; signed short int green; short int occupied; long int kmer_rank; long int kmer_sequence_value; short int kmer_length;};
    short int ybins = 100;
    short int xbins = 50;
    long int bins = xbins * ybins;
    struct hex_tile plot1[xbins+1][ybins+1];
    struct plot_kmer_position {long int kmer_sequence_value; short int kmer_length; long int kmer_rank; short int x; short int y;};
    struct plot_kmer_position plot1_kmer_position[bins];
    long int ybin;
    long int xbin;
    long int kmer_number;
    long int compared_to_kmer_number;
    double total_distance;
    double min_total_distance;
    short int x_distance;
    short int y_distance;
    double hex_distance;
    short int kmer_limit = 10;
    double ykat = 5;
    double xkat = 2.258;
    double hyp = pow(pow(ykat,2)+pow(xkat,2),0.5);    
    long int tile_value;
    short int edit_distance;
    struct bitiupac_structure seq1;
    bitiupac_structure_init(&seq1);
    struct bitiupac_structure seq2;
    bitiupac_structure_init(&seq2);
    struct alignment a;
    
    black_score = 2 * signal_counts * (Nlength - current_kmer_length) / number_of_kmers;
    
    for (ybin = 0; ybin < ybins; ybin++) for (xbin = 0; xbin < xbins; xbin++) 
    {
        (plot1[xbin][ybin]).red = 0;
        (plot1[xbin][ybin]).green = 0;
        (plot1[xbin][ybin]).occupied = 0;
        (plot1[xbin][ybin]).kmer_rank = 1000;
    }
    
    /* SETS FIRST KMER TO MIDDLE 
    (plot1_kmer_position[0]).x = xbins/2;
    (plot1_kmer_position[0]).y = ybins/2;
    (plot1_kmer_position[0]).kmer_sequence_value = (top_kmers[0]).kmer;
    (plot1_kmer_position[0]).kmer_length = (top_kmers[0]).kmer_length;
    (plot1[(plot1_kmer_position[0]).x][(plot1_kmer_position[0]).y]).red = 255;
    // ((top_kmers[0]).incidence / max_total_score) * 255;
    (plot1[(plot1_kmer_position[0]).x][(plot1_kmer_position[0]).y]).green = 0;
    (plot1[(plot1_kmer_position[0]).x][(plot1_kmer_position[0]).y]).kmer_rank = 0;
    (plot1[(plot1_kmer_position[0]).x][(plot1_kmer_position[0]).y]).occupied = 1;
    
    for(kmer_number = 1; kmer_number < bins; kmer_number++)
    {
    //if((top_kmers[kmer_number]).kmer_length != 6 | (top_kmers[kmer_number]).max_gap_length != 0) continue;
    
    min_total_distance = 10000;
    for (ybin = 0; ybin < ybins; ybin++) for (xbin = 0; xbin < xbins; xbin++) 
    {
    if((plot1[xbin][ybin]).occupied == 1) continue;
    for (total_distance = 0, compared_to_kmer_number = 0; compared_to_kmer_number < kmer_number; compared_to_kmer_number++)
    {
        if (compared_to_kmer_number > kmer_limit) break;
        x_distance = 2 * (xbin - (plot1_kmer_position[compared_to_kmer_number]).x);
        y_distance = ybin - (plot1_kmer_position[compared_to_kmer_number]).y;
        // if (x_distance > y_distance) hex_distance = x_distance;
        // else hex_distance = y_distance;
        hex_distance = pow(pow(y_distance * 2 * ykat, 2) + pow(x_distance * (xkat*2 + hyp), 2), 0.5);
        
        Bitseq_to_bitiupac((top_kmers[compared_to_kmer_number]).kmer, (top_kmers[compared_to_kmer_number]).kmer_length, (top_kmers[compared_to_kmer_number]).local_max_max_gap_length, ((top_kmers[compared_to_kmer_number]).local_max_max_gap_length == 0)*(top_kmers[compared_to_kmer_number]).kmer_length/2, &seq1);
        Bitseq_to_bitiupac((top_kmers[kmer_number]).kmer, (top_kmers[kmer_number]).kmer_length, (top_kmers[kmer_number]).local_max_max_gap_length, ((top_kmers[kmer_number]).local_max_max_gap_length == 0)*(top_kmers[kmer_number]).kmer_length/2, &seq2);       
        
        edit_distance = Align(&a, &seq1, &seq2);
        if (edit_distance == 0) continue;
        total_distance += abs (hex_distance / (6 * hyp) - edit_distance) / (log(edit_distance));
    }
        if (total_distance < min_total_distance)
            {
                min_total_distance = total_distance;
                (plot1_kmer_position[kmer_number]).x = xbin;
                (plot1_kmer_position[kmer_number]).y = ybin;
                (plot1_kmer_position[kmer_number]).kmer_rank = kmer_number;
                (plot1_kmer_position[kmer_number]).kmer_sequence_value = (top_kmers[kmer_number]).kmer;
                (plot1_kmer_position[kmer_number]).kmer_length = (top_kmers[kmer_number]).kmer_length;
            }
    }
        (plot1[(plot1_kmer_position[kmer_number]).x][(plot1_kmer_position[kmer_number]).y]).occupied = 1;
        
        tile_value = (top_kmers[kmer_number]).incidence;
        if (tile_value > black_score) (plot1[(plot1_kmer_position[kmer_number]).x][(plot1_kmer_position[kmer_number]).y]).red = 255 * (tile_value - black_score) / (max_total_score-black_score);
        else (plot1[(plot1_kmer_position[kmer_number]).x][(plot1_kmer_position[kmer_number]).y]).green = 255 * (black_score - tile_value) / (0-black_score);
        (plot1[(plot1_kmer_position[kmer_number]).x][(plot1_kmer_position[kmer_number]).y]).kmer_rank = (plot1_kmer_position[kmer_number]).kmer_rank;
        
    }

    
    /* PLOTS HEXBIN GRAPH 
    
    for (ybin = 0; ybin < ybins; ybin++, y_position += ykat)  
    {
        for (x_position = 0 + (ybin % 2) * (xkat+hyp), xbin = 0; xbin < xbins; xbin++, x_position += xkat+hyp+hyp+xkat)
        {
            fprintf(outfile, "<polygon points=\"%.1f,%.1f %.1f,%.1f %.1f,%.1f %.1f,%.1f %.1f,%.1f %.1f,%.1f\" style=\"fill:rgb(%i,%i,%i);stroke:black;stroke-width:0.3\"/>", 
                    x_position+xkat, (double) y_position, 
                    x_position+xkat+hyp, (double) y_position, 
                    x_position+xkat+hyp+xkat, y_position+ykat, 
                    x_position+xkat+hyp, y_position+ykat+ykat,
                    x_position+xkat, y_position+ykat+ykat,
                    (double) x_position, y_position+ykat,
                    (plot1[xbin][ybin]).red, (plot1[xbin][ybin]).green, 0);
            if ((plot1[xbin][ybin]).kmer_rank < 9) fprintf(outfile, "<text x=\"%li\" y=\"%li\" fill = \"yellow\" stroke = \"yellow\" font-size=\"%i\" font-family = \"%s\" >%li</text>\n", x_position+1, y_position+8, 10, font, (plot1[xbin][ybin]).kmer_rank + 1);
        }
    }
    */
    
    /* SORTS ACCORDING TO PREFERRED */
    if (preferred_length != 0) qsort (top_kmers, number_of_kmers, sizeof(struct kmer_incidence_table), Sort_according_to_preferred);
    
    x_position = 0;
    for (heatmap = 0; heatmap < number_of_heatmaps; heatmap++)
    {
        
    if (heatmap != 0 | preferred_length != 0) for(current_heatmap_row = 0; current_heatmap_row < number_of_heatmap_rows; current_heatmap_row++)
    {
        current_kmer = (top_kmers[current_heatmap_row]).kmer;
        current_kmer_length = (top_kmers[current_heatmap_row]).kmer_length;
        true_gap_position = current_kmer_length / 2;
        length_correction = pow(4,(too_long_kmer-current_kmer_length-1));
        
        /* NORMALIZES HEATMAP ROW */
        for(total_score = 0, current_gap_position = 0, current_gap_length = 0, min_score = 1E10, max_score = -1E10, min_enrichment = 10000000; current_gap_length < Nlength - current_kmer_length; current_gap_length++, current_gap_position = true_gap_position)
        {
            current_count = results[file_number][current_kmer_length][current_gap_position][current_gap_length][current_kmer];
            current_background_count = results[file_number-1][current_kmer_length][current_gap_position][current_gap_length][current_kmer];
            current_enrichment = (double) ((float) current_count * size_factor + kmer_pseudocount) / ((float) current_background_count + kmer_pseudocount);
            if (current_enrichment < min_enrichment) min_enrichment = current_enrichment;
            if (current_count > max_score) max_score = current_count;
            if (current_count < min_score) min_score = current_count;
            total_score += current_count;
        }
        total_score /= length_correction;
        max_score /= length_correction;
        if(max_score == 0) max_score = 0.001;
  
        x_position = 5;
        /* PRINTS OUT LABEL */
        if (current_heatmap_row == 0)
        {
        fprintf(outfile, "<text  x=\"");
        fprintf(outfile, "%li", x_position);
        fprintf(outfile, "\" y=\"");
        fprintf(outfile, "%li", y_position + tile_height - 10);
        fprintf(outfile, "\" fill = \"black\" stroke = \"black\" font-size=\"%i\" font-family = \"", font_size);
        fprintf(outfile, "%s", font);
        fprintf(outfile, "\" >");
        if(heatmap == 0) strcpy(heatmap_name, "Halfsite_spacing");
        if(heatmap == 1) strcpy(heatmap_name, "Top_local_max_count");
        if(heatmap == 2) strcpy(heatmap_name, "Local_max_enrichment");
        if(heatmap == 3) strcpy(heatmap_name, "Enrichment");
        if(heatmap == 4) strcpy(heatmap_name, "Local_max_containing");
        if(heatmap == 5) strcpy(heatmap_name, "Incidence_ranked");
        if(heatmap == 6) strcpy(heatmap_name, "Local_max_ranked");
        fprintf(outfile, "%s</text>\n", heatmap_name);     
        }
        
        y_position += tile_height;
        
        /* PRINTS OUT KMER */
        fprintf(outfile, "<text  x=\"");
        fprintf(outfile, "%li", x_position);
        fprintf(outfile, "\" y=\"");
        fprintf(outfile, "%li", y_position + tile_height - 7);
        fprintf(outfile, "\" fill = \"black\" stroke = \"none\" font-size=\"%i\" font-family = \"", font_size);
        fprintf(outfile, "%s", font);
        fprintf(outfile, "\" >");
        for(position = current_kmer_length-1; position > -1 ; position--) {if(current_kmer_length - position - 1  == current_gap_position)fprintf(outfile,"."); fprintf(outfile, "%c", forward[(current_kmer & mask_ULL[1][position]) >> (position * 2)]);}
        for(; position < too_long_kmer-current_kmer_length-1 ;position++) fprintf(outfile," "); /* FILLS SPACES IF KMER SHORTER THAN MAX LENGTH */
        fprintf(outfile, "</text>\n");

        
        /* PRINTS KMER LOGO FOR THE HIGHEST INCIDENCE LOCAL MAX GAP WIDTH (OR HIGHEST INCIDENCE GAP LENGTH IF NO LOCAL MAX) */
        x_position += (font_size-9) * (too_long_kmer) + 5;
        
        /* IF PREFERRED, PRINTS ORIENTATION */
        if (heatmap == 0)
        {
            if (current_heatmap_row == 2 || current_heatmap_row == 5 || current_heatmap_row == 8) mid1 = -3;  else mid1 = 3;
            if (current_heatmap_row == 1) {ac1 = 1; ac2 = 0;} else {ac1 = 0; ac2 = 1;}
            if (current_heatmap_row > 3 && current_heatmap_row < 7) {ac1 = 0; ac2 = 0;}
            if (current_heatmap_row > 6) {ac1 = 1; ac2 = 1;}
            fprintf(outfile, "<polyline points =\"");
            fprintf(outfile, "%li,%li %li,%li %li,%li \" ", x_position-xshift-mid1,y_position+yshift, x_position-xshift+mid1,y_position+yshift+3, x_position-xshift-mid1,y_position+yshift+6);
            fprintf(outfile, "fill = \"%s\" stroke = \"%s\" stroke-width = \"2\"/>", arrow_color[ac1], arrow_color[ac1]);
            x_position += 10;
            if (current_heatmap_row == 3 || current_heatmap_row == 6 || current_heatmap_row == 9) mid2 = -3;  else mid2 = 3;
            fprintf(outfile, "<polyline points =\"");
            fprintf(outfile, "%li,%li %li,%li %li,%li \" ", x_position-xshift-mid2,y_position+yshift, x_position-xshift+mid2,y_position+yshift+3, x_position-xshift-mid2,y_position+yshift+6);
            fprintf(outfile, "fill = \"%s\" stroke = \"%s\" stroke-width = \"2\"/>", arrow_color[ac2], arrow_color[ac2]);
        }
        
        
        current_gap_length = (top_kmers[current_heatmap_row]).local_max_max_gap_length;
        if(current_gap_length < 0) current_gap_length = (top_kmers[current_heatmap_row]).max_gap_length;
        if (current_gap_length == 0) current_gap_position = 0;
        else current_gap_position = true_gap_position;
            
        for (position = current_kmer_length-1; position >= 0; position--)
        {
        /* IDENTIFIES THE FOUR KMER COUNTS FOR THIS POSITION */
            for(total_count = 0, nucleotide_value = 0; nucleotide_value < 4; nucleotide_value++)
            {
            mutant_kmer = ((current_kmer & (~(mask_ULL[1][position]))) | (nucleotide_value << (position * 2))) & lowmask_ULL[current_kmer_length-1];
            current_count = results[file_number][current_kmer_length][current_gap_position][current_gap_length][mutant_kmer];
            total_count += current_count;
            }
            for(nucleotide_value = 0; nucleotide_value < 4; nucleotide_value++)
            {
            mutant_kmer = ((current_kmer & (~(mask_ULL[1][position]))) | (nucleotide_value << (position * 2))) & lowmask_ULL[current_kmer_length-1];
            order[0][nucleotide_value] = nucleotide_value;
            order[1][nucleotide_value] = ((double) results[file_number][current_kmer_length][current_gap_position][current_gap_length][mutant_kmer]) / total_count;
            }
        
        /* BUBBLE SORT */
        for (counter = 0; counter < 3; counter++) 
        {
            for(nucleotide_value = counter; nucleotide_value < 4; nucleotide_value++)
            {
                if (order[1][counter] < order[1][nucleotide_value]) 
                {
                    swap = order[0][counter];
                    order[0][counter] = order[0][nucleotide_value];
                    order[0][nucleotide_value] = swap;
                    swap = order[1][counter];
                    order[1][counter] = order[1][nucleotide_value];
                    order[1][nucleotide_value] = swap;
                }
            }
        }
        for(font_position = y_position + 3.5, nucleotide_value = 0; nucleotide_value < 4; nucleotide_value++)
        {
        /* PRINTS OUT SCALED PATH NUCLEOTIDES */
        fprintf(outfile, "<use xlink:href=\"#%c\" ", forward[(int) order[0][nucleotide_value]]);
        fprintf(outfile, " transform=\"translate(%li,%f) scale(%f,%f)\" visibility=\"visible\" />\n", x_position, font_position + (order[1][nucleotide_value] * pwm_position_height), ((double) pwm_position_width) / 10, order[1][nucleotide_value] * ((double) pwm_position_height) / 10);
        font_position += (order[1][nucleotide_value] * pwm_position_height);
        }
        x_position += pwm_position_width;
        if(position == current_gap_position + (current_kmer_length % 2) * (current_gap_position != 0)) 
        {
        /* PRINTS OUT Ns */
        fprintf(outfile, "<text  x=\"");
        fprintf(outfile, "%li", x_position);
        fprintf(outfile, "\" y=\"");
        fprintf(outfile, "%li", y_position + tile_height - 10);
        fprintf(outfile, "\" fill = \"gray\" stroke = \"gray\" font-size=\"%i\" font-family = \"", font_size/2);
        fprintf(outfile, "%s", font);
        fprintf(outfile, "\" >");
        if(current_gap_length < 10) fprintf(outfile, "-");
        fprintf(outfile,"-%iN-", current_gap_length);
        fprintf(outfile, "</text>\n");
        x_position += (font_size/2-5.6) * 5;
        }
        }

        x_position += pwm_position_width * (too_long_kmer - current_kmer_length - 1); /* ADDS SPACE IF SHORTER THAN MAX KMER */
        
        /* PRINTS BAR */
        x_position += 5;
        blue = 250;
        red = 250;
        green = 250;
        fprintf(outfile, " <rect x=\"%li\" y=\"%li\" width=\"%li", x_position, y_position + tile_height/3, max_bar_width - (max_bar_width * total_score) / max_total_score);
        fprintf(outfile, "\" height=\"");
        fprintf(outfile, "%i", tile_height/3);
        fprintf(outfile, "\" style=\"");
        fprintf(outfile, "fill:rgb(");
        fprintf(outfile, "%i,%i,%i);stroke-width:%i;stroke:rgb(0,0,0)\"/>\n", red, green, blue, stroke);  
        x_position += max_bar_width - (max_bar_width * total_score) / max_total_score;
        
        blue = 256;
        red = 0;
        green = 0;
        fprintf(outfile, " <rect x=\"%li\" y=\"%li\" width=\"%li", x_position, y_position + tile_height/3, (max_bar_width * total_score) / max_total_score);
        fprintf(outfile, "\" height=\"");
        fprintf(outfile, "%i", tile_height/3);
        fprintf(outfile, "\" style=\"");
        fprintf(outfile, "fill:rgb(");
        fprintf(outfile, "%i,%i,%i);stroke-width:%i;stroke:rgb(0,0,0)\"/>\n", red, green, blue, stroke);       
        blue = 0; 
        x_position += (max_bar_width * total_score) / max_total_score;
        
        /* PRINTS SCALE */
        if (current_heatmap_row == 0) 
        {
        for(counter = 0; counter < Nlength - shortest_kmer; counter++) 
        {
            fprintf(outfile, " <text  x=\"%li\" y=\"%li\" fill = \"black\" stroke = \"black\" font-size=\"14\" font-family = \"%s\">%li</text>\n", x_position + counter * tile_width + tile_width / 2 - (counter > 9) * 5, y_position-3, font, counter);
        }
        fprintf(outfile, " <text  xml:space=\"preserve\" x=\"%li\" y=\"%li\" fill = \"black\" stroke = \"black\" font-size=\"28\" font-family = \"%s\"> C  D  F     S      Localmax  Fold      Max     Fold    Total</text>\n", x_position + counter * tile_width + tile_width / 2 - (counter > 9) * 5, y_position-10, font, counter);
        }
        
        x_position += 5;
        
        /* PRINTS GUIDE LINE */
        fprintf(outfile, "<polyline points =\"");
        fprintf(outfile, "%i,%li %li,%li \" ", 0,y_position,x_position, y_position);
        fprintf(outfile, "fill = \"none\" stroke = \"gray\" stroke-width = \"0.2\"/>");
        
        /* PRINTS ONE LINE OF HEATMAP TILES */
        for(current_gap_length = 0, current_gap_position = 0; current_gap_length < Nlength - current_kmer_length; current_gap_length++, current_gap_position = true_gap_position)
        {
            /* COUNT TILE */
            /* DETERMINES COLOR */
            black_score = (max_score+min_score) / 2;
            red = 0;
            green = 0;
            current_count = results[file_number][current_kmer_length][current_gap_position][current_gap_length][current_kmer];
            if (current_count > black_score) red = (current_count - black_score) * 256 / (max_score-black_score);
            else green = (black_score - current_count) * 256 / (black_score-min_score);
            
            /* PRINTS OUT TILE */
            fprintf(outfile, " <rect x=\"%li\" y=\"%li\" width=\"%i", x_position, y_position, tile_width);
            fprintf(outfile, "\" height=\"");
            fprintf(outfile, "%i", tile_height);
            fprintf(outfile, "\" style=\"");
            fprintf(outfile, "fill:rgb(");
            fprintf(outfile, "%i,%i,%i);stroke-width:%i;stroke:rgb(0,0,0)\"/>\n", red, green, blue, stroke);
            
            /* ENRICHMENT TILE */
            /* DETERMINES COLOR */
            red = 0;
            green = 0;
            black_score = 1;
            current_background_count = results[file_number-1][current_kmer_length][current_gap_position][current_gap_length][current_kmer];
            current_enrichment = (double) ((float) current_count * size_factor + kmer_pseudocount) / ((float) current_background_count + kmer_pseudocount);
            if (current_enrichment > black_score) red = (current_enrichment-1) * 256 / ((top_kmers[current_heatmap_row]).max_enrichment-1);
            else green = (1-current_enrichment) * 256 / (1-min_enrichment);
            
            /* PRINTS OUT TOP RIGHT TRIANGLE */
            fprintf(outfile, "<polygon points=\"%.0f,%.2f %.0f,%.2f %.0f,%.2f\" style=\"fill:rgb(%i,%i,%i);stroke:black;stroke-width:0.3\"/>", x_position+3* ((double) tile_width)/4, y_position+0.45, (double) x_position+tile_width, y_position+0.45, (double) x_position+tile_width, y_position+0.45+ ((double) tile_height)/4, red, green, blue);
            
            /* BOX TILE fprintf(outfile, " <rect x=\"%.0f\" y=\"%.2f\" width=\"%.0f", x_position+3* ((double) tile_width)/4, y_position+0.45,((double) tile_width) /4);
            fprintf(outfile, "\" height=\"");
            fprintf(outfile, "%i", tile_height/4);
            fprintf(outfile, "\" style=\"");
            fprintf(outfile, "fill:rgb(");
            fprintf(outfile, "%i,%i,%i);stroke-width:0.5;stroke:rgb(0,0,0)\"/>\n", red, green, blue, stroke);
            */
            
            /* PRINTS YELLOW CIRCLE IF LOCAL MAX */
            if(Localmax(results, file_number, shortest_kmer, too_long_kmer, current_kmer_length, current_gap_position, current_gap_length, current_kmer, count_also_spaced_kmers, kmer_length_difference_cutoff, kmerstring))
            {
                fprintf(outfile, " <circle cx=\"%li\" cy=\"%li\" r=\"%i\" fill=\"yellow\" stroke-width=\"%i\"/>\n", x_position + tile_width/2, y_position + tile_width/2, tile_width/4, 0);
            }
            x_position += tile_width;
        }
        
        /* SAVES ONE LINE OF OUTPUT VALUES TO TSV_FILE */
        fprintf(tsv_outfile, "%s\t", heatmap_name);
        for(position = current_kmer_length-1; position > -1 ; position--) {if(current_kmer_length - position - 1  == current_gap_position) fprintf(tsv_outfile,"."); fprintf(tsv_outfile, "%c", forward[(current_kmer & mask_ULL[1][position]) >> (position * 2)]);}
        fprintf(tsv_outfile, "\t");
        for(current_gap_length = 0, current_gap_position = 0; current_gap_length < Nlength - current_kmer_length; current_gap_length++, current_gap_position = true_gap_position) fprintf(tsv_outfile, "%li\t", results[file_number][current_kmer_length][current_gap_position][current_gap_length][current_kmer]);
        for(current_gap_length = 0, current_gap_position = 0; current_gap_length < Nlength - current_kmer_length; current_gap_length++, current_gap_position = true_gap_position) fprintf(tsv_outfile, "%i\t", Localmax(results, file_number, shortest_kmer, too_long_kmer, current_kmer_length, current_gap_position, current_gap_length, current_kmer, count_also_spaced_kmers, kmer_length_difference_cutoff, kmerstring));
        for(current_gap_length = 0, current_gap_position = 0; current_gap_length < Nlength - current_kmer_length; current_gap_length++, current_gap_position = true_gap_position) 
        {
        current_count = results[file_number][current_kmer_length][current_gap_position][current_gap_length][current_kmer];
        current_background_count = results[file_number-1][current_kmer_length][current_gap_position][current_gap_length][current_kmer];
        current_enrichment = (double) ((float) current_count * size_factor + kmer_pseudocount) / ((float) current_background_count + kmer_pseudocount);
        fprintf(tsv_outfile, "%.2f\t", current_enrichment);
        }
        
        peak_features = Peak_to_peak_distance(results[file_number][current_kmer_length], current_gap_position, current_kmer, Nlength - current_kmer_length, 0.1, peak_features);
        shape = Peak_shape(results[file_number][current_kmer_length],current_gap_position, current_kmer, Nlength - current_kmer_length);
        
        /* GENERATES PROTEIN / DNA CALL */
        call = '-';
        if(peak_features[0] > 8 && peak_features[0] < 12 && shape < 1) call = 'D';
        if((peak_features[0] < 9 || peak_features[0] > 11) && shape > 2) call = 'P';
        
        fprintf(tsv_outfile,"%c\t%.0f\t%.1f\t%.2f\t%li\t%.1f\t%.0f\t%.1f\t%li\n", call,peak_features[0], peak_features[1], shape, (top_kmers[current_heatmap_row]).local_max_max_incidence, (top_kmers[current_heatmap_row]).local_max_enrichment, max_score, (top_kmers[current_heatmap_row]).max_enrichment, total_score);
        
        x_position += tile_width * (current_kmer_length - shortest_kmer); /* ADDS SPACE IF LONGER THAN MIN KMER */
        /* PRINTS MAX COUNT*/
        
        /* PRINTS PEAK TO PEAK DISTANCE AND RATIO */
        fprintf(outfile, "<text xml:space=\"preserve\" x=\"");
        fprintf(outfile, "%li", x_position);
        fprintf(outfile, "\" y=\"");
        fprintf(outfile, "%li", y_position + tile_height - 7);
        fprintf(outfile, "\" fill = \"black\" stroke = \"none\" font-size=\"%i\" font-family = \"", font_size);
        fprintf(outfile, "%s", font);
        fprintf(outfile, "\" >");
        fprintf(outfile,"%2c%3.0f%4.1f%7.2f%9li%9.1f%9.0f%9.1f%9li", call, peak_features[0], peak_features[1], shape, (top_kmers[current_heatmap_row]).local_max_max_incidence, (top_kmers[current_heatmap_row]).local_max_enrichment, max_score, (top_kmers[current_heatmap_row]).max_enrichment, total_score);
        fprintf(outfile, "</text>\n"); 
        
        /* SPECIAL EFFECTS FOR PREFERRED KMER HEATMAP */
        if (current_kmer_length == preferred_length && heatmap == 0)
        {
            if ((top_kmers[current_heatmap_row+1]).preferred == 0) break;
            if ((top_kmers[current_heatmap_row+1]).preferred <= 6 && first_shift == 0) {y_position += 20; first_shift = 1;} 
            if ((top_kmers[current_heatmap_row+1]).preferred <= 3 && second_shift == 0) {y_position += 20; second_shift = 1;} 
        }
    }
        
    if (heatmap == 0)
    {
    /* FINDS ALL KMERS WITH AT LEAST ONE HIGH SCORING LOCAL MAX AND SORTS THEM ACCORDING TO INCIDENCE */
    qsort (top_kmers, number_of_kmers, sizeof(struct kmer_incidence_table), Sort_according_to_incidence); /* IN CASE THERE IS NO LOCAL MAX */
    qsort (top_kmers, number_of_kmers, sizeof(struct kmer_incidence_table), Sort_according_to_local_max_max_incidence);
    }
        
    /* SORTS ACCORDING TO MAX LOCAL MAX ENRICHMENT */
    if (heatmap == 1) qsort (top_kmers, number_of_kmers, sizeof(struct kmer_incidence_table), Sort_according_to_local_max_enrichment);
        
    /* SORTS ACCORDING TO MAX ENRICHMENT */
    if (heatmap == 2) qsort (top_kmers, number_of_kmers, sizeof(struct kmer_incidence_table), Sort_according_to_enrichment);
        
    /* SORTS ALL KMERS WITH AT LEAST ONE LOCAL MAX ACCORDING TO INCIDENCE */
    if (heatmap == 3) 
    {
    qsort (top_kmers, number_of_kmers, sizeof(struct kmer_incidence_table), Sort_according_to_local_max);
    for(counter = 0; counter < number_of_kmers; counter++) if ((top_kmers[counter]).local_max == 0) break;
    qsort (top_kmers, counter, sizeof(struct kmer_incidence_table), Sort_according_to_incidence);
    }
 
    /* SORTS ACCORDING TO INCIDENCE */
    if (heatmap == 4) qsort (top_kmers, number_of_kmers, sizeof(struct kmer_incidence_table), Sort_according_to_incidence);
        
    /* SORTS KMERS BY NUMBER OF ANY LOCAL MAXES */
    if (heatmap == 5) qsort (top_kmers, number_of_kmers, sizeof(struct kmer_incidence_table), Sort_according_to_any_local_max);
        
    y_position += 50;
        
    }
    
    
    fprintf(outfile, "</svg>\n");
    fclose(outfile);
    fprintf(tsv_outfile, "END\n");
    fclose(tsv_outfile);
    free(top_kmers);
    free(tsv_filename);
    free(heatmap_name);
    free(arrow_color[0]);
    free(arrow_color[1]);
    free(arrow_color);
    free(peak_features);
    for (counter = 0; counter < 3; counter++) free(order[counter]);
    free(order);
    return(0);
}



/* GENERATES ADJACENT DINUC MARKOV MODEL a FROM BASE DEPENDENCY MATRIX d */
short int Generate_fADM (struct adjacent_dinucleotide_model *a, struct base_dependency_matrix *d, short int length)
{
    short int dinucleotide;
    short int mononucleotide1;
    short int mononucleotide2;    
    short int first;
    double total;
    short int seed_start = Nlength - 1;
    
    /* NORMALIZES COUNTS FOR FIRST DINUCLEOTIDE INSIDE SEED */
    for(dinucleotide = 0, total = 0; dinucleotide < 16; dinucleotide++) total += (*d).incidence[seed_start+1][seed_start][dinucleotide];
    for(dinucleotide = 0; dinucleotide < 16; dinucleotide++) (*a).fraction[dinucleotide][0] = ((double) (*d).incidence[seed_start+1][seed_start][dinucleotide]) / total;
    
    /* NORMALIZES COUNTS FOR OTHER CONSECUTIVE DINUCLEOTIDES INSIDE SEED (CONDITIONAL PROBABILITIES) */
    for(first = 1; first < length-1; first++)
    {
        for(mononucleotide1=0; mononucleotide1 < 4; mononucleotide1++)
        {
            for(total = 0, mononucleotide2=0, dinucleotide = mononucleotide1 * 4; mononucleotide2 < 4; mononucleotide2++, dinucleotide++) total += (*d).incidence[seed_start+first+1][seed_start+first][dinucleotide];        
            for(mononucleotide2=0, dinucleotide = mononucleotide1 * 4; mononucleotide2 < 4; mononucleotide2++, dinucleotide++) (*a).fraction[dinucleotide][first] = ((double) (*d).incidence[seed_start+first+1][seed_start+first][dinucleotide]) / total;      
        }
    }
    
    return(0);    
}

/* GENERATES ADJACENT DINUC MARKOV MODEL a FROM BASE DEPENDENCY MATRIX d */
short int Generate_ADM (struct adjacent_dinucleotide_model *a, struct base_dependency_matrix *d, short int length)
{
    short int dinucleotide;
    short int mononucleotide1;
    short int mononucleotide2;    
    short int first;
    double total;
    double subtotal;
    short int seed_start = Nlength - 1;
    (*a).width = length;
    length--;

    
    /* COUNTS MONONUCLEOTIDES AND NORMALIZES COUNTS FOR CONSECUTIVE DINUCLEOTIDES INSIDE SEED (CONDITIONAL PROBABILITIES) */
    for(first = 0; first < length; first++)
    {
        for(total = 0, dinucleotide = 0; dinucleotide < 16; dinucleotide++) total += (*d).incidence[seed_start+first+1][seed_start+first][dinucleotide];   
        for(mononucleotide1=0; mononucleotide1 < 4; mononucleotide1++)
        {
            for(subtotal = 0, mononucleotide2=0, dinucleotide = mononucleotide1 * 4; mononucleotide2 < 4; mononucleotide2++, dinucleotide++) subtotal += (*d).incidence[seed_start+first+1][seed_start+first][dinucleotide];        
            for(mononucleotide2=0, dinucleotide = mononucleotide1 * 4; mononucleotide2 < 4; mononucleotide2++, dinucleotide++) (*a).fraction[dinucleotide][first] = ((double) (*d).incidence[seed_start+first+1][seed_start+first][dinucleotide]) / subtotal;
            (*a).mononuc_fraction[mononucleotide1][first] = ((double) subtotal) / total;
        }
    }
    
    /* LAST POSITION OF MONONUCLEOTIDES */
    for(mononucleotide1=0; mononucleotide1 < 4; mononucleotide1++)
    {
        for(subtotal = 0, dinucleotide = mononucleotide1; dinucleotide < 16; dinucleotide += 4) subtotal += (*d).incidence[seed_start+first][seed_start+first-1][dinucleotide]; 
        (*a).mononuc_fraction[mononucleotide1][first] = ((double) subtotal) / total;
    }
    
    return(0);    
}

/* PRINTS ADJACENT DINUC MARKOV MODEL a */
short int Print_ADM (struct adjacent_dinucleotide_model *a)
{
    char *forward;
    if (rna == 0) forward = dnaforward;
    else forward = rnaforward;
    short int dinucleotide;   
    short int first;
    for(dinucleotide = 0; dinucleotide < 16; dinucleotide++)
    {
        printf("\n");
        for(first = 0; first < (*a).width-1; first++) printf("%.3f\t", (*a).fraction[dinucleotide][first]);
        printf("ADM_DI\t%c%c",forward[(dinucleotide & 12)>>2], forward[dinucleotide & 3]);
    }
    for(dinucleotide = 0; dinucleotide < 4; dinucleotide++)
    {
        printf("\n");
        for(first = 0; first < (*a).width; first++) printf("%.3f\t", (*a).mononuc_fraction[dinucleotide][first]);
        printf("ADM_MONO_%c", forward[dinucleotide]);
    }
    return(0);
}


/* LOG TRANSFORMS ADM */
short int Log_fold_affinity_ADM (struct adjacent_dinucleotide_model *a)
{
    short int block = 0;
    short int line = 0;
    short int pwm_position = 0;
    double column_maximum;
    for (pwm_position = 0; pwm_position < (*a).width-1; pwm_position++)
    {
        for(block = 0; block < 13; block +=4)
        {
        for(column_maximum = 0, line = 0; line <= 3; line++) if ((*a).fraction[line+block][pwm_position] > column_maximum) column_maximum = (*a).fraction[line+block][pwm_position];
        for(line = 0; line <= 3; line++) 
        {
            if ((*a).fraction[line+block][pwm_position] != 0) (*a).fraction[line+block][pwm_position] = log10((*a).fraction[line+block][pwm_position] / column_maximum);
            else (*a).fraction[line+block][pwm_position] = log10(pseudocount);
        }
        }
        
    }
    
    for (pwm_position = 0; pwm_position < (*a).width; pwm_position++)  
    {
    for(column_maximum = 0, line = 0; line <= 3; line++) if ((*a).mononuc_fraction[line][pwm_position] > column_maximum) column_maximum = (*a).mononuc_fraction[line][pwm_position];
    for(line = 0; line <= 3; line++) 
    {
        if ((*a).mononuc_fraction[line][pwm_position] != 0) (*a).mononuc_fraction[line][pwm_position] = log10((*a).mononuc_fraction[line][pwm_position] / column_maximum);
        else (*a).mononuc_fraction[line][pwm_position] = log10(pseudocount);
    }
    }
    Print_ADM(a);
    return(0);
}

/* SUBROUTINE THAT CONVERTS A FOLD AFFINITY ADM TO A LOG RATIO ADM (DEFINES SPECIFICITY OF EACH NUCLEOTIDE AT EACH POSITION) */
short int Log_ratio_ADM (struct adjacent_dinucleotide_model *a)
{
    short int block = 0;
    short int line = 0;
    short int pwm_position = 0;
    for (pwm_position = 0; pwm_position < (*a).width-1; pwm_position++)
    {
        for(block = 0; block < 13; block +=4)
        {
        for(line = 0; line < 4; line++) 
        {
            (*a).fraction[block+line][pwm_position] = log10(((*a).fraction[block+line][pwm_position] + pseudocount) / (pseudocount + 1 - (*a).fraction[block+line][pwm_position]));
        }
        }
    }
    for (pwm_position = 0; pwm_position < (*a).width; pwm_position++)  
    {
        for(line = 0; line < 4; line++) 
        {
            (*a).mononuc_fraction[line][pwm_position] = log10(((*a).mononuc_fraction[line][pwm_position] + pseudocount) / (pseudocount + 1 - (*a).mononuc_fraction[line][pwm_position]));
        }
    }
    
    Print_ADM(a);
    return(0);
}


/* SUBROUTINE THAT EXTRACTS A PWM FROM AN ADM */
short int PWM_from_ADM (struct adjacent_dinucleotide_model *a, struct normalized_pwm *p)
{
    short int line = 0;
    short int pwm_position = 0;
    for(line = 0; line < 4; line++) for (pwm_position = 0; pwm_position < (*a).width; pwm_position++)
    {
        (*p).fraction[line][pwm_position] = (*a).mononuc_fraction[line][pwm_position];    
    }
    (*p).width = (*a).width;
    return(0);
}


/* SUBROUTINE THAT NORMALIZES ADM a TO YIELD PROBABILITIES */
short int normalize_ADM (struct adjacent_dinucleotide_model *a)
{
    short int line = 0;
    short int block = 0;
    short int pwm_position = 0;
    double subtotal;
    for (pwm_position = 0; pwm_position < (*a).width; pwm_position++) 
    {
    for(subtotal = 0, line = 0; line < 4; line++) subtotal += (*a).mononuc_fraction[line][pwm_position];
    for(line = 0; line < 4; line++) (*a).mononuc_fraction[line][pwm_position] = (*a).mononuc_fraction[line][pwm_position] / subtotal;    
    if (pwm_position < (*a).width - 1) for(block=0; block < 13; block += 4)
    {
    for(subtotal = 0, line = 0; line < 4; line++) subtotal += (*a).fraction[block+line][pwm_position];
    for(line = 0; line < 4; line++) (*a).fraction[block+line][pwm_position] = (*a).fraction[block+line][pwm_position] / subtotal;
    }
    }
    return(0);
}

/* SUBROUTINE THAT GENERATES ADM FROM A NORMALIZED PWM */
short int PWM_to_ADM (struct adjacent_dinucleotide_model *a, struct normalized_pwm *p)
{
    short int line = 0;
    short int block = 0;
    short int pwm_position = 0;
    for(line = 0; line < 4; line++) for (pwm_position = 0; pwm_position < (*a).width - 1; pwm_position++)
    {
    (*a).mononuc_fraction[line][pwm_position] = (*p).fraction[line][pwm_position]; 
    for(block=0; block < 4; block++) (*a).fraction[block*4+line][pwm_position] = (*p).fraction[block][pwm_position] * (*p).fraction[line][pwm_position+1]; 
    }
    (*a).mononuc_fraction[line][pwm_position] = (*p).fraction[line][pwm_position]; 
    (*a).width = (*p).width;
    normalize_ADM(a);
    return(0);
}


/* SUBROUTINE THAT LOADS AN ADM */
short int Load_ADM (struct adjacent_dinucleotide_model *a, char *filename)
{
    long int counter;
    char text1;
    short int line = 0;
    short int pwm_position = 0;
    char *current_string;
    current_string = malloc(200);
    FILE *pwmfile;
    if ((pwmfile = fopen(filename, "r")) == '\0') {printf("\nNo File: %s", filename); exit (2);}
    for(line = 0; line < 16;)
    {
        for(counter = 0; counter < 30; counter++)
        {
            text1 = getc(pwmfile);
            if (text1 == EOF || text1 == '\n' || text1 == '\t')
            {
                current_string[counter] = '\0'; 
                if (counter > 0 && (current_string[0] == '0' || current_string[0] == '1' || current_string[0] == '2' || current_string[0] == '3' || current_string[0] == '4' || current_string[0] == '5' || current_string[0] == '6' || current_string[0] == '7' || current_string[0] == '8' || current_string[0] == '9' || current_string[0] == ' ' || current_string[0] == '-')) 
                {
                    (*a).fraction[line][pwm_position] = atof(current_string);
                    printf("\t%f", (*a).fraction[line][pwm_position]);
                    fflush(stdout);
                    pwm_position++; 
                }
                if (text1 == '\n' || text1 == EOF) {(*a).width = pwm_position; line++; pwm_position = 0;}
                break;
            }
            current_string[counter]=text1;
            /* printf ("%c", text1); */
        }
    }
    pwm_position = 0;
    for(line = 0; line < 4;)
    {
        for(counter = 0; counter < 30; counter++)
        {
            text1 = getc(pwmfile);
            if (text1 == EOF || text1 == '\n' || text1 == '\t')
            {
                current_string[counter] = '\0'; 
                if (counter > 0 && (current_string[0] == '0' || current_string[0] == '1' || current_string[0] == '2' || current_string[0] == '3' || current_string[0] == '4' || current_string[0] == '5' || current_string[0] == '6' || current_string[0] == '7' || current_string[0] == '8' || current_string[0] == '9' || current_string[0] == ' ' || current_string[0] == '-')) 
                {
                    (*a).mononuc_fraction[line][pwm_position] = atof(current_string);
                    printf("\t%f", (*a).mononuc_fraction[line][pwm_position]);
                    fflush(stdout);
                    pwm_position++; 
                }
                if (text1 == '\n' || text1 == EOF) {(*a).width = pwm_position; line++; pwm_position = 0;}
                break;
            }
            current_string[counter]=text1;
            /* printf ("%c", text1); */
        }
    }
    
    free (current_string);
    Print_ADM(a);
    if (text1 == EOF && line != 19) return(1);
    else return (0);
}

/* SUBROUTINE THAT FINDS BEST SCORE AND POSITION FOR KMER IN A WIDER ADM */
double Kmerscore_ADM (struct adjacent_dinucleotide_model *a, __uint128_t kmer_sequence_value, short int kmer_length, short int *kmermatch_position)
{
    short int nucleotide;
    short int rev_nucleotide;
    short int kmer_position;
    double score = 0;
    double revscore = 0;
    double best_score = -100;
    signed short int pwm_position;
    __uint128_t reverse_complement_kmer_sequence_value;
    __uint128_t current_kmer_sequence_value;
    __uint128_t current_reverse_complement_kmer_sequence_value;
    
    reverse_complement_kmer_sequence_value = Reverse_complement_sequence_value(kmer_sequence_value, kmer_length);
    
    /* printf("\nFORWARD_SEQ: "); Kmerprint(kmer_sequence_value, kmer_length); */
    
    for (pwm_position = 0; pwm_position < (*a).width - kmer_length; pwm_position++)
    {
        /* GENERATES KMER SCORE */
        
      
        /* CONDITIONAL PROBABILITIES */
        for(kmer_position = 0, score= 0, revscore = 0, current_kmer_sequence_value = kmer_sequence_value, current_reverse_complement_kmer_sequence_value = reverse_complement_kmer_sequence_value; ; kmer_position++, current_kmer_sequence_value >>= 2, current_reverse_complement_kmer_sequence_value >>= 2) 
        {
            nucleotide = (current_kmer_sequence_value & 15);
 //           nucleotide = ((nucleotide & 3) << 2) + ((nucleotide & 12) >> 2);
            rev_nucleotide = (current_reverse_complement_kmer_sequence_value & 15);
   //         rev_nucleotide = ((rev_nucleotide & 3) << 2) + ((rev_nucleotide & 12) >> 2);
            // printf("\nDinucleotide %c%c at kmer_position %i and pwm_position %i, score addition %.2f", forward[(nucleotide & 12) >> 2], forward[nucleotide & 3], kmer_position, pwm_position, (*a).fraction[nucleotide][(*a).width-1-pwm_position-kmer_position]);
            score += (*a).fraction[nucleotide][(*a).width-1-pwm_position-kmer_position];
            revscore += (*a).fraction[rev_nucleotide][(*a).width-1-pwm_position-kmer_position];

        
        /* ENTRY POSITION */
        if(kmer_position > kmer_length-3)
        {
        nucleotide = (current_kmer_sequence_value & 12) >> 2;
        score += (*a).mononuc_fraction[nucleotide][(*a).width-1-pwm_position-kmer_position];
        rev_nucleotide = (current_reverse_complement_kmer_sequence_value & 12) >> 2;
        revscore += (*a).mononuc_fraction[rev_nucleotide][(*a).width-1-pwm_position-kmer_position];   
        // printf("\nMononucleotide %c at kmer_position %i and pwm_position %i, score addition %.2f", forward[nucleotide], kmer_position, pwm_position, (*a).mononuc_fraction[nucleotide][(*a).width-1-pwm_position-kmer_position]); 
        break;
        }
        }
        
        if (revscore > score) score = revscore;
        /* printf ("\nScore %f at position %i", score, (*p).width-1-pwm_position-kmer_position); */
        if (score >= best_score) /* SCORES PWM AT THIS POSITION */
        {
            // printf("\nBEST THUS FAR %f", score);
            best_score = score;
            *kmermatch_position = kmer_position;
            /* printf ("\nMATCH at %i: score %f higher than cut-off %f", match.position[number_of_matches], score, cut_off); */
        }
        /* else printf ("\nNo match at %i: %i is not %i", Nlength - current_position - query_sequence_length, query_sequence_ULL, test_sequence_ULL & query_sequence_ULL); */
    }
            // printf("\n*** BEST ***  %f\n", best_score);
    return (best_score);
}

/* SUBROUTINE THAT ADDS FLANKS TO AN ADM */
short int Add_flanks_to_ADM(struct adjacent_dinucleotide_model *b, struct adjacent_dinucleotide_model *a, short int length, double fill_value)
{
    short int nucleotide;
    short int pwm_position = (*a).width + 2 * length;

    /* DINUCLEOTIDES */
    for ( ; pwm_position >= (*a).width + length - 1; pwm_position--) for (nucleotide = 0; nucleotide < 16; nucleotide++) (*b).fraction[nucleotide][pwm_position] = fill_value;
    for ( ; pwm_position >= length; pwm_position--) for (nucleotide = 0; nucleotide < 16; nucleotide++) (*b).fraction[nucleotide][pwm_position] = (*a).fraction[nucleotide][pwm_position-length];
    for ( ; pwm_position >= 0; pwm_position--) for (nucleotide = 0; nucleotide < 16; nucleotide++) (*b).fraction[nucleotide][pwm_position] = fill_value;
    
    /* MONONUCLEOTIDES */
    pwm_position = (*a).width + 2 * length;
    for ( ; pwm_position >= (*a).width + length; pwm_position--) for (nucleotide = 0; nucleotide < 4; nucleotide++) (*b).mononuc_fraction[nucleotide][pwm_position] = fill_value;
    for ( ; pwm_position >= length; pwm_position--) for (nucleotide = 0; nucleotide < 4; nucleotide++) (*b).mononuc_fraction[nucleotide][pwm_position] = (*a).mononuc_fraction[nucleotide][pwm_position-length];
    for ( ; pwm_position >= 0; pwm_position--) for (nucleotide = 0; nucleotide < 4; nucleotide++) (*b).mononuc_fraction[nucleotide][pwm_position] = fill_value;
    
    (*b).width = (*a).width + 2 * length;
        Print_ADM(b);
    return(0);
}

/* CONVERTS SEQUENCE VALUE TO IUPAC SEQUENCE VALUE */
sequence_value_to_iupac_sequence_value (__uint128_t sequence_value_ULL, short int sequence_length, short int gap_position, short int gap_length, short int add_flanks, struct bitiupac_structure iupac_sequence_value)
{
    short int position;
    
}

/* ALIGNS TWO PWMs, returns structure that contains the orientation of PWM1 that aligns better */
struct normalized_pwm *Align_pwm(struct normalized_pwm *qp, struct normalized_pwm *rc_qp, struct normalized_pwm *mfp)
{
short int strand;
short int position;
short int current_position;
double *position_score;
double sum_x_squared;
double sum_y_squared;
double sum_xy;
position_score = malloc(sizeof(double)*(*qp).width+5);
for (position = 0; position < (*qp).width; position++) position_score[position] = 0;
double current_difference = 0;
(*qp).match.score = -10000000;
short int nucleotide;
double euclidian_distance;
double x_ic;
double y_ic;
double max_ic;
double kl_divergence;
struct normalized_pwm *fw_qp;
fw_qp = qp;
    
    
    for (strand = 1; strand > -2; strand -= 2)
    {
    for (position = 0; position < (*mfp).width - (*qp).width; position++)
    {
        for (current_difference = 0, current_position = 0; current_position < (*qp).width; current_position++)
        {
            for (x_ic = 2, y_ic = 2, euclidian_distance = 0, kl_divergence = 0, sum_x_squared = 0, sum_y_squared = 0, sum_xy = 0, position_score[current_position] = 0, nucleotide = 0; nucleotide < 4; nucleotide++)
            {
                /* K-L DIVERGENCE */
                kl_divergence += log(((*mfp).fraction[nucleotide][position+current_position]+pseudocount)/((*qp).fraction[nucleotide][current_position]+pseudocount)) * (*mfp).fraction[nucleotide][position+current_position];
                
                /* INFORMATION CONTENTS */
                x_ic += (*mfp).fraction[nucleotide][position+current_position] * log2 ((*mfp).fraction[nucleotide][position+current_position] + pseudocount);
                y_ic += (*qp).fraction[nucleotide][current_position] * log2 ((*qp).fraction[nucleotide][current_position] + pseudocount);
               
                /* EUCLIDIAN DISTANCE */
                euclidian_distance += fabs ((*mfp).fraction[nucleotide][position+current_position] - (*qp).fraction[nucleotide][current_position]);
                
                /* UNCENTERED CORRELATION */
                sum_x_squared += (double) pow((*mfp).fraction[nucleotide][position+current_position], 2);
                sum_y_squared += (double) pow((*qp).fraction[nucleotide][current_position], 2);
                sum_xy += (double) (*mfp).fraction[nucleotide][position+current_position] * (*qp).fraction[nucleotide][current_position];
                
                /* printf("\nComparing pwm1 position %i to pwm2 position %i at nucleotide %i, %.2f-%.2f current kl, euc, xic and yic are %.2f, %.2f, %.2f and %.2f", current_position, position+current_position, nucleotide, (*qp).fraction[nucleotide][current_position], (*mfp).fraction[nucleotide][position+current_position], kl_divergence, euclidian_distance, x_ic, y_ic); */
            }
            /* position_score[current_position] = sum_xy/(sqrt(sum_y_squared)*sqrt(sum_x_squared)); */
            max_ic = x_ic; if (y_ic > max_ic) max_ic = y_ic;
            position_score[current_position] = max_ic * (1 - euclidian_distance);
            current_difference += position_score[current_position];
        }
        if (current_difference > (*qp).match.score) 
        {
        /* printf("\n\nCURRENT BEST SCORE AT %i IS %.2f", position, current_difference); */
        (*qp).match.score = current_difference;   
        (*qp).match.position = position;
        (*qp).match.strand = strand;
        for (current_position = 0; current_position < (*qp).width; current_position++) (*qp).position_score[current_position] = position_score[current_position];
        }
        /* printf("\n\nPOSITION %i score is %.2f", position, current_difference); */
    }

/* CHANGE PWM TO REVCOMP */
qp = rc_qp;
/* CARRY BEST SCORE TO REVCOMP PWM */
(*qp).match.score = (*fw_qp).match.score;
/* (*qp).match.position = (*fw_qp).match.position;
(*qp).match.strand = (*fw_qp).match.strand;
for (current_position = 0; current_position < (*qp).width; current_position++) (*qp).position_score[current_position] = (*qp_fw).position_score[current_position]; */
        
}
if ((*qp).match.strand == -1) return (rc_qp);
else return (fw_qp);
}

/* ########################## MAIN PROGRAM ############################# */
/* ########################## MAIN PROGRAM ############################# */
/* ########################## MAIN PROGRAM ############################# */
/* ########################## MAIN PROGRAM ############################# */
/* ########################## MAIN PROGRAM ############################# */
int main (int argc, char *argv[])

{

	Nlength++;


	time_t  t0, t1; 
	t0 = time('\0');

	double *error_values;
	error_values = malloc(sizeof(double) * 10 + 5);
	error_values[0] = 0;
	
	/* STRINGS FOR DIFFERENT ORIENTATIONS */
	char **orientation_string;
	orientation_string = malloc(sizeof(char *) * 3 + 5); 
	orientation_string[0] = "HT";
	orientation_string[1] = "HH";
	orientation_string[2] = "TT";
	
	char **repeat_report;
	repeat_report = malloc(sizeof(char *) * 6 + 5);
	repeat_report[0] = "--";
	repeat_report[1] = "Odd_Repeat";
	repeat_report[2] = "Even_Repeat";
	repeat_report[3] = "Odd_Palindrome";
	repeat_report[4] = "Even_Palindrome";
	
	long int *last_kmer;
	last_kmer = malloc(sizeof(long int) * 100 + 5);	

	double pvalue;
	short int shortcounter;
	long int counter;
	long int counter2;
	long int counter3;
	
	/* clears pvalue cache */
	for (counter = 0; counter < 1000; counter++) for (counter2 = 0; counter2 < 1000; counter2++) pvalue_cache[counter][counter2] = 2;
	
	FILE *open_file;
	char **file_name;
	file_name = malloc (sizeof (char *) * (number_of_files + 1) + 5);
	for (counter = 0; counter < number_of_files; counter++) {file_name[counter] = malloc(1000); strcpy(file_name[counter], "no file");}
	short int file_number = 0;
	double *number_of_sequences_analyzed;
    number_of_sequences_analyzed = malloc(sizeof(double) * (number_of_files + 1) + 5);
	double *number_of_sequences_with_hits;
    number_of_sequences_with_hits = malloc(sizeof(double) * (number_of_files + 1) + 5);
	double *number_of_sequences_with_no_hits;
    number_of_sequences_with_no_hits = malloc(sizeof(double) * (number_of_files + 1) + 5);
	for (counter = 0; counter < number_of_files; counter++) 
	{
	number_of_sequences_analyzed[counter] = 0;
	number_of_sequences_with_hits[counter] = 0;
	number_of_sequences_with_no_hits[counter] = 0;
	}
	short int firstnoncommandposition = 0;
	short int print_counts = 0;
	short int extendedoutput = 0;
	short int remove_non_unique = 0;
	short int print_input_sequences = 0;
	short int print_p_values = 0;
	short int output_all_gap_lengths_in_one_line = 0;
	char *linecommand;
	short int too_long_kmer = 8;
	short int shortest_kmer = 8;
	
	short int kmer_length;	
	long int line = 0;
	long int current_kmer;
	short int current_kmer_length;
	short int current_gap_position;
	short int current_gap_length;
	
	char *current_sequence = malloc(1000); 
	strcpy(current_sequence, "INITIALVALUE");
	char *tempstring = malloc(10000);
	__uint128_t current_sequence_value_ULL;
	__uint128_t forward_sequence_value_ULL;
	__uint128_t deleted_sequence_value_ULL;
	__uint128_t position_value;
	char *forward_lc = "acgt";
	char *dnareverse = "TGCA";
	signed short int position;
	short int start_position;
	short int end_position;
	short int nucleotide_value;
	short int kmer_count = 1;
	long int kmer_incidence;
	long int kmer2_incidence;
	double max_incidence;
	long int total_incidence;
	float current_fold_expected;
    short int nucleotide;
    short int max_nucleotide;
    long int max_kmer_count;
    long int test_kmer;
    long int test_kmer_count;
    short int iupac_bits;
    short int represents_n_nucleotides;
    long int represents_n_kmers;
    long int top_normalized_count;
    double kmer_length_difference_cutoff = 0.35;
    double iupac_cutoff = 0.25;
	short int position2;
    short int localmaxes = 0;
    
    __uint128_t query_sequence_value = 2*64 + 3*4; /* SEQUENCE: GATA */
	short int query_sequence_length = 4; 
    
	/* FORWARD AND REVERSE NUCLEOTIDE COUNT PWM STRUCTURES */
	struct count_pwm *nc;
	nc = malloc(sizeof(struct count_pwm) * 3 + 5);
	count_pwm_init(&nc[0], "EMPTY", max_Nlength, 0);
	count_pwm_init(&nc[1], "EMPTY", max_Nlength, 0);
	short int print_nucleotides = 0;
	short int print_frequencies = 0;
	char text1;
	long int charcounter;
	short int deletion_size = 1;
	short int Nmer_position;
	short int too_long_nmer_position;
	
	double minimum_kmer_count = 1;
	short int strand;
	
	short int non_unique_flag = 0;
	short int current_sorted_list = 0;
	short number_of_unordered_sequences = 0;
	struct sequence_incidence_table *index_position;
	struct sequence_incidence_table **sorted_list;
	sorted_list = malloc (sizeof(struct sequence_incidence_table *) * 2 + 5);
	struct sequence_incidence_table *unordered_list;
	unordered_list = malloc(4000);
	struct sequence_incidence_table **sorted_index;
	sorted_index = malloc(sizeof(struct sequence_incidence_table *) * 4098 + 10);
	short int current_sequence_contains_match = 0;
	
	long int current_number_of_different_sequences = 0;
	
    long int matches_to_filter = 0;
    
	short int eof_reached;
    
    double **flank_kmer_expected_count;
    long int *palindromic_hits;
    palindromic_hits = malloc(number_of_files * sizeof(long int) + 5);    
    flank_kmer_count = malloc(number_of_files * sizeof(long int **) + 5);
    flank_kmer_expected_count = malloc(number_of_files * sizeof(double **) + 5);
    for(file_number = 0; file_number < number_of_files; file_number++)
    {
    palindromic_hits[file_number] = 0;
    flank_kmer_count[file_number] = malloc(256 * sizeof(long int) + 5);
    flank_kmer_expected_count[file_number] = malloc(256 * sizeof(double) + 5);
    for(counter = 0; counter < 256; counter++) 
    {
    flank_kmer_count[file_number][counter] = 0;
    flank_kmer_expected_count[file_number][counter] = 0;
    }
    }
    long int *match_orientations = '\0';
    short int match_length = 0;
    
    /* VARIABLES AND STRUCTURES FOR HIT POSITION DATA GENERATION */
    struct hit_position_matrix hit_position;
    hit_position_matrix_init(&hit_position, "hit positions", max_Nlength, 0);
	short int correct_position = 1;
    char limit_hits_to_strand = 'N';
    short int limit_hits_to_position = 0;
    char *exclude_strand;
    exclude_strand = malloc(sizeof(char *) * max_Nlength + 5);
    short int *exclude_position;
    exclude_position = malloc(sizeof(short int *) * max_Nlength + 5);
    for (counter = 0; counter < max_Nlength; counter++) {exclude_position[counter] = 0; exclude_strand[counter] = 'N';}
    
    /* FLAGS */
    short int xyplot = 0;
	long int max_counts = 0;
	short int max_spacing = 0;
	short int max_orientation = 0;
    short int flank_with_pwm = 0;
	short int even_background = 0;
	short int only_palindromes = 0;
    short int palindrome_correction = 1;
    short int information_content_output = 0;
    short int count_also_spaced_kmers = 1;
    short int dinucleotide_properties = 0;
    short int searchstring_is_palindrome = 0;
	short int flank = 0;
	short int kmer_table = 0;
    short int multinomial = 0;
    short int complex_background = 1;
    short int center_spaced_kmer_table = 0;
    short int svg_only = 0;
    short int match_filter = 0;
    short int modeldistance = 0;
    short int difference_logo = 0;
    short int negative_values_allowed = 0;
    short int offset = 0;
    short int lowpaircount = 0;
    short int number_of_heatmap_rows = 20;
    char *user_specified_output_file;
    user_specified_output_file = malloc(1000);
    user_specified_output_file[0] = '\0';
    short int **max_align_score;
    max_align_score = malloc(sizeof(short int *) * 3 + 5);
    max_align_score[0] = malloc(sizeof(short int) * 3 + 5);
    max_align_score[1] = malloc(sizeof(short int) * 3 + 5);
    max_align_score[0][0] = 0;
    max_align_score[0][1] = 0;
    max_align_score[1][0] = 0;
    max_align_score[1][1] = 0;
    char *forward;
    forward = dnaforward;
    char *nucleotide_bitiupac = dna_bitiupac;
    char *nucleotide_iupac = dna_iupac;
    
    COMMAND = malloc(5000);
    strcpy(COMMAND, "");
    for(counter = 0; counter < argc; counter++) {strcat(COMMAND, argv[counter]);strcat(COMMAND, " ");}

/* COMMAND PARSER */	
/* PARSES OPTIONS */
if (argc >= 2 + firstnoncommandposition)
{
for (; ; firstnoncommandposition++)
{
linecommand = argv[1 + firstnoncommandposition];

if (linecommand[0] == '-' && linecommand[1] != '\0')
	{
        if (strcmp(linecommand, "--pwm") == 0) {flank_with_pwm = 1; flank = 1;}
		if (strcmp(linecommand, "--f") == 0) {flank = 1; kmer_count = 0;}
        if (strcmp(linecommand, "--logo") == 0) svg_only = 1;
        if (strcmp(linecommand, "--difflogo") == 0) {svg_only = 1; difference_logo = 1; paths = 1;}
        if (strcmp(linecommand, "--dist") == 0) modeldistance = 1;
        if (strcmp(linecommand, "--pwmalign") == 0) pwm_align = 1;
        if (strcmp(linecommand, "-editdist") == 0) {lowpaircount = 1; count_also_spaced_kmers = 1;};
        if (strcmp(linecommand, "-paths") == 0) paths = 1;
		if (strcmp(linecommand, "-barcodelogo") == 0) barcodelogo = 1;
        if (strcmp(linecommand, "-heightscaledbars") == 0) scale_bars = 1;
        if (strcmp(linecommand, "-maxheightscaledbars") == 0) {scale_bars = 1; max_scale_bar = 1;}
        if (strcmp(linecommand, "-colorscaledbars") == 0) gray_bars = 1;
        if (strcmp(linecommand, "-label") == 0) barcodelogolabels = 1;
        if (strcmp(linecommand, "-neg") == 0) {negative_values_allowed = 1; paths = 1;}
        if (strcmp(linecommand, "-noname") == 0) noname = 1;
        if (strcmp(linecommand, "-rna") == 0) {rna = 1; nucleotide_bitiupac = rna_bitiupac; nucleotide_iupac = rna_iupac; forward = rnaforward;}
        if (strcmp(linecommand, "-x") == 0) extendedoutput = 1;
		if (strcmp(linecommand, "-u") == 0) remove_non_unique = 1;
		if (strcmp(linecommand, "-c") == 0) print_counts = 1;
		if (strcmp(linecommand, "-p") == 0)	print_p_values = 1;
		if (strcmp(linecommand, "-i") == 0) print_input_sequences = 1;
		if (strcmp(linecommand, "-s") == 0)	output_all_gap_lengths_in_one_line = 1;
		if (strcmp(linecommand, "-n") == 0) print_nucleotides = 1;
		if (strcmp(linecommand, "-q") == 0) print_frequencies = 1;
        if (strcmp(linecommand, "-a") == 0) align_matches = 1;
        if (strcmp(linecommand, "-xyplot") == 0) xyplot = 1; 
        if (strcmp(linecommand, "-contacts") == 0) contacts = 1;
        if (strcmp(linecommand, "-printlocalmax") == 0) print_local_max = 1;
		if (strcmp(linecommand, "-nogaps") == 0) count_also_spaced_kmers = 0;
		if (strcmp(linecommand, "-allgaps") == 0) count_also_spaced_kmers = 2;
        if (strcmp(linecommand, "-both") == 0) count_both_instances_of_palindromic_hit = 1;
        if (strcmp(linecommand, "-nocall") == 0) nocall = 1;
        if (strcmp(linecommand, "-CpG") == 0) {methylCGcompare = 1; dinucleotide_properties = 1;}
        if (strcmp(linecommand, "-bothifnotequal") == 0) count_unequal_hits_only = 1;
        if (strcmp(linecommand, "-forwardonly") == 0) count_only_forward_instance_of_palindromic_hit = 1;
        if (strcmp(linecommand, "-reverseonly") == 0) count_only_reverse_instance_of_palindromic_hit = 1;
        if (strcmp(linecommand, "-forwardifequal") == 0) prefer_forward_strand = 1;
        if (strcmp(linecommand, "-reverseifequal") == 0) prefer_reverse_strand = 1;
		if (strcmp(linecommand, "-dimer") == 0) only_palindromes = 1;
		if (strcmp(linecommand, "-dinuc") == 0) dinucleotide_properties = 1;
		if (strcmp(linecommand, "-kmer") == 0) kmer_table = 1;
		if (strcmp(linecommand, "-ic") == 0) information_content_output = 1;
		if (strcmp(linecommand, "-e") == 0) {even_background = 1; file_number++;}
        if (strcmp(linecommand, "-14N") == 0) Nlength = 15;
		if (strcmp(linecommand, "-26N") == 0) Nlength = 27;
        if (strcmp(linecommand, "-30N") == 0) Nlength = 31;
        if (strcmp(linecommand, "-40N") == 0) Nlength = 41;
		if (strcmp(linecommand, "-mono") == 0) complex_background = 0;

        if (linecommand[0] == '-' && linecommand[1] == 'o' && linecommand[2] == '=') strcpy(user_specified_output_file, linecommand+3);
        
		if (linecommand[0] == '-' && linecommand[1] == 'm' && linecommand[2] == '=') multinomial = atoi(linecommand+3);
        
		if (linecommand[0] == '-' && linecommand[1] == 'f' && linecommand[2] == 'k' && linecommand[3] == '=') flank_kmer_pos = atoi(linecommand+4);

        if (linecommand[0] == '-' && linecommand[1] == 'm' && linecommand[2] == 'f' && linecommand[3] == '=') match_filter = atoi(linecommand+4);

        if (linecommand[0] == '-' && linecommand[1] == 'k' && linecommand[2] == 'l' && linecommand[3] == '=') 
        {
        shortest_kmer = atoi(linecommand+4);
        too_long_kmer = shortest_kmer+1;
        }
        
        if (linecommand[0] == '-' && linecommand[1] == 'l' && linecommand[2] == 'o' && linecommand[7] == '=')
		{
			kmer_length_difference_cutoff = atof(linecommand+8);
		}
        
        if (linecommand[0] == '-' && linecommand[1] == 'i' && linecommand[2] == 'u' && linecommand[6] == '=')
		{
			iupac_cutoff = atof(linecommand+7);
		}

        if (linecommand[0] == '-' && linecommand[1] == 'h' && linecommand[2] == 'r' && linecommand[6] == '=')
		{
			number_of_heatmap_rows = atof(linecommand+7);
		}
        
        if (linecommand[0] == '-' && linecommand[1] == 'e' && linecommand[2] == 'd' && linecommand[9] == '=') 
        {
            lowpaircount = 1;
            local_max_min_percent = ((double) atoi(linecommand+10))/100;
        }
        
        if (linecommand[0] == '-' && linecommand[1] == 'l' && linecommand[2] == 'i' && linecommand[3] == 'm' && linecommand[4] == '=')
		{
			limit_hits_to_position = atoi(linecommand+5);
			limit_hits_to_strand = linecommand[strlen(linecommand)-1];
		}

if (linecommand[0] == '-' && linecommand[1] == 'e' && linecommand[2] == 'x' && linecommand[3] == 'c' && linecommand[4] == '=')
{
for(counter = 0, linecommand += 5 ; counter < max_Nlength ; counter++, linecommand++)
{
exclude_position[counter] = atoi(linecommand);
for (; *linecommand != ',' && *linecommand != '\0'; linecommand++);
exclude_strand[counter] = *(linecommand - 1);
if (*linecommand == '\0') break; 
}
}

if (linecommand[0] == '-' && linecommand[1] == 'm' && linecommand[2] == 'a' && linecommand[3] == 't' && linecommand[6] == '=')
{
    char *halfsite1;
    char *halfsite2;
    halfsite1=malloc(1000);
    halfsite2=malloc(1000);
    for (counter = 7; counter < strlen(linecommand), linecommand[counter] != ','; counter++);
    strcpy(halfsite2, linecommand+counter+1);
    strcpy(halfsite1, linecommand+7);
    halfsite1[counter-7] = '\0';
    //printf("\nHalfsites %s,%s, %i\n", halfsite1, halfsite2, counter);
    for(query_sequence_value = 0, position = 0, position_value = pow(4,strlen(halfsite1)-1); position < strlen(halfsite1); position++, position_value /= 4)
	{
		for (nucleotide_value = 0; nucleotide_value < 4 && halfsite1[position] != forward[nucleotide_value]; nucleotide_value++);
		if(nucleotide_value == 4) {printf("\nERROR IN MATCH SEQUENCE\n"); exit (1);}
		query_sequence_value += position_value * nucleotide_value;
	}
    long int halfsite1_value = query_sequence_value;
    for(query_sequence_value = 0, position = 0, position_value = pow(4,strlen(halfsite2)-1); position < strlen(halfsite2); position++, position_value /= 4)
	{
		for (nucleotide_value = 0; nucleotide_value < 4 && halfsite2[position] != forward[nucleotide_value]; nucleotide_value++);
		if(nucleotide_value == 4) {printf("\nERROR IN MATCH SEQUENCE\n"); exit (1);}
		query_sequence_value += position_value * nucleotide_value;
	}
    long int halfsite2_value = query_sequence_value;    
    
    match_orientations = malloc (sizeof(long int) * 10 + 5);
    short int shift1 = strlen(halfsite1) * 2;
    short int shift2 = strlen(halfsite2) * 2;
    
    /* HETERODIMER */
    match_orientations[0] = halfsite1_value << shift1 | halfsite2_value;
    match_orientations[1] = halfsite2_value << shift2 | halfsite1_value;
    match_orientations[2] = Reverse_complement_sequence_value_li(halfsite1_value, strlen(halfsite1)) << shift1 | halfsite2_value;
    match_orientations[3] = halfsite1_value << shift1 | Reverse_complement_sequence_value(halfsite2_value, strlen(halfsite2));
    
    /* HOMODIMER1 */
    match_orientations[4] = halfsite1_value << shift1 | halfsite1_value;
    match_orientations[5] = Reverse_complement_sequence_value(halfsite1_value, strlen(halfsite1)) << shift1 | halfsite1_value;
    match_orientations[6] = halfsite1_value << shift1 | Reverse_complement_sequence_value(halfsite1_value, strlen(halfsite1));

    /* HOMODIMER2 */
    match_orientations[7] = halfsite2_value << shift2 | halfsite2_value;
    match_orientations[8] = Reverse_complement_sequence_value(halfsite2_value, strlen(halfsite2)) << shift2 | halfsite2_value;
    match_orientations[9] = halfsite2_value << shift2 | Reverse_complement_sequence_value(halfsite2_value, strlen(halfsite2));
    match_length = strlen(halfsite1) + strlen(halfsite2);
}

}
else break;
}
}

/* INITIALIZES VARIABLES WHOSE SIZE OR USE IS DEPENDENT ON ARGUMENTS */
__uint128_t left_position_value = 1;
left_position_value <<= ((Nlength-2) * 2);
    
long int ***align_score_histogram;
align_score_histogram = malloc(sizeof(long int *) * 3 + 5);
align_score_histogram[0] = malloc(sizeof(long int *) * 3 + 5);
align_score_histogram[1] = malloc(sizeof(long int *) * 3 + 5);
align_score_histogram[0][0] = malloc(sizeof(long int) * Nlength + 5);
align_score_histogram[0][1] = malloc(sizeof(long int) * Nlength + 5);
align_score_histogram[1][0] = malloc(sizeof(long int) * Nlength + 5);
align_score_histogram[1][1] = malloc(sizeof(long int) * Nlength + 5);    
for (counter = 0; counter < Nlength; counter++) {align_score_histogram[0][0][counter] = 0; align_score_histogram[1][0][counter] = 0; align_score_histogram[0][1][counter] = 0; align_score_histogram[1][1][counter] = 0;}
    
/* PRIMARY mask_ULL FOR EACH SEPARATE NUCLEOTIDE */
for (mask_ULL[1][0] = 3, counter = 0; counter < Nlength; counter++) mask_ULL[1][counter+1] = mask_ULL[1][counter] << 2;

    
	/* FLANK TOOL VARIABLES */
	char *searchstring;
	searchstring = malloc(1000);
	char *valuestring;
	valuestring = malloc(1000);
    

	short int *number_of_matches;
	number_of_matches = malloc(sizeof(short int) * 200 + 5);
	number_of_matches[0] = 0;
	number_of_matches[1] = 0;
	struct match *match;
	match = malloc(sizeof(struct match) * 3 + 5);
	match_init(&match[0], Nlength);
	match_init(&match[1], Nlength);
	match_init(&match[2], Nlength);
	struct match *dinucmatch;
	dinucmatch = malloc(sizeof(struct match) * 3 + 5);
	match_init(&dinucmatch[0], Nlength);
	match_init(&dinucmatch[1], Nlength);
	match_init(&dinucmatch[2], Nlength);
	struct count_connecting_matrix cm;
	count_connecting_matrix_init(&cm, "all hits connecting matrix", Nlength, 0);
    struct count_connecting_matrix **cms_to_heatmap;
    cms_to_heatmap = malloc(sizeof(struct count_connecting_matrix *) * 4 + 5);
    cms_to_heatmap[0] = &cm;
    
	short int iupac_query = 0;
	
	struct count_connecting_matrix two_hits_connecting_matrix;
	count_connecting_matrix_init(&two_hits_connecting_matrix, "two hit connecting matrix", Nlength, 0);	
	struct count_connecting_matrix two_hits_fold_connecting_matrix;
	count_connecting_matrix_init(&two_hits_fold_connecting_matrix, "fold two hit connecting matrix", Nlength, 0);
    
	struct count_pwm *all_hits_pwm;
	all_hits_pwm = malloc (sizeof(struct count_pwm) * 3 + 5);
	struct count_pwm *one_hit_pwm;
	one_hit_pwm = malloc (sizeof(struct count_pwm) * 3 + 5);
	struct count_pwm *one_hit_exponential;
	one_hit_exponential = malloc (sizeof(struct count_pwm) * 3 + 5);
    struct count_pwm *all_hits_exponential;
	all_hits_exponential = malloc (sizeof(struct count_pwm) * 3 + 5);
	short int reverse_strand_position;
	struct count_pwm ***two_hits_pwm;
	two_hits_pwm = malloc(sizeof(struct count_pwm *) * 3 + 5);
	two_hits_pwm[0] = malloc(sizeof(struct count_pwm *) * 5 + 5);
	two_hits_pwm[1] = malloc(sizeof(struct count_pwm *) * 5 + 5);
	for (counter = 0; counter < 5; counter++) 
	{
        two_hits_pwm[0][counter] = malloc (sizeof(struct count_pwm) * (Nlength * 2) + 5);
        two_hits_pwm[1][counter] = malloc (sizeof(struct count_pwm) * (Nlength * 2) + 5);	
	}
	
	
	for(file_number = 0; file_number < number_of_files; file_number++)
	{
        count_pwm_init(&all_hits_pwm[file_number], "empty", Nlength * 2, 0);
        count_pwm_init(&one_hit_pwm[file_number], "empty", Nlength * 2, 0);
        count_pwm_init(&one_hit_exponential[file_number], "empty", Nlength * 2, 0);
        count_pwm_init(&all_hits_exponential[file_number], "empty", Nlength * 2, 0);
        for (counter2 = 0; counter2 < 4; counter2++) for (counter3 = 0; counter3 < Nlength * 2; counter3++) count_pwm_init(&two_hits_pwm[file_number][counter2][counter3], "empty", Nlength * 2, 0);
	}
	
	/* QUERY PWM STRUCTURE */
	struct normalized_pwm qp;
	normalized_pwm_init(&qp, "empty", Nlength * 2, 0);

    /* MATCH FILTER PWM STRUCTURE */
	struct normalized_pwm mfp;
	normalized_pwm_init(&mfp, "empty", Nlength * 2, 0);

    /* DINUC MULTINOMIAL 2 PWM STRUCTURE */
	struct normalized_pwm mmp;
	normalized_pwm_init(&mmp, "empty", Nlength * 2, 0);
    
	/* BACKGROUND QUERY PWM STRUCTURE */
	struct normalized_pwm bqp;
	normalized_pwm_init(&bqp, "empty", Nlength * 2, 0);
	
	/* PRINT PWM STRUCTURE */
	struct normalized_pwm p;
	normalized_pwm_init(&p, "empty", Nlength * 2, 0);
    
	/* NORMALIZED BACKGROUND PWM STRUCTURES */
	struct normalized_pwm background_pwm[2];
	normalized_pwm_init(&background_pwm[0], "empty", Nlength * 2, 0.25);
	normalized_pwm_init(&background_pwm[1], "empty", Nlength * 2, 0.25);
    
	/* NORMALIZED CONNECTING MATRIX */
	struct normalized_connecting_matrix cp;
	normalized_connecting_matrix_init(&cp, "normalized connecting matrix", Nlength * 2, 0);
    
    /* ADM STRUCTURES */
    struct adjacent_dinucleotide_model unflanked_adm;
    struct adjacent_dinucleotide_model flanked_adm;
    struct adjacent_dinucleotide_model *flanked_adm_p;
    
	short int first; 
	short int last;
	short int first_match_position;
	short int spacing;
	short int orientation;
	__uint128_t first_sequence_value_ULL;
	char *pwm_name;
	double cut_off = 9;
	double filter_cut_off = 9;
	double multi_2_cut_off = 9;
    
	/* VARIABLES AND STRUCTURES FOR LOGO GENERATION */
	short int current_logo = 0;
	struct normalized_pwm *np;
	np = malloc(sizeof(struct normalized_pwm) * 200 + 5);
	struct normalized_pwm **np_p;
	np_p = malloc(sizeof(struct normalized_pwm *) * 200 + 5);
	for (counter = 0; counter < 200; counter++) 
	{
        normalized_pwm_init(&np[counter], "empty", Nlength * 2, 0);
        np_p[counter] = &np[counter];
	}
	long int total_number_of_two_hits;
	long int total_possible_spacings;
    
	
/* VARIABLES FOR MULTINOMIAL */
double swap = 0;
file_number = 0;
double lambda = 1;
short int no_insertions = 0;
short int multinomial_2_dinuc = 0;
long int background_kmer_count;
long int signal_kmer_count;
long int number_of_kmers = 0;
    
/* VARIABLES FOR INFORMATION CONTENT */
double background_info = 16;
double signal_info = 16;
double background_nonhit_info = -100;
double signal_nonhit_info = -100;
double fraction;	
double total_exp;
double total_obs;
    
double current_ic_score;
double max_ic_score = 0;
double ic_score = 0;
double delta_ic = 0;
short int max_first = 0;
short int max_second = 0;
short int end_trim = 15;
short int max_dinucleotide;
short int min_dinucleotide;
short int warning = 0;
double forward_pos_ic;
double reverse_pos_ic;
double min_fold;
double max_fold;
double current_fold;

short int load_pwm = 0;
short int load_adm = 0;
short int loaded_pwm_width = 0;
    
/* double **total_relative_deviation;
total_relative_deviation = malloc(sizeof(double *) * Nlength * 2 + 5);
for(counter = 0; counter < Nlength * 2; counter++) total_relative_deviation[counter] = malloc(sizeof(double) * Nlength * 2 + 5); */
double **uncentered_correlation;
uncentered_correlation = malloc(sizeof(double *) * Nlength * 2 + 5);
for(counter = 0; counter < Nlength * 2; counter++) uncentered_correlation[counter] = malloc(sizeof(double) * Nlength * 2 + 5);
long double sum_x_squared;
long double sum_y_squared;
long double sum_xy;

/* PWM ALIGN OPTION ENTIRE MAIN PROGRAM CODE */
if (pwm_align == 1 && argc >= 3 + firstnoncommandposition)
{
short int number_of_pwms = 2;
short int width = 0;
    
strcpy(searchstring, argv[1 + firstnoncommandposition]);
printf("\n\nMATRIX 1:\t%s", searchstring);
Load_pwm (&qp, searchstring, 0);
strcpy(qp.name, searchstring);    
Normalize_pwm(&qp);
// Log_ratio_pwm(&qp);

strcpy(tempstring, argv[2 + firstnoncommandposition]);
printf("\n\nMATRIX 2:\t%s", tempstring); 
Load_pwm (&mfp, tempstring, 0);
strcpy(mfp.name, tempstring); 
Normalize_pwm(&mfp);
// Log_ratio_pwm(&mfp);

if (argc >= 4 + firstnoncommandposition)
{
number_of_pwms = 3;
strcpy(valuestring, argv[3 + firstnoncommandposition]);
printf("\n\nMATRIX 3:\t%s", valuestring); 
Load_pwm (&p, valuestring, 0);
strcpy(p.name, valuestring); 
Normalize_pwm(&p);
// Log_ratio_pwm(&p);        
}

width = qp.width;
if (number_of_pwms == 3 && p.width > qp.width) width = p.width;
    
Add_flanks(&mfp, width, 0.25);

struct normalized_pwm rc_qp; 
normalized_pwm_init (&rc_qp, qp.name, qp.width, 0);
Reverse_complement_normalized_pwm(&rc_qp, &qp);    
 
struct normalized_pwm *qp_p;
qp_p = Align_pwm(&qp, &rc_qp, &mfp);
mfp.match.position = mfp.match.position + width; 


struct normalized_pwm *p_p;
if (number_of_pwms == 3)
{
struct normalized_pwm rc_p; 
normalized_pwm_init (&rc_p, p.name, p.width, 0);
Reverse_complement_normalized_pwm(&rc_p, &p);  

p_p = Align_pwm(&p, &rc_p, &mfp);    

}
    
strcat (searchstring, "_");
strcat (searchstring, tempstring);
if (number_of_pwms == 3) 
{
strcat (searchstring, "_");
strcat (searchstring, valuestring);        
}
strcat (searchstring, ".svg");
 
Load_pwm (&mfp, tempstring, 1);

np_p[0] = qp_p;
np_p[1] = &mfp;
if (number_of_pwms == 3) 
{
    if (contacts == 0) np_p[2] = p_p;
    else  
    {
    np_p[1] = p_p;
    np_p[2] = &mfp;
    }
    nocall = 2;
}
else nocall = 1;

noname = 1;
   
Svg_logo(searchstring, number_of_pwms, np_p, '\0', '\0', '\0', '\0', '\0', '\0', 0, 0);
    
printf("\n\nPWM1 best match score to PWM2 is %.2f position is %i strand %i", (*qp_p).match.score, (*qp_p).match.position - width + 1, (*qp_p).match.strand);
if (number_of_pwms == 3)  {
        printf("\nPWM3 best match score to PWM2 is %.2f position is %i strand %i\n", (*p_p).match.score, (*p_p).match.position - width + 1, (*p_p).match.strand);

short int pwm_position;
struct normalized_pwm *first_match_pwm_p;
struct normalized_pwm *second_match_pwm_p;
if ((*p_p).match.position < (*qp_p).match.position) {first_match_pwm_p = p_p; second_match_pwm_p = qp_p;} else {first_match_pwm_p = qp_p; second_match_pwm_p = p_p;}

/* PRINTS ALIGNMENT SCORES BY BASE */
char *overlapcall[] = {"none", "OVERLAP"};
char *strengthcall[] = {"strong", "weak"};
printf("\nSIDE\tNAME\tINPUT_PWM_POS\tSCORE\tOVERLAP\tALIGNMENT");
for(counter = 0; counter < (*first_match_pwm_p).width; counter++) printf("\nleft\t%s\t%i\t%.2f\t%s\t%s\talign_score", (*first_match_pwm_p).name, (*first_match_pwm_p).original_position[counter] + 1, (*first_match_pwm_p).position_score[counter], overlapcall[counter + (*first_match_pwm_p).match.position >= (*second_match_pwm_p).match.position], strengthcall[(*first_match_pwm_p).position_score[counter] < 0.2]);
for(counter = 0; counter < (*second_match_pwm_p).width; counter++) printf("\nright\t%s\t%i\t%.2f\t%s\t%s\talign_score", (*second_match_pwm_p).name, (*second_match_pwm_p).original_position[counter] + 1, (*second_match_pwm_p).position_score[counter], overlapcall[(*second_match_pwm_p).match.position + counter < (*first_match_pwm_p).match.position + (*first_match_pwm_p).width], strengthcall[(*first_match_pwm_p).position_score[counter] < 0.25]);
    
/* PRINTS TOTAL GAP WIDTH */
printf("\n\nComposite site total gap width: \t%3i", abs((*p_p).match.position - (*qp_p).match.position) - (*first_match_pwm_p).width);
    
if (contacts == 1)
{
/* PRINTS CORE GAP WIDTH */
for(counter = 0; counter < (*second_match_pwm_p).width; counter++) if ((*second_match_pwm_p).fraction[8][counter]+(*second_match_pwm_p).fraction[9][counter]+(*second_match_pwm_p).fraction[10][counter]+(*second_match_pwm_p).fraction[11][counter] != 0) break;
for(counter2 = (*first_match_pwm_p).width-1; counter2 > 0; counter2--, counter++) if ((*first_match_pwm_p).fraction[8][counter2]+(*first_match_pwm_p).fraction[9][counter2]+(*first_match_pwm_p).fraction[10][counter2]+(*first_match_pwm_p).fraction[11][counter2] != 0) break;
printf("\nComposite site core gap width:  \t%3li", abs((*p_p).match.position - (*qp_p).match.position) - (*first_match_pwm_p).width + counter);
    
/* PRINTS BACKBONE GAP WIDTH */
for(counter = 0; counter < (*second_match_pwm_p).width; counter++) if ((*second_match_pwm_p).fraction[4][counter]+(*second_match_pwm_p).fraction[5][counter]+(*second_match_pwm_p).fraction[6][counter]+(*second_match_pwm_p).fraction[7][counter] + (*second_match_pwm_p).fraction[12][counter]+(*second_match_pwm_p).fraction[13][counter]+(*second_match_pwm_p).fraction[14][counter]+(*second_match_pwm_p).fraction[15][counter] != 0) break;
for(counter2 = (*first_match_pwm_p).width-1; counter2 > 0; counter2--, counter++) if ((*first_match_pwm_p).fraction[4][counter2]+(*first_match_pwm_p).fraction[5][counter2]+(*first_match_pwm_p).fraction[6][counter2]+(*first_match_pwm_p).fraction[7][counter2] + (*first_match_pwm_p).fraction[12][counter2]+(*first_match_pwm_p).fraction[13][counter2]+(*first_match_pwm_p).fraction[14][counter2]+(*first_match_pwm_p).fraction[15][counter2] != 0) break;
printf("\nComposite site backbone gap width:  \t%3li", abs((*p_p).match.position - (*qp_p).match.position) - (*first_match_pwm_p).width + counter);

/* RELOADS ORIGINAL MIDDLE (COMPOSITE) PWM */
Load_pwm (&mfp, tempstring, 0);
/* ADDS UP CONTACTS */
short int position_on_left_pwm;
short int position_on_right_pwm; 
for (nucleotide = 4; nucleotide < 16; nucleotide++) for(counter = 0; counter < mfp.width; counter++) 
{
mfp.fraction[nucleotide][counter] = 0;
position_on_left_pwm = counter - (*first_match_pwm_p).match.position + width + 1;
position_on_right_pwm = counter - (*second_match_pwm_p).match.position + width + 1;
if (position_on_left_pwm >= 0 && position_on_left_pwm < (*first_match_pwm_p).width) mfp.fraction[nucleotide][counter] += (*first_match_pwm_p).fraction[nucleotide][position_on_left_pwm];
if (position_on_right_pwm >= 0 && position_on_right_pwm < (*second_match_pwm_p).width) mfp.fraction[nucleotide][counter] += (*second_match_pwm_p).fraction[nucleotide][position_on_right_pwm];
}
/* PRINTS COMPOSITE PWM WITH TOTAL CONTACTS */
printf("\n\nOriginal composite PWM (2) with contacts from the individual PWMs (1,3)");
for (nucleotide = 0; nucleotide < 16; nucleotide++)
    {
    printf ("\n%s_Composite_PWM\t%s", mfp.name, pwm_row_ids[nucleotide]);
    for(pwm_position = 0; pwm_position < mfp.width; pwm_position++) 
    printf("\t%.0f", mfp.fraction[nucleotide][pwm_position]);
    }
    
/* GENERATES CORE-TO-CORE PWM */
for(counter = 0; counter < (*first_match_pwm_p).width; counter++) if ((*first_match_pwm_p).fraction[8][counter]+(*first_match_pwm_p).fraction[9][counter]+(*first_match_pwm_p).fraction[10][counter]+(*first_match_pwm_p).fraction[11][counter] != 0) break;
for(counter2 = (*second_match_pwm_p).width-1; counter2 > 0; counter2--) if ((*second_match_pwm_p).fraction[8][counter2]+(*second_match_pwm_p).fraction[9][counter2]+(*second_match_pwm_p).fraction[10][counter2]+(*second_match_pwm_p).fraction[11][counter2] != 0) break;
printf("\n\nCore-to-core PWM start position:\t%3li\n", (*first_match_pwm_p).match.position + counter - width + 1);
printf("Core-to-core PWM end position:   \t%3li\n", (*second_match_pwm_p).match.position + (*second_match_pwm_p).width - 1 - counter2 - width);
printf("\nCore-to-core PWM");

for (nucleotide = 0; nucleotide < 16; nucleotide++)
    {
    printf ("\n%s_CTC_PWM\t%s", mfp.name, pwm_row_ids[nucleotide]);
    for(pwm_position = (*first_match_pwm_p).match.position + counter - width; pwm_position < (*second_match_pwm_p).match.position + (*second_match_pwm_p).width - 1 - counter2 - width; pwm_position++) printf("\t%.0f", mfp.fraction[nucleotide][pwm_position]);
    }
printf("\n");

}
}
       exit(0);
}
    
    
/* LOGO OPTION ENTIRE MAIN PROGRAM CODE */
if (svg_only == 1 && argc >= 2 + firstnoncommandposition)
{
strcpy(searchstring, argv[1 + firstnoncommandposition]);
qp.negative_values_allowed = negative_values_allowed;
    
Load_pwm (&qp, searchstring, 1);
strcpy(qp.name, searchstring);
    
    if (difference_logo == 1 && argc >= 3 + firstnoncommandposition)
    {
        offset = 0;
        if(argc >= 4 + firstnoncommandposition) offset = atoi(argv[3 + firstnoncommandposition]);
        strcpy(tempstring, argv[2 + firstnoncommandposition]);
        Load_pwm (&mfp, tempstring, 1);
        Subtract_normalized_pwm(&qp,&mfp, offset);
        qp.negative_values_allowed = 1; 
        strcat(qp.name, "_minus_");
        strcat(qp.name, tempstring);
    }

np_p[0] = &qp;
strcpy(tempstring, qp.name);
strcpy(searchstring, qp.name);
strcat(tempstring, ".svg");
    if (argc >= 3 + firstnoncommandposition && difference_logo == 0) {strcpy(tempstring, argv[2 + firstnoncommandposition]); strcpy(searchstring, argv[2 + firstnoncommandposition]);}
Svg_logo(tempstring, 1, np_p, '\0', '\0', '\0', '\0', '\0', '\0', 0, 0);
char *system_command;
system_command = malloc (2000);
strcpy(system_command, "convert ");
strcat(system_command, tempstring);
strcat(system_command, " ");
strcat(system_command, searchstring);
strcat(system_command, ".png");
printf("\n%s\n", system_command);
system(system_command);
exit(0);
} 
/* END OF LOGO OPTION */

/* DIST OPTION ENTIRE MAIN PROGRAM CODE */
    if (modeldistance == 1 && argc >= 3 + firstnoncommandposition)
    {
        short int gapped = 0;
        short int gappedkmerlength = 4;
        short int gap;
        long double onetwo;
        long double twoone;
        long double gapped_correlation;
        long double firstscore;
        long double secondscore;
        double **firstpwm;
        double **secondpwm;
        signed short int firstlength;
        signed short int secondlength;
        signed short int secondpwmwidth;
        signed short int firstpwmwidth;
        signed short int max_pwm_width;
        long long int kmervalue;
        long long int kmervalue2;

        
        strcpy(searchstring, argv[1 + firstnoncommandposition]);
        printf("\n\nNORMALIZED MATRIX 1:\t%s", searchstring);
        Load_pwm (&qp, searchstring, 0);
        firstlength = qp.width;
        Normalize_pwm(&qp);
        printf("\n\nLOG MATRIX 1:\t%s", searchstring);
        Log_ratio_pwm(&qp);
        
        strcpy(tempstring, argv[2 + firstnoncommandposition]);
        printf("\n\nNORMALIZED MATRIX 2:\t%s", tempstring); 
        Load_pwm (&mfp, tempstring, 0);
        secondlength = mfp.width; 
        Normalize_pwm(&mfp);
        printf("\n\nLOG MATRIX 2:\t%s", tempstring);
        Log_ratio_pwm(&mfp);

        if(firstlength > secondlength) too_long_kmer = firstlength;
        else too_long_kmer = secondlength;
        max_pwm_width = too_long_kmer;
        
        if (argc >= 4 + firstnoncommandposition) 
        {
        if (strcmp(argv[3 + firstnoncommandposition], "gapped") == 0) gapped = 1;
        else firstlength = atoi(argv[3 + firstnoncommandposition]);   
        }
        if (argc >= 5 + firstnoncommandposition) 
        {
        secondlength = atoi(argv[4 + firstnoncommandposition]);
        secondlength = secondlength - secondlength % 2;
        gappedkmerlength = secondlength / 2;
        }
        if (gapped == 1) 
        {
        firstlength = secondlength / 2;
        if (max_pwm_width < gappedkmerlength * 2) gappedkmerlength = max_pwm_width / 2;
        }

        if(firstlength > secondlength) too_long_kmer = firstlength;
        else too_long_kmer = secondlength;
        
        Add_flanks(&qp, too_long_kmer-2, -0.60206);
        Add_flanks(&mfp, too_long_kmer-2, -0.60206);

        firstpwm = qp.fraction;
        firstpwmwidth = qp.width;
        
        secondpwm = mfp.fraction;
        secondpwmwidth = mfp.width;
        
        long long int last_kmer_value;
        
        sum_x_squared = 0;
        sum_y_squared = 0;
        sum_xy = 0;
        long int different = 0;
        long int same = 0;
        
        last_kmer_value = pow(4,firstlength);
        #pragma omp parallel for schedule(static) reduction(+:sum_xy, sum_x_squared, sum_y_squared) default(shared) private(kmervalue,firstscore,secondscore)
        for(kmervalue = last_kmer_value; kmervalue > 0; kmervalue--)
        {
        firstscore = fastKmerscore(firstpwm, firstpwmwidth, kmervalue, firstlength);
        secondscore = fastKmerscore(secondpwm, secondpwmwidth, kmervalue, firstlength);
        sum_x_squared += pow(10,firstscore+firstscore);
        sum_y_squared += pow(10,secondscore+secondscore);
        sum_xy += pow(10,firstscore+secondscore);
       /*     if(firstscore != secondscore) different++;
            else same++; */
        }
        onetwo = sum_xy/(sqrt(sum_y_squared)*sqrt(sum_x_squared));

        /* printf("\ndifferent %li, same %li\n", different, same); */
               
        if(firstlength != secondlength)
        {
        sum_x_squared = 0;
        sum_y_squared = 0;
        sum_xy = 0;
        
        last_kmer_value = pow(4,secondlength);
        #pragma omp parallel for schedule(static) reduction(+:sum_xy, sum_x_squared, sum_y_squared) default(shared) private(kmervalue,firstscore,secondscore)
        for(kmervalue = 0; kmervalue < last_kmer_value; kmervalue++)
        {
        firstscore = fastKmerscore(firstpwm, firstpwmwidth, kmervalue, secondlength);
        secondscore = fastKmerscore(secondpwm, secondpwmwidth, kmervalue, secondlength);
        sum_x_squared += pow(10,firstscore+firstscore);
        sum_y_squared += pow(10,secondscore+secondscore);
        sum_xy += pow(10,firstscore+secondscore);
        }
        twoone = sum_xy/(sqrt(sum_y_squared)*sqrt(sum_x_squared));
        }
        else twoone = onetwo;
        

        /* GAPPED kmers */
        if (gapped == 1)
        {
        
        sum_x_squared = 0;
        sum_y_squared = 0;
        sum_xy = 0;
        
        last_kmer_value = pow(4,gappedkmerlength);
        
#pragma omp parallel for schedule(static) reduction(+:sum_xy, sum_x_squared, sum_y_squared) default(shared) private(kmervalue,kmervalue2, gap, firstscore,secondscore)
        for(gap = max_pwm_width - 2 * gappedkmerlength; gap >= 0; gap--)
        {
        for(kmervalue = last_kmer_value; kmervalue > 0; kmervalue--)
        {
        for(kmervalue2 = last_kmer_value; kmervalue2 > 0; kmervalue2--)
        {

            firstscore = fastgappedKmerscore(firstpwm, firstpwmwidth, kmervalue, kmervalue2, gappedkmerlength, gap);
            secondscore = fastgappedKmerscore(secondpwm, secondpwmwidth, kmervalue, kmervalue2, gappedkmerlength, gap);
            sum_x_squared += pow(10,firstscore+firstscore);
            sum_y_squared += pow(10,secondscore+secondscore);
            sum_xy += pow(10,firstscore+secondscore);
        }
        }
        }
        gapped_correlation = sum_xy/(sqrt(sum_y_squared)*sqrt(sum_x_squared));
            
        printf("\n\nSimilarity (cosine angle correlation) \n%s to %s (%i-mer)\t%Lg\n%s to %s (%i-mer)\t%Lg\n%s to %s (gapped %i-mer)\t%Lg\n", searchstring, tempstring, firstlength, onetwo, tempstring, searchstring, secondlength, twoone, searchstring, tempstring, gappedkmerlength * 2, gapped_correlation);
        }
        else
        {
        printf("\n\nSimilarity (cosine angle correlation) \n%s to %s (%i-mer)\t%Lg\n%s to %s (%i-mer)\t%Lg\n", searchstring, tempstring, firstlength, onetwo, tempstring, searchstring, secondlength, twoone);
        }
        
        t1 = time('\0');
        printf ("\nTime: %ld seconds\n", (long) (t1 - t0));
        exit(0);
    } 
  
/* END OF DIST OPTION */
    
/* CHECKS IF WINFLAT IS INSTALLED */
if (print_p_values == 1) if (Winflat(10, 10, 100, 100) == 2) {printf ("\nERROR: Winflat not installed in path\n"); argc = 1; fflush(stdout);}

/* PARSES ARGUMENTS */
if ((argc >= 5 + firstnoncommandposition || (flank == 1 && argc >= 3 + firstnoncommandposition)) && strcmp("HELP", argv[1]) != 0 && strcmp("help", argv[1]) != 0 && strcmp("-h", argv[1]) != 0 && strcmp("-H", argv[1]) != 0)
{
/* FLANK TOOL */
if (flank == 1)
{
strcpy(file_name[0], argv[1 + firstnoncommandposition]);
if (strcmp(file_name[0], "-") == 0) {even_background = 1; file_number++;}
strcpy(file_name[1], argv[2 + firstnoncommandposition]);
strcpy(searchstring, argv[3 + firstnoncommandposition]);
/* CHECKS IF SEARCHSTRING IS A NUCLEOTIDE SEQUENCE, IF NOT, ASSUMES IT IS A FILENAME FOR PWM */
for (counter = 0; counter < strlen(searchstring); counter++) if (searchstring[counter] != 'A' & searchstring[counter] != 'C' & searchstring[counter] != 'G' & searchstring[counter] != 'T') flank_with_pwm = 1;
if (flank_with_pwm == 1)
{
iupac_query = 1;
flank_with_pwm = 1;
for (counter = 0; counter < strlen(searchstring); counter++) 
{
for (counter2 = 0; (counter2 < strlen(nucleotide_iupac)) & (searchstring[counter] != nucleotide_iupac[counter2]); counter2++);
if (counter2 == strlen(nucleotide_iupac)) {iupac_query = 0; flank_with_pwm = 1; break;}
}
}

if (argc >= 5 + firstnoncommandposition)
{
strcpy(valuestring, argv[4 + firstnoncommandposition]);
if (valuestring[strlen(valuestring)-1] == '%') {valuestring[strlen(valuestring)-1] = '\0'; minimum_kmer_count = atof(valuestring) / 100;}
else minimum_kmer_count = atof(valuestring);
}
if (cut_off != 0) cut_off = log10(minimum_kmer_count);
else cut_off = log10(pseudocount);
if (argc >= 6 + firstnoncommandposition && flank_with_pwm == 1) minimum_kmer_count = atof(argv[5 + firstnoncommandposition]);
pwm_name = searchstring;
}
/* NMER COUNT TOOL */
else
{
strcpy(file_name[0], argv[1 + firstnoncommandposition]);
strcpy(file_name[1], argv[2 + firstnoncommandposition]);
shortest_kmer = atoi(argv[3 + firstnoncommandposition]);
too_long_kmer = atoi(argv[4 + firstnoncommandposition])+1;
minimum_kmer_count = atof(argv[5 + firstnoncommandposition]);
if (kmer_table == 1 && argc >= 6 + firstnoncommandposition) 
{
strcpy (searchstring, argv[6 + firstnoncommandposition]);

/* CHECKS IF SEARCHSTRING IS A NUCLEOTIDE SEQUENCE, IF NOT, ASSUMES IT IS A FILENAME FOR PWM */
for (counter = 0; counter < strlen(searchstring); counter++) if (searchstring[counter] != 'A' & searchstring[counter] != 'C' & searchstring[counter] != 'G' & searchstring[counter] != 'T') load_pwm = 1;
if (load_pwm == 1)
{
iupac_query = 1;
load_pwm = 0;
for (counter = 0; counter < strlen(searchstring); counter++) 
{
for (counter2 = 0; (counter2 < strlen(nucleotide_iupac)) & (searchstring[counter] != nucleotide_iupac[counter2]); counter2++);
if (counter2 == strlen(nucleotide_iupac)) {iupac_query = 0; load_pwm = 1; break;}
}
}

if (load_pwm == 1) 
{
if (strstr(searchstring+strlen(searchstring)-4, ".adm") == 0)
{
Load_pwm (&qp, searchstring, 1);
Log_ratio_pwm(&qp);
}
else
{
    load_adm = 1;
    printf ("\nloads ADM");
    adjacent_dinucleotide_model_init(&unflanked_adm,"unflanked_adm", Nlength);
    adjacent_dinucleotide_model_init(&flanked_adm,"unflanked_adm", Nlength+2*(too_long_kmer-2));
    Load_ADM (&unflanked_adm, searchstring);
    Log_ratio_ADM (&unflanked_adm);
    PWM_from_ADM (&unflanked_adm, &qp);
    Add_flanks_to_ADM(&flanked_adm, &unflanked_adm, too_long_kmer-2, -0.60206);
    /* short int *jj;
    Kmerscore_ADM(&flanked_adm, 5762,8, jj);
    Kmerscore_ADM(&flanked_adm, 24763,8, jj);
    Kmerscore_ADM(&flanked_adm, 24762,8, jj);
    exit(0); */
}
    }
else 
{
Multinomial_pwm(&qp, searchstring);
if (match_filter != 0) Multinomial_pwm(&mfp, searchstring);
}
loaded_pwm_width = qp.width;
Add_flanks(&qp, too_long_kmer-2, -0.60206);
}
}
}
else
/* SOMETHING WRONG, PRINT INSTRUCTIONS */
{
printf ("Usage: \n./spacek --tool -options [background file name] [sample file name] [shortest kmer] [longest kmer] [minimum count or fold (with option -s) for printing] [PWM file (only with option -kmer)] \n\tCount number of occurrences of subsequences (kmers) of the indicated size. Indicates if count is higher than sequences within 1 Huddinge distance.\n\tWrites also a graphical summary to an svg file. \n\n  Other tools:\n  --f [background file name or - for none] [sample file name] [kmer sequence (consensus or IUPAC)] [minimum incidence] or\n  --pwm [background file name or - for none] [sample file name] [pwm file name] [score cut-off e.g. 0.1 or 10%%] [minimum incidence] \n\tFind occurrences of input sequence / pwm and print flanking preference, spacing and orientation information. \n\tGraphical summary of data is also written to an svg file\n\n  --logo [PWM file name] [output file name (opt)]\n\tGenerate SVG logo file from input PWM, also generates .png if convert is installed in path\n  --difflogo [PWM1 file name] [PWM2 file name] [offset (opt)] \n\tGenerate difference logo (PWM1-PWM2)\n  --dist [PWM1 file name] [PWM2 file name] ['gapped' or kmer1 length (opt)] [length of gapped kmer or kmer2 (opt)] \n\tCalculate uncentered correlation between maximum scores from PWM1 and PWM2 for all kmers. Default kmers are ungapped, default lengths are PWM widths (8 for gapped)\n  --pwmalign [PWM1 file name] [PWM2 file name] [PWM3 file name (optional)]\n\tAlign PWMs and generate a svg with logos aligned to PWM2, with boxes indicating similar base positions.\n\tOption -contacts shows also base contacts of PWM1 and PWM3 on PWM2.\n\nOptions:\n\nLogo and sequence formatting:\n  -contacts\tShow hydrogen bond contacts (PWMs where contact information is included must be used)\n  -paths\tGenerate PWM logos using paths and not fonts\n  -neg\tAllow negative values in PWM (using paths)\n  -nocall\t do not include kmer-based prediction calls into logos\n  -noname\tDo not include filename in the logo\n  -barcodelogo\tGenerate PWM barcode logos using colored rectangles\n  -heightscaledbars\tScale height of barcode logos based on nucleotide frequency\n  -maxheightscaledbars\tScale height of barcode logos based on maximum nucleotide frequency at position\n  -colorscaledbars\tScale color of barcode logos based on nucleotide frequency\n  -label\tInclude white IUPAC label to barcodelogo (cutoff = more than 50%% of max)\n  -rna\tInput sequence is from RNA, use uracil (U) instead of thymidine (T) in sequences and logos\n  -CpG\tHighlight CpGs in logos if they occur with lower (opaque) or higher (black outline) frequency in sample file than in control file \n\nInput sequences:\n  -14N\tSet sequence length to 14N (default is 20N)\n  -30N\tSet sequence length to 30N\n  -40N\tSet sequence length to 40N\n  -u\tUse only unique input sequences\n\nKmer counting:\n  -nogaps\tCount only kmers without gaps\n  -allgaps\tCount kmers with gaps in any position (default is only middle 1 or 2 positions)\n  -longer≈[kmer_length _difference_cutoff for local max] (default = 0.4)]\n  -iupac≈[cutoff for making a base to an iupac] (default = 0.25 of maximum base)]\n\nMatching and PWM generation:\n  -m=[number]\tGenerate PWM using multinomial [number] distribution\n  -lim=[position followed by strand (one or both must be given)]\t Show only hits at indicated strand and/or position (e.g. -lim=F, -lim=5 or -lim=4F)\n  -exc=[position][strand],[position][strand]‚..,[position][strand]\tExclude indicated positions and strands (e.g. -exc=4,3F,6R)\n  -both\tCount both instances of palindromic hits\n  -bothifnotequal\tCount both instances of palindromic hits if hit scores are not equal\n  -forwardonly\tCount only forward instance of palindromic hits (default counts hit from strand with better score, if score equal alternates between strands)\n  -reverseonly\tCount only reverse instance of palindromic hits\n  -forwardifequal\tCount only forward instance of palindromic hits if scores equal (default alternates between strands)\n  -reverseifequal\tCount only reverse instance of palindromic hits if scores equal\n  -e\tUse even background instead of background from file\n  -mono\tUse mononucleotide background instead of multinomial background\n  -fk=[position] count 4-mers at given position relative to match (e.g. -fk=-4 counts 4-mers that precede search sequence)\n  -mf=[multinomial] match filter: consider only sequences with one hit to indicated multinomial as true one hits\n  -kl=[kmer length] kmer length used to calculate lambda (default 8)\n\nText output:\n  -q\tPrint frequencies instead of counts\n  -n\tPrint nucleotide counts for full sequences\n  -s\tPrint values for different gap lengths on same line\n  -c\tPrint raw counts from both input files\n  -p\tPrint p-values (Winflat program needs to be installed in path)\n  -i\tPrint incidence of all input sequences\n  -kmer\tPrint kmer count table with scores for each kmer against the PWM indicated\n  -editdist[≈min_local_max_percent_cutoff (default = 10)] (for --f option)\tPrint table with max count kmer pairs, all local maxima, and all cloud counts for each edit distance\n  -dimer\tPrint counts for dimeric sequences only\n  -ic\tOutput information content for all spacings\n  -dinuc\tPrint dinucleotide data to output and svg (counting dinucleotides uses multinomial 2)\n  -x\tExtended output for debugging\n\nSvg output:\n  -o=[output file name] base name of svg output file(s)\n  -match=[kmer1],[kmer2]\tInclude spacing and orientation heatmap for these kmers in kmer summary svg\n  -hrows≈[number of rows for kmer svg summary heatmaps (default = 20)]\n  -xyplot\tGenerate scatterplot of kmer counts for the shortest kmer, with CpG-containing kmers indicated in red\n\n ");

    
exit(0);
} 
 
	/* struct normalized_pwm rc; */
 
    /* INITIALIZES ALIGN SCORES */
    struct alignscore *all_hits_align_scores;
    all_hits_align_scores = malloc(sizeof(struct alignscore) * 3 + 5);
    alignscore_structure_init(&all_hits_align_scores[0]);
    alignscore_structure_init(&all_hits_align_scores[1]);
    struct alignscore current_align_score;
    alignscore_structure_init(&current_align_score);
    
	/* INITIALIZES DINUCLEOTIDE PROPERTIES STRUCTURES IF DINUCLEOTIDE PROPERTIES ARE ANALYZED */
	struct dinucleotide_properties di;
	struct dinucleotide_properties_matrix *one_hit_di;
	struct base_dependency_matrix *dep;
	struct base_dependency_matrix expected_dinucleotides;
		
	if (dinucleotide_properties == 1)
	{
	dinucleotide_properties_init(&di);
	one_hit_di = malloc(sizeof(struct dinucleotide_properties_matrix) * (number_of_files + 1) + 5);
	dep = malloc(sizeof(struct base_dependency_matrix) * (number_of_files + 1) + 5);
	dinucleotide_properties_matrix_init(&one_hit_di[0], &di, "dinucleotide1", Nlength * 2, 0, query_sequence_length);
	dinucleotide_properties_matrix_init(&one_hit_di[1], &di, "dinucleotide2", Nlength * 2, 0, query_sequence_length);
	base_dependency_matrix_init(&dep[0], "Background_dinucleotide_dependencies", Nlength * 2);
	base_dependency_matrix_init(&dep[1], "Dinucleotide_dependencies", Nlength * 2);
	base_dependency_matrix_init(&expected_dinucleotides, "Expected_dinucleotides", Nlength * 2);
	}
	
	
	/* CALCULATES MULTINOMIAL PWM IF MULTINOMIAL INDICATED OR LOADS PWM IF FLANK WITH PWM = 1 */
	if (flank_with_pwm == 1 || (multinomial != 0 && kmer_table == 0))
	{
	/* shortest_kmer = 8; */
	too_long_kmer = shortest_kmer+1;
	if (query_sequence_length < 12) kmer_count = 1;
	no_insertions = 1;
	flank_with_pwm = 1;
	/* printf("\nLoads Pwm: %s", pwm_name); */
	if (multinomial == 0 && iupac_query == 0) 
	{
	Load_pwm (&qp, pwm_name, 1);
	}
	else 
	{
    if (multinomial < 2 && dinucleotide_properties == 1) multinomial_2_dinuc = 1;
	if (iupac_query == 1) 
    {
    Iupac_to_pwm(&qp, searchstring); /* CALCULATES PWM FOR IUPAC */
    if (match_filter != 0) Iupac_to_pwm(&mfp, searchstring); /* CALCULATES PWM FOR MATCH FILTER */
    if (multinomial_2_dinuc == 1) Iupac_to_pwm(&mmp, searchstring); /* CALCULATES PWM FOR DINUC MULTINOMIAL 2 FILTER */ 
    }
    else 
    {
    Multinomial_pwm(&qp, searchstring); /* CALCULATES PWM FOR MULTINOMIAL */
    if (match_filter != 0) Multinomial_pwm(&mfp, searchstring); /* CALCULATES PWM FOR MATCH FILTER */
    if (multinomial_2_dinuc == 1) Multinomial_pwm(&mmp, searchstring); /* CALCULATES PWM FOR DINUC MULTINOMIAL 2 FILTER */ 
	}
    }
	query_sequence_length = qp.width;
	/* ADDS PWMs TO LOGO LIST */
	Copy_normalized_pwm(&np[current_logo], &qp);
	strcpy (np[current_logo].name, pwm_name);
	
		
	if (multinomial == 0 && iupac_query == 0) Log_fold_affinity_pwm(&qp);
	else 
	{
	cut_off = 0 - multinomial - 0.0001;
    multi_2_cut_off = 0 - 2 - 0.0001;
	filter_cut_off = 0 - match_filter - 0.0001;
    sprintf(tempstring, " - multinomial %i", multinomial);
	strcat(np[current_logo].name, tempstring);
	/* Reverse_complement_normalized_pwm(&rc, &qp);
	for (counter = 0; counter < 4; counter++) 
	{
	printf("\n");
 	for (counter2 = 0; counter2 < rc.width; counter2++) printf ("\t%.1f", rc.fraction[counter][counter2]); 
	} */
	}
	current_logo++;
	Invert_pwm(&bqp, &qp, query_sequence_length / 3);
	}
	
	/* CALCULATES QUERY SEQUENCE VALUE IF FLANK IS USED */
	if (flank == 1)
	{
	if ((multinomial ^ flank_with_pwm) == 0 && iupac_query == 0) 
	{
	query_sequence_length = strlen(searchstring);
	for(query_sequence_value = 0, position = 0, position_value = pow(4,query_sequence_length-1); position < query_sequence_length; position++, position_value /= 4)
	{
		for (nucleotide_value = 0; nucleotide_value < 4 && searchstring[position] != forward[nucleotide_value]; nucleotide_value++);
		if(nucleotide_value == 4) {printf("\nERROR IN QUERY SEQUENCE\n"); exit (1);}
		query_sequence_value += position_value * nucleotide_value;
	}
	}
	/* TESTS IF QUERY SEQUENCE IS A PALINDROME */
        searchstring_is_palindrome = Is_this_string_iupac_palindrome(searchstring, -1 , -1); 
        /* (Is_this_sequence_dimer(query_sequence_value, query_sequence_length) == 4); */
	if (searchstring_is_palindrome == 1) printf("\nPalindrome");
	if (extendedoutput == 1) printf("\nsequence_value_ULL = %ju", (uintmax_t) query_sequence_value); 
	}
	

	if (dinucleotide_properties == 1) one_hit_di[1].query_sequence_length = query_sequence_length;

	
			
	/* ALLOCATES MEMORY AND CLEARS INDEX IF ONLY UNIQUE SEQUENCES USED */
	if(remove_non_unique == 1)
	{
		sorted_list[0] = malloc(max_number_of_sequences * sizeof(struct sequence_incidence_table) + 10);
		sorted_list[1] = malloc(max_number_of_sequences * sizeof(struct sequence_incidence_table) + 10);
	}
	
	/* GENERATES mask_ULLS FOR EXTRACTION OF EACH KMER STRING */
	/* GENERATES KMER mask_ULLS */
	for(current_kmer_length = 2; current_kmer_length <= Nlength; current_kmer_length++)
	{
		for (start_position = 0; start_position < Nlength-current_kmer_length; start_position++)
		{
			for (mask_ULL[current_kmer_length][start_position]=mask_ULL[1][start_position], position = start_position+1; position < current_kmer_length+start_position; position++) {mask_ULL[current_kmer_length][start_position] += mask_ULL[1][position];
			/* printf("\n%i\t%i\t%i", current_kmer_length, start_position, mask_ULL[current_kmer_length][start_position]); */}
		}
	}
	
	/* GENERATES HIGH AND LOW mask_ULLS FOR DELETIONS */
	for (lowmask_ULL[0] = mask_ULL[1][0], position = 1; position < Nlength-2; position++) lowmask_ULL[position] = lowmask_ULL[position-1]+mask_ULL[1][position];
	for (highmask_ULL[Nlength-2] = mask_ULL[1][Nlength-2], position = Nlength-3; position > 0; position--) highmask_ULL[position] = highmask_ULL[position+1]+mask_ULL[1][position];

	/* INITIALIZES RESULT TABLE */
	long int *****results;
    long long int variable_length;
	results = malloc(sizeof(long int *) * 3 + 5);

	/* long int *results[number_of_files][too_long_kmer+1][too_long_kmer+2][Nlength-shortest_kmer+3]; */
	long int *current_value;
	short int kmer_length_size;
	short int gap_position_size;
	/* printf("\n\nsizeof int %i\t%i", sizeof(current_value), sizeof(current_value[0])); */
	/* clears results */
	for(counter = 0; counter < number_of_files; counter++)
	{
	results[counter] = malloc(sizeof(long int *) * too_long_kmer + 5);
	for(current_kmer_length = shortest_kmer; current_kmer_length < too_long_kmer; current_kmer_length++)
	{
	kmer_length_size = ((current_kmer_length-2) * (count_also_spaced_kmers != 0)) + 2;
	results[counter][current_kmer_length] = malloc(sizeof(long int *) * kmer_length_size + 5);
		for(current_gap_position = 0; current_gap_position < kmer_length_size; current_gap_position++)
		{
        /* DO NOT ALLOCATE MEMORY FOR ALL GAPS IF ONLY CENTER GAPS COUNTED */
        if (count_also_spaced_kmers == 0 && current_gap_position > 1) continue;
        if (count_also_spaced_kmers == 1 && (current_gap_position != current_kmer_length / 2 && current_gap_position != current_kmer_length / 2 + current_kmer_length % 2 && current_gap_position > 1)) continue;
            
            
		gap_position_size = ((Nlength-current_kmer_length+1) * (count_also_spaced_kmers != 0)) + 1;
		results[counter][current_kmer_length][current_gap_position] = malloc(sizeof(long int *) * gap_position_size + 5);
			for(current_gap_length = 0; current_gap_length < gap_position_size; current_gap_length++)
			{
        /* DO NOT ALLOCATE MEMORY FOR ALL GAPS IF ONLY CENTER GAPS COUNTED */        
                if (current_gap_position > 1 && current_gap_length == 0) continue;
                if (current_gap_position == 0 && current_gap_length > 1) break;
                if (count_also_spaced_kmers != 2 && current_gap_position == 1 && current_gap_length == 1) break;
                
                variable_length = pow(4, current_kmer_length);
                results[counter][current_kmer_length][current_gap_position][current_gap_length] = malloc(variable_length * sizeof(long int *) + 5);
                // memset (results[counter][current_kmer_length][current_gap_position][current_gap_length], 0, pow(4, current_kmer_length) * sizeof(long int *) + 5);
				current_value = results[counter][current_kmer_length][current_gap_position][current_gap_length];
				for(current_kmer = 0; current_kmer < variable_length; current_kmer++)
				{
				current_value[current_kmer] = 0; 
				}
			}
		}
	}
	}
	
	/* FILE MAIN LOOP */
	for (; file_number < 2; fclose(open_file), file_number++)
	{
	open_file = fopen(file_name[file_number], "r");
	if (open_file == '\0') 
	{
		printf("File: %s not found\n\n", file_name[file_number]); 
		if(print_counts == 1) number_of_files = 1;
		else exit(1);
	}

	/* CLEARS DIFFERENT SEQUENCES LIST */
	if(remove_non_unique == 1) 
	{
		if (extendedoutput == 1) printf ("Clears sequence list\n");
		for (counter = 0; counter < 101; counter++) 
		{
			/* printf("\t%i", counter); */
			unordered_list[counter].sequence_value_ULL = 1E9;
			unordered_list[counter].incidence = 0;
		}
		
		for (counter = 0; counter < max_number_of_sequences; counter++) 
		{
				sorted_list[0][counter].sequence_value_ULL = 0;
				sorted_list[0][counter].incidence = 0;
				sorted_list[1][counter].sequence_value_ULL = 0;
				sorted_list[1][counter].incidence = 0;
		}

		number_of_unordered_sequences = 0;
		for(counter = 0; counter < 4097; counter++) sorted_index[counter] = '\0';
		sorted_index[counter] = sorted_list[0];
	}
		if (extendedoutput == 1) {printf ("Sequence list cleared\n");fflush(stdout);}
		
	/* SEQUENCE LINE LOADING LOOP */
	for(eof_reached = 0, line = 1; ; line++)
	{
	
	start_of_line_loop:
	
	
	
		/* TAKES ONE SEQUENCE FROM FILE */
		for(charcounter = 0; ; ) 
		{
			text1 = getc(open_file);
			/* printf("%c2", text1); */
			if (text1 == EOF) {eof_reached = 1; break;}
			if (text1 == '\n') break;
			if (text1 == 'A' || text1 == 'C' || text1 == 'G' || text1 == 'T' || text1 == 'N') /* ONLY ACCEPTS NUCLEOTIDES IN CAPS, N ONLY ACCEPTED SO THAT ERROR WILL BE DIFFERENT */
			{
			current_sequence[charcounter] = text1;
			charcounter++;
			}
		}
		current_sequence[charcounter] = '\0';
		/* printf ("\nSequence: %s", current_sequence); */
		if(eof_reached == 0 && strlen(current_sequence) != Nlength-1) {printf("\nWrong sequence length on line %li", line); goto start_of_line_loop;}
	
	/* CHECKS IF AT END OF FILE */
	if (eof_reached == 1)
	{
	printf("\nEOF encountered in file %i on line %li\n", file_number, line);
	break;
	}
		
	/* STRAND LOOP */
	for(non_unique_flag = 0, strand = 0; strand < 2 && non_unique_flag == 0; strand++)
	{
	
	/* printf ("\tstrand %i", strand); */
	
	/* CALCULATES INTEGER VALUE CORRESPONDING TO SEQUENCE N-mer */
	if (strand == 0)
	{	
	/* FORWARD STRAND */
	for(current_sequence_value_ULL = 0, position = 0, position_value = left_position_value; position < Nlength-1; position++, position_value /= 4)
	{
		for (nucleotide_value = 0; nucleotide_value < 4 && current_sequence[position] != dnaforward[nucleotide_value]; nucleotide_value++);
		if(nucleotide_value == 4) {/* printf("\nSEQUENCE ERROR AT POSITION %i, %i \n", line, position); */ goto start_of_line_loop;}
		else if (file_number == 0 && (multinomial && complex_background) == 0) nc[strand].incidence[nucleotide_value][position]++; /* COUNTS NUCLEOTIDE */
		current_sequence_value_ULL += position_value * nucleotide_value;
	}
	forward_sequence_value_ULL = current_sequence_value_ULL;
	/* printf("\nsequence_value_ULL = %llu", current_sequence_value_ULL); */
	}
	else
	{
	/* REVERSE STRAND */
	for(current_sequence_value_ULL = 0, position = Nlength-2, position_value = left_position_value; position > -1; position--, position_value /= 4)
	{
		for (nucleotide_value = 0; nucleotide_value < 4 && current_sequence[position] != dnareverse[nucleotide_value]; nucleotide_value++);
		if(nucleotide_value == 4) printf("\nSEQUENCE ERROR AT POSITION %li, %i \n", line, position);
		else if (file_number == 0 && (multinomial && complex_background) == 0) nc[strand].incidence[nucleotide_value][Nlength - position - 2]++; /* COUNTS NUCLEOTIDE */
		current_sequence_value_ULL += position_value * nucleotide_value;
	}	
	}
	
        /* printf("\nsequence, value, sequence : %s\t%ju\t", current_sequence, current_sequence_value_ULL);
        for(position = Nlength-2; position > -1 ; position--) printf("%c", forward[(current_sequence_value_ULL & mask_ULL[1][position]) >> (position * 2)]); */

	/* IF REMOVE NON UNIQUE IS ONE, INDEXES FULL-LENGTH SEQUENCES TO FIND IF SAME SEQUENCE REOCCURS */
	if(strand == 0 && remove_non_unique == 1)
	{
		for (shortcounter = 0; shortcounter < number_of_unordered_sequences; shortcounter++)
		{
			/* printf ("\nComparing %i to %i", current_sequence_value_ULL, unordered_list[shortcounter].sequence_value_ULL); */
			if (current_sequence_value_ULL == unordered_list[shortcounter].sequence_value_ULL) {non_unique_flag = 1; break;}
		}
		if (non_unique_flag == 1)
		{
		unordered_list[shortcounter].incidence++; 
		/* printf("\nUNORDERED LIST HIT: Sequence %i found %i times", unordered_list[shortcounter].sequence_value_ULL, unordered_list[shortcounter].incidence); */ 
		continue;
		}
		unordered_list[number_of_unordered_sequences].sequence_value_ULL = current_sequence_value_ULL;
		unordered_list[number_of_unordered_sequences].incidence = 1;
		index_position = sorted_index[current_sequence_value_ULL / 65536];
		if (index_position != '\0')
		{
		/* printf("\nLOOKING IF %i IS NON-UNIQUE IN SORTED LIST, STARTING FROM %i", current_sequence_value_ULL, (*sorted_index[current_sequence_value_ULL / 65536]).sequence_value_ULL);
		printf("\nRANGE %i", sorted_index[4097] - index_position); */
		for (; index_position < sorted_index[4097]; index_position++)
		{
			/* printf ("\nzzComparing %i to %i", current_sequence_value_ULL, (*index_position).sequence_value_ULL); */
			if (current_sequence_value_ULL == (*index_position).sequence_value_ULL) {non_unique_flag = 1; break;}
			if (current_sequence_value_ULL < (*index_position).sequence_value_ULL) break;
		}
		if (non_unique_flag == 1) 
		{
		(*index_position).incidence++; 
		/* printf("\nSORTED LIST HIT: Sequence %i found %i times", (*index_position).sequence_value_ULL, (*index_position).incidence); */
		continue;
		}
		}
		if (number_of_unordered_sequences > 99) 
			{
				/* printf("100 found-"); */
				number_of_unordered_sequences = Sort_and_index (sorted_list, unordered_list, current_sorted_list, number_of_unordered_sequences, sorted_index);
				if (current_sorted_list == 1) current_sorted_list = 0;
				else current_sorted_list = 1;
				/* printf("\nCurrent_number_of_different_sequences: %i", sorted_index[4097] - sorted_list[current_sorted_list]); */
			}
		else number_of_unordered_sequences++;
	}
		
		
		if (extendedoutput == 1)
		{
			printf("\nLINE %li, STRAND %i, integer_value %ju, sequence: ", line, strand, (uintmax_t) current_sequence_value_ULL);
			/* SEQUENCE PRINT */
			for(position = Nlength-2; position > -1 ; position--) printf("%c", forward[(current_sequence_value_ULL & mask_ULL[1][position]) >> (position * 2)]); 
		}
	
	/* only counts nucleotides from first file if flank is used 
	if (flank == 1 && file_number == 0) continue; */

	
	if (flank == 1)
	{
	/* FLANK TOOL */
	
    /* FINDS IF THIS SEQUENCE HAS MATCH IN EITHER STRAND */
    if (strand == 0)
    {
        /* COUNTS MATCHES TO FILTER IF FILTER SET */
        if (match_filter != 0)
        {
            Findpwmmatch (&mfp, filter_cut_off, current_sequence_value_ULL, &match[0]);
            Findpwmmatch (&mfp, filter_cut_off, Reverse_complement_sequence_value(current_sequence_value_ULL, Nlength-1), &match[1]);
            Remove_palindromic_matches(&match,query_sequence_length);
            matches_to_filter = match[0].position[0] + match[1].position[0];
        }
        
        /* COUNTS MATCHES TO MULTINOMIAL 2 IF DINUC SET AND MULTINOMIAL IS LESS THAN 2 */
        if (multinomial_2_dinuc != 0)
        {
            Findpwmmatch (&mmp, multi_2_cut_off, current_sequence_value_ULL, &dinucmatch[0]);
            Findpwmmatch (&mmp, multi_2_cut_off, Reverse_complement_sequence_value(current_sequence_value_ULL, Nlength-1), &dinucmatch[1]);
            Remove_palindromic_matches(&dinucmatch,query_sequence_length);
            if (dinucmatch[0].position[0] == 1 && dinucmatch[1].position[0] == 0) Multinomial_add_to_dinucleotide_matrix(&one_hit_di[file_number], &di, &dep[file_number],  &qp, dinucmatch[0].position[1], dinucmatch[0].score[1], multi_2_cut_off, forward_sequence_value_ULL);
            if (dinucmatch[0].position[0] == 0 && dinucmatch[1].position[0] == 1) Multinomial_add_to_dinucleotide_matrix(&one_hit_di[file_number], &di, &dep[file_number], &qp, dinucmatch[1].position[1], dinucmatch[1].score[1], multi_2_cut_off, Reverse_complement_sequence_value(current_sequence_value_ULL, Nlength-1));
        }
        
        /* FINDS number_of_matches of QUERY SEQUENCE OR PWM in CURRENT SEQUENCE */
        if (flank_with_pwm == 0) 
        {
        number_of_matches[0] = Findexactmatch (query_sequence_value, query_sequence_length, current_sequence_value_ULL, &match[0]);
        number_of_matches[1] = Findexactmatch (query_sequence_value, query_sequence_length, Reverse_complement_sequence_value(current_sequence_value_ULL, Nlength-1), &match[1]);
        }
        else 
        {
        number_of_matches[0] = Findpwmmatch (&qp, cut_off, current_sequence_value_ULL, &match[0]);
        number_of_matches[1] = Findpwmmatch (&qp, cut_off, Reverse_complement_sequence_value(current_sequence_value_ULL, Nlength-1), &match[1]);
        }
    if (number_of_matches[0] + number_of_matches[1] > 0) 
    {
    current_sequence_contains_match = 1;
    }
    else 
    {
    current_sequence_contains_match = 0;
    number_of_sequences_with_no_hits[file_number]++;
    }
    }

    /* EXCLUDES ENTIRE STRAND IF INDICATED */
    if ((limit_hits_to_strand == 'F' || limit_hits_to_strand == 'f') && limit_hits_to_position == 0) number_of_matches[1] = 0;
    if ((limit_hits_to_strand == 'R' || limit_hits_to_strand == 'r') && limit_hits_to_position == 0) number_of_matches[0] = 0;
        
        
	/* REMOVES MATCHES FROM REVERSE STRAND IF QUERY SEQUENCE IS A PALINDROME
	if (searchstring_is_palindrome == 1) number_of_matches[1] = 0; */
	

	
	/* printf("\n%i number_of_matches in strand %i at positions: ", number_of_matches[strand], strand);
	for (counter = 0; counter < number_of_matches[strand]; counter++) 
	{
	if (strand == 0) printf ("%i\t", match.position[strand][counter]);
	else printf ("%i\t", Nlength - match.position[strand][counter] - query_sequence_length + 1);
	} */
	
	
	/* CHECKS IF HIT POSITION & STRAND LIMITED BY USER AND FINDS IF ANY OF THE HITS number_of_matches MATCH USER DEFINED POSITION(S) */
	if (limit_hits_to_position != 0)
	{
	correct_position = 0;
	if (limit_hits_to_strand != 'R' && limit_hits_to_strand != 'r') for (counter = 1; counter <= number_of_matches[0]; counter++) if (match[0].position[counter] == limit_hits_to_position) correct_position = 1;
	if (limit_hits_to_strand != 'F' && limit_hits_to_strand != 'f') for (counter = 1; counter <= number_of_matches[1]; counter++) if (Nlength - query_sequence_length + 1 - match[1].position[counter] == limit_hits_to_position) correct_position = 1;
	}
        
	/* CHECKS IF HIT POSITION & STRAND EXCLUDED BY USER AND FINDS IF ANY OF THE HITS IN EXCLUDED POSITIONS */
	if (exclude_position[0] != 0)
	{
	correct_position = 1;
	for (counter2 = 0; exclude_position[counter2] != 0; counter2++)
	{
	if (exclude_strand[counter2] != 'R' && exclude_strand[counter2] != 'r') for (counter = 1; counter <= number_of_matches[0]; counter++) if (match[0].position[counter] == exclude_position[counter2]) correct_position = 0;
	if (exclude_strand[counter2] != 'F' && exclude_strand[counter2] != 'f') for (counter = 1; counter <= number_of_matches[1]; counter++) if (Nlength - query_sequence_length + 1 - match[1].position[counter] == exclude_position[counter2]) correct_position = 0;
	}

	}
	
    /* REMOVES PALINDROMIC MATCHES */
    if (strand == 0 && number_of_matches[0] > 0 && number_of_matches[1] > 0) 
    {
    Remove_palindromic_matches(&match,query_sequence_length);
    palindromic_hits[file_number] += number_of_matches[0] + number_of_matches[1] - match[0].position[0] - match[1].position[0];
    number_of_matches[0] = match[0].position[0];
    number_of_matches[1] = match[1].position[0];
    }
    if (match_filter == 0) matches_to_filter = number_of_matches[0] + number_of_matches[1]; /* SEPARATE FILTER IS NOT SET */
    
	if (strand == 1 && correct_position == 1)
	{
	if (number_of_matches[0] + number_of_matches[1] > 0) number_of_sequences_with_hits[file_number]++;
	/* FLANK (ONLY SEQUENCES WITH ONE OCCURRENCE) */
	/* ADDS NUCLEOTIDES FROM number_of_matches TO ONE HIT RESULT PWM */
	
    /* ONE HIT ONLY (MATCHES TO FILTER IS 1) */
    if (matches_to_filter == 1)
    {
    if (number_of_matches[1] == 0 && number_of_matches[0] == 1) 
	{
    /* ADDS FLANK KMER */
    if(flank_kmer_pos != -100 && (Nlength - match[0].position[1])-flank_kmer_pos < Nlength && (Nlength - match[0].position[1])-flank_kmer_pos > 3) flank_kmer_count[file_number][(forward_sequence_value_ULL >> (((Nlength - match[0].position[1])-flank_kmer_pos-4)*2)) & 255]++;
        
	if (multinomial == 0) 
    {
        Add_to_pwm (&one_hit_pwm[file_number], match[0].position[1], forward_sequence_value_ULL, &background_pwm[0]);
    	if (dinucleotide_properties == 1 && multinomial_2_dinuc == 0) Add_to_dinucleotide_matrix(&one_hit_di[file_number], &di, &dep[file_number], match[0].position[1], forward_sequence_value_ULL);
    }
	else 
    {
        Multinomial_add_to_pwm (&one_hit_pwm[file_number], &qp, match[0].position[1], match[0].score[1], cut_off, forward_sequence_value_ULL, &background_pwm[0]);
        if (dinucleotide_properties == 1 && multinomial_2_dinuc == 0) Multinomial_add_to_dinucleotide_matrix(&one_hit_di[file_number], &di, &dep[file_number],  &qp, match[0].position[1], match[0].score[1], cut_off, forward_sequence_value_ULL);
	}
    }
	
    if (number_of_matches[0] == 0 && number_of_matches[1] == 1) 
	{
        /* ADDS FLANK KMER */
        if(flank_kmer_pos != -100 && (Nlength - match[1].position[1])-flank_kmer_pos < Nlength && (Nlength - match[1].position[1])-flank_kmer_pos > 3) flank_kmer_count[file_number][(current_sequence_value_ULL >> (((Nlength - match[1].position[1])-flank_kmer_pos-4)*2)) & 255]++;
        
        if (multinomial == 0)
        {
        Add_to_pwm (&one_hit_pwm[file_number], match[1].position[1], current_sequence_value_ULL, &background_pwm[1]);
        if (dinucleotide_properties == 1 && multinomial_2_dinuc == 0) Add_to_dinucleotide_matrix(&one_hit_di[file_number], &di, &dep[file_number], match[1].position[1], current_sequence_value_ULL);
        }
        else 
        {
        Multinomial_add_to_pwm (&one_hit_pwm[file_number], &qp, match[1].position[1], match[1].score[1], cut_off, current_sequence_value_ULL, &background_pwm[1]);
        if (dinucleotide_properties == 1 && multinomial_2_dinuc == 0) Multinomial_add_to_dinucleotide_matrix(&one_hit_di[file_number], &di, &dep[file_number], &qp, match[1].position[1], match[1].score[1], cut_off, current_sequence_value_ULL);
        }
	}
	}
    /* END OF ONE HIT ONLY */
        
	/* FLANK (SEQUENCES WITH TWO OCCURRENCES) */
	if ((number_of_matches[0] == 2 || number_of_matches[1] == 2) && number_of_matches[0] + number_of_matches[1] == 2 && file_number == 1) 
	{
		orientation = head_to_tail;
		if (number_of_matches[0] == 2) 
		{
			spacing = match[0].position[1] - match[0].position[2];
			if (spacing > 49) Exit_with_error("spacing too long", error_values);
			first_match_position = match[0].position[2];
			first_sequence_value_ULL = forward_sequence_value_ULL;
			/* printf("\nSame orientation in forward strand "); */
			Add_to_pwm (&two_hits_pwm[file_number][orientation][spacing], first_match_position, forward_sequence_value_ULL, &background_pwm[0]);
		}
		else 
		{
			spacing = match[1].position[1] - match[1].position[2];
			first_match_position = match[1].position[2];
			first_sequence_value_ULL = current_sequence_value_ULL;
			/* printf("\nSame orientation in reverse strand "); */
			Add_to_pwm (&two_hits_pwm[file_number][orientation][spacing], first_match_position, first_sequence_value_ULL, &background_pwm[1]);
		}
		two_hits_connecting_matrix.incidence[orientation][spacing]++;
		/* printf ("%s, spacing %i, %i:th occurrence, position %i", orientation_string[orientation], spacing, two_hits_pwm[file_number][orientation][spacing].max_counts+1, first_match_position);*/
		
	}
	
	if (number_of_matches[0] == 1 && number_of_matches[1] == 1 && file_number == 1)
	{
		first_match_position = Nlength - match[1].position[1] - query_sequence_length + 1;
		/* printf("\nOpposite orientation, values forward %i, reverse %i", match[0].position[1], first_match_position); */
		spacing = match[0].position[1] - first_match_position;
		if (spacing > 0) 
		{
			orientation = tail_to_tail;
			first_sequence_value_ULL = forward_sequence_value_ULL;
		}
		else
		{
			orientation = head_to_head;
			spacing = first_match_position - match[0].position[1];
			first_match_position = match[0].position[1];
			first_sequence_value_ULL = forward_sequence_value_ULL;
		}
	/* printf ("\nOpposite orientation %s, spacing %i, %i:th occurrence, position %i", orientation_string[orientation], spacing, two_hits_pwm[file_number][orientation][spacing].max_counts+1, first_match_position); */
	Add_to_pwm (&two_hits_pwm[file_number][orientation][spacing], first_match_position, first_sequence_value_ULL, &background_pwm[0]);
	two_hits_connecting_matrix.incidence[orientation][spacing]++;
	}
	
	/* SPACING AND ORIENTATION */
	/* FILLS CONNECTING MATRIX */
	/* HEAD TO TAIL (SAME) SPACINGS AND PWM FILL */
	/* FORWARD STRAND */
    if (align_matches == 1) {max_align_score[0][0] = 0; max_align_score[1][0] = 0; max_align_score[0][1] = 0; max_align_score[1][1] = 0;}
	for (counter = 1; counter <= number_of_matches[0]; counter++) 
	{
    if (align_matches == 1) Add_to_alignscore (forward_sequence_value_ULL, &all_hits_align_scores[file_number], &current_align_score, match[0].position[counter], max_align_score[file_number]);
	if (multinomial == 0) Add_to_pwm (&all_hits_pwm[file_number], match[0].position[counter], forward_sequence_value_ULL, &background_pwm[0]);
	else Multinomial_add_to_pwm (&all_hits_pwm[file_number], &qp, match[0].position[counter], match[0].score[counter], cut_off, forward_sequence_value_ULL, &background_pwm[0]);
	}
	if (file_number == 1) for (counter = 1; counter <= number_of_matches[0]; counter++) 
	{
	hit_position.incidence[0][match[0].position[counter]]++;
	for (counter2 = counter; counter2 <= number_of_matches[0]; counter2++) cm.incidence[head_to_tail][(match[0].position[counter]-match[0].position[counter2])]++;
	}
	/* REVERSE STRAND */
	for (counter = 1; counter <= number_of_matches[1]; counter++)
	{
    if (align_matches == 1) Add_to_alignscore (current_sequence_value_ULL, &all_hits_align_scores[file_number], &current_align_score, match[1].position[counter], max_align_score[file_number]);
	if (multinomial == 0) Add_to_pwm (&all_hits_pwm[file_number], match[1].position[counter], current_sequence_value_ULL, &background_pwm[1]);
	else Multinomial_add_to_pwm (&all_hits_pwm[file_number], &qp, match[1].position[counter], match[1].score[counter], cut_off, current_sequence_value_ULL, &background_pwm[1]);
	}
	if (file_number == 1) for (counter = 1; counter <= number_of_matches[1]; counter++)
	{
	hit_position.incidence[1][match[1].position[counter]-1]++;
	for (counter2 = counter; counter2 <= number_of_matches[1]; counter2++) cm.incidence[head_to_tail][(match[1].position[counter]-match[1].position[counter2])]++;
	}
	/* HEAD TO HEAD and TAIL TO TAIL (OPPOSITE) */
	if (file_number == 1) for (counter = 1; counter <= number_of_matches[0]; counter++) for (counter2 = 1; counter2 <= number_of_matches[1]; counter2++)
	{
	reverse_strand_position = Nlength - match[1].position[counter2] - query_sequence_length + 1;
	if (match[0].position[counter] >= reverse_strand_position) cm.incidence[tail_to_tail][(match[0].position[counter]-reverse_strand_position)]++;
	if (match[0].position[counter] <= reverse_strand_position) cm.incidence[head_to_head][(reverse_strand_position-match[0].position[counter])]++;
	}
    if (align_matches == 1 && max_align_score[file_number][0] > 0) align_score_histogram[file_number][0][max_align_score[file_number][0]]++;
    if (align_matches == 1 && max_align_score[file_number][1] > 0) align_score_histogram[file_number][1][max_align_score[file_number][1]]++;        
	}
	
	/* END OF FLANK TOOL */
	}
	
	
	/* KMER COUNT TOOL */
	if (kmer_count == 1)
	{
	/* COUNTS KMERS */
	
	/* KMERS WITH NO GAPS */

	for (current_kmer_length = shortest_kmer; current_kmer_length < too_long_kmer; current_kmer_length++)
	{
		/* printf("\n\nSequence length: %i", current_kmer_length); */
		end_position = Nlength-current_kmer_length;
		for (position = 0; position < end_position; position++)
		{
			/* current_value = results[current_kmer_length][0][0]; */
            
			results[file_number][current_kmer_length][0][0][(current_sequence_value_ULL & mask_ULL[current_kmer_length][position]) >> (position * 2)]++;
			/* printf("\n%i", (current_sequence_value_ULL & mask_ULL[current_kmer_length][position]) >> (position * 2)); */
			/* fflush(stdout); */
		}
		if (strand == 1 && current_sequence_contains_match == 0) for (position = 0; position < end_position; position++)
		{
		/* COUNTS KMERS FROM NON-MATCHING SEQUENCE */
			results[file_number][current_kmer_length][1][0][(forward_sequence_value_ULL & mask_ULL[current_kmer_length][position]) >> (position * 2)]++;
			results[file_number][current_kmer_length][1][0][(current_sequence_value_ULL & mask_ULL[current_kmer_length][position]) >> (position * 2)]++;
		}
	}
	}
	
	if (kmer_count == 1 && no_insertions == 0 && count_also_spaced_kmers != 0)
	{
	/* KMERS WITH GAPS */
	/* GENERATES DELETIONS INTO CURRENT SEQUENCE */
	for (deletion_size = 1; deletion_size < Nlength - shortest_kmer; deletion_size++)
	{
		Nmer_position = 1; 
        too_long_nmer_position = Nlength-1-deletion_size;
		for (; Nmer_position < too_long_nmer_position; Nmer_position++)
		{
			deleted_sequence_value_ULL = (current_sequence_value_ULL & lowmask_ULL[Nmer_position-1]) ^ ((current_sequence_value_ULL & highmask_ULL[Nmer_position+deletion_size]) >> (deletion_size * 2));
			/* printf ("\n\n%i\t%i\t%i", deletion_size, Nmer_position, deleted_sequence_value_ULL); */
			/* SEQUENCE PRINT */
			/* for(position = Nlength-2; position > -1 ; position--) printf("%c", forward[(deleted_sequence_value_ULL & mask_ULL[1][position]) >> (position * 2)]);*/
			
			/* FINDS KMERS WITH GAPS FROM DELETED SEQUENCE */
			for (current_kmer_length = shortest_kmer; current_kmer_length < too_long_kmer; current_kmer_length++)
			{
                if (count_also_spaced_kmers == 1) 
                {
                position = Nmer_position - current_kmer_length + current_kmer_length / 2;
                end_position = position + 1 + current_kmer_length % 2;
                }
                else 
                {
                position = Nmer_position - current_kmer_length + 1;
                end_position = Nmer_position;
				}
                if (position < 0) position = 0;
				if (end_position > Nlength - current_kmer_length - deletion_size) end_position = Nlength - current_kmer_length - deletion_size;
                for(; position < end_position; position++)
				{

					results[file_number][current_kmer_length][current_kmer_length-Nmer_position+position][deletion_size][((deleted_sequence_value_ULL & mask_ULL[current_kmer_length][position])) >> (position * 2)]++;
				}
																		   
			}
			
		}
	}
	
	/* END OF KMER COUNT TOOL */
	}
	
	
	/* END OF STRAND LOOP */
	}
	
	
	/* END OF LINE LOOP */
	}


	if (remove_non_unique == 1) 
		{
		number_of_unordered_sequences = Sort_and_index (sorted_list, unordered_list, current_sorted_list, number_of_unordered_sequences, sorted_index);
		if (current_sorted_list == 1) current_sorted_list = 0;
		else current_sorted_list = 1;
		/* printf("\nCurrent_number_of_different_sequences , unordered: %i , %i", sorted_index[4097] - sorted_list[current_sorted_list], number_of_unordered_sequences);*/
		number_of_sequences_analyzed[file_number] = (double) ((sorted_index[4097] - sorted_list[current_sorted_list]) + number_of_unordered_sequences);
		/* (sorted_index[4097] - sorted_list[current_sorted_list]) + number_of_unordered_sequences; */
		}
		else number_of_sequences_analyzed[file_number] = (double) (line-1);
		printf("\nNumber of sequences in file %i: %.0f", file_number, number_of_sequences_analyzed[file_number]);
		
		/* PRINTS INPUT SEQUENCES AND THEIR INCIDENCE */ 
		if (print_input_sequences == 1 && remove_non_unique == 1)
		{
			printf("\n\nUnique_sequences: %.2f\nSequence", number_of_sequences_analyzed[file_number]);
			if (remove_non_unique == 1) printf("\tIncidence");
			for(counter = 0; counter < number_of_sequences_analyzed[file_number]; counter++)
			{
				printf("\n");
				for(position = Nlength-2; position > -1 ; position--) printf("%c", forward[(sorted_list[1-current_sorted_list][counter].sequence_value_ULL & mask_ULL[1][position]) >> (position * 2)]);
				if (remove_non_unique == 1) printf("\t%li", sorted_list[1-current_sorted_list][counter].incidence);
			}
		
		} 
	
	/* printf ("\n1-%i", np[0].width); */
		
		
	/* CALCULATES BACKGROUND NUCLEOTIDE DISTRIBUTION FROM FIRST FILE */
	if (file_number == 0 && even_background == 0 && (multinomial && complex_background) == 0) 
	{
	Count_to_normalized_pwm(&background_pwm[0], &nc[0]);
	Count_to_normalized_pwm(&background_pwm[1], &nc[1]);
	}
	
	/* printf ("\n2-%i", np[0].width); */
	/* END OF MAIN (FILE) LOOP */
	}
	
exit_with_style:
	

		
	
	if (flank == 1)
	{
	
        
	/* ESTIMATES LAMBDA FROM KMER COUNT DATA (TAKES TOTAL COUNT OF LOWER HALF OF 8-MERS) */
	if (kmer_count == 1)
	{
	/* SORT KMERS ACCORDING TO INCIDENCE */
	number_of_kmers = pow(4,shortest_kmer);
	background_kmer_count = 0;
	signal_kmer_count = 0;
	current_kmer_length = shortest_kmer;
    if(lowpaircount == 1)
    {
    struct sumtable sums;
    sumtable_init(&sums, current_kmer_length);    
    Lowpaircounts (&sums, results[1][current_kmer_length][0][0], 0, current_kmer_length);
    sumtable_init(&sums, current_kmer_length);
    }
	qsort (results[0][current_kmer_length][0][0], number_of_kmers, sizeof(long int), Numeric_sort_long_ints);
	qsort (results[1][current_kmer_length][0][0], number_of_kmers, sizeof(long int), Numeric_sort_long_ints); 
	for (counter = number_of_kmers * 0.25; counter < number_of_kmers * 0.75; counter++)
	{
	signal_kmer_count += results[1][current_kmer_length][0][0][counter]; 
	background_kmer_count += results[0][current_kmer_length][0][0][counter];
	/* printf("\n%i", counter); 
	fflush(stdout); */
	}
	lambda =  (signal_kmer_count * number_of_sequences_analyzed[0]) / (background_kmer_count * number_of_sequences_analyzed[1]);
	}
	else lambda =  (number_of_sequences_with_hits[0] * number_of_sequences_analyzed[1]) / (number_of_sequences_with_hits[1] * number_of_sequences_analyzed[0]);
	
	Count_to_normalized_pwm(&np[current_logo], &all_hits_pwm[0]);
	strcpy (np[current_logo].name,"All hits background ");
	sprintf(tempstring, " : %li", all_hits_pwm[0].max_counts);
	strcat(np[current_logo].name, tempstring);
	current_logo++;

	Count_to_normalized_pwm(&np[current_logo], &all_hits_pwm[1]);
	strcpy (np[current_logo].name,"All hits uncorrected ");
	sprintf(tempstring, " : %li", all_hits_pwm[1].max_counts);
	strcat(np[current_logo].name, tempstring);
	current_logo++;
	
        /* GENERATES EXPONENTIAL PWM */
        double pseudocount;
        double pseudocount2;
        double swap2;
        double sizefactor;
        sizefactor = (number_of_sequences_analyzed[0] / number_of_sequences_analyzed[1]);
        if (multinomial != 0 && even_background == 0 && complex_background == 1)
        {
            pseudocount = all_hits_pwm[1].max_counts * 0.01;
            pseudocount2 = one_hit_pwm[1].max_counts * 0.01;
            for (counter = 0; counter < all_hits_pwm[0].width; counter++) for (counter2 = 0; counter2 < 4; counter2++)
            {
            swap2 = 10000 * (((double) all_hits_pwm[1].incidence[counter2][counter]+pseudocount) * sizefactor / ((double) all_hits_pwm[0].incidence[counter2][counter] + pseudocount * sizefactor)-lambda);
            all_hits_exponential[1].incidence[counter2][counter] = (int) swap2;
            swap2 = 10000 * (((double) one_hit_pwm[1].incidence[counter2][counter]+pseudocount2) * sizefactor / (((double) one_hit_pwm[0].incidence[counter2][counter]) + pseudocount2 * sizefactor)-lambda);
            one_hit_exponential[1].incidence[counter2][counter] = (int) swap2;
            }
        }

        
	/* SUBTRACTS BACKGROUND MULTINOMIAL DISTRIBUTION FROM SIGNAL */
	if (multinomial != 0 && even_background == 0 && complex_background == 1)
	{
	for (counter = 0; counter < all_hits_pwm[0].width; counter++) for (counter2 = 0; counter2 < 4; counter2++)
	{
	swap = (all_hits_pwm[1].incidence[counter2][counter]) * sizefactor - lambda * all_hits_pwm[0].incidence[counter2][counter];
	/* printf("\nPosition %i %i: background counts %.3f, signal counts %.3f, corrected value %.3f", counter2, counter,  all_hits_pwm[0].incidence[counter2][counter],  all_hits_pwm[1].incidence[counter2][counter], swap); */
	all_hits_pwm[1].incidence[counter2][counter] = swap;
	swap = (one_hit_pwm[1].incidence[counter2][counter]) * sizefactor - lambda * one_hit_pwm[0].incidence[counter2][counter];
	one_hit_pwm[1].incidence[counter2][counter] = swap;
	}
	}
	
        
	
	file_number = 1;
	Count_to_normalized_pwm(&p, &all_hits_pwm[file_number]);

	/* PRINTS PWMs */
	/* ALL HITS */
    short int left_flank_length = 10;
    short int right_flank_length = 12;
	first = Nlength - left_flank_length - 1;
	last = Nlength + right_flank_length + query_sequence_length - 1;
	long int position;
	short int nucleotide;
	fflush(stdout);
	printf ("\nAll Hits \tPosition");
	printf ("\nAll Hits               ");
	for (position = first; position < last; position++)
	{
	printf ("\t%li", (position - first - left_flank_length + 1));
	}
	for (counter2 = 0; counter2 < 4; counter2++)
	{
	printf ("\nAll Hits");
	printf ("\t%c", forward[counter2]);
	for (counter = first; counter < last;  counter++) 
	{
	if (print_frequencies == 0) printf ("\t%.0f", all_hits_pwm[file_number].incidence[counter2][counter]);
	else printf ("\t%0.3f", p.fraction[counter2][counter]);
	}
	}
	
	fflush(stdout);
	/* ONLY ONE HIT */
	Count_to_normalized_pwm(&p, &one_hit_pwm[file_number]);
	printf ("\n\nOne Hit   \tPosition");
	printf ("\nOne Hit              ");
	for (counter = first; counter < last;  counter++) printf ("\t%li", counter - first - left_flank_length + 1);
	for (counter2 = 0; counter2 < 4; counter2++)
	{
	printf ("\nOne Hit     ");
	printf ("\t%c", forward[counter2]);
	for (counter = first; counter < last; counter++) 
	{
	if (print_frequencies == 0) printf ("\t%.0f", one_hit_pwm[file_number].incidence[counter2][counter]);
	else printf ("\t%0.3f", p.fraction[counter2][counter]);
	}
	}

    /* ONLY ONE HIT, EXPONENTIAL */
    Count_to_normalized_pwm(&p, &one_hit_exponential[file_number]);
    printf ("\n\nOneHit_Exp\tPosition");
    printf ("\nOneHit_Exp          ");
    for (counter = first; counter < last;  counter++) printf ("\t%li", counter - first - left_flank_length + 1);
    for (counter2 = 0; counter2 < 4; counter2++)
    {
    printf ("\nOneHit_Exp ");
            printf ("\t%c", forward[counter2]);
            for (counter = first; counter < last; counter++) 
            {
                if (print_frequencies == 0) printf ("\t%.0f", one_hit_exponential[file_number].incidence[counter2][counter]);
                else printf ("\t%0.3f", p.fraction[counter2][counter]);
            }
        }
        
	/* ADDS PWMs TO LOGO LIST */
	Count_to_normalized_pwm(&np[current_logo], &all_hits_pwm[file_number]);
	strcpy (np[current_logo].name,"All hits");
	sprintf(tempstring, " : %li", all_hits_pwm[file_number].max_counts);
	strcat(np[current_logo].name, tempstring);
	current_logo++;	
		
	Count_to_normalized_pwm(&np[current_logo], &one_hit_pwm[file_number]);
	strcpy (np[current_logo].name,"One hit  ");
	sprintf(tempstring, " : %li", one_hit_pwm[file_number].max_counts);
	strcat(np[current_logo].name, tempstring);
	current_logo++;
	
    Count_to_normalized_pwm(&np[current_logo], &one_hit_exponential[file_number]);
    strcpy (np[current_logo].name,"One hit exponential");
    sprintf(tempstring, " : %li", one_hit_pwm[file_number].max_counts); /* ACTUAL MAX COUNTS THAT WERE USED, EXPONENTIAL IS NORMALIZED TO 10000 */
    strcat(np[current_logo].name, tempstring);
    current_logo++;
        
		
	/* TWO HITS */
	for (total_possible_spacings = 0, total_number_of_two_hits = 0, orientation = 0; orientation < 3; orientation++) for (spacing = query_sequence_length; spacing < Nlength - query_sequence_length; spacing++) 
	{
	total_number_of_two_hits += two_hits_pwm[file_number][orientation][spacing].max_counts;
	total_possible_spacings += (Nlength - query_sequence_length - spacing) * 2;
	}
	printf ("\nTotal number of two hits: %li", total_number_of_two_hits);
	printf ("\nTotal possible spacings: %li", total_possible_spacings);
	for (orientation = 0; orientation < 3; orientation++) for (spacing = 0; spacing < 30; spacing++) 
	{
	Count_to_normalized_pwm(&p, &two_hits_pwm[file_number][orientation][spacing]);
	if (two_hits_pwm[file_number][orientation][spacing].max_counts != 0) current_fold_expected = (two_hits_pwm[file_number][orientation][spacing].max_counts * total_possible_spacings) / (all_hits_pwm[file_number].max_counts * (Nlength - query_sequence_length - spacing) * (pow(all_hits_pwm[file_number].max_counts / number_of_sequences_analyzed[1],2) * exp(-all_hits_pwm[file_number].max_counts / number_of_sequences_analyzed[1]) / 2) );
	else current_fold_expected = 0;
	two_hits_fold_connecting_matrix.incidence[orientation][spacing] = (long int) (current_fold_expected * 100);
	/* printf("\n%i,%i,%i,%f, %i",orientation, spacing, two_hits_pwm[file_number][orientation][spacing].max_counts, current_fold_expected, two_hits_fold_connecting_matrix.incidence[orientation][spacing]); */
	if (two_hits_pwm[file_number][orientation][spacing].max_counts > minimum_kmer_count)
	{
	if (two_hits_pwm[file_number][orientation][spacing].max_counts > max_counts) {max_counts = two_hits_pwm[file_number][orientation][spacing].max_counts; max_orientation = orientation; max_spacing = spacing;}
	Count_to_normalized_pwm(&np[current_logo], &two_hits_pwm[file_number][orientation][spacing]);
	strcpy(np[current_logo].name, orientation_string[orientation]);
	sprintf(tempstring, "-%i", spacing);
	strcat(np[current_logo].name, tempstring);
	sprintf(tempstring, " : %li", two_hits_pwm[file_number][orientation][spacing].max_counts);
	strcat(np[current_logo].name, tempstring);
	sprintf(tempstring, " (%.1f)", current_fold_expected);
	strcat(np[current_logo].name, tempstring);
	current_logo++;
	printf ("\n%s-%i", orientation_string[orientation], spacing);
	printf ("\n%s-%i", orientation_string[orientation], spacing);
	printf ("\tTwo hits with spacing %i, orientation %s, max_counts %li", spacing, orientation_string[orientation], two_hits_pwm[file_number][orientation][spacing].max_counts);
	printf ("\n%s-%i", orientation_string[orientation], spacing);
	for (counter = 20; counter <= 43;  counter++) printf ("\t%li", counter - 26);
	for (counter2 = 0; counter2 < 4; counter2++)
	{
	printf ("\n%s-%i", orientation_string[orientation], spacing);
	printf ("\t%c", forward[counter2]);
	for (counter = 20; counter < 43;  counter++) 
	{
	if (print_frequencies == 0) printf ("\t%.0f", two_hits_pwm[file_number][orientation][spacing].incidence[counter2][counter]);
	else printf ("\t%0.3f", p.fraction[counter2][counter]);
	}
	}
	}
	}
	
	
	/* PRINTS MULTIPLE HIT CONNECTING MATRIX */
	Count_to_normalized_connecting_matrix (&cp, &cm);
	printf ("\n\n! Two or more hits connecting matrix\n! Orientation            \tSpacing\n!                  ");
	for (counter = 0; counter < Nlength - query_sequence_length; counter++) printf ("\t%li", counter);
	printf ("\n! ------------------------------------------------------------------------------------------------------------------------\n! Head to tail > > :");
	for (counter = 0; counter < Nlength - query_sequence_length; counter++) 
	{
	if (print_frequencies == 0) printf ("\t%li", cm.incidence[head_to_tail][counter]);
	else printf ("\t%0.3f", cp.fraction[head_to_tail][counter]);
	}
	printf ("\n! Head to head > < :");
	for (counter = 0; counter < Nlength - query_sequence_length; counter++) 
	{
	if (print_frequencies == 0) printf ("\t%li", cm.incidence[head_to_head][counter]);
	else printf ("\t%0.3f", cp.fraction[head_to_head][counter]);
	}
	printf ("\n! Tail to tail < > :");
	for (counter = 0; counter < Nlength - query_sequence_length; counter++) 
	{
	if (print_frequencies == 0) printf ("\t%li", cm.incidence[tail_to_tail][counter]);
	else printf ("\t%0.3f", cp.fraction[tail_to_tail][counter]);
	}
	
	/* PRINTS TWO HIT CONNECTING MATRIX */
	Count_to_normalized_connecting_matrix (&cp, &two_hits_connecting_matrix);
	printf ("\n\n2! Exactly two hits connecting matrix\n2! Orientation            \tSpacing\n2!                  ");
	for (counter = 0; counter < Nlength - query_sequence_length; counter++) printf ("\t%li", counter);
	printf ("\n2! ------------------------------------------------------------------------------------------------------------------------\n2! Head to tail > > :");
	for (counter = 0; counter < Nlength - query_sequence_length; counter++) 
	{
	if (print_frequencies == 0) printf ("\t%li", two_hits_connecting_matrix.incidence[head_to_tail][counter]);
	else printf ("\t%0.3f", cp.fraction[head_to_tail][counter]);
	}
	printf ("\n2! Head to head > < :");
	for (counter = 0; counter < Nlength - query_sequence_length; counter++) 
	{
	if (print_frequencies == 0) printf ("\t%li", two_hits_connecting_matrix.incidence[head_to_head][counter]);
	else printf ("\t%0.3f", cp.fraction[head_to_head][counter]);
	}
	printf ("\n2! Tail to tail < > :");
	for (counter = 0; counter < Nlength - query_sequence_length; counter++) 
	{
	if (print_frequencies == 0) printf ("\t%li", two_hits_connecting_matrix.incidence[tail_to_tail][counter]);
	else printf ("\t%0.3f", cp.fraction[tail_to_tail][counter]);
	}
	
	
	
	}
	else
	/* PRINTS COUNTS OR FOLD CHANGES */
	for (file_number = 0; file_number < number_of_files; file_number++)
	{
	if (print_counts == 0 && file_number > 0) break;
        
	if (output_all_gap_lengths_in_one_line == 0)
	{
	/* ALL LINES SEPARATELY */
	if (center_spaced_kmer_table == 0) printf("\nBackground\tSignal\tSequence_length\tGap_position\tGap length\tKmer\tIUPAC\tredundancy\tBackground count\tSignal count\tLength_normalized_count\tFold_change\tLocal_max\tRepeat\n");
    else printf("\nSequence\tBackground_count\tSignal_count\tNormalized_incidence\tScore\tRepeat\t:SK\n");
	for(current_kmer_length = shortest_kmer; current_kmer_length < too_long_kmer; current_kmer_length++)
	{
		current_gap_position = 0; 
        kmer_length_size = current_kmer_length;
		gap_position_size = (Nlength-current_kmer_length-1) * (count_also_spaced_kmers != 0) + 1;
		for(; current_gap_position < kmer_length_size; current_gap_position++)
		{
			for(current_gap_length = 0; current_gap_length < gap_position_size; current_gap_length++)
			{
                if (current_gap_length == 0 && current_gap_position == 0);
                else 
                {
                if (count_also_spaced_kmers == 1 && current_gap_position != current_kmer_length / 2 && current_gap_position != current_kmer_length / 2 + current_kmer_length % 2) continue;
                if(current_gap_length == 0 && current_gap_position > 0) continue;
                if(current_gap_position == 0 && current_gap_length > 0) break;
                }

				for(current_kmer = 0; current_kmer < pow(4, current_kmer_length); current_kmer++)
				{
					if(only_palindromes == 1 && Is_this_sequence_dimer(current_kmer, current_kmer_length) == 0) continue;
                    top_normalized_count = results[file_number+1-print_counts][current_kmer_length][current_gap_position][current_gap_length][current_kmer] >> ((too_long_kmer-current_kmer_length-1) << 1);
					if(top_normalized_count >= minimum_kmer_count) 
					{
                        if (center_spaced_kmer_table == 0)
                        {
						printf("%s\t", file_name[file_number]);
						if (print_counts == 0) printf("%s\t", file_name[file_number+1]);
						printf("%i\t%i\t%i\t", current_kmer_length, current_gap_position, current_gap_length);
						}
                        for(position = current_kmer_length-1; position > -1 ; position--)
							{
								if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");
								printf("%c", forward[(current_kmer & mask_ULL[1][position]) >> (position * 2)]);
							}
						kmer_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][current_kmer];
                        
                        /* PRINTS IUPAC */
                        printf("\t");
                        
                        for(represents_n_kmers = 1, position = current_kmer_length-1; position > -1 ; position--) 
                        {
                            if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");
                            for(iupac_bits = 0, max_nucleotide = 0, max_kmer_count = 0, nucleotide = 0; nucleotide < 4; nucleotide++)
                            {
                            test_kmer = ((current_kmer & (~(mask_ULL[1][position]))) | (nucleotide << (position * 2))) & lowmask_ULL[current_kmer_length-1];
// printf("  compared to: "); for(position2 = current_kmer_length-1; position2 > -1 ; position2--) {if(current_kmer_length - position2 - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(test_kmer & mask_ULL[1][position2]) >> (position2 * 2)]);}
                            test_kmer_count = results[file_number+1][current_kmer_length][current_gap_position][current_gap_length][test_kmer];
                                         if (test_kmer_count > max_kmer_count) 
                                         {
                                         max_nucleotide = nucleotide;
                                         max_kmer_count = test_kmer_count;
                                             // printf("\tMAX");
                                         }
                            }
                            iupac_bits |= (1 << max_nucleotide);
                            for(represents_n_nucleotides = 1, nucleotide = 0; nucleotide < 4; nucleotide++)
                            {
                            if(nucleotide == max_nucleotide) continue;
                            test_kmer = (current_kmer & (~(mask_ULL[1][position]))) | (nucleotide << (position * 2)) & lowmask_ULL[current_kmer_length-1];
                            test_kmer_count = results[file_number+1][current_kmer_length][current_gap_position][current_gap_length][test_kmer];
                            /* IUPACs BASE IF IT IS MORE THAN CUTOFF FRACTION OF MAX */
                                if(test_kmer_count > max_kmer_count * iupac_cutoff) 
                            {
                            iupac_bits |= (1 << nucleotide); 
                            represents_n_nucleotides++;
                            }
                            }
                            printf("%c", nucleotide_bitiupac[iupac_bits]);
                            represents_n_kmers *= represents_n_nucleotides;
                        }
                        printf("\t%li", represents_n_kmers);
						if (print_counts == 1) printf("\t%li\t%.1f\t%s\n", kmer_incidence , ((double) kmer_incidence ) * pow(4,current_kmer_length) / (Nlength-current_kmer_length-current_gap_length), repeat_report[Is_this_sequence_dimer(current_kmer, current_kmer_length)]);
						else 
						{
							kmer2_incidence = results[file_number+1][current_kmer_length][current_gap_position][current_gap_length][current_kmer];
							printf("\t%li\t%li\t%li\t%.1f\t%s", kmer_incidence, kmer2_incidence, top_normalized_count, (double) (number_of_sequences_analyzed[file_number]/number_of_sequences_analyzed[file_number+1]) * (double) kmer2_incidence / (double) kmer_incidence, 
                                   repeat_report[Is_this_sequence_dimer(current_kmer, current_kmer_length)]);
							if (print_p_values == 1) printf("\t%f", Winflat(kmer_incidence, kmer2_incidence, number_of_sequences_analyzed[file_number], number_of_sequences_analyzed[file_number+1]));
                            if (Localmax(results, file_number+1, shortest_kmer, too_long_kmer, current_kmer_length, current_gap_position, current_gap_length, current_kmer, count_also_spaced_kmers, kmer_length_difference_cutoff, tempstring)) printf("\tlocal_max"); 
                            else printf("\t");
							if (center_spaced_kmer_table == 0) printf("\n");
                            else printf("\t:SK\n");
						}
					}
				}
			}
		}
	}
	}
	else
	{
		/* ALL GAP LENGTHS IN ONE LINE */
		printf("\n__%s\t\t\t\tGAP LENGTHS\nFile            \tLength\tSequence\tTotal\tMax\t", file_name[file_number]);
		for (current_gap_length = 0; current_gap_length < Nlength-shortest_kmer; current_gap_length++) printf("\t%i", current_gap_length);
		printf("\tRepeat\n");
		for(current_kmer_length = shortest_kmer; current_kmer_length < too_long_kmer; current_kmer_length++)
		{
			for(current_gap_position = 0; current_gap_position < current_kmer_length; current_gap_position++)
                {
                if (current_gap_length == 0 && current_gap_position == 0);
                else if (count_also_spaced_kmers == 1 && current_gap_position != current_kmer_length / 2 && current_gap_position != current_kmer_length / 2 + current_kmer_length % 2) continue;
			
				for(current_kmer = 0; current_kmer < pow(4, current_kmer_length); current_kmer++)
				{
                    localmaxes = 0;
					max_incidence = 0;
					total_incidence = 0;
					if(only_palindromes == 1 && Is_this_sequence_dimer(current_kmer, current_kmer_length) == 0) continue;
					for(current_gap_length = 0; current_gap_length < Nlength-current_kmer_length; current_gap_length++)
					{
						/* printf("\n%i\t%i\t%i\t%i", current_kmer_length, current_gap_position, current_gap_length, current_kmer);
						fflush(stdout); */
						if (current_gap_length == 0) kmer_incidence = results[file_number][current_kmer_length][0][0][current_kmer];
						else kmer_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][current_kmer];
						if (print_counts == 0)
						{
							if (current_gap_length == 0) kmer2_incidence = results[file_number+1][current_kmer_length][0][0][current_kmer];
							else kmer2_incidence = results[file_number+1][current_kmer_length][current_gap_position][current_gap_length][current_kmer];
							total_incidence += kmer2_incidence;
							if (((number_of_sequences_analyzed[file_number]/number_of_sequences_analyzed[file_number+1]) * (double) (kmer2_incidence+1) / (double) (kmer_incidence+1)) > max_incidence) max_incidence = ((number_of_sequences_analyzed[file_number]/number_of_sequences_analyzed[file_number+1]) * (double) (kmer2_incidence+1) / (double) (kmer_incidence+1));
						}
						else 
						{
							total_incidence += kmer_incidence;
							if (kmer_incidence > max_incidence) max_incidence = kmer_incidence;
						}
					}
					if (max_incidence >= minimum_kmer_count) 
					{
						
						printf("__%s\t", file_name[file_number]);
						if (print_counts == 0) printf("%s\t", file_name[file_number+1]);
						printf("%i\t", current_kmer_length);
						for(position = current_kmer_length-1; position > -1 ; position--)
						{
							if(current_kmer_length - position - 1  == current_gap_position) printf("n");
							printf("%c", forward[(current_kmer & mask_ULL[1][position]) >> (position * 2)]);
						}
                        palindrome_correction = 1 + (Is_this_sequence_dimer(current_kmer, current_kmer_length) == 4);
						printf("\t%li\t%.0f\t", total_incidence/palindrome_correction, max_incidence/palindrome_correction);
						for(current_gap_length = 0; current_gap_length < Nlength-current_kmer_length; current_gap_length++)
						{
							if (current_gap_length == 0) kmer_incidence = results[file_number][current_kmer_length][0][0][current_kmer];
							else kmer_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][current_kmer];
							if (print_counts == 1) printf("\t%li", kmer_incidence/palindrome_correction);
							else 
							{
								if (current_gap_length == 0) kmer2_incidence = results[file_number+1][current_kmer_length][0][0][current_kmer];
								else kmer2_incidence = results[file_number+1][current_kmer_length][current_gap_position][current_gap_length][current_kmer];
								if (print_p_values == 1) printf("\t%f", Winflat(kmer_incidence, kmer2_incidence, number_of_sequences_analyzed[file_number], number_of_sequences_analyzed[file_number+1]));
								else printf("\t%f", (number_of_sequences_analyzed[file_number]/number_of_sequences_analyzed[file_number+1]) * (double) kmer2_incidence / (double) kmer_incidence);
							}
                            if (kmer2_incidence >= minimum_kmer_count && print_counts == 0) if (Localmax(results, file_number+1, shortest_kmer, too_long_kmer, current_kmer_length, current_gap_position, current_gap_length, current_kmer, count_also_spaced_kmers, kmer_length_difference_cutoff, tempstring)) localmaxes++;
						}
						printf("\t%s", repeat_report[Is_this_sequence_dimer(current_kmer, current_kmer_length)]);
                        if(localmaxes != 0) printf("\tlocal_max");
                        printf("\n");
					}
				}
			}
		}
	}
        if (count_also_spaced_kmers != 0) Kmer_svg(user_specified_output_file, results, file_number+1, shortest_kmer, too_long_kmer, count_also_spaced_kmers, kmer_length_difference_cutoff, minimum_kmer_count, number_of_sequences_analyzed[file_number], number_of_sequences_analyzed[file_number+1], number_of_heatmap_rows, tempstring, match_orientations, match_length);
	}

/* PRINTS NUCLEOTIDE COUNTS */
if (print_nucleotides == 1)
{
printf("\n\n# Nucleotide counts");
for(strand = 0; strand < 2; strand++) 
{
if (strand == 0) printf ("\n# Forward strand\n#\t"); else printf ("\n# Reverse strand\n#\t");
for (counter = 1; counter < Nlength; counter++) printf ("\t%li", counter);
for (nucleotide_value = 0; nucleotide_value < 4; nucleotide_value++)
{
printf ("\n#\t%c", forward[nucleotide_value]);
for (position = 0; position < Nlength-1; position++) printf ("\t%.0f", nc[strand].incidence[nucleotide_value][position]);
}
}
}

if (print_nucleotides == 1)
{
printf("\n\n# Nucleotide frequencies");
for(strand = 0; strand < 2; strand++) 
{
if (strand == 0) printf ("\n# Forward strand\n#\t"); else printf ("\n# Reverse strand\n#\t");
for (counter = 1; counter < Nlength; counter++) printf ("\t%li", counter);
for (nucleotide_value = 0; nucleotide_value < 4; nucleotide_value++)
{
printf ("\n#\t%c", forward[nucleotide_value]);
for (position = 0; position < Nlength-1; position++) printf ("\t%.3f", background_pwm[strand].fraction[nucleotide_value][position]);
}
}
}

fflush(stdout);

if (flank == 1)
{
strcpy (background_pwm[0].name, "Forward backgr");
strcpy (background_pwm[1].name, "Reverse backgr");
np_p[current_logo] = &background_pwm[0];
current_logo++;
np_p[current_logo] = &background_pwm[1];
current_logo++;
if(user_specified_output_file[0] == '\0')
{
strcpy(tempstring, file_name[1]);
strcat(tempstring, "_logo.svg");
}
else strcpy(tempstring, user_specified_output_file);
    
cm.number_of_total_matches = all_hits_pwm[1].max_counts;
two_hits_connecting_matrix.number_of_total_matches = all_hits_pwm[1].max_counts;
two_hits_fold_connecting_matrix.number_of_total_matches = all_hits_pwm[1].max_counts;
cms_to_heatmap[1] = &two_hits_connecting_matrix;
cms_to_heatmap[2] = &two_hits_fold_connecting_matrix;

    
    /* NORMALIZES HIT POSITIONS AND PRINTS THEM */
    for (max_counts = 0, counter2 = 0; counter2 < 2; counter2++) for(counter = 1; counter < Nlength - query_sequence_length; counter++) max_counts += hit_position.incidence[counter2][counter];
    for (counter2 = 0; counter2 < 2; counter2++) for(counter = 1; counter < Nlength - query_sequence_length; counter++) hit_position.fraction[counter2][counter] = ((Nlength - query_sequence_length) * 2 - 2) * (double) hit_position.incidence[counter2][counter] / (double) max_counts;
    
    printf ("\n\n+\n+ Hit positions");
    for(counter = 1; counter <= Nlength - query_sequence_length; counter++) printf ("\t%li", counter);
    printf ("\tInformation content\n+ Forward: ");
    for(counter = 1; counter <= Nlength - query_sequence_length; counter++) 
    {
        if (print_frequencies == 1) printf ("\t%.1f", hit_position.fraction[0][counter]);
        else printf ("\t%li", hit_position.incidence[0][counter]);
    }
    forward_pos_ic = Information_content(hit_position.incidence[0], 1, Nlength - query_sequence_length + 1);
    if (forward_pos_ic > 1) warning++;
    if (forward_pos_ic > 2) warning++;
    printf ("\t%.2f",  forward_pos_ic);
    printf ("\n+ Reverse: ");
    for(counter = 1; counter <= Nlength - query_sequence_length; counter++) 
    {
        if (print_frequencies == 1)  printf ("\t%.1f", hit_position.fraction[1][Nlength - query_sequence_length - counter]);
        else printf ("\t%li", hit_position.incidence[1][Nlength - query_sequence_length - counter]);
    }
    reverse_pos_ic = Information_content(hit_position.incidence[1], 0, Nlength - query_sequence_length - 1);
    if (reverse_pos_ic > 1) warning++;
    if (reverse_pos_ic > 2) warning++;
    printf ("\t%.2f", reverse_pos_ic);
    
    /* PRINTS HIT STATISTICS AND LAMBDA */
    printf ("\n\n\tTotal sequences\tHits\tFraction\tLower half 8mer total\tMax 8mer count");
    if (kmer_count == 1) printf ("\nBackground\t%.0f\t%.0f\t%.0f\t%li\t%li", number_of_sequences_analyzed[0], number_of_sequences_with_hits[0], 100 * number_of_sequences_with_hits[0] / number_of_sequences_analyzed[0], background_kmer_count, results[0][current_kmer_length][0][0][number_of_kmers-1]);
    if (kmer_count == 1) printf ("\nSignal    \t%.0f\t%.0f\t%.0f\t%li\t%li", number_of_sequences_analyzed[1], number_of_sequences_with_hits[1], 100 * number_of_sequences_with_hits[1] / number_of_sequences_analyzed[1], signal_kmer_count, results[1][current_kmer_length][0][0][number_of_kmers-1]);
    /* printf ("\nLower half 8mer counts: background %i, signal %i", background_kmer_count,  signal_kmer_count); */
    
    printf ("\nLambda %.3f", lambda);
    /* printf("\nshortest kmer = %i and too long kmer = %i bak %.0f sig %.0f", shortest_kmer, too_long_kmer, number_of_sequences_analyzed[0], number_of_sequences_analyzed[1]); */

    
/* DINUCLEOTIDE PROPERTIES */
if (dinucleotide_properties == 1) 
{
Expected_dinucleotides (&dep[1], &expected_dinucleotides, searchstring, 1000000000);

    
for (counter = end_trim; counter < dep[1].width - end_trim; counter++) 
{
for (counter2 = end_trim; counter2 < dep[1].width - end_trim; counter2++)
{
ic_score = dep[1].total_relative_deviation[counter][counter2];
    
if (ic_score > max_ic_score) 
{
for (total_exp = 0, total_obs = 0, counter3 = 0; counter3 < 16; counter3++) 
{
total_obs += (double) dep[1].incidence[counter][counter2][counter3];
total_exp += (double) expected_dinucleotides.incidence[counter][counter2][counter3];
}
if (total_obs > minimum_kmer_count)
{
max_ic_score = ic_score;
delta_ic = ic_score;
max_first = counter;
max_second = counter2;
for (min_fold = 1000, max_fold = -1000, counter3 = 0; counter3 < 16; counter3++) 
{
current_fold = (total_exp * (double) dep[1].incidence[counter][counter2][counter3]) / (total_obs * (double) expected_dinucleotides.incidence[counter][counter2][counter3]);
if (current_fold < min_fold) {min_fold = current_fold; min_dinucleotide = counter3;}
if (current_fold > max_fold) {max_fold = current_fold; max_dinucleotide = counter3;}
}
}
}
}
}

(*(np_p[3])).pairwise_correlation[1].first_base = max_first;
(*(np_p[3])).pairwise_correlation[1].second_base = max_second;
(*(np_p[3])).pairwise_correlation[1].max_dinucleotide = max_dinucleotide;
(*(np_p[3])).pairwise_correlation[1].min_dinucleotide = min_dinucleotide;
(*(np_p[3])).pairwise_correlation[1].max_fold_change = max_fold;
(*(np_p[3])).pairwise_correlation[1].min_fold_change = min_fold;
(*(np_p[3])).pairwise_correlation[1].delta_ic = delta_ic;
strcpy((*(np_p[3])).seed, searchstring);
Svg_logo(tempstring, current_logo, np_p, cms_to_heatmap, &one_hit_di[1], &dep[0], &dep[1], &expected_dinucleotides, all_hits_align_scores, lambda, warning);
}
else 
{
Svg_logo(tempstring, current_logo, np_p, cms_to_heatmap, '\0', '\0', '\0', '\0', all_hits_align_scores, lambda, warning);
}


}

   
	/* OUTPUTS KMER TABLE */
	if (kmer_table == 1 && flank != 1)
	{
	short int *kmermatch_position;
	kmermatch_position = malloc(sizeof(short int) * 2 + 5);
	current_gap_position = current_kmer_length / 2;
	current_gap_length = 0;
	last_kmer[0] = pow(4, (too_long_kmer-1));
	for (file_number = 0; file_number < number_of_files; file_number++) 
	{
        
    for(current_gap_length = 0; current_gap_length < loaded_pwm_width-current_kmer_length + 2; current_gap_length++)
    {
	for (current_kmer_length = shortest_kmer; current_kmer_length < too_long_kmer; current_kmer_length++) 
	{
	last_kmer[current_kmer_length] = pow(4, current_kmer_length);
    if (file_number == 1) printf("\nKMER\tsequence");
    else printf("\nBACKGROUND KMER\tsequence");
    printf("\t%imer\t%imer_expected", current_kmer_length, current_kmer_length);
	}
	if (file_number == 1) printf("\t::: %i:", current_gap_length);
	else printf("\t::B %i:", current_gap_length); 
	for (current_kmer = 0; current_kmer < last_kmer[0]; current_kmer++)
	{
	printf("\n%li\t", current_kmer);
	for(position = too_long_kmer - 2; position >= shortest_kmer & current_kmer <= last_kmer[position]; position--) 
    {
    printf("%c", forward_lc[(current_kmer & mask_ULL[1][position]) >> (position * 2)]);
    if(position == current_kmer_length / 2) for(counter=0;counter < current_gap_length;counter++) printf("n");
    }
    for(; position > -1 ; position--) 
    {
    printf("%c", forward[(current_kmer & mask_ULL[1][position]) >> (position * 2)]);
    if(position == current_kmer_length / 2) for(counter=0;counter < current_gap_length;counter++) printf("n");
	}
    for (current_kmer_length = shortest_kmer; current_kmer_length < too_long_kmer; current_kmer_length++) 
	{
	if (current_kmer < last_kmer[current_kmer_length])
	{
    /* TESTS IF CURRENT SEQUENCE IS A PALINDROME */
    palindrome_correction = 1 + (Is_this_sequence_dimer(current_kmer, current_kmer_length) == 4);
    if (current_gap_length == 0) kmer_incidence = results[file_number][current_kmer_length][0][0][current_kmer];
    else kmer_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][current_kmer];
	printf("\t%li", kmer_incidence/palindrome_correction);
	printf("\t%.2f", Kmerscore (&qp, current_kmer, current_kmer_length, kmermatch_position));
    printf("\t%.2f", fastgappedKmerscore (qp.fraction, qp.width, current_kmer & lowmask_ULL[current_kmer_length/2-1], (current_kmer >> current_kmer_length) & lowmask_ULL[current_kmer_length/2-1],current_kmer_length/2, current_gap_length));
    if(load_adm == 1) printf("\t%.2f", Kmerscore_ADM (&flanked_adm, current_kmer, current_kmer_length, kmermatch_position));
	}
	else printf("\t-\t-");
	}
    if (file_number == 1) printf("\t::: %i:", current_gap_length);
    else printf("\t::B %i:", current_gap_length);  
	}
    }
	}
	}


    
    
	/* COUNTS INFORMATION CONTENT OF HIT AND NON-HIT K-MER DISTRIBUTIONS */
	current_gap_position = 0;
	current_gap_length = 0;

    for (current_kmer_length = shortest_kmer; current_kmer_length < too_long_kmer; current_kmer_length++)
	{
	background_info = 2 * current_kmer_length;
	signal_info = 2 * current_kmer_length;
	if (flank == 1) 
	{
	background_nonhit_info = 2 * current_kmer_length;
	signal_nonhit_info = 2 * current_kmer_length;
	}
	number_of_kmers = pow(4,current_kmer_length);
	for (counter = 0; counter < number_of_kmers; counter++) 
	{
	fraction = results[1][current_kmer_length][0][0][counter] / (number_of_sequences_analyzed[1] * 2 * (Nlength - current_kmer_length));
	if (results[1][current_kmer_length][0][0][counter] != 0) signal_info += fraction * log2 (fraction);
	fraction = results[0][current_kmer_length][0][0][counter] / (number_of_sequences_analyzed[0] * 2 * (Nlength - current_kmer_length));
	if (results[0][current_kmer_length][0][0][counter] != 0) background_info += fraction * log2 (fraction);
	if (flank == 1)
	{
	fraction = results[1][current_kmer_length][1][0][counter] / (number_of_sequences_with_no_hits[1] * 2 * (Nlength - current_kmer_length));
	if (results[1][current_kmer_length][1][0][counter] != 0) signal_nonhit_info += fraction * log2 (fraction);
	fraction = results[0][current_kmer_length][1][0][counter] / (number_of_sequences_with_no_hits[0] * 2 * (Nlength - current_kmer_length));
	if (results[0][current_kmer_length][1][0][counter] != 0) background_nonhit_info += fraction * log2 (fraction);
	}
	}
	printf ("\n\nInformation content for kmer length %i\tAll sequences\tSequences not hit\nSignal          \t%.2f\t%.2f\tbits\nBackground      \t%.2f\t%.2f\tbits", current_kmer_length, signal_info, signal_nonhit_info, background_info, background_nonhit_info);
	}


/* CALCULATES AND PRINTS ALL INFORMATION CONTENTS */
if (information_content_output == 1 && flank != 1)
{
	double signal_information_content[too_long_kmer][too_long_kmer][too_long_kmer];
	double background_information_content[too_long_kmer][too_long_kmer][too_long_kmer];
	long int signal_count_sum[too_long_kmer][too_long_kmer][too_long_kmer];
	long int background_count_sum[too_long_kmer][too_long_kmer][too_long_kmer];
	printf("\n\nINFORMATION CONTENT\tKmer\tGap_position\tGap_length\tSignal_bits\tBackground_bits\tSignal_counts\tBackground_counts\tIC");
	
	for(current_kmer_length = shortest_kmer; current_kmer_length < too_long_kmer; current_kmer_length++)
	{
		for(current_gap_position = 0; current_gap_position < current_kmer_length; current_gap_position++)
		{
			for(current_gap_length = (current_gap_position > 0); current_gap_length < Nlength-current_kmer_length; current_gap_length++)
			{
			signal_information_content[current_kmer_length][current_gap_position][current_gap_length] = 2 * current_kmer_length;
			background_information_content[current_kmer_length][current_gap_position][current_gap_length] = 2 * current_kmer_length;
			signal_count_sum[current_kmer_length][current_gap_position][current_gap_length] = 0;
			background_count_sum[current_kmer_length][current_gap_position][current_gap_length] = 0;
			if (current_gap_position == 0 & current_gap_length != 0) continue;
			
				for(current_kmer = 0; current_kmer < pow(4, current_kmer_length); current_kmer++)
				{
					signal_count_sum[current_kmer_length][current_gap_position][current_gap_length] += results[1][current_kmer_length][current_gap_position][current_gap_length][current_kmer];
					background_count_sum[current_kmer_length][current_gap_position][current_gap_length] += results[0][current_kmer_length][current_gap_position][current_gap_length][current_kmer];
					fraction = results[1][current_kmer_length][current_gap_position][current_gap_length][current_kmer] / (number_of_sequences_analyzed[1] * 2 * (Nlength - current_kmer_length - current_gap_length));
					/* printf("\nkmer length %i,position %i,gap length %i, current_kmer %i, fraction %f", current_kmer_length, current_gap_position, current_gap_length, current_kmer, fraction); */
					if (results[1][current_kmer_length][current_gap_position][current_gap_length][current_kmer] != 0) signal_information_content[current_kmer_length][current_gap_position][current_gap_length] += fraction * log2 (fraction);
					fraction = results[0][current_kmer_length][current_gap_position][current_gap_length][current_kmer] / (number_of_sequences_analyzed[0] * 2 * (Nlength - current_kmer_length - current_gap_length));
					if (results[0][current_kmer_length][current_gap_position][current_gap_length][current_kmer] != 0) background_information_content[current_kmer_length][current_gap_position][current_gap_length] += fraction * log2 (fraction);
				}
				printf("\n\t%i\t%i\t%i\t%.3f\t%.3f\t%li\t%li\tIC", current_kmer_length, current_gap_position, current_gap_length, signal_information_content[current_kmer_length][current_gap_position][current_gap_length], background_information_content[current_kmer_length][current_gap_position][current_gap_length], signal_count_sum[current_kmer_length][current_gap_position][current_gap_length], background_count_sum[current_kmer_length][current_gap_position][current_gap_length]);

			}
		}
	}
}



/* PRINTS DINUCLEOTIDE MATRIX */
if (dinucleotide_properties == 1)
{
for (counter = 0; counter < (one_hit_di[1]).number_of_dinucleotide_properties; counter++) 
{
printf ("\n%s", (one_hit_di[1]).dinucleotide_property_string[counter]);
for (counter2 = 0; counter2 < (one_hit_di[1]).width; counter2++) printf ("\t%.2f", (one_hit_di[1]).score[counter][counter2]);
}



/*
for (counter = 0; counter < np[2].width; counter++)
{
printf("\n%i", counter);
for (counter2 = 0; counter2 < np[2].width; counter2++) printf("\t%.2f", dep.information_content[counter][counter2] - expected_dinucleotides.information_content[counter][counter2] );
}
*/
char *in_out_string = "oi";
short int first_is_in = 0;    
short int second_is_in = 0; 
short int first_pos_relative_to_seed = 0;
short int second_pos_relative_to_seed = 0;
long int max_obs;
long int current_total_dinuc_count;
long int current_total_dinuc_background_count;
double difference;
    
/* PRINTS DINUCLEOTIDE DATA */
    printf("\ndinucleotide\tpos1\tpos2\tmononuc1_fraction\tmononuc2_fraction\tdinuc_fraction\tdinuc_background_fraction\tcount\tbackground_count\texpected(mono)\texpected(h1)\tobsIC\texpIC\tdeltaIC\tscore\tfold_dev\ttotal_counts\tmax_count\ttotal_rel_dev\texp_min_deviation\teo_correlation\tpermutated_correlation\tic\tin/out\tCpG");
for (counter = 0; counter < dep[1].width; counter++) 
{
for (counter2 = 0; counter2 < dep[1].width; counter2++)
{
/* CALCULATES SUMMARY STATISTICS */
for (/*total_relative_deviation[counter][counter2] = 0, uncentered_correlation[counter][counter2] = 0, sum_x_squared = 0, sum_y_squared = 0, sum_xy = 0, */ max_obs = 0, total_exp = 0, total_obs = 0, counter3 = 0; counter3 < 16; counter3++) 
{
total_obs += (double) dep[1].incidence[counter][counter2][counter3];
total_exp += (double) expected_dinucleotides.incidence[counter][counter2][counter3]; 
if(dep[1].incidence[counter][counter2][counter3] > max_obs) max_obs = dep[1].incidence[counter][counter2][counter3];
 /*   printf("\n%f\t%f", total_obs, total_exp); */
/* CALCULATES TOTAL RELATIVE DEVIATION
total_relative_deviation[counter][counter2] += abs(total_obs * ((double) expected_dinucleotides.incidence[counter][counter2][counter3] / total_exp) - dep[1].incidence[counter][counter2][counter3])/total_obs; */
}
/*for (counter3 = 0; counter3 < 16; counter3++) total_relative_deviation[counter][counter2] += abs((double) expected_dinucleotides.incidence[counter][counter2][counter3] - (double) dep[1].incidence[counter][counter2][counter3]) / total_obs;*/
    
/* PRINTS LARGE TABLE */
for (current_total_dinuc_count = 0, current_total_dinuc_background_count = 0, counter3 = 0; counter3 < 16; counter3++) 
{
current_total_dinuc_count += dep[1].incidence[counter][counter2][counter3];
current_total_dinuc_background_count += dep[0].incidence[counter][counter2][counter3]; 
}
for (counter3 = 0; counter3 < 16; counter3++) 
{
first_pos_relative_to_seed = counter2-Nlength+2;
second_pos_relative_to_seed = counter-Nlength+2;
    
/* DETERMINES IF POSITIONS ARE INSIDE SEED */
first_is_in = 0;
second_is_in = 0;
if (first_pos_relative_to_seed > 0 && first_pos_relative_to_seed <= query_sequence_length) first_is_in = 1; 
if (second_pos_relative_to_seed > 0 && second_pos_relative_to_seed <= query_sequence_length) second_is_in = 1;
    if (counter2 < counter) 
    {
    printf("\n%c%c\t%i\t%i\t_%i:%i_\t%.2f\t%.2f\t%.2f\t%.2f\t%li\t%li\t%.1f\t%.1f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.0f\t%li\t%.3f\t%.3f\t%.3f\t%.3f\tic\t%c%c\t", forward[(counter3 & 12) >> 2], forward[counter3 & 3], first_pos_relative_to_seed, second_pos_relative_to_seed, first_pos_relative_to_seed, second_pos_relative_to_seed, np[2].fraction[(counter3 & 12) >> 2][counter2], np[2].fraction[(counter3 & 3)][counter], ((double) dep[1].incidence[counter][counter2][counter3] / current_total_dinuc_count), ((double) dep[0].incidence[counter][counter2][counter3] / current_total_dinuc_background_count), dep[1].incidence[counter][counter2][counter3], dep[0].incidence[counter][counter2][counter3], total_obs * ((double) expected_dinucleotides.incidence[counter][counter2][counter3] / total_exp), Expected_mismatched_dinucleotides(counter3, dep[1].incidence[counter][counter2]), dep[1].information_content[counter][counter2], expected_dinucleotides.information_content[counter][counter2], dep[1].information_content[counter][counter2] - expected_dinucleotides.information_content[counter][counter2], dep[1].information_content[counter][counter2] * (dep[1].information_content[counter][counter2] - expected_dinucleotides.information_content[counter][counter2]), (total_exp * (double) dep[1].incidence[counter][counter2][counter3]) / (total_obs * (double) expected_dinucleotides.incidence[counter][counter2][counter3]), total_obs, max_obs, dep[1].total_relative_deviation[counter][counter2], dep[1].count_statistic_expected_total_relative_deviation[counter][counter2], dep[1].eo_correlation[counter][counter2], dep[1].permutated_correlation[counter][counter2], in_out_string[first_is_in], in_out_string[second_is_in]);
        
        difference = (double) dep[1].incidence[counter][counter2][counter3] / current_total_dinuc_count - (double) dep[0].incidence[counter][counter2][counter3] /current_total_dinuc_background_count;
        if (difference < -0.2) printf ("low");
        if (difference > 0.2) printf ("high");

    }
if ((counter3 == 6) && counter - counter2 == 1) printf ("\tCpG");

}
    
}
}

   
printf("\n\nSum of relative deviations between expected and observed dinucleotide counts\npos1\tpos2\tdeviation\texp_min_deviation\tuncentered correlation\texcess information\tscore\trdv");
for (counter2 = Nlength - 1; counter2 - Nlength <= query_sequence_length; counter2++) 
{
for (counter = counter2 + 1; counter - Nlength + 2 <= query_sequence_length; counter++)
{
first_pos_relative_to_seed = counter2-Nlength+2;
second_pos_relative_to_seed = counter-Nlength+2;
    printf("\n%i\t%i\t%.3f\t%.3f\t%.3f\t%.3f\t%.2f\t%.2f\trdv", first_pos_relative_to_seed, second_pos_relative_to_seed, dep[1].total_relative_deviation[counter][counter2], dep[1].count_statistic_expected_total_relative_deviation[counter][counter2], dep[1].eo_correlation[counter][counter2], dep[1].permutated_correlation[counter][counter2], dep[1].information_content[counter][counter2] - expected_dinucleotides.information_content[counter][counter2], dep[1].information_content[counter][counter2] * (dep[1].information_content[counter][counter2] - expected_dinucleotides.information_content[counter][counter2]));
    /* printf("\t%f",dep[1].total_relative_deviation[counter][counter2]); */
}
}

    struct adjacent_dinucleotide_model adm;
    adjacent_dinucleotide_model_init(&adm, "adm", query_sequence_length);
    Generate_ADM(&adm, &dep[1], query_sequence_length);
    Print_ADM(&adm);
    strcpy(tempstring, file_name[1]);
    strcat(tempstring, "_ADM");
    strcpy(adm.name, tempstring);
    strcat(tempstring, "_logo.svg");
    Svg_logo_ADM (tempstring, 0, 0, &adm, 0.05, 0.1, 1000, 0.1);
    
}
                   
/*
for (counter = 0; counter < (one_hit_di[1]).number_of_dinucleotide_properties; counter++) 
{
printf ("\n%s", (one_hit_di[1]).dinucleotide_property_string[counter]);
for (counter2 = 0; counter2 < (one_hit_di[1]).width; counter2++) printf ("\t%i", (one_hit_di[1]).count[counter][counter2]);
}
*/
    double bak_correl;
    double sig_correl;
    /* PRINTS FLANK KMER DATA */
    if(flank_kmer_pos != -100)
    {
    bak_correl = Kmer_mono_correlation(flank_kmer_count[0], flank_kmer_expected_count[0], 4);
    sig_correl = Kmer_mono_correlation(flank_kmer_count[1], flank_kmer_expected_count[1], 4);
    printf("\n\nFLANK KMERS AT POSITION %i", flank_kmer_pos);
    printf("\nKmer\tback_obs\tback_exp\tsignal_obs\tsignal_exp");
    for(counter=0;counter<256;counter++) printf("\n%c%c%c%c\t%li\t%.1f\t%li\t%.1f\tfK", forward[counter >> 6], forward[(counter >> 4) & 3], forward[(counter >> 2) & 3], forward[counter & 3], flank_kmer_count[0][counter], flank_kmer_expected_count[0][counter], flank_kmer_count[1][counter], flank_kmer_expected_count[1][counter]);
    printf("\nCorrelation with mononucleotide-based expected kmers:\t%.3f\t%.3f", bak_correl, sig_correl);
    }

    short int gap;
    /* PRINTS ALIGN SCORES */
    if(align_matches == 1)
    {
        /* for(counter = 0; counter < Nlength; counter++) for (counter2 = 0; counter2 < Nlength; counter2++)
        {
            printf("\n%i,%i\t", counter, counter2);
            Seqprint(mask_ULL[counter][counter2], Nlength-1);    
        } */
        printf("\n\nALIGN SCORES ::a\n");
        for(counter = 0; counter < Nlength*2-2; counter++) printf("%li\t",counter-Nlength+1);
        printf("::a\n");
        for(strand = 0; strand < 2; strand++)
        {
            for(counter = 0; counter < Nlength*2-2; counter++) printf("\t %.1f", ((double) all_hits_align_scores[1].score[counter][strand])/((double) all_hits_align_scores[1].count[counter][strand]+0.01));
            printf("\t::a\n");
        }
        printf("\n\nDirect Repeat abundance");
        for(counter = 0; counter < Nlength; counter++) 
        {
        printf ("\n%li", counter); 
        for(gap = 0; gap < Nlength; gap++) printf("\t%li", all_hits_align_scores[1].direct_repeat[counter][gap]);
        }
        printf("\n\nInverted Repeat abundance");
        for(counter = 0; counter < Nlength; counter++) 
        {
        printf ("\n%li", counter); 
        for(gap = 0; gap < Nlength; gap++) printf("\t%li", all_hits_align_scores[1].inverted_repeat[counter][gap]);
        }
        
printf("\n\nDirect repeat length histogram\tSignal\tBackground\tRatio");      
        for(counter = 0; counter < Nlength; counter++) 
        {
            printf("\n%li\t%li\t%li\t%.1f", counter, align_score_histogram[1][0][counter], align_score_histogram[0][0][counter], (double) align_score_histogram[1][0][counter]/((double) align_score_histogram[0][0][counter] + 0.001));
        }

printf("\n\nInverted repeat length histogram\tSignal\tBackground\tRatio");      
        for(counter = 0; counter < Nlength; counter++) 
        {
            printf("\n%li\t%li\t%li\t%.1f", counter, align_score_histogram[1][1][counter], align_score_histogram[0][1][counter], (double) align_score_histogram[1][1][counter]/((double) align_score_histogram[0][1][counter] + 0.001));
        }
        
    }
    
    if (xyplot == 1) 
    {
    if (user_specified_output_file[0] == '\0') 
    {
    strcpy(searchstring, file_name[0]);
    strcat(searchstring, "_");
    strcat(searchstring, file_name[1]);
    }
    else strcpy(searchstring, user_specified_output_file);
    strcat(searchstring, "_");
    sprintf(valuestring, "%i", shortest_kmer);
    strcat(searchstring, valuestring);   
    strcat(searchstring, "mer_xyplot.svg");
    XYplot(searchstring, results, shortest_kmer, count_also_spaced_kmers, kmer_length_difference_cutoff, file_name[0], file_name[1], tempstring);
    }
    printf ("\n\nPalindromic hits\nBackground:\t%li\nSignal    :\t%li", palindromic_hits[0],palindromic_hits[1]);
	printf ("\n");
	fflush(stdout);

	t1 = time('\0');
	printf ("\nTime: %ld seconds\n", (long) (t1 - t0));
}
