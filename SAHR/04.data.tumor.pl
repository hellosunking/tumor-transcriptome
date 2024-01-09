#!/usr/bin/perl
##
## Date  :
#
use strict;
use warnings;

my $cancer = "$ARGV[1]";
my $data_dir  = "$ARGV[2]";

my @tumor;
my @normal;
open IN, "$ARGV[0]" or die( "$!" );

while( <IN> ) {
        chomp;
        my @l = split /\t/;
#case_id  
###TCGA-A2-A0CT 01A/11A
       # if( -s "/lustre/home/dxhu/project/project1_1_GeneCT_test/data/RNA-seq.new.pipeline.2022/$cancer/$data_dir/$l[0]-11A.gz" ) {
        # push @normal, "$l[0]-11A";
	#	}
	if( -s "$data_dir/$cancer/Count/$l[0]-01A.gz" ) {
	push @tumor, "$l[0]-01A";		
		}        

}

close IN;

my @all=@tumor;
my @exp;
my %cnt;
foreach my $f ( @all ) {
        my %e;
        open COUNT, "zcat $data_dir/$cancer/Count/$f.gz |" or die( "$!" );
        while( <COUNT> ) {
                chomp;                  ##ENSG00000000003.15      TSPAN6  protein_coding  5452    2731    2721    97.9460 29.4016 33.2659
                my @l = split /\t/;     ##gene_id gene_name       gene_type       unstranded      stranded_first  stranded_second tpm_unstranded  fpkm_unstranded fpkm_uq_unstranded
                next unless $l[0]=~/^ENSG/;
		next unless $l[2]=~/^protein_coding/ || $l[2]=~/^lncRNA/;
		next if $l[1]=~/^MT-/ || $l[1]=~/^RP[SL]/ || $l[1]=~/^MRP[SL]/;
                $e{$l[0]} = $l[7];
                $cnt{$l[0]} ++;
       		
	}
        close COUNT;

        push @exp, \%e;
}

print join("\t", "Gene", @all), "\n";
foreach my $g (sort keys %cnt) {
        if( $cnt{$g} != $#all+1 ) {
                print STDERR "WARNING: $g is invalid!\n";
        } else {
                my @here = map { $_->{$g} } @exp;
                print join("\t", $g, @here), "\n";
        }
}

