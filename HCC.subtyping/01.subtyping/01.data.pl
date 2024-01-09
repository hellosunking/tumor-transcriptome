#!/usr/bin/perl
##
## Date  :
#
use strict;
use warnings;

my @case;
my $tNum;
open IN, "$ARGV[0]" or die( "$!" );
while( <IN> ) {
        chomp; # 
        my @l = split /\t/; # SID     Vital   Days    Age     Stage   Gender  Race 
	if( -s "~/TCGA/RNA-seq.new.pipeline.2022/HCC/Count/$l[0]-01A.gz" ) {
               push @case, "$l[0]";
                ++ $tNum;
        }
}

close IN;

my @exp;
my %cnt;
foreach my $f ( @case ) {
        my %e;
        open COUNT, "zcat ~/TCGA/RNA-seq.new.pipeline.2022/HCC/Count/$f-01A.gz |" or die( "$!" );
        while( <COUNT> ) {
                chomp;                  ##ENSG00000000003.15      TSPAN6  protein_coding  5452    2731    2721    97.9460 29.4016 33.2659
                my @l = split /\t/;     ##gene_id gene_name       gene_type       unstranded      stranded_first  stranded_second tpm_unstranded  fpkm_unstranded fpkm_uq_unstranded
                next unless $l[0]=~/^ENSG/;
		next unless $l[2]=~/^protein_coding/ || $l[2]=~/^lncRNA/;
		next if $l[1]=~/^MT-/ || $l[1]=~/^RP[SL]/ || $l[1]=~/^MRP[SL]/;
                $e{$l[0]} = $l[3];
                $cnt{$l[0]} ++;
       		
	}
        close COUNT;

        push @exp, \%e;
}

print join("\t", "Gene", @case), "\n";
foreach my $g (sort keys %cnt) {
        if( $cnt{$g} != $#case+1 ) {
                print STDERR "WARNING: $g is invalid!\n";
        } else {
                my @here = map { $_->{$g} } @exp;
                print join("\t", $g, @here), "\n";
        }
}

