#!/usr/bin/perl
##
## Date  :
#
use strict;
use warnings;

my @case;
my @control;
my $tNum;
my $nNum;
open IN, "$ARGV[0]" or die( "$!" );

while( <IN> ) {
        chomp;
        my @l = split /\t/; # SID     Vital   Days    Age     Stage   Gender  Race    Subtype
        		#TCGA-CC-5258    1       129     17586   Stage_II        male    asian   C0
	if( -s "~/TCGA/RNA-seq.new.pipeline.2022/HCC/Count/$l[0]-01A.gz" ) {
                if($l[7]=~/C[01]/){
                push @control, "$l[0]";
                 ++ $nNum;
                  }
                elsif($l[7]=~/C2/) {
                push @case, "$l[0]";
                ++ $tNum;
                }
          }
}

open OUT,">$ARGV[1]" or die ( "$!");

print OUT "$nNum\t$tNum\n";


close IN;

my @all=(@control,@case);
my @exp;
my %cnt;
foreach my $f ( @all ) {
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

print join("\t", "Gene", @all), "\n";
foreach my $g (sort keys %cnt) {
        if( $cnt{$g} != $#all+1 ) {
                print STDERR "WARNING: $g is invalid!\n";
        } else {
                my @here = map { $_->{$g} } @exp;
                print join("\t", $g, @here), "\n";
        }
}

