#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;
#use KSLIB::loadGenome qw/loadGenome/;
#use KSLIB::cigarUtil qw/fix_seq_from_CIGAR extract_size_from_CIGAR/;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <in.gene><cancer.type> [dir=.] [filter.zero=0|1]\n\n";
#	print STDERR "\nThis program is designed to \n\n";
	exit 2;
}

our $target  = $ARGV[0];
our $dataDIR = $ARGV[2] || '.';
our $filter0 = $ARGV[3] || 0;

our $cutoff;

## load survival info
our %surv;
open IN, "$dataDIR/clinical.info.training" or die( "$!" );
<IN>;
while( <IN> ) {
	chomp;
	my @l = split /\t/, $_;#case_id vital_status days_for_survival age_at_diagnosis ajcc_pathologic_t primary_diagnosis
	next if $l[2] eq "NA";
	if( $l[1] eq "Dead" ) {
		$surv{$l[0]} = "1\t$l[2]\t$l[3]\t$l[4]";
	} elsif( $l[1] eq "Alive" ) {
		next unless $l[2] >0;	## discard records with 0 following-up days
		$surv{$l[0]} = "0\t$l[2]\t$l[3]\t$l[4]";
	}
}
close IN;

print "SampleID\tType\tFPKM\tCategory\tVital\tDays\tAge\tStage\n";

## load adjacent normal
my @adj;
open IN, "$dataDIR/normal.matrix" or die( "$!" );
my $header = <IN>;
chomp( $header );
my @samples = split /\t/, $header;
while( <IN> ) {
	chomp;
	my @l = split /\t/, $_;
	next unless $l[0] eq $target;
	@adj = @l;
}
close IN;

foreach my $i ( 1..$#adj ) {
	print join("\t", $samples[$i], "Adj", $adj[$i], "L", 0, 0, 0, "Adj"), "\n"
}

## load tumors
process_FPKM( "$dataDIR/tumor.matrix", $ARGV[1] );
#process_FPKM( "$dataDIR/non-asian.HTseq.FPKM.matrix", "White" );

sub process_FPKM {
	my $file = shift;
	my $code = shift;

	my %fpkm;
	my $zero = 0;
	my $matched = 0;
	open IN, "$file" or die( "$!" );
	my $header = <IN>;
	chomp( $header );
	my @samples = split /\t/, $header;
	while( <IN> ) {
		chomp;
		my @l = split /\t/, $_;
		next unless $l[0] eq $target;
		foreach my $i ( 1..$#samples ) {
			my $id = $samples[$i];
			$id =~ s/-\d\d\S$//;
			if( exists $surv{$id} ) {
				$fpkm{$id} = $l[$i];
				++ $zero if $l[$i] < 1e-10;
				++ $matched;
			}
		}
		last;
	}
	close IN;

	#my $cutoff;
	if( $zero > $matched / 3 ) {	## if more than 1/3 of the values are zero, then use zero as cutoff
		$cutoff = 1e-10;
	} else {
		my @sorted = sort {$a<=>$b} values %fpkm;
		$cutoff = $sorted[ $matched/2-1 ];
	}
	foreach my $id ( keys %fpkm ) {
		my $cat = ($fpkm{$id} > $cutoff) ? 'H' : 'L';
		print join("\t", $id, $code, $fpkm{$id}, $cat, $surv{$id}), "\n";
	}
}

open OUT, ">>../Gene.exp.cutoff" or die ("$!");
print OUT $target,"\t",$cutoff,"\n";
