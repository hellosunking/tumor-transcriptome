#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;
#use KSLIB::loadGenome qw/loadGenome/;
#use KSLIB::cigarUtil qw/fix_seq_from_CIGAR extract_size_from_CIGAR/;

## ENSEMBL gene info

my %info;
open IN, "/lustre/home/dxhu/project/project4_desurv/paper_figure/all.program/01.SAHR/data/ENSEMBL.v101.info" or die( "$!" );
while( <IN> ) {
	chomp;
	my @l = split /\t/;     ##ENSG00000237613       FAM138A lncRNA
	$info{$l[0]} = "$l[1]\t$l[2]";
}
close IN;

my $gen=$ARGV[0];
my $degs = load_DESeq2_tsv( $gen );

print "#Gene\tSymbol\tBiotype\tDEG\tfc\tp\n";

foreach my $g ( sort keys %$degs ) {
	my ($fc, $p) = split /\t/, $degs->{$g};
	next unless (abs($fc)>1 && $p<0.01);	
        my $deg = ($fc > 0) ? "Up" : "Down";
	$g =~ /^(ENSG\d+)/;
	my $eid = $1;
	my ( $sym, $bio ) = ("NA", "NA");
	($sym, $bio) = split /\t/, $info{$eid} if exists $info{$eid};

	print join("\t", $g, $sym, $bio, $deg, $fc, $p), "\n";
}

sub load_DESeq2_tsv {
	my $file = shift;

	my %DEG;
	open IN, "$file" or die( "$!" );
	<IN>;	## header;
	while( <IN> ) {
		chomp;
		my @l = split /\t/;	##Gene	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj
		next if $l[-1] eq 'NA';
		$DEG{$l[0]} = "$l[2]\t$l[-1]";
	}
	close IN;

	return \%DEG;
}

