#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;
#use KSLIB::loadGenome qw/loadGenome/;
#use KSLIB::cigarUtil qw/fix_seq_from_CIGAR extract_size_from_CIGAR/;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <in.DESeq2.filtered> [anno=hg38|mm10]\n\n";
	print STDERR "\nThis step is NOT NEEDED for newer versions of RNAseq pipeline using RefSeq annotations.\n\n";
	exit 2;
}

my $anno = $ARGV[1] || 'hg38';

my %id2sym;
open IN, "/lustre/home/ksun/Genomes/$anno/ENSEMBL.v101.info" or die( "$!" );
while( <IN> ) {
	chomp;
	my @l = split /\t/;	##ENSG00000223972	DDX11L1	transcribed_unprocessed_pseudogene
	$id2sym{$l[0]} = $l[1]."\t".$l[2];
}
close IN;

open IN, "$ARGV[0]" or die( "$!" );
<IN>;	## header
print "#Symbol\tProtein\tENSID\tlog2FC\tP-adj\tDEG\n";
while( <IN> ) {
	chomp;
	my @l = split /\t/;	##Gene	log2FC	P-adj
	$l[0] =~ s/\.\d+$//;
	my $sym = $id2sym{$l[0]} || 'NA';
	my $deg = ($l[1] > 0) ? "Up" : "Down";
	print "$sym\t$l[0]\t$l[1]\t$l[2]\t$deg\n";
}
close IN;


