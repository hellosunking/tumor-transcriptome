#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <in.DEseq2.tsv> [Criteria:std|loose|qvalue|Q-val=0.01] [log2FC=1]\n\n";
	print STDERR "Log(1.5) = 0.5849625\n\n";
	exit 2;
}

my ($qval, $fc) = ( 0.01, 1 );

if( $#ARGV > 0 ) {
	if( $ARGV[1] eq 'loose' ) {	## loose criteria
		$qval = 0.05;
		$fc   = log(1.5)/log(2);
	} elsif( $ARGV[1] eq 'standard' || $ARGV[1] eq 'std' ){
		$qval = 0.01;
		$fc   = 1;
	} elsif( $ARGV[1] eq 'qvalue' ) {	## ultra-relax mode
		$qval = 0.05;
		$fc   = 0;
	}else {
		$qval = $ARGV[1];
		$fc   = $ARGV[2];
	}
}

open IN, "$ARGV[0]" or die( "$!" );
my $header = <IN>;
#print "Gene$header";
print "Gene\tlog2FC\tP-adj\n";
while( <IN> ) {
	chomp;
	my @l = split /\t/;	##Gene	baseMean	log2FoldChange	lfcSE	stat	qvalue	padj
	next if $l[-1]=~/NA/i;
#	print "$_\n" if abs($l[2])>$fc && $l[-1]<$qval;
	print "$l[0]\t$l[2]\t$l[-1]\n" if abs($l[2])>$fc && $l[-1]<$qval;
}
close IN;


