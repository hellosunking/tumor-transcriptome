#!/usr/bin/perl

#
# Author: Ahfyth
#

use strict;
use warnings;


my (%surv, %HR);
open IN, "$ARGV[0]" or die( "$!" );
while( <IN> )
{
	chomp;
	my @l = split /\t/;	##KIF21A	Oncogene	pval-hr	HR
	$surv{$l[0]} = $l[1];
	$HR{$l[0]}   = $l[3];
}
close IN;

my @genes = sort keys %surv;

#my %cutoff;	## expression cutoff for each gene, which is used to update the model
my %anno;	## expression category of all the genes
my %sampleinfo;
foreach my $g ( @genes )
{
#	my $type = $surv{$g};
#	print STDERR "Loading $g: $type\n";
	open IN, "$ARGV[1]/$g" or die( "$!:$ARGV[1]" );
	<IN>;	## header
	while( <IN> )
	{
		chomp;			#SampleID        Type    FPKM    Category        Vital   Days    Age     Stage
		my @l = split /\t/;	##barcode	vital	days	expression
		next if $l[5] eq "NA";
		next if $l[5] <= 0;
		$anno{$l[0]}->{$g} = $l[3];
		$sampleinfo{$l[0]} = "$l[4]\t$l[5]";
	}
	close IN;
}

open OUT, ">$ARGV[2]" or die( "$!" );
print OUT join("\t", "#Sample", "Vital", "Days", @genes), "\n";
foreach my $s ( sort keys %sampleinfo )
{
	my @expinfo = map { $anno{$s}->{$_} || 'X' } @genes;
	print OUT join("\t", $s, $sampleinfo{$s}, @expinfo), "\n";
}
close OUT;

