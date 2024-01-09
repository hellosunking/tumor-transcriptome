#!/usr/bin/perl

#
# Author: Ahfyth
#

use strict;
use warnings;

if( $#ARGV < 0 )
{
	print STDERR "\nUsage: $0 <survival.info> [Cancer=UnknowCancer] [pval=0.05]\n";
	print STDERR "\nThis program is designed to \n\n";
	exit 2;
}

my $cancer = $ARGV[1] || 'UnknowCancer';
my $cutoff = $ARGV[2] || 0.05;

my %stat;

print "Gene\tSymbol\tDEG.Status\tSurv.Pval\tSurv.Status\tHR\n";
open IN, "$ARGV[0]" or die( "$!" );
<IN>;	## header
while( <IN> )
{
	chomp;
	next unless /\S+/;	
	my @l = split /\t/;    #Gene.Symbol DEG     pval    obs1    obs2    exp1    exp2    HR
				#ENSG00000000005.6.TNMD  Down    0.131639188168189       5       10      7.87849874557594        7.12150125442406        2.25497492697869

	$l[0]=~/(ENSG\d+\.\d+)\.(\w+)/;
	my $gid=$1;
	my $sym=$2;
	my $flag = '';
	if( $l[2] > $cutoff )
	{
		$flag = 'NA';
	}
	else
	{
		if( $l[3]>$l[5] && $l[4]<$l[6] )	## the higher the expression, the easier to die
		{
			if( $l[1] =~ /up/i )
			{
				$flag = 'Oncogene';
			}
			else
			{
				$flag = 'DownSaver';
			}
		}
		elsif( $l[3]<$l[5] && $l[4]>$l[6] )	## the lower the expression, the easier to die
		{
			if( $l[1] =~ /down/i )
			{
				$flag = 'TumorSuppressor';
			}
			else
			{
				$flag = 'UpSaver';
			}
		}
		else	## inconsistent data
		{
			$flag = 'BADDATA';
		}
	}
	print "$1\t$2\t$l[1]\t$l[2]\t$flag\t$l[7]\n" unless $flag eq 'NA';
	$stat{$flag} ++;
}
close IN;

my @type = sort keys %stat;
my @cnt  = map {$stat{$_}} @type;

print STDERR join("\t", "Type", @type), "\n",
			join("\t", $cancer, @cnt),"\n";

