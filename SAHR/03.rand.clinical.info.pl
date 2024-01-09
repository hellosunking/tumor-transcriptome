#!/usr/bin/perl

#
# Author: Ahfyth
#

use strict;
use warnings;

if( $#ARGV < 1 )
{
	print STDERR "\nUsage: $0 <in.info> <out.prefix> [training.proportion=0.75]\n";
	print STDERR "\nThis program is designed to \n\n";
	exit 2;
}

my $training;
my $testing;

my $training_proportion = $ARGV[2] || 0.75;


open IN, "$ARGV[0]" or die( "$!" );
my $header = <IN>;
$training  = $header;
$testing   = $header;
while( <IN> )
{
	if( rand() < $training_proportion )
	{
		$training .= $_;
	}
	else
	{
		$testing  .= $_;
	}
}
close IN;

open OUT, ">$ARGV[1].training" or die( "$!" );
print OUT $training;
close OUT;

open OUT, ">$ARGV[1].testing" or die( "$!" );
print OUT $testing;
close OUT;


