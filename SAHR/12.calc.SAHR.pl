#!/usr/bin/perl

#
# Author: Ahfyth
#

use strict;
use warnings;

if( $#ARGV < 3 )
{
	print STDERR "\nUsage: $0 <topgene.anno> <matrix > <clinical.info> <out.file=LUSC.SAHR>\n\n";
	exit 2;
}

my (%surv, %HR);
open IN, "$ARGV[0]" or die( "$!" );
while( <IN> )
{
	chomp;
	my @l = split /\t/;	##AL590428.1      0.00176827868963492     UpSaver 3.99675486871364
	$surv{$l[0]} = $l[2];
	$HR{$l[0]}   = $l[3];
}
close IN;

my %other;
open IN, "$ARGV[2]" or die( "$!" );
<IN>;
while( <IN> )
{
	chomp;
	my @l = split /\t/; ##submitter_id project_id tumor_stage age_at_diagnosis
#case_submitter_id       vital_status    days_for_survival       age_at_diagnosis        ajcc_pathologic_t
#TCGA-CV-7255    Dead    64      11931   T4a
	$other{$l[0]} = "$l[4]\t$l[3]";
}
close IN;

open OUT, ">$ARGV[3]" or die( "$!" );

open IN, "$ARGV[1]" or die( "$!" );
my $header = <IN>;
chomp( $header );
my @gene = split /\t/, $header;
my %cnt;
print OUT "Sample\tVital\tDays\tSAHR\tStage\tAge\n";
while( <IN> )
{
	chomp;
	my @l = split /\t/;	##Sample	Vital	Days	Gene1	Gene2 ...
	next unless $l[0] =~ /^TCGA/;
	my $SAHR = 0;	## score of aggregated hazard ratio
	foreach my $i (3..$#l)
	{
		my $g = $gene[$i];
		my $type = $surv{$g};

		my $direction;
		if( $l[$i] eq 'H' )
		{
			if( $type eq 'Oncogene' || $type eq 'DownSaver' )
			{
				$direction = 1;
			}
			else
			{
				$direction = -1;
			}
		}
		elsif( $l[$i] eq 'L' )
		{
			if( $type eq 'TumorSuppressor' || $type eq 'UpSaver' )
			{
				$direction = 1;
			}
			else
			{
				$direction = -1;
			}
		}
		else	## no data for this gene, skip
		{
			print STDERR "No data for $l[0]: $g!\n";
			next;
		}
		
		$SAHR += $direction * $HR{$g};
	}
#	my $code;
#	if(    $SAHR < -50 ){ $code = 'Q1' }
#	elsif( $SAHR <   0 ){ $code = 'Q2' }
#	elsif( $SAHR <  50 ){ $code = 'Q3' }
#	else                { $code = 'Q4' }
#	$cnt{$code} ++;
	print OUT join("\t", $l[0], $l[1], $l[2], $SAHR, $other{$l[0]}||"NA\tNA"), "\n";
}
close IN;
close OUT;

#foreach my $k ( sort keys %cnt )
#{
#	print STDERR "$k\t$cnt{$k}\n";
#}

