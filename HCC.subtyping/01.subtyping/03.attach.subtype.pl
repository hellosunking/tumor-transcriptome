#!/usr/bin/perl


use strict;
use warnings;


if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <selected.sample.info> <cluster.info>\n\n";
#	print STDERR "\nThis program is designed to \n\n";
	exit 2;
}

my %cluster;

open IN, "$ARGV[1]" or die( "$!" );
while( <IN> ) {
	chomp;
	my @l = split /\t/;	##TCGA.2V.A95S.01A	0
	#next unless $l[0]=~s/\.01A.gz$//;
	$l[0] =~ s/\./-/g;
	$cluster{$l[0]} = "C$l[1]";
}
close IN;

open IN, "$ARGV[0]" or die( "$!" );
#my $h = <IN>;
#chomp($h);
print "SID\tVital\tDays\tAge\tStage\tGender\tRace\tSubtype\n";
while( <IN> ) {
	chomp;
	my @l = split /\t/; ##case_submitter_id vital_status days_to_survival age_at_diagnosis ajcc_pathologic_stage gender race
	next unless $l[0]=~/^TCGA/;
	next if $l[1] eq 'NA' || $l[2] eq 'NA';
	next unless exists $cluster{$l[0]};
	if( $l[1] eq "Alive" ) {
		$l[1] = 0;
	} elsif ( $l[1] eq "Dead" ) {
		$l[1] = 1;
	}
	print join("\t", $l[0],$l[1],$l[2],$l[3],$l[4],$l[5],$l[6], $cluster{$l[0]}), "\n";
}
close IN;

