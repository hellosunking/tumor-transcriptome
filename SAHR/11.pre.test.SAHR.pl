#!/usr/bin/perl
##
## Date: 20221103
#
use strict;
use warnings;

if( $#ARGV < 5 )
{
        print STDERR "\nUsage: $0 <topgene.annno> <cutoff.file> <clinical.info> <output.file> <cancer type > <directory> \n";
        print STDERR "\nThis program is designed to \n\n";
        exit 2;
}

my $cancer = "$ARGV[4]";
my $data_dir  = "$ARGV[5]";

my (%surv, %HR);
open IN, "$ARGV[0]" or die( "$!" );
while( <IN> )
{
        chomp;
        my @l = split /\t/;     ##AL590428.1      0.00176827868963492     UpSaver 3.99675486871364
        $surv{$l[0]} = $l[2];
#       $HR{$l[0]}   = $l[2];
}
close IN;
#
my @genes = sort keys %surv;

## load cutoff
my %cutoff;
open IN, "$ARGV[1]" or die( "$!" );
while( <IN> )
{
      chomp;
      my @l = split /\t/;     ##DKK1 337.6176
      $cutoff{$l[0]}=$l[1] if exists $surv{$l[0]};
}
close IN;
#
open OUT, ">$ARGV[3]" or die( "$!" );
print OUT join("\t", "#Sample", "Vital", "Days", @genes), "\n";

## Loading sample code and corresponding RNAseq data
my %clinical;
#my @info;
open IN, "$ARGV[2]" or die( "$!" );
my $header=<IN>;
while( <IN> )
{
       chomp;	
       my  @info = split /\t/;     
#case_submitter_id       vital_status    days_for_survival       age_at_diagnosis        ajcc_pathologic_t
#TCGA-CR-7372    Alive   759     16610   T1
	next if $info[2] eq "NA";
#	next if $info[2] <= 0;
        if( $info[1] eq "Dead" ) {
                $clinical{$info[0]} = "1\t$info[2]";
        } elsif( $info[1] eq "Alive" ) {
        #        next unless $info[2] > 0;        ## discard records with 0 following-up days
                $clinical{$info[0]} = "0\t$info[2]";
        }
}
close IN;

open IN, "$ARGV[2]" or die( "$!" );
while( <IN> )
{
        chomp;
        my @info = split /\t/;
	next unless (exists $clinical{$info[0]});
	if ( -s "$data_dir/$cancer/Count/$info[0]-01A.gz" ) 
	{  
		my %exp;
		open COUNT, "zcat $data_dir/$cancer/Count/$info[0]-01A.gz |" or die ("$!:");
		while( <COUNT> ) 
		{
			chomp;
			my @l = split /\t/;     ##gene_id gene_name       gene_type       unstranded      stranded_first  stranded_second tpm_unstranded  fpkm_unstranded fpkm_uq_unstranded
			next unless $l[0]=~/^ENSG/;
			if( exists $cutoff{$l[0]} )
			{
				if( $l[7] > $cutoff{$l[0]} )
				{
					$exp{$l[0]} = 'H';
				}else
				{
					$exp{$l[0]} = 'L';
				}
			}	
		}		
		close COUNT;
		print OUT join("\t", $info[0], $clinical{$info[0]}, map {$exp{$_}||'X'} @genes ), "\n";
	}
}
	
close IN;

close OUT;

