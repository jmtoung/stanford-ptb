#!/usr/bin/perl
#$ -cwd ### use current directory
#$ -M jmtoung@stanford.edu ### email address
##$ -m ea ### mail is to be sent at abort and end time
#$ -j y ### combine stdout and stderr
#$ -t 1-1 ### array job #
#$ -V ### use current environment variables
##$ -pe DJ 12 ### parallel threads
#$ -l h_vmem=14G ### request memory
#$ -l h_rt=12:00:00

use lib '/home/jmtoung/dev';
use GetLineCount;
use strict;
use File::Basename;

my $ID = $ENV{SGE_TASK_ID} - 1;

my $cov = "/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Data/allG1G2-1000G-HRS/step4-AddPhenotype/step4-OmniDataFromSuyash-Phase3/1000G-as-controls/exclude-CG-AT-snps/exclude-141/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_filtered_ld.cov";


############
open(COV,$cov) or die "[STDERR]: can't open $cov: $!\n";

my %files;
while(<COV>) {
	chomp;
	my @split = split('\t');
	next if $split[0] eq 'FID';

	my $dataset = $split[14];

	next if $dataset eq '1000g';

	if ($dataset eq 'allg1g2') {
		my $ptbSubtype = $split[28]; 
		
		if ($ptbSubtype != 1) {
			print STDERR "filtered out allg1g2 individual b/c not spontaneous: $split[1]\n";
			next;
		}
	}
	

	# superCode is predictedPopulation
	my $superCode = $split[26];

	my $fh;
	unless (exists $files{$superCode}) {
		open($files{$superCode},">"."IndividualLists-SPONTANEOUS-".$superCode . ".txt") or die "[STDERR]: can't open $superCode file\n";
	}

	$fh = $files{$superCode};

	print {$fh} join("\t",@split[0..1]), "\n";

}

print STDERR "completed\t$ID\n";
