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

my $cov = "/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Data/allG1G2-1000G-HRS/step3-PCA/step4-OmniDataFromSuyash-Phase3/1000G-as-controls/exclude-CG-AT-snps/exclude-141/20160602_NADAV/step13-SNPRelate/all_merged_mode_5_filtered_ld.cov";

my ($covName,$covDir,$covExt) = fileparse($cov,'\.cov');

my $out = $covName . ".cov";

open(OUT,">".$out) or die "[STDERR]: can't open $out: $!\n";

open(COV,$cov) or die "[STDERR]: can't open $cov: $!\n";

my %stats;
while(<COV>) {
	chomp;
	my @split = split('\t');

	if ($split[0] eq 'FID') {
		push(@split,'pheno_ptb');
		print OUT join("\t",@split), "\n";
		next;
	}

	my $dataset = $split[15];

	my $ptb;

	if ($dataset eq 'allg1g2') {
		$ptb = 2;
	} else {
		$ptb = 1;
	}

	$stats{$ptb}++;

	push(@split,$ptb);
	print OUT join("\t",@split), "\n";
}

foreach my $stat (keys %stats) {
	print STDERR $stat, "\t", $stats{$stat}, "\n";
}

print STDERR "completed\t$ID\n";
