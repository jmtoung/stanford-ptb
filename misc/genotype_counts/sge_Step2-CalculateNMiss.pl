#!/usr/bin/perl
#$ -cwd ### use current directory
#$ -M jmtoung@stanford.edu ### email address
##$ -m ea ### mail is to be sent at abort and end time
#$ -j y ### combine stdout and stderr
#$ -t 1-1 ### array job #
#$ -V ### use current environment variables
##$ -pe DJ 12 ### parallel threads
##$ -l h_vmem=50G ### request memory
#$ -l h_rt=100:00:00

use lib '/home/jmtoung/dev';
use GetLineCount;
use strict;
use File::Basename;

my $ID = $ENV{SGE_TASK_ID} - 1;

my $file = "/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Data/allG1G2-1000G-HRS/step2-MergeFiles/step4-OmniDataFromSuyash-Phase3/step4-ExcludeRelatedsMissings/exclude-CG-AT-snps/exclude-141/step9-ExtractGenotypes-ByDataset-PhenoPTBSubtype-BY-COLUMN/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141.genotypecounts";

open(FILE,$file) or die "[STDERR]: can't open $file: $!\n";

my $OUT;

my $out = $file . ".nmiss";

open($OUT,">".$out) or die "[STDERR]: can't open $out: $!\n";

my %counts;
my $currentSnp;

while(<FILE>) {
	chomp;
	my @split = split('\t');
	
	my $rsId = $split[1];
	my $pop = $split[6];

	my $dataset = $split[4];

	my $subtype = $split[5];

	if (!$currentSnp or $currentSnp ne $rsId) {
		printData(\%counts, $currentSnp, $OUT);
		%counts = ();		
		$currentSnp = $rsId;
	}

	my @counts = split(',',$split[8]);

	foreach my $count (@counts) {
		my @c = split(':',$count);

		next if $c[0] eq '00';

		$counts{$pop} += $c[1];	


	}
}

printData(\%counts, $currentSnp, $OUT);


sub printData {
	my ($counts, $currentSnp, $FH) = @_;

	return if !defined $currentSnp;

	foreach my $pop (sort keys %{$counts}) {
		print {$FH} join("\t",$currentSnp, $pop, $counts->{$pop}), "\n";
	}
}


print STDERR "completed\t$ID\n";
