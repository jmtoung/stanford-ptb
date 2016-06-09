#!/usr/bin/perl
#$ -cwd ### use current directory
#$ -M jmtoung@stanford.edu ### email address
##$ -m ea ### mail is to be sent at abort and end time
#$ -j y ### combine stdout and stderr
#$ -t 1-1 ### array job #
#$ -V ### use current environment variables
##$ -pe DJ 12 ### parallel threads
#$ -l h_vmem=8G ### request memory
##$ -l h_rt=100:00:00

use lib '/home/jmtoung/dev';
use GetLineCount;
use strict;
use File::Basename;

my $ID = $ENV{SGE_TASK_ID} - 1;

my $file = "/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Analysis/RunPlink/step4-OmniDataFromSuyash-Phase3/allG1G2-1000G-HRS-ccGWAS-covar/pc10ev-sex/2014-11-22_SPONTANEOUS-PTB-WITHOUT-1000G/Step7-AddAdjustedPValues/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_ALL_ADD.refBase.annotation.exonicVariant.adjusted.txt";

my $counts = "/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Data/allG1G2-1000G-HRS/step2-MergeFiles/step4-OmniDataFromSuyash-Phase3/step4-ExcludeRelatedsMissings/exclude-CG-AT-snps/exclude-141/step9-ExtractGenotypes-ByDataset-PhenoPTBSubtype-BY-COLUMN/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141.genotypecounts.case-control-counts_SPONTANEOUS-ONLY-EXCLUDE-1000G.txt";

open(COUNTS,$counts) or die "[STDERR]: can't open $counts: $!\n";

my ($fileName,$fileDir,$fileExt) = fileparse($file,'\.txt');

my $out = $fileName . ".counts" . $fileExt;

open(OUT,">".$out) or die "[STDERR]: can't open $out: $!\n";

my %data;
while(<COUNTS>) {
	chomp;
	my @split = split('\t');
	my $rs = $split[0];
	my $pop = $split[1];
	
	$data{$pop}{$rs}{'case'} = $split[2];
	$data{$pop}{$rs}{'control'} = $split[3];

}

open(FILE,$file) or die "[STDERR]: can't open $file: $!\n";

my $lineNo = 1;
while(<FILE>) {
	chomp;
	my @split = split('\t');

	my $rs = $split[1];
	my $pop = $split[9];
	
	if ($lineNo == 1) {
		push(@split,'case_genotype_counts','control_genotype_counts', 'case_allele_counts', 'control_allele_counts');
	} else {

		my @alleleCounts;	
		if (exists $data{$pop}{$rs}{'case'}) {
			my $genoCounts = $data{$pop}{$rs}{'case'};
			push(@split,$genoCounts);
			my $alleleCounts = calculateAlleleCounts($genoCounts);
			push(@alleleCounts,$alleleCounts);
		} else {
			print STDERR "missingCounts:\t$rs\t$pop\n";
			push(@split,'NA');
			push(@alleleCounts,'NA');
		} 

		if (exists $data{$pop}{$rs}{'control'}) {
			my $genoCounts = $data{$pop}{$rs}{'control'};
			push(@split,$genoCounts);
			my $alleleCounts = calculateAlleleCounts($genoCounts);
			push(@alleleCounts,$alleleCounts);
		} else {
			print STDERR "missingCounts:\t$rs\t$pop\n";
			push(@split,'NA');
			push(@alleleCounts,'NA');
		}

		push(@split,@alleleCounts);

	}

	print OUT join("\t",@split), "\n";

	$lineNo++;

}

sub calculateAlleleCounts {
	my $counts = shift;

	my %alleles;
	foreach my $c (split(',',$counts)) {
		my @c = split(':',$c);

		my @alleles = split('',$c[0]);
		
		foreach my $a (@alleles) {
			$alleles{$a} += $c[1];

		}
	}

	my @counts;
	foreach my $a (sort keys %alleles) {
		push(@counts,join(":",$a, $alleles{$a}));
	}

	return join(",",@counts);

}

print STDERR "completed\t$ID\n";
