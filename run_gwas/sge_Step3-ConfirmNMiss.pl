#!/usr/bin/perl
#$ -cwd ### use current directory
#$ -M jmtoung@stanford.edu ### email address
##$ -m ea ### mail is to be sent at abort and end time
#$ -j y ### combine stdout and stderr
#$ -t 1-5 ### array job #
#$ -V ### use current environment variables
##$ -pe DJ 12 ### parallel threads
##$ -l h_vmem=12G ### request memory
##$ -l h_rt=100:00:00

use strict;
use File::Basename;

my $ID = $ENV{SGE_TASK_ID} - 1;

my $genotypeCounts = "/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Data/allG1G2-1000G-HRS/step2-MergeFiles/step4-OmniDataFromSuyash-Phase3/step4-ExcludeRelatedsMissings/exclude-CG-AT-snps/exclude-141/step9-ExtractGenotypes-ByDataset-PhenoPTBSubtype-BY-COLUMN/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141.genotypecounts.nmiss-SPONTANEOUS-ONLY";

my @files = qw(
/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Analysis/RunPlink/step4-OmniDataFromSuyash-Phase3/allG1G2-1000G-HRS-ccGWAS-covar/pc10ev-sex/2014-11-05_SPONTANEOUS-PTB/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_AFR.assoc.logistic
/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Analysis/RunPlink/step4-OmniDataFromSuyash-Phase3/allG1G2-1000G-HRS-ccGWAS-covar/pc10ev-sex/2014-11-05_SPONTANEOUS-PTB/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_AMR.assoc.logistic
/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Analysis/RunPlink/step4-OmniDataFromSuyash-Phase3/allG1G2-1000G-HRS-ccGWAS-covar/pc10ev-sex/2014-11-05_SPONTANEOUS-PTB/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_EAS.assoc.logistic
/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Analysis/RunPlink/step4-OmniDataFromSuyash-Phase3/allG1G2-1000G-HRS-ccGWAS-covar/pc10ev-sex/2014-11-05_SPONTANEOUS-PTB/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_EUR.assoc.logistic
/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Analysis/RunPlink/step4-OmniDataFromSuyash-Phase3/allG1G2-1000G-HRS-ccGWAS-covar/pc10ev-sex/2014-11-05_SPONTANEOUS-PTB/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_SAS.assoc.logistic
);

my ($fileName,$fileDir,$fileExt) = fileparse($files[$ID],'\.assoc\.logistic');

my @fileName = split('_',$fileName);

my $pop = $fileName[-1];

print STDERR "pop:\t$pop\n";

open(COUNTS,$genotypeCounts) or die "[STDERR]: Can't open $genotypeCounts: $!\n";

my %counts;
while(<COUNTS>) {
	chomp;
	my @split = split('\t');
	my $pop_this = $split[1];
	my $nmiss = $split[2];
	my $rsId = $split[0];

	next unless $pop_this eq $pop;

	$counts{$rsId} = $nmiss;
}

open(FILE,$files[$ID]) or die "[STDERR]: can't open $files[$ID]: $!\n";

my $confirmed = 0;
my $na_confirmed = 0;
while(<FILE>) {
	chomp;
	$_ =~ s/^\s+//;
	my @split = split('\s+');

	next if $split[0] eq 'CHR';

	my $rsId = $split[1];
	my $nmiss = $split[5];
	my $pvalue = $split[-1];

	if (exists $counts{$rsId}) {
		if ($nmiss != $counts{$rsId}) {
			if ($pvalue eq 'NA') {
				$na_confirmed++;
			} else {
				print STDERR "wrong_counts\t$rsId\n";
			}
		} else {
			$confirmed++;
		}
	} else {
		die "[STDERR]: missing $rsId in counts\n";
	}
}

print STDERR "number_confirmed:\t$confirmed\n";
print STDERR "number_na_confirmed:\t$na_confirmed\n";
##########

print STDERR "completed\t$ID\n";
