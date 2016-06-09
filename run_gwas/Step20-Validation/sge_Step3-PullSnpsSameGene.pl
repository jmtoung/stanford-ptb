#!/usr/bin/perl
#$ -cwd ### use current directory
#$ -M jmtoung@stanford.edu ### email address
##$ -m ea ### mail is to be sent at abort and end time
#$ -j y ### combine stdout and stderr
#$ -t 1-1 ### array job #
#$ -V ### use current environment variables
##$ -pe DJ 12 ### parallel threads
#$ -l h_vmem=12G ### request memory
#$ -l h_rt=100:00:00

use strict;
use File::Basename;
use lib '/home/jmtoung/dev';
use GetLineCount;

my $ID = $ENV{SGE_TASK_ID} - 1;

my $sig = "/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Analysis/RunPlink/step4-OmniDataFromSuyash-Phase3/allG1G2-1000G-HRS-ccGWAS-covar/pc10ev-sex/2014-11-22_SPONTANEOUS-PTB-WITHOUT-1000G/Step15-AddWith1000GResults/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_ALL_ADD.refBase.annotation.exonicVariant.adjusted.counts.otherPValue.confirmNmiss.rsID.otherResults.fdrLePt05.txt";

my $full = "/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Analysis/RunPlink/step4-OmniDataFromSuyash-Phase3/allG1G2-1000G-HRS-ccGWAS-covar/pc10ev-sex/2014-11-22_SPONTANEOUS-PTB-WITHOUT-1000G/Step15-AddWith1000GResults/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_ALL_ADD.refBase.annotation.exonicVariant.adjusted.counts.otherPValue.confirmNmiss.rsID.otherResults.txt";



# 
open(SIG,$sig) or die "[STDERR]:can't open $sig: $!\n";
my %genes;
my %snps;
while(<SIG>) {
	chomp;
	my @split = split('\t');

	my $gene = $split[11];
	my $kpg = $split[32];
	my $pop = $split[9];

	next if $gene eq 'gencodeGene';

	next if $gene =~ /dist/;

	$genes{$pop}{$gene}++;
	$snps{$pop}{$kpg}++;
}

my %files;
open(FULL,$full) or die "[STDERR]: can't open $full: $!\n";
while(<FULL>) {
	chomp;
	my @split = split('\t');

	my $gene = $split[11];
	my $kpg = $split[32];
	my $pop = $split[9];

	next if (exists $snps{$pop}{$kpg});

	next unless (exists $genes{$pop}{$gene});

	unless (exists $files{$pop}) {
		my $file = "snpsToValidate-" . $pop . "_sameGene.txt";
		open($files{$pop},">".$file) or die "[STDERR]: can't open $file: $!\n";
	}

	print {$files{$pop}} join("\t",$kpg, $split[0], $split[2], @split), "\n";
}


print STDERR "completed\t$ID\n";
