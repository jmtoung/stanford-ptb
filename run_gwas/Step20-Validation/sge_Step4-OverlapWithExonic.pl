#!/usr/bin/perl
#$ -cwd ### use current directory
#$ -M jmtoung@stanford.edu ### email address
##$ -m ea ### mail is to be sent at abort and end time
#$ -j y ### combine stdout and stderr
#$ -t 1-6 ### array job #
#$ -V ### use current environment variables
##$ -pe DJ 12 ### parallel threads
#$ -l h_vmem=12G ### request memory
#$ -l h_rt=100:00:00

use strict;
use File::Basename;
use lib '/home/jmtoung/dev';
use GetLineCount;

my $ID = $ENV{SGE_TASK_ID} - 1;

my $map = "/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Analysis/RunPlink/step4-OmniDataFromSuyash-Phase3/allG1G2-1000G-HRS-ccGWAS-covar/pc10ev-sex/2014-11-22_SPONTANEOUS-PTB-WITHOUT-1000G/Step20-Validation/GWASRep.map";

my @files = qw(
/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Analysis/RunPlink/step4-OmniDataFromSuyash-Phase3/allG1G2-1000G-HRS-ccGWAS-covar/pc10ev-sex/2014-11-22_SPONTANEOUS-PTB-WITHOUT-1000G/Step20-Validation/snpsToValidate-AFR.txt
/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Analysis/RunPlink/step4-OmniDataFromSuyash-Phase3/allG1G2-1000G-HRS-ccGWAS-covar/pc10ev-sex/2014-11-22_SPONTANEOUS-PTB-WITHOUT-1000G/Step20-Validation/snpsToValidate-EUR.txt
/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Analysis/RunPlink/step4-OmniDataFromSuyash-Phase3/allG1G2-1000G-HRS-ccGWAS-covar/pc10ev-sex/2014-11-22_SPONTANEOUS-PTB-WITHOUT-1000G/Step20-Validation/snpsToValidate-AFR_proxySnps.txt
/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Analysis/RunPlink/step4-OmniDataFromSuyash-Phase3/allG1G2-1000G-HRS-ccGWAS-covar/pc10ev-sex/2014-11-22_SPONTANEOUS-PTB-WITHOUT-1000G/Step20-Validation/snpsToValidate-EUR_proxySnps.txt
/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Analysis/RunPlink/step4-OmniDataFromSuyash-Phase3/allG1G2-1000G-HRS-ccGWAS-covar/pc10ev-sex/2014-11-22_SPONTANEOUS-PTB-WITHOUT-1000G/Step20-Validation/snpsToValidate-AFR_sameGene.txt
/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Analysis/RunPlink/step4-OmniDataFromSuyash-Phase3/allG1G2-1000G-HRS-ccGWAS-covar/pc10ev-sex/2014-11-22_SPONTANEOUS-PTB-WITHOUT-1000G/Step20-Validation/snpsToValidate-EUR_sameGene.txt
);

my ($fileName,$fileDir,$fileExt) = fileparse($files[$ID],'\.txt');

my $out = $fileName . "_overlapExome" . $fileExt;

open(OUT,">".$out) or die "[STDERR]: can't open $out: $!\n";

open(MAP,$map ) or die "[STDERR]: can't open $map: $!\n";

my %map;
while(<MAP>) {
	chomp;
	my @split = split('\t');
	my $key = join("_",$split[0],$split[3]);
	$map{$key}++;
}

open(FILE,$files[$ID]) or die "[STDERR]: can't open $files[$ID]: $!\n";

while(<FILE>) {
	chomp;
	my @split = split('\t');
	my $key = join("_",$split[1],$split[2]);

	if (exists $map{$key}) {
		push(@split,"overlap-exome");
	} else {
		push(@split,"no-overlap");
	}

	print OUT join("\t",@split), "\n";
}


print STDERR "completed\t$ID\n";
