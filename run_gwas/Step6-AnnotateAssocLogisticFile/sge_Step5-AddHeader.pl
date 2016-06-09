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

my $originalFile = "/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Analysis/RunPlink/step4-OmniDataFromSuyash-Phase3/allG1G2-1000G-HRS-ccGWAS-covar/pc10ev-sex/2014-11-22_SPONTANEOUS-PTB-WITHOUT-1000G/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_ALL_ADD.assoc.logistic";

my $file = "/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Analysis/RunPlink/step4-OmniDataFromSuyash-Phase3/allG1G2-1000G-HRS-ccGWAS-covar/pc10ev-sex/2014-11-22_SPONTANEOUS-PTB-WITHOUT-1000G/Step6-AnnotateAssocLogisticFile/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_ALL_ADD.refBase.annotation.exonicVariant.txt";

my $temp = $file . ".temp";
open(ORIGINAL, $originalFile) or die "[STDERR]: can't open $originalFile: $!\n";

my @header;
while(<ORIGINAL>) {
	chomp;
	my @split = split('\t');
	@header=@split;
	$split[0] eq 'CHR' or die "[STDERR]: not what expected\n";
	last;
}

open(TEMP,">".$temp) or die "[STDERR]: can't open $temp: $!\n";

push(@header,'gencodeFeature','gencodeGene','exonic_variant_type','exonic_variant_transcript');

print TEMP join("\t",@header), "\n";

open(FILE,$file) or die "[STDERR]: can't open $file: $!\n";

while(<FILE>) {
	chomp;
	my @split = split('\t');

	scalar(@split) == scalar(@header) or die "[STDERR]: wrong # columns\n";

	print TEMP join("\t",@split), "\n";
}
close(FILE);

runCommand("mv $temp $file");

print STDERR "completed\t$ID\n";
