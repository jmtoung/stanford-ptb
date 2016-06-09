#!/usr/bin/perl
#$ -cwd ### use current directory
#$ -M jmtoung@stanford.edu ### email address
##$ -m ea ### mail is to be sent at abort and end time
#$ -j y ### combine stdout and stderr
#$ -t 1-1 ### array job #
#$ -V ### use current environment variables
##$ -pe DJ 12 ### parallel threads
#$ -l h_vmem=8G ### request memory
#$ -l h_rt=100:00:00

use lib '/home/jmtoung/dev';
use GetLineCount;
use strict;
use File::Basename;

my $ID = $ENV{SGE_TASK_ID} - 1;

my $file = "/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Analysis/RunPlink/step4-OmniDataFromSuyash-Phase3/allG1G2-1000G-HRS-ccGWAS-covar/pc10ev-sex/2014-11-22_SPONTANEOUS-PTB-WITHOUT-1000G/Step6-AnnotateAssocLogisticFile/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_ALL_ADD.refBase.annotation.txt";

my $exonic = "/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Analysis/RunPlink/step4-OmniDataFromSuyash-Phase3/allG1G2-1000G-HRS-ccGWAS-covar/pc10ev-sex/2014-11-22_SPONTANEOUS-PTB-WITHOUT-1000G/Step6-AnnotateAssocLogisticFile/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_ALL_ADD.refBase.annovarInput.exonic_variant_function";

my ($fileName,$fileDir,$fileExt) = fileparse($file,'\.txt');

my $out = $fileName . ".exonicVariant.txt";

open(OUT,">".$out) or die "[STDERR]: can't open $out: $!\n";

open(EXONIC,$exonic) or die "[STDERR]: Can't open $exonic: $!\n";
my %exonic;
while(<EXONIC>) {
	chomp;
	my @split = split('\t');
	
	my $lineNo = $split[0];
	$exonic{$lineNo} = \@split;
}

open(FILE,$file) or die "[STDERR]: can't open $file: $!\n";

# start at lineNo 2 because first was header
my $lineNo = 2;
while(<FILE>) {
	chomp;
	my @split = split('\t');

	my $key = "line" . $lineNo;

	if (exists $exonic{$key}) {
		my $exonic = $exonic{$key};
	
		# type of variant (non synon/synon)
		push(@split,$exonic->[1]);
	
		push(@split,$exonic->[2]);

	} else {
		push(@split,'NA','NA');
	}

	print OUT join("\t",@split), "\n";

	$lineNo++;
}


print STDERR "completed\t$ID\n";
