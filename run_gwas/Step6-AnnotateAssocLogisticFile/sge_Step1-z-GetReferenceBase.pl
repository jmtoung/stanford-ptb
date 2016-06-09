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

my $ref = "/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Data/allG1G2-1000G-HRS/step1-PrepFiles/allG1G2/step1-ExtractReferenceBases/allG1G2.refbases";

my $file = "/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Analysis/RunPlink/step4-OmniDataFromSuyash-Phase3/allG1G2-1000G-HRS-ccGWAS-covar/pc10ev-sex/2014-11-22_SPONTANEOUS-PTB-WITHOUT-1000G/Step6-AnnotateAssocLogisticFile/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_ALL_ADD.annovarInput";

my ($fileName,$fileDir,$fileExt) = fileparse($file,'\.annovarInput');

my $out = $fileName . ".refBase" . $fileExt;

open(FILE,$file) or die "[STDERR]: can't open $file: $!\n";
open(OUT,">".$out) or die "[STDERR]: can't open $out: $!\n";

open(REF,$ref) or die "[STDERR]: can't open $ref: $!\n";
my %ref;
while(<REF>) {
	chomp;
	my @split = split('\t');
	my $key = join("_",$split[0],$split[3],$split[1]);

	if (exists $ref{$key}) {
		die "[STDERR]: exists $key\n";
	} else {
		$ref{$key} = uc($split[4]);
	}
}

my $lineNo = 1;
while(<FILE>) {
	chomp;

	my @split = split('\t');

	if ($lineNo == 1){ 
		print OUT join("\t",@split), "\n";
		$lineNo++;
		next;
	}

	$split[4] = $split[8];

	# look up ref
	my $chrom = $split[5];
	$chrom = 'X' if $chrom eq '23';
	$chrom = 'Y' if $chrom eq '24';
	$chrom = 'XY' if $chrom eq '25';
	$chrom = 'MT' if $chrom eq '26';

	my $key = join("_",$chrom,$split[7],$split[6]);

	if (exists $ref{$key}) {
		$split[3] = $ref{$key};
	} else {
		die "[STDERR]: doesn't exist ref for $key\n";
	}

	print OUT join("\t",@split), "\n";
}


print STDERR "completed\t$ID\n";
