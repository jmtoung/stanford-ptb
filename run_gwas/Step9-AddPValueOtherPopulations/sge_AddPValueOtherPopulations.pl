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

my $file = "/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Analysis/RunPlink/step4-OmniDataFromSuyash-Phase3/allG1G2-1000G-HRS-ccGWAS-covar/pc10ev-sex/2014-11-22_SPONTANEOUS-PTB-WITHOUT-1000G/Step8-AddAlleleGenotypeCounts/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_ALL_ADD.refBase.annotation.exonicVariant.adjusted.counts.txt";

my ($fileName,$fileDir,$fileExt) = fileparse($file,'\.txt');

my $out = $fileName . ".otherPValue" . $fileExt;

open(OUT,">".$out) or die "[STDERR]: can't open $out: $!\n";

open(FILE,$file) or die "[STDERR]: can't open $file: $!\n";
my %data;
my %pops;

while(<FILE>) {
	chomp;
	my @split = split('\t');

	next if $split[0] eq 'CHR';

	my $chrom = $split[0];
	my $rs = $split[1];
	my $bp = $split[2];

	my $pop = $split[9];
	my $pvalue = $split[8];

	my $key = join("_",$chrom,$rs,$bp);

	$pops{$pop}++;

	if (exists $data{$pop}{$key}) {
		die "[STDERR]: Conflict\n";
	} else {
		$data{$pop}{$key} = $pvalue;
	}

}

open(FILE,$file) or die "[STDERR]: can't open $file: $!\n";

my $lineNo = 1;
my @pops = sort keys %pops;

while(<FILE>) {
	chomp;
	my @split = split('\t');

	my $chrom = $split[0];
	my $rs = $split[1];
	my $bp = $split[2];

	my $pop = $split[9];
	my $pvalue = $split[8];

	my $key = join("_",$chrom,$rs,$bp);
	
	if ($lineNo == 1) {
		foreach my $pop (@pops) {
			push(@split,"p_value_" . $pop);
		}
	} else {
		foreach my $pop (@pops) {
			if (exists $data{$pop}{$key}) {
				push(@split,$data{$pop}{$key});
			} else {
				push(@split,'NA');
			}
		}
	}

	print OUT join("\t",@split), "\n";

	$lineNo++;

}


print STDERR "completed\t$ID\n";
