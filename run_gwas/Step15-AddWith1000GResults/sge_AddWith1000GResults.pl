#!/usr/bin/perl
#$ -cwd ### use current directory
#$ -M jmtoung@stanford.edu ### email address
##$ -m ea ### mail is to be sent at abort and end time
#$ -j y ### combine stdout and stderr
#$ -t 1-1 #5 ### array job #
#$ -V ### use current environment variables
##$ -pe DJ 12 ### parallel threads
#$ -l h_vmem=16G ### request memory
#$ -l h_rt=100:00:00

use strict;
use File::Basename;

my $ID = $ENV{SGE_TASK_ID} - 1;

my $file = "/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Analysis/RunPlink/step4-OmniDataFromSuyash-Phase3/allG1G2-1000G-HRS-ccGWAS-covar/pc10ev-sex/2014-11-22_SPONTANEOUS-PTB-WITHOUT-1000G/Step11-ReplaceRSIds/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_ALL_ADD.refBase.annotation.exonicVariant.adjusted.counts.otherPValue.confirmNmiss.rsID.txt";

my $otherFile = "/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Analysis/RunPlink/step4-OmniDataFromSuyash-Phase3/allG1G2-1000G-HRS-ccGWAS-covar/pc10ev-sex/2014-11-05_SPONTANEOUS-PTB/Step11-ReplaceRSIds/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_ALL_ADD.refBase.annotation.exonicVariant.adjusted.counts.otherPValue.confirmNmiss.rsID.txt";

my $otherTag = "with1000G-GWAS";
######################

open(OTHER,$otherFile) or die "[STDERR]: can't open $otherFile: $!\n";
my %other;
while(<OTHER>) {
	chomp;
	my @split = split('\t');
	next if $split[0] eq 'CHR';

	## CHROM, BP, POP
	my $key = join("_",$split[0], $split[2], $split[9]);
	## P VALUE, FDR
	my @data = ($split[8], $split[20]);
	
	if (exists $other{$key}) {
		die "[STDERR]: already exists for $key\n";
	} else {
		$other{$key} = \@data;
	}
}

my ($fileName,$fileDir,$fileExt) = fileparse($file,'\.txt');
my $out = $fileName . ".otherResults" . $fileExt;

open(OUT,">".$out) or die "[STDERR]: Can't open $out: $!\n";

open(FILE,$file) or die "[STDERR]: can't open $file: $!\n";
my $notFound = 0;
while(<FILE>) {
	chomp;
	my @split = split('\t');

	my $snp = $split[1];
	
	unless ($split[0] eq 'CHR') {
		my $key = join("_",$split[0], $split[2], $split[9]);
		if (exists $other{$key}) {

			push(@split,$other{$key}->[0], $other{$key}->[1]);

			my $pvalue = $split[8];

			if ($pvalue ne 'NA' && $other{$key}->[0] ne 'NA') {
				my $logP = log($pvalue/$other{$key}->[0]);
				push(@split,$logP);
			} else {
				push(@split,'NA');
			}
		} else {
			push(@split,'NA','NA','NA');
			$notFound++;
		}
	} else {
		my @header = ('p-value', 'fdr', '-log(P/Pold)');
		foreach my $h (@header) {
			push(@split,$otherTag . "_" . $h);
		}
	}

	print OUT join("\t",@split), "\n";
}

print STDERR "notFound:\t$notFound\n";

print STDERR "completed\t$ID\n";


sub log10 {
	my $n = shift;
	return log($n)/log(10);
}
