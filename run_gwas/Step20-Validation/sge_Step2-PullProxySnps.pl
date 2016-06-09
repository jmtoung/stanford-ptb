#!/usr/bin/perl
#$ -cwd ### use current directory
#$ -M jmtoung@stanford.edu ### email address
##$ -m ea ### mail is to be sent at abort and end time
#$ -j y ### combine stdout and stderr
#$ -t 1-2 ### array job #
#$ -V ### use current environment variables
##$ -pe DJ 12 ### parallel threads
#$ -l h_vmem=12G ### request memory
#$ -l h_rt=100:00:00

use strict;
use File::Basename;
use lib '/home/jmtoung/dev';
use GetLineCount;

my $ID = $ENV{SGE_TASK_ID} - 1;

my $map = "/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Data/allG1G2-1000G-HRS/step2-MergeFiles/step4-OmniDataFromSuyash-Phase3/step4-ExcludeRelatedsMissings/exclude-CG-AT-snps/exclude-141/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141.map";

my @snps = qw(
/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Analysis/RunPlink/step4-OmniDataFromSuyash-Phase3/allG1G2-1000G-HRS-ccGWAS-covar/pc10ev-sex/2014-11-22_SPONTANEOUS-PTB-WITHOUT-1000G/Step20-Validation/snpsToValidate-AFR.txt
/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Analysis/RunPlink/step4-OmniDataFromSuyash-Phase3/allG1G2-1000G-HRS-ccGWAS-covar/pc10ev-sex/2014-11-22_SPONTANEOUS-PTB-WITHOUT-1000G/Step20-Validation/snpsToValidate-EUR.txt
);

my $threshold = 0.8;

my @r2 = qw(
/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Data/allG1G2-1000G-HRS/step2-MergeFiles/step4-OmniDataFromSuyash-Phase3/step4-ExcludeRelatedsMissings/exclude-CG-AT-snps/exclude-141/step20-CalculateR2/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_AFR.ld
/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Data/allG1G2-1000G-HRS/step2-MergeFiles/step4-OmniDataFromSuyash-Phase3/step4-ExcludeRelatedsMissings/exclude-CG-AT-snps/exclude-141/step20-CalculateR2/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_EUR.ld
);


my ($snpName,$snpDir,$snpExt) = fileparse($snps[$ID],'\.txt');

my $out = $snpDir . $snpName . "_proxySnps.txt";

open(OUT,">".$out) or die "[STDERR]: can't open $out: $!\n";

open(MAP,$map) or die "[STDERR]: can't open $map: $!\n";
my %map;
while(<MAP>) {
	chomp;
	my @split = split('\t');
	my $rs = $split[1];
	if (exists $map{$rs}) {
		die "[STDERR]: $rs already exists\n";
	}
	$map{$rs} = \@split;
}

open(SNPS,$snps[$ID]) or die "[STDERR]: can't open $snps[$ID]: $!\n";
my %snps;
while(<SNPS>) {
	chomp;
	my @split = split('\t');
#	if (exists $snps{$split[3]}) {
#		die "[STDERR]: already saw $split[3]\n";
#	}

	$snps{$split[3]} = \@split;
}

open(R2,$r2[$ID]) or die "[STDERR]: can't open $r2[$ID]: $!\n";

while(<R2>) {
	chomp;
	$_ =~ s/^\s+//;
	my @split = split('\s+');
	
	my $snp1 = $split[2];
	my $snp2 = $split[5];
	my $r2 = $split[6];

	next unless $r2 >= $threshold;
	my $m;

	if (exists $snps{$snp1}) {
		if (exists $map{$snp2}) {
			$m = $map{$snp2};
		} else {
			die "[STDERR]: can't find map for $snp1\n";
		}
		print OUT join("\t",$snp2, $m->[0], $m->[3], $snp1, @{$snps{$snp1}}), "\n";		
	}

	if (exists $snps{$snp2}) {
		if (exists $map{$snp1}) {
			$m = $map{$snp1};
		} else {
			die "[STDERR]: can't find map for $snp2\n";
		}
		print OUT join("\t",$snp1, $m->[0], $m->[3], $snp2, @{$snps{$snp2}}), "\n";
	}

}


print STDERR "completed\t$ID\n";
