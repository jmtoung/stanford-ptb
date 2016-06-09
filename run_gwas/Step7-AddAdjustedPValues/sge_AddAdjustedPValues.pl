#!/usr/bin/perl
#$ -cwd ### use current directory
#$ -M jmtoung@stanford.edu ### email address
##$ -m ea ### mail is to be sent at abort and end time
#$ -j y ### combine stdout and stderr
#$ -t 1-1 ### array job #
#$ -V ### use current environment variables
##$ -pe DJ 12 ### parallel threads
#$ -l h_vmem=20G ### request memory
#$ -l h_rt=100:00:00

use lib '/home/jmtoung/dev';
use GetLineCount;
use strict;
use File::Basename;

my $ID = $ENV{SGE_TASK_ID} - 1;

my $file = "/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Analysis/RunPlink/step4-OmniDataFromSuyash-Phase3/allG1G2-1000G-HRS-ccGWAS-covar/pc10ev-sex/2014-11-22_SPONTANEOUS-PTB-WITHOUT-1000G/Step6-AnnotateAssocLogisticFile/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_ALL_ADD.refBase.annotation.exonicVariant.txt";

my ($fileName,$fileDir,$fileExt) = fileparse($file,'\.txt');

my $out = $fileName . ".adjusted" . $fileExt;

my @adjusted = qw(
/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Analysis/RunPlink/step4-OmniDataFromSuyash-Phase3/allG1G2-1000G-HRS-ccGWAS-covar/pc10ev-sex/2014-11-22_SPONTANEOUS-PTB-WITHOUT-1000G/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_AFR.assoc.logistic.adjusted
/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Analysis/RunPlink/step4-OmniDataFromSuyash-Phase3/allG1G2-1000G-HRS-ccGWAS-covar/pc10ev-sex/2014-11-22_SPONTANEOUS-PTB-WITHOUT-1000G/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_AMR.assoc.logistic.adjusted
/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Analysis/RunPlink/step4-OmniDataFromSuyash-Phase3/allG1G2-1000G-HRS-ccGWAS-covar/pc10ev-sex/2014-11-22_SPONTANEOUS-PTB-WITHOUT-1000G/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_EAS.assoc.logistic.adjusted
/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Analysis/RunPlink/step4-OmniDataFromSuyash-Phase3/allG1G2-1000G-HRS-ccGWAS-covar/pc10ev-sex/2014-11-22_SPONTANEOUS-PTB-WITHOUT-1000G/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_EUR.assoc.logistic.adjusted
/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Analysis/RunPlink/step4-OmniDataFromSuyash-Phase3/allG1G2-1000G-HRS-ccGWAS-covar/pc10ev-sex/2014-11-22_SPONTANEOUS-PTB-WITHOUT-1000G/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_SAS.assoc.logistic.adjusted
);

# load adjusted values

my @header;
my %data;
foreach $a (@adjusted) {
	open(ADJUSTED,$a) or die "[STDERR]: Can't open $a: $!\n";

	my ($aName,$aDir,$aExt) = fileparse($a,'\.assoc\.logistic\.adjusted');
	my @aName = split('_',$aName);
	my $pop = $aName[-1];
	
	while(<ADJUSTED>) {
		chomp;
		$_ =~ s/^\s//;
		my @split = split('\s+|\t+');
		shift @split if $split[0] eq '';

		if ($split[0] eq 'CHR') {
			unless(@header) {
				@header = @split[2..(scalar(@split) - 1)];
			}
		} else {
			my $key = join("_",@split[0..1]);
			if (exists $data{$pop}{$key}) {
				die "[STDERR]: already exists for $pop in $key\n";
			} else {
				$data{$pop}{$key} = \@split;
			}
		}
	}
}

open(FILE,$file) or die "[STDERR]: can't open $file: $!\n";
open(OUT,">".$out) or die "[STDERR]: can't open $out: $!\n";

while(<FILE>) {
	chomp;
	my @split = split('\t');

	if ($split[0] eq 'CHR') {
		@header = (@split,@header);
		print OUT join("\t",@header), "\n";
	} else {
		my $pop = $split[-5];
		my $key = join("_",@split[0..1]);
		if (exists $data{$pop}{$key}) {
			push(@split,@{$data{$pop}{$key}}[2..(scalar(@{$data{$pop}{$key}}) - 1)]);
			scalar(@split) == scalar(@header) or die "[STDERR]: wah\n";
			print OUT join("\t",@split), "\n";
		} else {
			die "[STDERR]: Can't find $key for $pop\n";
		}
	}

}

print STDERR "completed\t$ID\n";
