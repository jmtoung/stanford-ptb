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

my $ID = $ENV{SGE_TASK_ID} - 1;

my @files = qw(
/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Analysis/RunPlink/step4-OmniDataFromSuyash-Phase3/allG1G2-1000G-HRS-ccGWAS-covar/pc10ev-sex/2014-11-22_SPONTANEOUS-PTB-WITHOUT-1000G/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_AFR.assoc.logistic
/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Analysis/RunPlink/step4-OmniDataFromSuyash-Phase3/allG1G2-1000G-HRS-ccGWAS-covar/pc10ev-sex/2014-11-22_SPONTANEOUS-PTB-WITHOUT-1000G/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_AMR.assoc.logistic
/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Analysis/RunPlink/step4-OmniDataFromSuyash-Phase3/allG1G2-1000G-HRS-ccGWAS-covar/pc10ev-sex/2014-11-22_SPONTANEOUS-PTB-WITHOUT-1000G/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_EAS.assoc.logistic
/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Analysis/RunPlink/step4-OmniDataFromSuyash-Phase3/allG1G2-1000G-HRS-ccGWAS-covar/pc10ev-sex/2014-11-22_SPONTANEOUS-PTB-WITHOUT-1000G/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_EUR.assoc.logistic
/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Analysis/RunPlink/step4-OmniDataFromSuyash-Phase3/allG1G2-1000G-HRS-ccGWAS-covar/pc10ev-sex/2014-11-22_SPONTANEOUS-PTB-WITHOUT-1000G/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_SAS.assoc.logistic
);

my ($fileName,$fileDir,$fileExt) = fileparse($files[0],'\.assoc\.logistic');

my @fileName = split('_',$fileName);

pop @fileName;

my $out = $fileDir . join("_",@fileName,"ALL") . $fileExt;

open(OUT, ">".$out) or die "[STDERR]: can't open $out: $!\n";

my $headerPrinted = 0;
foreach my $file (@files) {
	open(FILE,$file) or die "[STDERR]: can't open $file: $!\n";

	my ($fName,$fDir,$fExt) = fileparse($file,'\.assoc\.logistic');
	my @fName = split('_', $fName);
	my $pop = pop @fName;

	my $lineNo = 1;

	while(<FILE>) {
		chomp;
		$_ =~ s/^(\s+)//;

		my @split = split('\s+|\t+');

		if ($lineNo == 1) {
			unless ($headerPrinted) {
				print OUT join("\t",@split,"POP"), "\n";
			}
		} else {
			print OUT join("\t",@split,$pop), "\n";
		}

		$lineNo++;
	}
}

print STDERR "completed\t$ID\n";
