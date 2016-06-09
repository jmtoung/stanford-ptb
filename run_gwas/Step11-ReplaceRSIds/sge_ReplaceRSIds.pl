#!/usr/bin/perl
#$ -cwd ### use current directory
#$ -M jmtoung@stanford.edu ### email address
##$ -m ea ### mail is to be sent at abort and end time
#$ -j y ### combine stdout and stderr
#$ -t 1-1 #5 ### array job #
#$ -V ### use current environment variables
##$ -pe DJ 12 ### parallel threads
##$ -l h_vmem=12G ### request memory
##$ -l h_rt=100:00:00

use strict;
use File::Basename;

my $ID = $ENV{SGE_TASK_ID} - 1;


my $oldMap = "/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Data/allG1G2-fixed/allG1G2-fixed-rmmiss.map";
my $newMap = "/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Data/allG1G2-fixed/allG1G2-fixed-rmmiss-rs.map";

my $file = "/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Analysis/RunPlink/step4-OmniDataFromSuyash-Phase3/allG1G2-1000G-HRS-ccGWAS-covar/pc10ev-sex/2014-11-22_SPONTANEOUS-PTB-WITHOUT-1000G/Step10-ConfirmNMiss/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_ALL_ADD.refBase.annotation.exonicVariant.adjusted.counts.otherPValue.confirmNmiss.txt";

my ($fileName,$fileDir,$fileExt) = fileparse($file,'\.txt');
my $out = $fileName . ".rsID" . $fileExt;
open(OUT,">".$out) or die "[STDERR]: Can't open $out: $!\n";

my $command = "paste $oldMap $newMap |";

open(COMMAND,$command) or die "[STDERR]: can't open $command: $!\n";

my %data;
while(<COMMAND>) {

	chomp;
	my @split = split('\t');

	my $oldID = $split[1];
	my $newID = $split[5];

	if (exists $data{$oldID}) {
		die "[SDTERR]: already exists $oldID\n";
	} else {
		$data{$oldID} = $newID;
	}

}
open(FILE,$file) or die "[STDERR]: can't open $file: $!\n";

my $notFound = 0;
while(<FILE>) {
	chomp;
	my @split = split('\t');

	my $snp = $split[1];
	
	unless ($snp eq 'SNP') {
		if (exists $data{$snp}) {
			push(@split,$split[1]);
			$split[1] = $data{$snp};
		} else {
			$notFound++;
			push(@split,'NA');
			print STDERR "could not find rsID for $snp\n";
#			die "[STDERR]: can't find $snp: $!\n";
		}
	} else {
		push(@split,'kgpID');
	}

	print OUT join("\t",@split), "\n";
}

print STDERR "notFound:\t$notFound\n";

print STDERR "completed\t$ID\n";
