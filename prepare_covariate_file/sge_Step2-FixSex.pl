#!/usr/bin/perl
#$ -cwd ### use current directory
#$ -M jmtoung@stanford.edu ### email address
##$ -m ea ### mail is to be sent at abort and end time
#$ -j y ### combine stdout and stderr
#$ -t 1-1 ### array job #
#$ -V ### use current environment variables
##$ -pe DJ 12 ### parallel threads
#$ -l h_vmem=14G ### request memory
#$ -l h_rt=12:00:00

use lib '/home/jmtoung/dev';
use GetLineCount;
use strict;
use File::Basename;

my $ID = $ENV{SGE_TASK_ID} - 1;

my $cov = "all_merged_mode_5_filtered_ld.cov";

my $fam = "/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Data/allG1G2-1000G-HRS/step2-MergeFiles/step4-OmniDataFromSuyash-Phase3/step4-ExcludeRelatedsMissings/exclude-CG-AT-snps/exclude-141/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141.fam";

my ($covName,$covDir,$covExt) = fileparse($cov,'\.cov');

my $out = $covName . ".covTemp";

open(OUT,">".$out) or die "[STDERR]: can't open $out: $!\n";

open(COV,$cov) or die "[STDERR]: can't open $cov: $!\n";

open(FAM,$fam ) or die "[STDERR]: can't open $fam: $!\n";
my %data;
while(<FAM>) {
	chomp;
	my @split = split('\t|\s');
	my $sex = $split[4];

	if ($sex ne '2' && $sex ne '1') {
		die "[STDERR]: weird sex $sex\n";
	}

	my $iid = $split[1];
	$data{$iid} = $sex;
}

while(<COV>) {
	chomp;
	my @split = split('\t');

	if ($split[0] eq 'FID') {
		
	} else {
		my $sex;
		my $iid = $split[1];

		if (exists $data{$iid}) {
			$sex = $data{$iid};
		} else {
			die "[STDERR]: can't find sex for $iid\n";
		}
		$split[4] = $sex;
	}

	# take out extra sex
	splice(@split,12,1);
	print OUT join("\t",@split), "\n";
}

runCommand("mv $out $cov");

print STDERR "completed\t$ID\n";
