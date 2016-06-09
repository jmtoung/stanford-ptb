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

my $file = "/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Analysis/RunPlink/step4-OmniDataFromSuyash-Phase3/allG1G2-1000G-HRS-ccGWAS-covar/pc10ev-sex/2014-11-22_SPONTANEOUS-PTB-WITHOUT-1000G/Step9-AddPValueOtherPopulations/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_ALL_ADD.refBase.annotation.exonicVariant.adjusted.counts.otherPValue.txt";

my ($fileName,$fileDir,$fileExt) = fileparse($file,'\.txt');

my $out = $fileName . ".confirmNmiss" . $fileExt;

open(OUT,">".$out) or die "[STDERR]: Can't open $out: $!\n";

open(FILE,$file) or die "[STDERR]: can't open $file: $!\n";

my %stats;
while(<FILE>) {
	chomp;
	my @split = split('\t');

	if ($split[0] eq 'CHR') {
		push(@split,'nmiss_confirmed');
	} else {

		my $counts = 0;

		foreach my $i (22,23) {
			my @c = split(',',$split[$i]);
			foreach my $c (@c) {
				my @c2 = split(':',$c);
				next if $c2[0] eq '00';
				$counts += $c2[1];
			}
		}

		my $rsId = $split[1];
		my $nmiss = $split[5];
		my $pvalue = $split[8];

		if ($nmiss != $counts) {
			if ($pvalue eq 'NA') {
				$stats{'NA'}++;
				push(@split,'NA');
			} else {
				$stats{'Wrong'}++;
				push(@split,'Wrong');
			}
		} else {
			$stats{'Confirmed'}++;
			push(@split,'Confirmed');
		}
	} 

	print OUT join("\t",@split), "\n";
}

foreach my $stat (sort keys %stats) {
	print STDERR $stat, "\t",$stats{$stat}, "\n";
}
print STDERR "completed\t$ID\n";
