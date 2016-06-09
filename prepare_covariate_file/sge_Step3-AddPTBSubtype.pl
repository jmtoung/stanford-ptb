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

my $cov = "all_merged_mode_5_filtered_ld.cov";
my $subtypes = "/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Data/allG1G2-PTBsubtypes/GWAS-G1G2-PTBsubtypes.txt";

my ($covName,$covDir,$covExt) = fileparse($cov,'\.cov');

my $out = $covName . ".covTemp";

open(OUT,">".$out) or die "[STDERR]: can't open $out: $!\n";

my %sub;
open(SUB,$subtypes) or die "[STDERR]: can't open $subtypes: $!\n";
while(<SUB>) {
	chomp;
	my @split = split('\t');
	my $type = $split[0];

	next if $type eq 'PTBsubtype';
	
	$type = '-9' if $type eq '';

	my $iid = $split[15];

	if (exists $sub{$iid}) {
		die "[STDERR]: Duplicate $iid\n";
	}

	$sub{$iid} = $type;

}

my ($covName,$covDir,$covExt) = fileparse($cov,'\.cov');

open(COV,$cov) or die "[STDERR]: can't open $cov: $!\n";
my $count=0;
while(<COV>) {
	chomp;
	my @split = split('\t');

	my $toAdd;
	if ($split[0] eq 'FID') {
		$toAdd = 'pheno_ptbSubtype';
	} else {

		my $iid = $split[1];
		my $dataset = $split[14];
		
		if ($dataset eq 'allg1g2') {
			$count++;
			if (exists $sub{$iid}) {
				$toAdd = $sub{$iid};
			} else {
				$toAdd = '-9';
				print STDERR "allg1g2_no_subtype\t$iid\n";
			}
		} else {
			$toAdd = '-9';
		}
	}

	push(@split,$toAdd);
	print OUT join("\t",@split), "\n";
}

runCommand("mv $out $cov");

print STDERR "total number of allg1g2 individuals: $count\n";
print STDERR "completed\t$ID\n";
