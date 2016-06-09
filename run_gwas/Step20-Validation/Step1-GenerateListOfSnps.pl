#!/usr/bin/perl -w

use strict;

my $file = "/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Analysis/RunPlink/step4-OmniDataFromSuyash-Phase3/allG1G2-1000G-HRS-ccGWAS-covar/pc10ev-sex/2014-11-22_SPONTANEOUS-PTB-WITHOUT-1000G/Step15-AddWith1000GResults/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_ALL_ADD.refBase.annotation.exonicVariant.adjusted.counts.otherPValue.confirmNmiss.rsID.otherResults.fdrLePt05.txt";
my $tag = "snpsToValidate";

my %out;

open(FILE,$file) or die "[STDERR]: can't open $file: $!\n";

while(<FILE>) {
	chomp;

	my @split = split('\t');
	next if $split[0] eq 'CHR';

	my $pop = $split[9];	
	my $rsId = $split[1];
	my $chrom = $split[0];
	my $pos = $split[2];
	my $kpgId = $split[32];

	my @rsIds = split(';',$rsId);

	foreach my $r ( @rsIds) {
		unless (exists $out{$pop}) {
			my $file = $tag . "-$pop" . ".txt";
			open($out{$pop},">".$file) or die "[STDERR]: can't open $file: $!\n";
		}
		my $fh;
		$fh = $out{$pop};

		print {$fh} join("\t",$r,$chrom,$pos,$kpgId,$split[10],$split[11]), "\n"; 
	}
}
