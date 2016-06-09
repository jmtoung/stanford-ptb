#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use List::Util qw(min max);
use POSIX;

my $type; # type is 'allele' or 'genotype'
my $ped;
my $cov;
my $datasetIdx; # optional
my $phenoIdx; # optional
my $pheno2Idx; # optional
my $map; 
my $start;
my $end;
my $out;

my $options = GetOptions(
	"type=s" => \$type,
	"ped=s" => \$ped,
	"cov=s" => \$cov,
	"datasetIdx=i" => \$datasetIdx,
	"phenoIdx=i" => \$phenoIdx,
	"pheno2Idx=i" => \$pheno2Idx,
	"map=s" => \$map,
	"start=i" => \$start,
	"end=i" => \$end,
	"out=s" => \$out
);

unless ($type eq 'allele' or $type eq 'genotype') {
	die "[STDERR]: type should be allele or genotype\n";
}
print STDERR "type:\t$type\n";
print STDERR "ped:\t$ped\n";
print STDERR "cov:\t$cov\n";
print STDERR "datasetIdx:\t$datasetIdx\n" if defined $phenoIdx;
print STDERR "phenoIdx:\t$phenoIdx\n" if defined $phenoIdx;
print STDERR "pheno2Idx:\t$pheno2Idx\n" if defined $pheno2Idx;
print STDERR "map:\t$map\n";
print STDERR "start:\t$start\n" if $start;
print STDERR "end:\t$end\n" if $end;
print STDERR "out:\t$out\n";

open(PED,$ped) or die "[STDERR]: can't open $ped: $!\n";
open(MAP,$map) or die "[STDERR]: can't open $map: $!\n";
open(OUT,">".$out) or die "[STDERR]: can't open $out: $!\n";

open(COV,$cov) or die "[STDERR]: can't open $cov: $!\n";
my %cov;
while(<COV>) {
	chomp;
	my @split = split('\t');
	next if $split[0] eq 'FID';
	my $iid = $split[1];

	my $dataset = $split[$datasetIdx] if defined $datasetIdx;
	if (defined $datasetIdx) {
		unless (defined $dataset) {
			die "[STDERR]: undefined dataset $dataset\n";
		}
	}

	my $pheno = $split[$phenoIdx] if defined $phenoIdx;
	if (defined $phenoIdx) {
		unless (defined $pheno) {
			die "[STDERR]: undefined pheno $pheno\n";
		}
	}

	my $pheno2 = $split[$pheno2Idx] if defined $pheno2Idx;
	if (defined $pheno2Idx) {
		unless (defined $pheno2) {
			die "[STDERR]: undefined pheno2 $pheno2\n";
		}
	}

	# store the datasetIdx
	$cov{'dataset'}{$iid} = $dataset if defined $datasetIdx;
	$cov{'pheno'}{$iid} = $pheno if defined $phenoIdx;
	$cov{'pheno2'}{$iid} = $pheno2 if defined $pheno2Idx;

}

my %counts;
my $j=0;
while(<PED>) {

	chomp;
	
	my @split = split('\t|\s');
	my $iid = $split[1];

	my $dataset;
	if (defined $datasetIdx) {
		if (exists $cov{'dataset'}{$iid}) {
			$dataset = $cov{'dataset'}{$iid};
		} else {
			die "[STDERR]: can't find dataset for $iid\n";
		}
	} else {
		$dataset = 'NA';
	}

	my $pheno;
	if (defined $phenoIdx) {
		if (exists $cov{'pheno'}{$iid}) {
			$pheno = $cov{'pheno'}{$iid};
		} else {
			die "[STDERR]: can't find pheno for $iid\n";
		}
	} else {
		$pheno = 'NA';
	}

	my $pheno2;
	if (defined $pheno2Idx) {
		if (exists $cov{'pheno2'}{$iid}) {
			$pheno2 = $cov{'pheno2'}{$iid};
		} else {
			die "[STDERR]: can't find pheno2 for $iid\n";
		}
	} else {
		$pheno2 = 'NA';
	}


	for (my $i = $start; $i < min($end, scalar(@split)); $i+=2) {
		my $map_row = floor($i/2) - 2;
	
		my @types;
		if ($type eq 'allele') {
			push(@types,$split[$i],$split[$i+1]);
		} elsif ($type eq 'genotype') {
			my @g = ($split[$i], $split[$i+1]);
			@g = sort @g;				
			push(@types, join("",@g));
		} else {
			die "[STDERR]: weird type\n";
		}

		foreach my $t (@types) {
			$counts{$map_row}{$dataset}{$pheno}{$pheno2}{$t}++;
		}

	}	
	$j++;
	print STDERR $j, "\n" if ($j % 1000) == 0;
	
}

my $map_row = 1;

while(<MAP>) {
	chomp;
	
	my @split = split('\t|\s');

	my $colIdx = ($map_row-1)*2 + 6;
	my $counts;

	if (exists $counts{$map_row}) {
		$counts = $counts{$map_row};
	} else {
		if ($colIdx >= $start && $colIdx <= $end) {
			die "[STDERR]: missing allele_counts for row $map_row\n";
		}
	}

	foreach my $dataset (sort keys %{$counts}) {

		foreach my $pheno (sort keys %{$counts->{$dataset}}) {

			foreach my $pheno2 (sort keys %{$counts->{$dataset}{$pheno}}) {
				my @types;
				my @counts;

				foreach my $t (sort keys %{$counts->{$dataset}{$pheno}{$pheno2}}) {
					push(@types,$t);
					push(@counts,join(":",$t,$counts->{$dataset}{$pheno}{$pheno2}{$t}));
				}

				my @output = @split;
				push(@output,$dataset);
				push(@output,$pheno);
				push(@output,$pheno2);
				push(@output,join(",",@types));
				push(@output,join(",",@counts));

				print OUT join("\t",@output), "\n";
			}

		}
	}

	$map_row++;
}

print STDERR "completed:ExtractAllelesFromPedFile.pl\n";
