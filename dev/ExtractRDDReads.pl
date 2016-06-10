#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use lib '/gpfs/fs121/h/toung/dev', '/gpfs/fs121/h/toung/lib/perl5/share/perl5';
use Database;
use Bio::DB::Sam;
use File::Basename;
use CalculateRDDStats_v5;
use Data::Dumper;

umask 0007;

my $rddFile;
my $region;
my $bam; ### don't give the home, start from aln
my $index;
my $alnDB; ### optional if you want to lookup alnID
my $alnID; ## optional
my $strand_specific;
my $unique_aln_only;
my $adapter_only;
my $minqual = 20;
my $combine_trim = 0;
my $edit_dist_filter = 0;

my $options = GetOptions(
	"rddFile=s" => \$rddFile,
	"region=s" => \$region,
	"bam=s" => \$bam,
	"index=s" => \$index,
	"alnDB=s" => \$alnDB,
	"alnID=s" => \$alnID,
	"strand_specific=i" => \$strand_specific,
	"unique_aln_only=i" => \$unique_aln_only,
	"adapter_only=i" => \$adapter_only,
	"minqual=i" => \$minqual,
	"combine_trim=i" => \$combine_trim,
	"edit_dist_filter=i" => \$edit_dist_filter
);

(print STDERR "rddFile:\t$rddFile\n") && -e $rddFile or die "[STDERR]: rddFile doesn't exist\n";
my ($regionPart,$regionStart,$regionEnd);
my $intervalUpdate = 1000;
if (defined $region) {
	($regionPart,$regionStart,$regionEnd) = split('-',$region);
	print STDERR "regionPart:\t$regionPart\n";
	print STDERR "regionStart:\t$regionStart\n";
	print STDERR "regionEnd:\t$regionEnd\n";
	$intervalUpdate = int($regionEnd - $regionStart + 1)/10;
	$intervalUpdate = 1 if $intervalUpdate < 1;
}

my @bam = split(',',$bam);
my @BAM; ### bam objects
foreach my $b (@bam) {
	($b =~ /unique|primary/ or die "[STDERR]: specified unique_aln_only but bam !~ /unique|primary/\n") if $unique_aln_only;
	my $BAM = Bio::DB::Sam->new(-bam => $b,-fasta => $index,-autoindex => 1);
        push(@BAM,$BAM);
}
print STDERR "numBamFiles:\t", scalar(@bam), "\n";

(print STDERR "index:\t$index\n") && -e $index or die "[STDERR]: index not defined\n";

(print STDERR "alnID:\t$alnID\n") && defined $alnID or die "[STDERR]: alnID not defined\n";

my @alnID = split(',',$alnID);
scalar(@alnID) == scalar(@bam) or die "[STDERR]: unequal number of alnID and bam\n";

(print STDERR "strand_specific:\t$strand_specific\n") && $strand_specific =~ /^(0|1)$/ or die "[STDERR]: strand specific not properly defined\n";
(print STDERR "unique_aln_only:\t$unique_aln_only\n") && $unique_aln_only =~ /^(0|1)$/ or die "[STDERR]: unique aln only not properly defined\n";
(print STDERR "adapter_only:\t$adapter_only\n") && $adapter_only =~ /^(0|1)$/ or die "[STDERR]: adapter only not defined \n";
(print STDERR "minqual:\t$minqual\n") && $minqual =~ /^[0-9]+$/ or die "[STDERR]: minqual not defined\n";
(print STDERR "combine_trim:\t$combine_trim\n") && $combine_trim =~ /^(0|1)$/ or die "[STDERR]: combine_trim not defined\n";
($edit_dist_filter == 0 || $edit_dist_filter == 1) && print STDERR "edit_dist_filter:\t$edit_dist_filter\n" or die "[STDERR]: edit_dist_filter not defined\n"; ### added 7/13/2012 2pm

### make output file
open(FILE,$rddFile) or die "[STDERR]: can't open $rddFile: $!\n";

my $counter = 1;

while(<FILE>) {
	chomp;

	next if (/^#/);
	
	my $lineNumber = $.;

	if (defined $region) {
		next if ($lineNumber < $regionStart);
		next if ($lineNumber > $regionEnd);
	}
	
	my @split = split('\t');
		
	my $chrom = $split[0];
	my $end = $split[2];
	my $strand = $split[3];
	my $ref_base = $split[4];
	my $rdd_base = $split[5];
	
	my @strands;
	if ($strand_specific) {
		$strand eq '+' or $strand eq '-' or die "[STDERR]: strand specific but strand is $strand\n";
		push(@strands,[$strand]);
	} else {
		push(@strands,['+','-']);
		if ($strand eq '-') {
			$ref_base =~ tr/ACGT/TGCA/;
			$rdd_base =~ tr/ACGT/TGCA/;
		}
	}
	my $RDD_OBJ = makeRddObject(\@BAM,$chrom,$end,$minqual,$combine_trim,$edit_dist_filter,\@bam,\@alnID); ###	\@full_bam,\@alnID);

	next if scalar(keys %{$RDD_OBJ}) == 0;

	foreach my $STRANDS (@strands) {
		makeFastq($RDD_OBJ,$STRANDS,$rdd_base,$adapter_only,$.);
	}

	print STDERR $counter . ".." if (($counter % $intervalUpdate) == 0); 

	$counter++;
}
close(FILE);

print STDERR "\n";
