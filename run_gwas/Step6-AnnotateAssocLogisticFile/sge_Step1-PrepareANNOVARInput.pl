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

my $script = "/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Analysis/RunPlink/step4-OmniDataFromSuyash-Phase3/allG1G2-1000G-HRS-ccGWAS-covar/pc-sex/2014-09-23/Step6-AnnotateAssocLogisticFile/PrepareANNOVARInput.pl";

my $file = "/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Analysis/RunPlink/step4-OmniDataFromSuyash-Phase3/allG1G2-1000G-HRS-ccGWAS-covar/pc10ev-sex/2014-11-22_SPONTANEOUS-PTB-WITHOUT-1000G/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_ALL_ADD.assoc.logistic";

my ($fileName,$fileDir,$fileExt) = fileparse($file,'\.assoc\.logistic');

my $out = $fileName . ".annovarInput";

my $command = "perl $script --file $file --out $out";

runCommand($command);

print STDERR "completed\t$ID\n";
