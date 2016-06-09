#!/usr/bin/perl
#$ -cwd ### use current directory
#$ -M jmtoung@gmail.com
#$ -m ea ### mail is to be sent at abort and end time
#$ -j y ### combine stdout and stderr
#$ -t 1-1 ### array job #
#$ -V ### use current environment variables
##$ -pe DJ 12 ### parallel threads
#$ -l h_vmem=30G ### request memory
#$ -l h_rt=12:00:00

use strict;
use lib '/home/jmtoung/dev';
use GetLineCount;
use File::Basename;

my $ID = $ENV{SGE_TASK_ID} - 1;

my $plink = "/srv/gsfs0/projects/butte/ptb_nadav/PLINK/plink";

my $bfile = "/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Data/allG1G2-1000G-HRS/step3-PCA/step4-OmniDataFromSuyash-Phase3/1000G-as-controls/exclude-CG-AT-snps/exclude-141/20160602_NADAV/all_merged_mode_5_filtered";

my ($bname,$bdir,$bext) = fileparse($bfile,'');

my $out = $bname . "_ldpt2";

my $command = "$plink --bfile $bfile --indep-pairwise 50 5 .2 --out $out";

runCommand($command);

print STDERR "completed\t$ID\n";
