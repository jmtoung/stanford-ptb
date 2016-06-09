#!/usr/bin/perl
#$ -cwd ### use current directory
#$ -M jmtoung@stanford.edu ### email address
##$ -m ea ### mail is to be sent at abort and end time
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

my $plink = "~/Bin/plink-1.07-x86_64/plink";

my $bfile = "/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Data/allG1G2-1000G-HRS/step2-MergeFiles/step4-OmniDataFromSuyash-Phase3/step4-ExcludeRelatedsMissings/exclude-CG-AT-snps/exclude-141/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141";

my $snpsToExclude = "/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Data/allG1G2-1000G-HRS/step3-PCA/step4-OmniDataFromSuyash-Phase3/1000G-as-controls/exclude-CG-AT-snps/exclude-141/step10-CombineFilters/combinedSnpsToRemove.txt";

my ($bname,$bdir,$bext) = fileparse($bfile,'');

my $out = $bname . "_filtered";

my $command = "$plink --bfile $bfile --exclude $snpsToExclude --recode --make-bed --out $out";

runCommand($command);

print STDERR "completed\t$ID\n";
