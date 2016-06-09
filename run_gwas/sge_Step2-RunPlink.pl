#!/usr/bin/perl
#$ -cwd ### use current directory
#$ -M jmtoung@gmail.com
#$ -m ea ### mail is to be sent at abort and end time
#$ -j y ### combine stdout and stderr
#$ -t 1-5 ### array job #
#$ -V ### use current environment variables
##$ -pe DJ 12 ### parallel threads
#$ -l h_vmem=12G ### request memory
#$ -l h_rt=100:00:00

use strict;
use File::Basename;

my $ID = $ENV{SGE_TASK_ID} - 1;

my $bfile = "/srv/gsfs0/projects/butte/ptb_nadav/merged_mode_5";

my $cov = "/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Data/allG1G2-1000G-HRS/step4-AddPhenotype/step4-OmniDataFromSuyash-Phase3/1000G-as-controls/exclude-CG-AT-snps/exclude-141/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_filtered_ld.cov";

my $pheno_name = "pheno_ptb";

my $covar = $cov;

my $covar_name = "SEX,ev1,ev2,ev3,ev4,ev5,ev6,ev7,ev8,ev9,ev10";

my $hwe = ".001";


my @lists = qw(
IndividualLists-SPONTANEOUS-AFR.txt
IndividualLists-SPONTANEOUS-AMR.txt
IndividualLists-SPONTANEOUS-EAS.txt
IndividualLists-SPONTANEOUS-EUR.txt
IndividualLists-SPONTANEOUS-SAS.txt
);

##########

my ($listName,$listDir,$listExt) = fileparse($lists[$ID],'\.txt');
my @listName = split('-',$listName);
my $superPopulation = $listName[-1];

my $plink = "/home/jmtoung/Bin/plink_linux_x86_64/plink";

my ($bname,$bdir,$bext) = fileparse($bfile,'');

my $out = $bname . "_" . $superPopulation;

my $command = "$plink --noweb --bfile $bfile --logistic --pheno $cov --pheno-name $pheno_name --covar $covar --covar-name $covar_name --keep $lists[$ID] --hwe $hwe --adjust --out $out";

print STDERR $command, "\n";

!system($command) or die "[STDERR]: can't run $command: \$!\n";

print STDERR "completed\t$ID\n";
