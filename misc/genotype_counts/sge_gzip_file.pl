#!/usr/bin/perl
#$ -cwd ### use current directory
#$ -M jmtoung@gmail.com ###jmtoung@stanford.edu ### email address
#$ -m ea ### mail is to be sent at abort and end time
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

my $file = "allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141.genotypecounts";

my $gzip = $file . ".gzip";

my $command = "gzip -c $file > $gzip";

!system($command) or die "[Error]: $!\n";
