#!/usr/bin/perl

use lib '/home/jmtoung/dev';
use GetLineCount;
use strict;
use File::Basename;
use Getopt::Long;
use List::Util qw[min max];

my $script = "/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Data/allG1G2-1000G-HRS/step2-MergeFiles/step4-OmniDataFromSuyash-Phase3/step4-ExcludeRelatedsMissings/exclude-CG-AT-snps/exclude-141/step9-ExtractGenotypes-ByDataset-PhenoPTB-BY-COLUMN/ExtractAllelesOrGenotypes-ByDataset-ByPheno-BY-COLUMN.pl";

my $ped = "/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Data/allG1G2-1000G-HRS/step2-MergeFiles/step4-OmniDataFromSuyash-Phase3/step4-ExcludeRelatedsMissings/exclude-CG-AT-snps/exclude-141/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141.ped";
my $map = "/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Data/allG1G2-1000G-HRS/step2-MergeFiles/step4-OmniDataFromSuyash-Phase3/step4-ExcludeRelatedsMissings/exclude-CG-AT-snps/exclude-141/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141.map";
my $cov = "/srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Data/allG1G2-1000G-HRS/step4-AddPhenotype/step4-OmniDataFromSuyash-Phase3/1000G-as-controls/exclude-CG-AT-snps/exclude-141/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_filtered_ld.cov";

# figure out number of columns in ped file
my $pedColCount = 4031506;

my $type = 'genotype';

my $datasetIdx = 14;

my $phenoIdx = 28; # ptb subtype

my $pheno2Idx = 26; # super population

my $numPerJob = 100000;

my ($pedName,$pedDir,$pedExt) = fileparse($ped,'\.ped');

my $outputDir = getPwd() . "/" . "sge";
runCommand("mkdir $outputDir") unless -d $outputDir;
my $outputTag = $pedName;
##############################################################



my @START;
my @END;
my @OUTPUT;
my $start = 6;
my $end;
while($start <= $pedColCount) {
	push(@START,$start);
	$end = min($pedColCount, $start + $numPerJob - 1);
	push(@END,$end);
	push(@OUTPUT,$outputDir . "/" . $outputTag . "." . formatInteger($start,6) . "-" . formatInteger($end,6) . "." . $type . "counts");
	$start = $end + 1;
}

my $START = join("\n",@START);
my $END = join("\n",@END);
my $OUTPUT = join("\n",@OUTPUT);

my $num_files = scalar(@START);

print STDOUT <<"END";
#!/usr/bin/perl 
#\$ -cwd ### use current directory
#\$ -M jmtoung\@stanford.edu
##\$ -m ea ### mail is to be sent at abort and end time
#\$ -j y ### combine stdout and stderr
#\$ -t 1-$num_files ### array job #
#\$ -V ### use current environment variables
##\$ -pe DJ 12 ### parallel threads
#\$ -l h_vmem=10G ### request memory
#\$ -l h_rt=100:00:00

use lib '/home/jmtoung/dev';
use GetLineCount;

my \$ID = \$ENV{SGE_TASK_ID} - 1;

my \$script = "$script";

my \$type = "$type";

my \$ped = "$ped";
my \$map = "$map";
my \$cov = "$cov";

my \$datasetIdx = $datasetIdx;

my \$phenoIdx = $phenoIdx;

my \$pheno2Idx = $pheno2Idx;

my \@start = qw(
$START
);

my \@end = qw(
$END
);

my \@output = qw(
$OUTPUT
);

my \$command = "perl \$script --type \$type --ped \$ped --map \$map --cov \$cov --datasetIdx \$datasetIdx --phenoIdx \$phenoIdx --pheno2Idx \$pheno2Idx --start \$start[\$ID] --end \$end[\$ID] --out \$output[\$ID]";

runCommand(\$command);

print "completed\\t\$ID\\n";

END


