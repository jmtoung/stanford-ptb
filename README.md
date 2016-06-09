Instructions:

1) Start with a bed file containing allg1g2, HRS, and 1000 genomes genotypes. 
example: /srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Data/allG1G2-1000G-HRS/step2-MergeFiles/step4-OmniDataFromSuyash-Phase3/step4-ExcludeRelatedsMissings/exclude-CG-AT-snps/exclude-141/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141.bed

/// ANCESTRY MATCHING PART
2) Create a filtered bed based HWE, allele freq, autosomal filter using the script 'pca_code/sge_Step11-MakeFilteredBed.pl'
Substitute $bfile with bed in 1) and use the $snpsToExclude that I have already prepared.

3) Get a list of SNPs that can be filtered based on LD using the script 'pca_code/step12-LD/sge_LD-Filter.pl'
Substitute $bfile with filtered bed from 2)

4) Create the filtered bed based on LD using the script 'pca_code/step12-LD/sge_Step2-MakeFilteredBed.pl'
Use $bfile from 2) and $pruneIn should be output from 3)

5) Next, you need to run this on your LOCAL machine. You need the R package 'SNPRelate', and I had issues installing it on the cluster. Pull down the bed file obtained from 4) down to your local machine. Install 'SNPRelate' package on R.

6) Run the 3 scripts in the folder 'pca_code/step13-SNPRelate'. You should get PCA plots from this. 
The code in the script 'Step2-MergePCs-with-CovariateFile-PredictPop.R' should give you a covariate file. The last few columns will contain 'SuperPopulationPrediction', which is the predicted ancestry group for HRS and allg1g2 based on 1000G ancestry matching.

7) Upload the covariate file from 6) (Step2 code) to the cluster.

8) 



