Couple of warnings:
- some code has hard-coded paths to $plink. Please double check all scripts and change to $plink of interest.

Instructions:

1) Start with a bed file containing allg1g2, HRS, and 1000 genomes genotypes. 
example: /srv/gsfs0/projects/butte/jmtoung/PTB-Genetics/Data/allG1G2-1000G-HRS/step2-MergeFiles/step4-OmniDataFromSuyash-Phase3/step4-ExcludeRelatedsMissings/exclude-CG-AT-snps/exclude-141/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141.bed

/// ANCESTRY MATCHING PART ///////////////////////////////////////////////////////////////////////////////////////

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

///////////////////////////////////////////////////////////////////////////////////////

/////// Update the master covariate file, fix a few things in this section ////////////////////////

8) Run the code in the 'prepare_covariate_file' folder:
- In 'sge_Step1...', the $cov should be the covariate file from 7). What we're doing here is creating the spontaneous preterm birth phenotype based on other covariates
- In 'sge_Step2...', we're fixing the sex covariate based on the original *.fam file. You should use the *.fam file associated with the bed from 1)
- In 'sge_Step3...', we are adding subtypes of preterm birth. Use the $subtypes file in the script, and the $cov from 'sge_Step2...'

9) Now we're done! ready to run the GWAS...

///// RUN GWAS PART //////////////////////////////////////////////////////////////////////

10) Run the scripts in 'run_gwas' in order
- 'sge_Step1...' -> creates lists of each population for GWAS
- 'sge_Step2...' -> runs the actual GWAS. Pay attention to line 50, that is the call to plink. Here, you should request --geno, --mind, --hwe, whatever PLINK options you want. 
- 'sge_Step3...' -> here, I'm doing a sanity check on whether the NMISS number plink reports corresponds to what I think it should be. Note that there is a $genotypeCounts file that is required in this script. You will need to generate this file yourself. To generate this file, please see the scripts in the 'misc/genotype_counts' folder.
-   * Start with running "MakeExtract..." to create an sge script. The sge script will call the perl 'Extract...' script (I think?).
-   * Then, do sge_Step2, sge_Step3, etc.
- 'sge_Step4...' -> Concatenate the results
- 'sge_Step5...' -> Extract the relevant ADD lines
- 'Step6...' -> In here, we annotate the results with the gene info. Note that I use ANNOVAR. You'll need to download your own version or call mine. Note that in these scripts there might be another file you'll need to generate much like in 'sge_Step3...'. You should be able to find my original code by digging around.
- 'Step7...' -> Add the adjusted P value from PLINK
- 'Step8...' -> Add genotype counts. Can get this from what you did in 'sge_Step3' with the scripts from 'misc/genotype_counts'
- 'Step9...' -> Add the p values from other populations
- 'Step10...' -> Sanity check
- 'Step11...' -> replace rs ids
- 'Step14...' -> extract top hits to plot them
- 'Step15...' -> can ignore this
- 'Step20...' -> extracted top hits for validation Sanger resequencing 



