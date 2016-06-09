library('SNPRelate')
library(tools)

### EDIT THESE
# bed file
bed.fn <- '/Users/jmtoung/Butte/Projects/PTB-Genetics/Data/allG1G2-1000G-HRS/step3-PCA/step4-OmniDataFromSuyash-Phase3/1000G-as-controls/exclude-CG-AT-snps/exclude-141/step12-LD/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_filtered_ld.bed'
fam.fn <- '/Users/jmtoung/Butte/Projects/PTB-Genetics/Data/allG1G2-1000G-HRS/step3-PCA/step4-OmniDataFromSuyash-Phase3/1000G-as-controls/exclude-CG-AT-snps/exclude-141/step12-LD/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_filtered_ld.fam'
bim.fn <- '/Users/jmtoung/Butte/Projects/PTB-Genetics/Data/allG1G2-1000G-HRS/step3-PCA/step4-OmniDataFromSuyash-Phase3/1000G-as-controls/exclude-CG-AT-snps/exclude-141/step12-LD/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_filtered_ld.bim'

##################################################################
# creates gds file
file.gds <- paste(c(file_path_sans_ext(basename(bed.fn)), '.gds'), collapse="")
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, file.gds)

### load gds file
genofile <- openfn.gds(file.gds)

# run PCA
pca <- snpgdsPCA(genofile)

#eigenval file
eigenval_file <- paste(c(file_path_sans_ext(basename(bed.fn)), '.eigenval'), collapse="")
write.table(pca$eigenval, eigenval_file, sep="\t")

#eigenvect file
eigenvect_file <-  paste(c(file_path_sans_ext(basename(bed.fn)), '.eigenvect'), collapse="")
write.table(pca$eigenvect, eigenvect_file, sep="\t")

# snp id file
snpid_file <-  paste(c(file_path_sans_ext(basename(bed.fn)), '.snpid'), collapse="")
write.table(pca$snp.id, snpid_file, sep="\t")

# sample id file
sampleid_file <-  paste(c(file_path_sans_ext(basename(bed.fn)), '.sampleid'), collapse="")
write.table(pca$sample.id, sampleid_file, sep="\t")

