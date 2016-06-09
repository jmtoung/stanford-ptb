library(tools)
library(class)

file_eigenval <- '/Users/jmtoung/Butte/Projects/PTB-Genetics/Data/allG1G2-1000G-HRS/step3-PCA/step4-OmniDataFromSuyash-Phase3/1000G-as-controls/exclude-CG-AT-snps/exclude-141/step13-SNPRelate/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_filtered_ld.eigenval'
file_eigenvect <- '/Users/jmtoung/Butte/Projects/PTB-Genetics/Data/allG1G2-1000G-HRS/step3-PCA/step4-OmniDataFromSuyash-Phase3/1000G-as-controls/exclude-CG-AT-snps/exclude-141/step13-SNPRelate/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_filtered_ld.eigenvect'
file_sampleid <- '/Users/jmtoung/Butte/Projects/PTB-Genetics/Data/allG1G2-1000G-HRS/step3-PCA/step4-OmniDataFromSuyash-Phase3/1000G-as-controls/exclude-CG-AT-snps/exclude-141/step13-SNPRelate/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_filtered_ld.sampleid'
file_snpid <- '/Users/jmtoung/Butte/Projects/PTB-Genetics/Data/allG1G2-1000G-HRS/step3-PCA/step4-OmniDataFromSuyash-Phase3/1000G-as-controls/exclude-CG-AT-snps/exclude-141/step13-SNPRelate/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_filtered_ld.snpid'

eigenval <- read.table(file_eigenval)
eigenvect <- read.table(file_eigenvect)
sampleid <- read.table(file_sampleid)
snpid <- read.table(file_snpid)

colnames(sampleid) <- 'sampleid'

# get covariate file
cov_file <- '/Users/jmtoung/Butte/Projects/PTB-Genetics/Data/allG1G2-1000G-HRS/step2-MergeFiles/step4-OmniDataFromSuyash-Phase3/step3-Merge1000G-HRS/allG1G2-HRS-1000G_ex141.cov'
cov_data <- read.table(cov_file,header=TRUE)

# make a data frame
pca_data <- data.frame(sampleid = sampleid, 
                   ev1 = eigenvect[,1],
                   ev2 = eigenvect[,2],
                   ev3 = eigenvect[,3],
                   ev4 = eigenvect[,4],
                   ev5 = eigenvect[,5],
                   ev6 = eigenvect[,6],
                   ev7 = eigenvect[,7],
                   ev8 = eigenvect[,8],
                   ev9 = eigenvect[,9],
                   ev10 = eigenvect[,10],
                   stringsAsFactors = FALSE)


# write the PCA data out
out_data <- paste(c(file_path_sans_ext(file_eigenval),'.pca'), collapse="")
write.table(pca_data, file=out_data, sep="\t", col.names=NA)

# merge PCA with cov file
merge_data <- merge(cov_data, pca_data, by.y="sampleid", by.x="IID")
merge_data <- merge_data[,c(2,1,3:ncol(merge_data))]

# predict population
pc_col_idx <- which(grepl('ev',colnames(merge_data)))
training <- merge_data[merge_data$dataset=="1000g",c(pc_col_idx,1:2)]
testing <- merge_data[merge_data$dataset!="1000g",c(pc_col_idx,1:2)]
super_pop_idx <- which(grepl('SuperPopulation',colnames(merge_data)))
training_classification <- merge_data[merge_data$dataset=="1000g",c(super_pop_idx,1:2)]

num_neighbors <- 3
knn_output <- knn(training[,1:length(pc_col_idx)], testing[,1:length(pc_col_idx)], training_classification[,1], num_neighbors, prob=TRUE)

# add knn-ouptut predictions to data
testing <- cbind(testing,knn_output)
non_pc_idx <- which(!grepl('ev',colnames(testing)))
testing2 <- testing[non_pc_idx]
testing2 <- testing2[,c(3,1,2)]
colnames(testing2) <- colnames(training_classification)
classes <- rbind(training_classification, testing2)
colnames(classes)[which(colnames(classes)=='SuperPopulation')] <- 'SuperPopulationPrediction'

data_out <- merge(merge_data, classes, by.x="IID", by.y="IID")

# check that FIDs match
if (!sum(data_out$FID.x == data_out$FID.y) == nrow(data_out)) {
  warning("FID x not equal to FID y")
}

# drop FID.y
data_out[,which(colnames(data_out)=="FID.y")] <- NULL

colnames(data_out)[which(colnames(data_out)=='FID.x')] <- 'FID'

data_out <- data_out[c(2,1,3:length(data_out))]

# write data
out_cov <- paste(c(file_path_sans_ext(file_eigenval),'.cov'), collapse="")
write.table(data_out, file=out_cov, sep="\t", col.names=TRUE, quote=FALSE,  row.names=FALSE)



