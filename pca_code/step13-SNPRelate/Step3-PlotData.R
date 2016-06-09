library('ggplot2')

# get covariate file
cov_file <- '/Users/jmtoung/Butte/Projects/PTB-Genetics/Data/allG1G2-1000G-HRS/step3-PCA/step4-OmniDataFromSuyash-Phase3/1000G-as-controls/exclude-CG-AT-snps/exclude-141/step13-SNPRelate/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_filtered_ld.cov'
setwd(dirname(cov_file))
cov_data <- read.table(cov_file,header=TRUE)

subtype_file <- '/Users/jmtoung/Desktop/PROJECT/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_filtered_ld.cov'
cov_data <- read.table(subtype_file,header=TRUE)
cov_data <- cov_data[!(cov_data$dataset == "allg1g2" & cov_data$pheno_ptbSubtype!=1),]

cov_data$dataset <- factor(cov_data$dataset)
levels(cov_data$dataset) <- c('1000 Genomes Project','Cases','Controls')
cov_data$dataset <- factor(cov_data$dataset, levels=c('Cases','Controls','1000 Genomes Project'))


### shuffle the data
cov_data <- cov_data[sample(1:nrow(cov_data), nrow(cov_data), replace=FALSE),]

### calculate number in each group
for (d in levels(cov_data$dataset)) {
  num <- sum(cov_data$dataset==d)
  levels(cov_data$dataset)[levels(cov_data$dataset)==d] <- paste(" ",d, " (n=",num,") ",sep="")
}


### PLOT PC1 vs PC2 for by datasets
p <- ggplot(cov_data, aes(ev1,ev2,colour=dataset)) + geom_point(alpha=I(1/5)) +
  xlab("\nPrincipal Component 1") + ylab("Principal Component 2\n") + labs(colour="") +
  theme_bw() +  
  theme(legend.direction="horizontal", legend.position="top") +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))
ggsave(p,filename="MANU_PCA_PLOT_1_2_BY_DATASET.jpeg",width=10,height=7,dpi=300)

### PLOT PC1 vs PC3 for 3 datasets
p <- ggplot(cov_data, aes(ev1,ev3,colour=dataset)) + geom_point(alpha=I(1/5)) +
  xlab("\nPrincipal Component 1") + ylab("Principal Component 3\n") + labs(colour="") +
  theme_bw() +  
  theme(legend.direction="horizontal", legend.position="top") + 
  guides(colour = guide_legend(override.aes = list(alpha = 1)))
ggsave(p,filename="MANU_PCA_PLOT_1_3_BY_DATASET.jpeg",width=10,height=7,dpi=300)

### PLOT PC1 vs PC2 for POPULATION

## remove 1000 Genomes
cov_data <- cov_data[!grepl("1000 Genomes", cov_data$dataset),]
for (pop in levels(cov_data$SuperPopulationPrediction)) {
  num <- sum(cov_data$SuperPopulationPrediction==pop)
  levels(cov_data$SuperPopulationPrediction)[levels(cov_data$SuperPopulationPrediction)==pop] <- paste(pop, " (n=",num,") ",sep="")
}
p <- ggplot(cov_data, aes(ev1,ev2,colour=SuperPopulationPrediction)) + geom_point(alpha=I(1/5)) +
  xlab("\nPrincipal Component 1") + ylab("Principal Component 2\n") + labs(colour="") +
  theme_bw() +  
  theme(legend.direction="horizontal", legend.position="top") + 
  guides(colour = guide_legend(override.aes = list(alpha = 1)), shape = guide_legend(override.aes = list(alpha = 1)))
ggsave(p,filename="MANU_PCA_PLOT_1_2_BY_POPULATION.jpeg",width=10,height=7,dpi=300)

p <- ggplot(cov_data, aes(ev1,ev3,colour=SuperPopulationPrediction)) + geom_point(alpha=I(1/5)) +
  xlab("\nPrincipal Component 1") + ylab("Principal Component 3\n") + labs(colour="") +
  theme_bw() +  
  theme(legend.direction="horizontal", legend.position="top") + 
  guides(colour = guide_legend(override.aes = list(alpha = 1)), shape = guide_legend(override.aes = list(alpha = 1)))
ggsave(p,filename="MANU_PCA_PLOT_1_3_BY_POPULATION.jpeg",width=10,height=7,dpi=300)


### plot birth weight
cases <- cov_data[grep("Cases",cov_data$dataset),]
p <- ggplot(cases, aes(bwt)) + geom_histogram() +
  theme_bw() + xlab("\nBirth Weight (grams)") + 
  ylab("Count\n") 
ggsave(p,filename="MANU_BIRTH_WEIGHT.jpeg",width=7,height=5,dpi=300)

# plot gestational age
cases$gage_weeks <- cases$gage/7
p <- ggplot(cases, aes(gage_weeks)) + geom_histogram() + 
  theme_bw() + xlab("\nGestational Age (weeks)") + 
  ylab("Count\n") + coord_cartesian(xlim = c(24.88, 30)) 
ggsave(p,filename="MANU_GAGE.jpeg",width=7,height=5,dpi=300)
# #### 
# 
# # plot old PC1 vs PC2 and PC1 vs PC3 (coloured by bpd status)
# p <- ggplot(cov_data[cov_data$dataset=="1000 Genomes Controls",], aes(ev1, ev2, colour=SuperPopulation)) + geom_point(alpha=I(1/3)) +
#   labs(colour="1000G Population") + 
#   xlab("\nPrincipal Component 1") + 
#   ylab("Principal Component 2\n") +
#   labs(title='Principal Components 1 vs 2 of Genetic Ancestry\n1000 Genomes Project Control Individuals') +
#   theme(legend.direction="horizontal", legend.position="top")
# ggsave(p, filename = '1000G-Controls-pc1-vs-pc2.pdf',width=8,height=7,dpi=300)
# 
# p <- ggplot(cov_data[cov_data$dataset=="1000 Genomes Controls",], aes(ev1, ev3, colour=SuperPopulation)) + geom_point(alpha=I(1/3)) +
#   labs(colour="1000G Population") + 
#   xlab("\nPrincipal Component 1") + 
#   ylab("Principal Component 3\n") +
#   labs(title='Principal Components 1 vs 3 of Genetic Ancestry\n1000 Genomes Project Control Individuals') +
#   theme(legend.direction="horizontal", legend.position="top")
# ggsave(p, filename = '1000G-Controls-pc1-vs-pc3.pdf',width=8,height=7,dpi=300)
# 
# 
# p <- ggplot(cov_data, aes(ev1, ev2, colour=dataset)) + geom_point(alpha=I(1/5)) +
#   labs(colour="Dataset") + 
#   xlab("\nPrincipal Component 1") + 
#   ylab("Principal Component 2\n") +
#   labs(title='Principal Components 1 vs 2 of Genetic Ancestry') +
#   theme(legend.direction="horizontal", legend.position="top")
# ggsave(p, filename = 'all-pc1-vs-pc2.pdf',width=8,height=7,dpi=300)
# 
# p <- ggplot(cov_data, aes(ev1, ev3, colour=dataset)) + geom_point(alpha=I(1/5)) +
#   labs(colour="Dataset") + 
#   xlab("\nPrincipal Component 1") + 
#   ylab("Principal Component 3\n") +
#   labs(title='Principal Components 1 vs 3 of Genetic Ancestry') +
#   theme(legend.direction="horizontal", legend.position="top")
# ggsave(p, filename = 'all-pc1-vs-pc3.pdf',width=8,height=7,dpi=300)
# 
