rm(list=ls())

library('ggplot2')
library('reshape2')

file <- '/Users/jtoung/Dropbox/2015-11-30_MakeGgplots/allG1G2-HRS-1000G_exRelatedsMissings_exCG-AT-SNPs_ex141_ALL_ADD.assoc.logistic'

data <- read.table(file, header=TRUE)

# sort by P
data <- data[order(data$POP, data$P),]
data$LOG10P <- -log10(data$P)

data$LEXP <- 0
# calculate the expected
cols <- c("CHR", "SNP", "BP", "A1", "TEST", "NMISS", "OR", "STAT")
for (POP in unique(data$POP)) {
  data.sub <- data[data$POP==POP,]
  
  EXP <- 1:nrow(data.sub)/(nrow(data.sub) + 1)
  LEXP <- -log10(EXP)
  
  data$LEXP[data$POP==POP] <- LEXP
  
  p <- ggplot(data[data$POP==POP,], aes(LEXP, LOG10P)) + geom_point() +
    geom_abline(intercept=0, slope=1, colour="red", linetype="dotted") +
    xlab("Expected P-Value (-log10)") + ylab("Observed P-Value (-log10)") + 
    theme_bw() 
  ggsave(p, filename=paste("qq_plot_", POP, ".jpeg", sep=""), width=8, height=7, dpi=300)
  
  
}

data.sample <- data[sample(1:nrow(data), nrow(data)/30),]

data <- data[sample(1:nrow(data)),]
p <- ggplot(data, aes(LEXP, LOG10P, colour=POP)) + geom_point() +
  geom_abline(intercept=0, slope=1, colour="black", linetype="dotted") +
  xlab("Expected P-Value (-log10)") + ylab("Observed P-Value (-log10)") + labs(colour="") +
  theme_bw() 
ggsave(p, filename=paste("qq_plot.jpeg", sep=""), width=8, height=7, dpi=300)

p <- p + coord_cartesian(xlim=c(0, 10), ylim=c(0, 10))
ggsave(p, filename=paste("qq_plot_10-10.jpeg", sep=""), width=8, height=7, dpi=300)
