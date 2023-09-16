# Make summary statistics, Manhattan plots, and QQ plots from regression association results 
## Load data and set working directory
setwd("Desktop")
require(data.table)
require(qqman)
require(ggplot2)
library(dplyr)
library(ggrepel)
library(cowplot)

line <- read.table("marker.txt")
infos <- fread("infos")
colnames(infos) <- c("SNP","REF","ALT","ALT_Frq","Rsq","Genotyped")
assoc <- fread("assoc")
colnames(assoc) <- c("CHR","SNP","POS","A1_OG","TEST","N","Beta","SE","L95","U95","STAT","Pvalue")

## Calculate lambda
assoc$CHISQ <- qchisq(assoc$Pvalue, 1, lower.tail=FALSE)
lambda <- data.frame(median(assoc$CHISQ) / qchisq(0.5, 1))
lambda_table <- data.frame(line, lambda)
write.table(lambda_table, file = "lambda.txt", col.names=FALSE, row.names=FALSE)

## Make QQ plot
tiff("QQ.tiff", width = 4.25, height = 4, unit = "in", res = 800, pointsize = 6)
qq(assoc$Pvalue)
dev.off()

## Prepare summary stat files
dat <- merge(infos, assoc, by = "SNP", all.y = T)

# Make columns
dat$markerID <- dat$SNP
dat$minorAllele <- dat$ALT
dat$majorAllele <- dat$REF
dat$beta <- ifelse(dat$ALT_Frq <= 0.5, dat$Beta, dat$Beta*(-1))
dat$se <- dat$SE
dat$maf <- dat$ALT_Frq
dat$p <- dat$Pvalue
dat$BP <- dat$POS
dat$OR <- exp(dat$beta)
dat0 <- dat[,c("CHR","BP","markerID","minorAllele","majorAllele","beta", "OR","se","maf", "p", "N")]

# To METAL
write.table(dat0, "Metal.tab", quote = F, sep = "\t", row.names = F)

# To COJO
dat$A1 <- dat$minorAllele
dat$A2 <- dat$majorAllele
dat$freq <- dat$maf
dat$b <- dat$beta
dat1 <- dat[,c("SNP","A1","A2","freq","b","se","p","N")]
write.table(dat1, "COJO.tab", quote = F, sep = "\t", row.names = F)

# Full STATS
dat2 = dat[,c("CHR","BP","SNP","minorAllele","majorAllele","beta","se","OR","maf", "p", "N", "Rsq", "Genotyped")]
write.table(dat2, "fullSTATS.tab", quote = F, sep = "\t", row.names = F)

# Reduce computational time by subsetting half of SNPs below the p threshold of 1e-4
dat$P = dat$p
dat$zscore = dat$beta/dat$se
gwasResults = dat[,c("SNP", "CHR", "BP", "P", "zscore")]

test1 <- subset(gwasResults, P < 1e-4)
test2 <- subset(gwasResults, P > 1e-4)
length <- (nrow(test2))/2
test3 <- sample_n(test2, length, replace = FALSE)
subdat <- rbind(test1, test3)

## Manhattan plot
tiff("ManH.tiff", width = 5, height = 3.5, unit = "in", res = 800, pointsize = 6)
manhattan(subdat)
dev.off()
