# This script lets you create a plot with multiple customized manhattans on top of each other 
## Set up libraries and load data 
setwd("Desktop")
library(data.table)
library(qqman)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(cowplot)

## swap out the files here with the summary statistics for your different cohorts/analyses
meta <- fread("fullSTATS_meta_allcov_noGBA.tab.gz")
colnames(meta) <- c("CHR","BP","SNP","minorAllele","majorAllele","maf","beta","se","OR","zscore","P")

NY <- fread("fullSTATS_NY_allcov_noGBA.tab.gz")
colnames(NY) <- c("CHR","BP","SNP","minorAllele","majorAllele","beta", "se","OR","maf", "P", "N","Rsq","Genotyped")

PPMI <- fread("fullSTATS_PPMI_allcov_noGBA.tab")
colnames(PPMI) <- c("CHR","BP","SNP","minorAllele","majorAllele","beta", "se","OR","maf", "P", "N","Rsq","Genotyped")

## Calculate lambdas 
meta$CHISQ <- qchisq(meta$P, 1, lower.tail=FALSE)
lambda_meta <- median(meta$CHISQ) / qchisq(0.5, 1)
lambda_meta

NY$CHISQ <- qchisq(NY$P, 1, lower.tail=FALSE)
lambda_NY <- median(NY$CHISQ) / qchisq(0.5, 1)
lambda_NY

PPMI$CHISQ <- qchisq(PPMI$P, 1, lower.tail=FALSE)
lambda_PPMI <- median(PPMI$CHISQ) / qchisq(0.5, 1)
lambda_PPMI

## QQ Plots 
### Can change to tiff or whatever you prefer for saving
png(file="QQ_meta_GCase_noGBA1.png", width = 5, height = 4, unit = "in", res = 800, pointsize = 6)
qq(meta$P, main = "Q-Q plot of Imputed GWAS p-values (META, GCase, No GBA1")
dev.off()

png(file="QQ_NY_GCase_noGBA1.png", width = 5, height = 4, unit = "in", res = 800, pointsize = 6)
qq(NY$p, main = "Q-Q plot of Imputed GWAS p-values (NY, GCase, No GBA1")
dev.off()

png(file="QQ_PPMI_GCase_noGBA1.png", width = 5, height = 4, unit = "in", res = 800, pointsize = 6)
qq(PPMI$p, main = "Q-Q plot of Imputed GWAS p-values (PPMI, GCase, No GBA1")
dev.off()

## Extract relevant data 
meta$zscore = meta$beta/meta$se
gwasResults_meta = meta[,c("SNP", "CHR", "BP", "P", "zscore")]

NY$zscore = NY$beta/NY$se
gwasResults_NY = NY[,c("SNP", "CHR", "BP", "P", "zscore")]

PPMI$zscore = PPMI$beta/PPMI$se
gwasResults_PPMI = PPMI[,c("SNP", "CHR", "BP", "P", "zscore")]

## Subset half of large p values to reduce computational time
test1 <- subset(gwasResults_meta, P < 1e-4) 
test2 <- subset(gwasResults_meta, P > 1e-4)
length <- (nrow(test2))/2
test3 <- sample_n(test2, length, replace = FALSE)
subdat_meta <- rbind(test1, test3)

test1 <- subset(gwasResults_NY, P < 1e-4) 
test2 <- subset(gwasResults_NY, P > 1e-4)
length <- (nrow(test2))/2
test3 <- sample_n(test2, length, replace = FALSE)
subdat_NY <- rbind(test1, test3)

test1 <- subset(gwasResults_PPMI, P < 1e-4) 
test2 <- subset(gwasResults_PPMI, P > 1e-4)
length <- (nrow(test2))/2
test3 <- sample_n(test2, length, replace = FALSE)
subdat_PPMI <- rbind(test1, test3)

## Set up to highlight peaks of interest
snpsOfInterest1 <- subset(subdat_meta, BP > 154400000 & BP < 159573111 & CHR == 1) # swap these out for whatever regions you want to highlight 
snpsOfInterest2 <- subset(subdat_meta, BP > 78020000 & BP < 78140000 & CHR == 17) # swap these out for whatever regions you want to highlight
snpsOfInterestAll <- rbind(snpsOfInterest1, snpsOfInterest2)
snpsOfInterest <- snpsOfInterestAll$SNP

dat_meta <- subdat_meta %>% 
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
  select(-chr_len) %>%
  left_join(subdat_meta, ., by=c("CHR"="CHR")) %>%
  arrange(CHR, BP) %>%
  mutate(BPcum=BP+tot) %>%
  mutate(is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no"))

dat_NY <- subdat_NY %>% 
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
  select(-chr_len) %>%
  left_join(subdat_NY, ., by=c("CHR"="CHR")) %>%
  arrange(CHR, BP) %>%
  mutate(BPcum=BP+tot) %>%
  mutate(is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no"))

dat_PPMI <- subdat_PPMI %>% 
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
  select(-chr_len) %>%
  left_join(subdat_PPMI, ., by=c("CHR"="CHR")) %>%
  arrange(CHR, BP) %>%
  mutate(BPcum=BP+tot) %>%
  mutate(is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no"))

axisdf_meta = dat_meta %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
axisdf_NY = dat_NY %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
axisdf_PPMI = dat_PPMI %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

## Check that for each analysis the SNP is actually present at the top of the signal for labelling purposes, aka make sure snpsToAnno isn't empty after the following commands
### Otherwise the label won't show up 
snpsToAnno_meta <- subset(dat_meta, SNP == "1:155205634:T:C" | SNP == "17:78056851:C:T")
  snpsToAnno_meta$Gene <- c("GBA1","GAA")
snpsToAnno_NY <- subset(dat_NY, SNP == "1:155205634:T:C" | SNP == "17:78056851:C:T")
  snpsToAnno_NY$Gene <- c("GBA1","GAA")
snpsToAnno_PPMI <- subset(dat_PPMI, SNP == "1:155205634:T:C" | SNP == "17:78096483:A:C")
  snpsToAnno_PPMI$Gene <- c("GBA1","GAA")
  
sig <- 5e-8

## Make plots for each cohort/analysis 
### You can adjust the height of the y scale for each to suit whatever your biggest -log10 p values are
meta_plot <- ggplot(dat_meta, aes(x=BPcum, y=-log10(P))) +
  geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("darkseagreen","aquamarine4"), 22 )) +
  scale_x_continuous( label = axisdf_meta$CHR, breaks= axisdf_meta$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,50)) +
  xlab("Chromsome") +
  geom_point(data=subset(dat_meta, is_highlight=="yes"), color="slateblue3", size=2) +
  theme_bw() +
  theme(legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()) +
  geom_label(data = snpsToAnno_meta, aes(x=BPcum+150000000, y=-log10(P)+2.25, label=Gene), color="black", 
             size=4 , angle=45, fontface="bold")


NY_plot <- ggplot(dat_NY, aes(x=BPcum, y=-log10(P))) +
  geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("darkseagreen","aquamarine4"), 22 )) +
  scale_x_continuous( label = axisdf_NY$CHR, breaks= axisdf_NY$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,40)) +
  xlab("Chromsome") +
  geom_point(data=subset(dat_NY, is_highlight=="yes"), color="slateblue3", size=2) +
  theme_bw() +
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  geom_label(data = snpsToAnno_NY, aes(x=BPcum+150000000, y=-log10(P)+2.25, label=Gene), color="black", 
             size=4 , angle=45, fontface="bold")

PPMI_plot <- ggplot(dat_PPMI, aes(x=BPcum, y=-log10(P))) +
  geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("darkseagreen","aquamarine4"), 22 )) +
  scale_x_continuous( label = axisdf_PPMI$CHR, breaks= axisdf_PPMI$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,15)) +
  xlab("Chromsome") +
  geom_point(data=subset(dat_PPMI, is_highlight=="yes"), color="slateblue3", size=2) +
  theme_bw() +
  theme(legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()) +
  geom_label(data = snpsToAnno_PPMI, aes(x=BPcum+150000000, y=-log10(P)+0.75, label=Gene), color="black", 
           size=4 , angle=45, fontface="bold")

# combine plots using cowplot 
tiff(file="GCase_multi_allcov_noGBA.tiff", width = 8, height = 12, unit = "in", res = 800, pointsize = 6)
cowplot::plot_grid(NY_plot + theme(axis.title.x = element_blank() ), 
                   PPMI_plot + theme(axis.title.x = element_blank() ), 
                   meta_plot,
                   nrow = 3,
                   labels = "auto")
dev.off()
