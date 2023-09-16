# Outliers, normality check, and interaction analyses

## Set your directory and load packages 
setwd("~/Desktop/GCase_working/new_plots_for_MS/")
library(data.table)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(ggpmisc)
library(plyr)

## Load data 
### These files have enzyme activity for GAA and GCase, as well as genotype status for top SNP in GAA locus
PPMI <- fread("../PPMI_enzyme-interaction_data.txt")
NY <- fread("../NY_enzyme-interaction_data.txt")
NY$`17:78061141:T:G` <- revalue(NY$`17:78061141:T:G_T`, c("0"="G/G", "2"="T/T", "1"="T/G"))
PPMI$`17:78061141:T:G` <- revalue(PPMI$`17:78061141:T:G_T`, c("0"="G/G", "2"="T/T", "1"="T/G"))

## Remove outliers from both cohorts 
Q_NY <- quantile(NY$GCase_Activity, probs=c(.25, .75), na.rm = FALSE)
Q_PPMI <- quantile(PPMI$GCase_Activity, probs=c(.25, .75), na.rm = FALSE)

iqr_NY <- IQR(NY$GCase_Activity)
iqr_PPMI <- IQR(PPMI$GCase_Activity)

up_NY <- Q_NY[2]+1.5*iqr_NY # Upper range  
low_NY <- Q_NY[1]-1.5*iqr_NY # Lower range

up_PPMI <- Q_PPMI[2]+1.5*iqr_PPMI # Upper range  
low_PPMI <- Q_PPMI[1]-1.5*iqr_PPMI # Lower range

NY_no_out <- subset(NY, NY$GCase_Activity > low_NY & NY$GCase_Activity < up_NY)
PPMI_no_out <- subset(PPMI, PPMI$GCase_Activity > low_PPMI & PPMI$GCase_Activity < up_PPMI)

### extract outliers to remove individuals from genetic data in GWAS script
NY_outliers <- subset(NY, NY$GCase_Activity < low_NY | NY$GCase_Activity > up_NY)
PPMI_outliers <- subset(PPMI, PPMI$GCase_Activity < low_PPMI | PPMI$GCase_Activity > up_PPMI)

write.table(NY_outliers, file = "NY_outliers.txt", sep = "\t",
            row.names = TRUE, col.names = NA)
write.table(PPMI_outliers, file = "PPMI_outliers.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

## Normality check 
### Make a normal dataset for comparison 
normal_data <- rnorm(200)

## Shapiro-Wilk test of normality 
shapiro.test(NY_no_out$GCase_Activity)
shapiro.test(PPMI_no_out$GCase_Activity)
shapiro.test(normal_data)

## Make the interaction plots
### All individuals 
all <- rbind(NY_no_out, PPMI_no_out, fill = TRUE)
all_cases <- all[ which(all$Status == 2),]
all_controls <- all[ which(all$Status == 1),]

all$`17:78061141:T:G` <- as.factor(all$`17:78061141:T:G`)
all_cases$`17:78061141:T:G` <- as.factor(all_cases$`17:78061141:T:G`)
all_controls$`17:78061141:T:G` <- as.factor(all_controls$`17:78061141:T:G`)

all_GAA_all <- ggplot(all, aes(x=GAA_activity, y=Gcase_activity, color=`17:78061141:T:G`)) + geom_point() + 
  scale_color_brewer(palette="Dark2") + 
  xlab("α-glucosidase activity (μmol/l/h)") + 
  ylab("GCase activity (μmol/l/h)") + 
  ggtitle("All samples") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  geom_smooth(method = "lm", se=F) +
  stat_cor(method = "pearson") +
  ylim(-1, 25) +
  theme_bw()

all_GAA_cases <- ggplot(all_cases, aes(x=GAA_activity, y=Gcase_activity, color=`17:78061141:T:G`)) + geom_point() +
  scale_color_brewer(palette="Dark2") + 
  xlab("α-glucosidase (μmol/l/h)") + 
  ylab("GCase activity (μmol/l/h)") +  
  ggtitle("Cases") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  geom_smooth(method = "lm", se=F) +
  stat_cor(method = "pearson") +
  ylim(-1, 25) +
  theme_bw()

all_GAA_controls <- ggplot(all_controls, aes(x=GAA_activity, y=Gcase_activity, color=`17:78061141:T:G`)) + geom_point() + 
  scale_color_brewer(palette="Dark2") + 
  xlab("α-glucosidase (μmol/l/h)") + 
  ylab("GCase activity (μmol/l/h)") + 
  ggtitle("Controls") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  geom_smooth(method = "lm", se=F) +
  stat_cor(method = "pearson") + 
  ylim(-1, 25) +
  theme_bw()

tiff(file="ALL_GCase_GAA_enzyme_corr_new.tiff", width = 10, height = 6, unit = "in", res = 800, pointsize = 6)
ggarrange(all_GAA_all, all_GAA_cases, all_GAA_controls, ncol = 3, nrow = 1, common.legend = TRUE, legend = "right")
dev.off()

### Columbia cohort 
NY_cases <- NY[ which(NY$Status == 2),]
NY_controls <- NY[ which(NY$Status == 1),]

NY_cases$`17:78061141:T:G` <- as.factor(NY_cases$`17:78061141:T:G`)
NY_controls$`17:78061141:T:G` <- as.factor(NY_controls$`17:78061141:T:G`)
NY$`17:78061141:T:G` <- as.factor(NY$`17:78061141:T:G`)

NY_GAA_all <- ggplot(NY, aes(x=GAA_activity, y=Gcase_activity, color=`17:78061141:T:G`)) + geom_point() + 
  scale_color_brewer(palette="Dark2") + 
  xlab("α-glucosidase (μmol/l/h)") + 
  ylab("GCase activity (μmol/l/h)") + 
  ggtitle("NY samples") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  geom_smooth(method = "lm", se=F) +
  stat_cor(method = "pearson") +
  ylim(-1, 25) +
  theme_bw()

NY_GAA_cases <- ggplot(NY_cases, aes(x=GAA_activity, y=Gcase_activity, color=`17:78061141:T:G`)) + geom_point() +
  scale_color_brewer(palette="Dark2") + 
  xlab("α-glucosidase (μmol/l/h)") + 
  ylab("GCase activity (μmol/l/h)") +  
  ggtitle("Cases") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  geom_smooth(method = "lm", se=F) +
  stat_cor(method = "pearson") +
  ylim(-1, 25) +
  theme_bw()

NY_GAA_controls <- ggplot(NY_controls, aes(x=GAA_activity, y=Gcase_activity, color=`17:78061141:T:G`)) + geom_point() + 
  scale_color_brewer(palette="Dark2") + 
  xlab("α-glucosidase (μmol/l/h)") + 
  ylab("GCase activity (μmol/l/h)") + 
  ggtitle("Controls") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  geom_smooth(method = "lm", se=F) +
  stat_cor(method = "pearson") + 
  ylim(-1, 25) +
  theme_bw()

tiff(file="NY_GCase_GAA_enzyme_corr_new.tiff", width = 10, height = 6, unit = "in", res = 800, pointsize = 6)
ggarrange(NY_GAA_all, NY_GAA_cases, NY_GAA_controls, ncol = 3, nrow = 1, common.legend = TRUE, legend = "right")
dev.off()

### PPMI cohort 
PPMI_cases <- PPMI[ which(PPMI$Status == 2),]
PPMI_controls <- PPMI[ which(PPMI$Status == 1),]

PPMI_cases$`17:78061141:T:G` <- as.factor(PPMI_cases$`17:78061141:T:G`)
PPMI_controls$`17:78061141:T:G` <- as.factor(PPMI_controls$`17:78061141:T:G`)

PPMI_GAA_all <- ggplot(PPMI, aes(x=GAA_Activity, y=GCase_Activity, color=`17:78061141:T:G`)) + geom_point() + 
  scale_color_brewer(palette="Dark2") + 
  xlab("α-glucosidase (μmol/l/h)") + 
  ylab("GCase activity (μmol/l/h)") + 
  ggtitle("PPMI samples") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  geom_smooth(method = "lm", se=F) +
  stat_cor(method = "pearson") +
  ylim(-1, 25) +
  theme_bw()

PPMI_GAA_cases <- ggplot(PPMI_cases, aes(x=GAA_Activity, y=GCase_Activity, color=`17:78061141:T:G`)) + geom_point() +
  scale_color_brewer(palette="Dark2") + 
  xlab("α-glucosidase (μmol/l/h)") + 
  ylab("GCase activity (μmol/l/h)") +  
  ggtitle("Cases") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  geom_smooth(method = "lm", se=F) +
  stat_cor(method = "pearson") +
  ylim(-1, 25) +
  theme_bw()

PPMI_GAA_controls <- ggplot(PPMI_controls, aes(x=GAA_Activity, y=GCase_Activity, color=`17:78061141:T:G`)) + geom_point() + 
  scale_color_brewer(palette="Dark2") + 
  xlab("α-glucosidase (μmol/l/h)") + 
  ylab("GCase activity (μmol/l/h)") + 
  ggtitle("Controls") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  geom_smooth(method = "lm", se=F) +
  stat_cor(method = "pearson") + 
  ylim(-1, 25) +
  theme_bw()

tiff(file="PPMI_GCase_GAA_enzyme_corr_new.tiff", width = 10, height = 6, unit = "in", res = 800, pointsize = 6)
ggarrange(PPMI_GAA_all, PPMI_GAA_cases, PPMI_GAA_controls, ncol = 3, nrow = 1, common.legend = TRUE, legend = "right")
dev.off()

## ANOVA
### Get column for absolute difference between GCase and GAA per person 
NY$dif_activity <- abs(NY$Gcase_activity-NY$GAA_activity)

### plots and stats 
ggplot(NY) +
  aes(x = `17:78061141:T:G_T`, y = dif_activity) +
  geom_boxplot()
aggregate(dif_activity ~ `17:78061141:T:G_T`,
          data = NY,
          function(x) round(c(mean = mean(x), sd = sd(x)), 2)
)

# ANOVA - two methods 
oneway.test(dif_activity ~ `17:78061141:T:G_T`, data = NY)

summary(aov(dif_activity ~ `17:78061141:T:G_T`, data = NY))

### Linear regression with interaction term 
## Basic covariates
fit_Gcase <- lm(Gcase_activity ~ `17:78061141:T:G_T` + GAA_activity + `17:78061141:T:G`*GAA_activity, data = NY) 
summary(fit_Gcase)

fit_GAA <- lm(GAA_activity ~ `17:78061141:T:G_T` + Gcase_activity + `17:78061141:T:G`*Gcase_activity, data = NY)
summary(fit_GAA)

## All covariates
fit_Gcase_allcov <- lm(Gcase_activity ~ `17:78061141:T:G` + GAA_activity + ASM_activity + GALC_activity + GLA_activity + Pheno + Age + Sex + G2019S + `17:78061141:T:G`*GAA_activity + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = NY)
summary(fit_Gcase_allcov)

fit_Gcase_PCs <- lm(Gcase_activity ~ `17:78061141:T:G` + GAA_activity + `17:78061141:T:G_T`*GAA_activity + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = NY)
summary(fit_Gcase_PCs)
