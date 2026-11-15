### David O'Donnell, adapted from Jinliang Yang
### Dec 1st, 2025
### data checking


## Install some required packages if have not installed

install.packages(c("Matrix","lme4"))
install.packages("devtools")
devtools::install_github("jyanglab/g3tools")

library("Matrix")
library("lme4")
library("devtools")
library("g3tools")

setwd("/Users/David_ODonnell/Documents/Github/n2qtl/")

pheno <- read.csv("profiling/Zea-Trial-1_Metadata.csv", na.strings = c("", "NA", "NULL", "missing"))
dim(pheno)
#784  22

pheno$NDFA_T1_B73_ref <- as.numeric(pheno$NDFA_T1_B73_ref)
pheno$NDFA_T2_B73_ref <- as.numeric(pheno$NDFA_T2_B73_ref)

#n2 <- subset(pheno, !is.na(X15NT1))

hist(pheno$d15N_T1, breaks=30)
hist(pheno$NDFA_T1_B73xMo17_ref, breaks=30)
hist(pheno$NDFA_T1_B73_ref, breaks=30)
hist(pheno$d15N_T2, breaks=30)
hist(pheno$NDFA_T2_B73xMo17_ref, breaks=30)
hist(pheno$NDFA_T2_B73_ref, breaks=30)
hist(pheno$HEIGHT_in, breaks=30)
hist(pheno$NODES, breaks=30)
hist(pheno$NODES_w_BR, breaks=30)
hist(pheno$NODES_w_AR, breaks=30)
hist(pheno$NODES_w_BRandAR, breaks=30)
hist(pheno$DAYS_TO_TASSEL, breaks=30)
hist(pheno$DAYS_TO_POLLEN, breaks=30)
hist(pheno$COBS, breaks=30)
hist(pheno$SILKS, breaks=30)


#BLUP and heritability estimation

fit1 <- get_BLUP(data = pheno, model = d15N_T1 ~ (1 | entry) + (1 | Longitude) 
         + (1 | Latitude), which.factor = "entry",
         outfile = "data/blup_evolution_d15N_T1.csv")

get_H2(fit1, numerator="entry", denominator=data.frame(f=c("entry", "Residual"),
                              df=c(1, 1)))
#0.6101074



fit2 <- get_BLUP(data = pheno, model = NDFA_T1_B73_ref ~ (1 | entry) + (1 | Longitude) 
                + (1 | Latitude), which.factor = "entry",
                outfile = "data/blup_evolution_NDFA_T1_B73_ref.csv")

get_H2(fit2, numerator="entry", denominator=data.frame(f=c("entry", "Residual"),
                                                      df=c(1, 1)))
#0.5831606



fit3 <- get_BLUP(data = pheno, model = d15N_T2 ~ (1 | entry) + (1 | Longitude) 
                 + (1 | Latitude), which.factor = "entry",
                 outfile = "data/blup_evolution_d15N_T2.csv")

get_H2(fit3, numerator="entry", denominator=data.frame(f=c("entry", "Residual"),
                                                       df=c(1, 1)))
#0.7252177



fit4 <- get_BLUP(data = pheno, model = NDFA_T2_B73_ref ~ (1 | entry) + (1 | Longitude) 
                 + (1 | Latitude), which.factor = "entry",
                 outfile = "data/blup_evolution_NDFA_T2_B73_ref.csv")

get_H2(fit4, numerator="entry", denominator=data.frame(f=c("entry", "Residual"),
                                                       df=c(1, 1)))
#0.7251031



fit5 <- get_BLUP(data = pheno, model = HEIGHT_in ~ (1 | entry) + (1 | Longitude) 
                 + (1 | Latitude), which.factor = "entry",
                 outfile = "data/blup_evolution_HEIGHT_in.csv")

get_H2(fit5, numerator="entry", denominator=data.frame(f=c("entry", "Residual"),
                                                       df=c(1, 1)))
#0.880308



fit6 <- get_BLUP(data = pheno, model = NODES ~ (1 | entry) + (1 | Longitude) 
                 + (1 | Latitude), which.factor = "entry",
                 outfile = "data/blup_evolution_NODES.csv")

get_H2(fit6, numerator="entry", denominator=data.frame(f=c("entry", "Residual"),
                                                       df=c(1, 1)))
#0.8669735



fit7 <- get_BLUP(data = pheno, model = NODES_w_BR ~ (1 | entry) + (1 | Longitude) 
                 + (1 | Latitude), which.factor = "entry",
                 outfile = "data/blup_evolution_NODES_w_BR.csv")

get_H2(fit7, numerator="entry", denominator=data.frame(f=c("entry", "Residual"),
                                                       df=c(1, 1)))
#0.6050732



fit8 <- get_BLUP(data = pheno, model = NODES_w_AR ~ (1 | entry) + (1 | Longitude) 
                 + (1 | Latitude), which.factor = "entry",
                 outfile = "data/blup_evolution_NODES_w_AR.csv")

get_H2(fit8, numerator="entry", denominator=data.frame(f=c("entry", "Residual"),
                                                       df=c(1, 1)))
#0.5898972



fit9 <- get_BLUP(data = pheno, model = NODES_w_BRandAR ~ (1 | entry) + (1 | Longitude) 
                 + (1 | Latitude), which.factor = "entry",
                 outfile = "data/blup_evolution_NODES_w_BRandAR.csv")

get_H2(fit9, numerator="entry", denominator=data.frame(f=c("entry", "Residual"),
                                                       df=c(1, 1)))
#0.6239602



fit10 <- get_BLUP(data = pheno, model = DAYS_TO_TASSEL ~ (1 | entry) + (1 | Longitude) 
                 + (1 | Latitude), which.factor = "entry",
                 outfile = "data/blup_evolution_DAYS_TO_TASSEL.csv")

get_H2(fit10, numerator="entry", denominator=data.frame(f=c("entry", "Residual"),
                                                       df=c(1, 1)))
#0.9030854



fit11 <- get_BLUP(data = pheno, model = DAYS_TO_POLLEN ~ (1 | entry) + (1 | Longitude) 
                 + (1 | Latitude), which.factor = "entry",
                 outfile = "data/blup_evolution_DAYS_TO_POLLEN.csv")

get_H2(fit11, numerator="entry", denominator=data.frame(f=c("entry", "Residual"),
                                                       df=c(1, 1)))
#0.9125291



fit12 <- get_BLUP(data = pheno, model = COBS ~ (1 | entry) + (1 | Longitude) 
                 + (1 | Latitude), which.factor = "entry",
                 outfile = "data/blup_evolution_COBS.csv")

get_H2(fit12, numerator="entry", denominator=data.frame(f=c("entry", "Residual"),
                                                       df=c(1, 1)))
#0.3539288



fit13 <- get_BLUP(data = pheno, model = SILKS ~ (1 | entry) + (1 | Longitude) 
                 + (1 | Latitude), which.factor = "entry",
                 outfile = "data/blup_evolution_SILKS.csv")

get_H2(fit13, numerator="entry", denominator=data.frame(f=c("entry", "Residual"),
                                                       df=c(1, 1)))
#0.1678926


