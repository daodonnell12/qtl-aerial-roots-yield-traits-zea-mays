### David O'Donnell, adapted from Jinliang Yang
### Dec 1st, 2025
### data checking
### Additional customization David O'Donnell


## Install some required packages if have not installed

install.packages("lme4")
install.packages("devtools")
devtools::install_github("jyanglab/g3tools")

library("lme4")
library("devtools")
library("g3tools")



pheno <- read.csv("profiling/Zea-Trial-1_Metadata_w-Long-Lat.csv")
dim(pheno)
#687  33

#n2 <- subset(pheno, !is.na(X15NT1))

hist(pheno$d15N_T1, breaks=30)
hist(pheno$NDFA_T1_B73.ref, breaks=30)
hist(pheno$d15N_T2, breaks=30)
hist(pheno$NDFA_T2_B73.ref, breaks=30)
hist(pheno$HEIGHT_in, breaks=30)
hist(pheno$NODES, breaks=30)
hist(pheno$NODES_w_BR, breaks=30)
hist(pheno$NODES_w_AR, breaks=30)
hist(pheno$NODES_w_BRandAR, breaks=30)
hist(pheno$DAYS_TO_TASSEL, breaks=30)
hist(pheno$DAYS_TO_POLLEN, breaks=30)
hist(pheno$COBS, breaks=30)
hist(pheno$SILKS, breaks=30)


fit <- get_BLUP(data = pheno, model = d15N_T1 ~ (1 | GROUP) + (1 | Longitude) 
         + (1 | Latitude), which.factor = "GROUP",
         outfile = "data/blup_evolution_d15N_T1.csv")

get_H2(fit, numerator="GROUP", denominator=data.frame(f=c("GROUP", "Residual"),
                              df=c(1, 1)))



fit2 <- get_BLUP(data = pheno, model = NDFA_T1_B73.ref ~ (1 | GROUP) + (1 | Longitude) 
                + (1 | Latitude), which.factor = "GROUP",
                outfile = "data/blup_evolution_NDFA_T1_B73.ref.csv")

get_H2(fit2, numerator="GROUP", denominator=data.frame(f=c("GROUP", "Residual"),
                                                      df=c(1, 1)))



fit3 <- get_BLUP(data = pheno, model = d15N_T2 ~ (1 | GROUP) + (1 | Longitude) 
                 + (1 | Latitude), which.factor = "GROUP",
                 outfile = "data/blup_evolution_d15N_T2.csv")

get_H2(fit3, numerator="GROUP", denominator=data.frame(f=c("GROUP", "Residual"),
                                                       df=c(1, 1)))



fit4 <- get_BLUP(data = pheno, model = NDFA_T2_B73.ref ~ (1 | GROUP) + (1 | Longitude) 
                 + (1 | Latitude), which.factor = "GROUP",
                 outfile = "data/blup_evolution_NDFA_T2_B73.ref.csv")

get_H2(fit4, numerator="GROUP", denominator=data.frame(f=c("GROUP", "Residual"),
                                                       df=c(1, 1)))



fit5 <- get_BLUP(data = pheno, model = HEIGHT_in ~ (1 | GROUP) + (1 | Longitude) 
                 + (1 | Latitude), which.factor = "GROUP",
                 outfile = "data/blup_evolution_HEIGHT_in.csv")

get_H2(fit5, numerator="GROUP", denominator=data.frame(f=c("GROUP", "Residual"),
                                                       df=c(1, 1)))



fit6 <- get_BLUP(data = pheno, model = NODES ~ (1 | GROUP) + (1 | Longitude) 
                 + (1 | Latitude), which.factor = "GROUP",
                 outfile = "data/blup_evolution_NODES.csv")

get_H2(fit6, numerator="GROUP", denominator=data.frame(f=c("GROUP", "Residual"),
                                                       df=c(1, 1)))



fit7 <- get_BLUP(data = pheno, model = NODES_w_BR ~ (1 | GROUP) + (1 | Longitude) 
                 + (1 | Latitude), which.factor = "GROUP",
                 outfile = "data/blup_evolution_NODES_w_BR.csv")

get_H2(fit7, numerator="GROUP", denominator=data.frame(f=c("GROUP", "Residual"),
                                                       df=c(1, 1)))



fit8 <- get_BLUP(data = pheno, model = NODES_w_AR ~ (1 | GROUP) + (1 | Longitude) 
                 + (1 | Latitude), which.factor = "GROUP",
                 outfile = "data/blup_evolution_NODES_w_AR.csv")

get_H2(fit8, numerator="GROUP", denominator=data.frame(f=c("GROUP", "Residual"),
                                                       df=c(1, 1)))



fit9 <- get_BLUP(data = pheno, model = NODES_w_BRandAR ~ (1 | GROUP) + (1 | Longitude) 
                 + (1 | Latitude), which.factor = "GROUP",
                 outfile = "data/blup_evolution_NODES_w_BRandAR.csv")

get_H2(fit9, numerator="GROUP", denominator=data.frame(f=c("GROUP", "Residual"),
                                                       df=c(1, 1)))



fit10 <- get_BLUP(data = pheno, model = DAYS_TO_TASSEL ~ (1 | GROUP) + (1 | Longitude) 
                 + (1 | Latitude), which.factor = "GROUP",
                 outfile = "data/blup_evolution_DAYS_TO_TASSEL.csv")

get_H2(fit10, numerator="GROUP", denominator=data.frame(f=c("GROUP", "Residual"),
                                                       df=c(1, 1)))



fit11 <- get_BLUP(data = pheno, model = DAYS_TO_POLLEN ~ (1 | GROUP) + (1 | Longitude) 
                 + (1 | Latitude), which.factor = "GROUP",
                 outfile = "data/blup_evolution_DAYS_TO_POLLEN.csv")

get_H2(fit11, numerator="GROUP", denominator=data.frame(f=c("GROUP", "Residual"),
                                                       df=c(1, 1)))



fit12 <- get_BLUP(data = pheno, model = COBS ~ (1 | GROUP) + (1 | Longitude) 
                 + (1 | Latitude), which.factor = "GROUP",
                 outfile = "data/blup_evolution_COBS.csv")

get_H2(fit12, numerator="GROUP", denominator=data.frame(f=c("GROUP", "Residual"),
                                                       df=c(1, 1)))



fit13 <- get_BLUP(data = pheno, model = SILKS ~ (1 | GROUP) + (1 | Longitude) 
                 + (1 | Latitude), which.factor = "GROUP",
                 outfile = "data/blup_evolution_SILKS.csv")

get_H2(fit13, numerator="GROUP", denominator=data.frame(f=c("GROUP", "Residual"),
                                                       df=c(1, 1)))



fit14 <- get_BLUP(data = pheno, model = NODES_w_AR ~ (1 | GROUP) + (1 | Longitude) 
                  + (1 | Latitude) + (1 | GROUP:Longitude) + (1 | GROUP:Latitude), which.factor = "GROUP",
                  outfile = "data/blup_evolution_NODES_w_AR_v2.csv")

get_H2(fit14, numerator="GROUP", denominator=data.frame(f=c("GROUP", "Residual"),
                                                        df=c(1, 1)))





fit21 <- get_BLUP(data = pheno, model = d15N_T1 ~ (1 | GROUP) + (1 | Longitude) 
                + (1 | Latitude) + (1 | GROUP:Longitude) + (1 | GROUP:Latitude), which.factor = "GROUP",
                outfile = "data/blup_evolution_d15N_T1_v2.csv")

get_H2(fit21, numerator="GROUP", denominator=data.frame(f=c("GROUP", "Residual"),
                                                      df=c(1, 1)))



fit22 <- get_BLUP(data = pheno, model = NDFA_T1_B73.ref ~ (1 | GROUP) + (1 | Longitude) 
                 + (1 | Latitude) + (1 | GROUP:Longitude) + (1 | GROUP:Latitude), which.factor = "GROUP",
                 outfile = "data/blup_evolution_NDFA_T1_B73.ref_v2.csv")

get_H2(fit22, numerator="GROUP", denominator=data.frame(f=c("GROUP", "Residual"),
                                                       df=c(1, 1)))



fit23 <- get_BLUP(data = pheno, model = d15N_T2 ~ (1 | GROUP) + (1 | Longitude) 
                 + (1 | Latitude) + (1 | GROUP:Longitude) + (1 | GROUP:Latitude), which.factor = "GROUP",
                 outfile = "data/blup_evolution_d15N_T2_v2.csv")

get_H2(fit23, numerator="GROUP", denominator=data.frame(f=c("GROUP", "Residual"),
                                                       df=c(1, 1)))



fit24 <- get_BLUP(data = pheno, model = NDFA_T2_B73.ref ~ (1 | GROUP) + (1 | Longitude) 
                 + (1 | Latitude) + (1 | GROUP:Longitude) + (1 | GROUP:Latitude), which.factor = "GROUP",
                 outfile = "data/blup_evolution_NDFA_T2_B73.ref_v2.csv")

get_H2(fit24, numerator="GROUP", denominator=data.frame(f=c("GROUP", "Residual"),
                                                       df=c(1, 1)))



fit25 <- get_BLUP(data = pheno, model = HEIGHT_in ~ (1 | GROUP) + (1 | Longitude) 
                 + (1 | Latitude) + (1 | GROUP:Longitude) + (1 | GROUP:Latitude), which.factor = "GROUP",
                 outfile = "data/blup_evolution_HEIGHT_in_v2.csv")

get_H2(fit25, numerator="GROUP", denominator=data.frame(f=c("GROUP", "Residual"),
                                                       df=c(1, 1)))



fit26 <- get_BLUP(data = pheno, model = NODES ~ (1 | GROUP) + (1 | Longitude) 
                 + (1 | Latitude) + (1 | GROUP:Longitude) + (1 | GROUP:Latitude), which.factor = "GROUP",
                 outfile = "data/blup_evolution_NODES_v2.csv")

get_H2(fit26, numerator="GROUP", denominator=data.frame(f=c("GROUP", "Residual"),
                                                       df=c(1, 1)))



fit27 <- get_BLUP(data = pheno, model = NODES_w_BR ~ (1 | GROUP) + (1 | Longitude) 
                 + (1 | Latitude) + (1 | GROUP:Longitude) + (1 | GROUP:Latitude), which.factor = "GROUP",
                 outfile = "data/blup_evolution_NODES_w_BR_v2.csv")

get_H2(fit27, numerator="GROUP", denominator=data.frame(f=c("GROUP", "Residual"),
                                                       df=c(1, 1)))



fit28 <- get_BLUP(data = pheno, model = NODES_w_AR ~ (1 | GROUP) + (1 | Longitude) 
                 + (1 | Latitude) + (1 | GROUP:Longitude) + (1 | GROUP:Latitude), which.factor = "GROUP",
                 outfile = "data/blup_evolution_NODES_w_AR_v2.csv")

get_H2(fit28, numerator="GROUP", denominator=data.frame(f=c("GROUP", "Residual"),
                                                       df=c(1, 1)))



fit29 <- get_BLUP(data = pheno, model = NODES_w_BRandAR ~ (1 | GROUP) + (1 | Longitude) 
                 + (1 | Latitude) + (1 | GROUP:Longitude) + (1 | GROUP:Latitude), which.factor = "GROUP",
                 outfile = "data/blup_evolution_NODES_w_BRandAR_v2.csv")

get_H2(fit29, numerator="GROUP", denominator=data.frame(f=c("GROUP", "Residual"),
                                                       df=c(1, 1)))



fit30 <- get_BLUP(data = pheno, model = DAYS_TO_TASSEL ~ (1 | GROUP) + (1 | Longitude) 
                  + (1 | Latitude) + (1 | GROUP:Longitude) + (1 | GROUP:Latitude), which.factor = "GROUP",
                  outfile = "data/blup_evolution_DAYS_TO_TASSEL_v2.csv")

get_H2(fit30, numerator="GROUP", denominator=data.frame(f=c("GROUP", "Residual"),
                                                        df=c(1, 1)))



fit31 <- get_BLUP(data = pheno, model = DAYS_TO_POLLEN ~ (1 | GROUP) + (1 | Longitude) 
                  + (1 | Latitude) + (1 | GROUP:Longitude) + (1 | GROUP:Latitude), which.factor = "GROUP",
                  outfile = "data/blup_evolution_DAYS_TO_POLLEN_v2.csv")

get_H2(fit31, numerator="GROUP", denominator=data.frame(f=c("GROUP", "Residual"),
                                                        df=c(1, 1)))



fit32 <- get_BLUP(data = pheno, model = COBS ~ (1 | GROUP) + (1 | Longitude) 
                  + (1 | Latitude) + (1 | GROUP:Longitude) + (1 | GROUP:Latitude), which.factor = "GROUP",
                  outfile = "data/blup_evolution_COBS_v2.csv")

get_H2(fit32, numerator="GROUP", denominator=data.frame(f=c("GROUP", "Residual"),
                                                        df=c(1, 1)))



fit33 <- get_BLUP(data = pheno, model = SILKS ~ (1 | GROUP) + (1 | Longitude) 
                  + (1 | Latitude) + (1 | GROUP:Longitude) + (1 | GROUP:Latitude), which.factor = "GROUP",
                  outfile = "data/blup_evolution_SILKS_v2.csv")

get_H2(fit33, numerator="GROUP", denominator=data.frame(f=c("GROUP", "Residual"),
                                                        df=c(1, 1)))

