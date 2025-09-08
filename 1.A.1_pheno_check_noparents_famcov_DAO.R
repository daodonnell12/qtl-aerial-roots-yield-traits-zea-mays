### Jinliang Yang
### Feb 23rd, 2017
### data checking
### Additional customization David O'Donnell


## Install some required packages if have not installed

install.packages("lme4")
install.packages(devtools)
#devtools::install_github("jyanglab/g3tools")

library("lme4")
library("g3tools")




pheno <- read.csv("data/1_QTL_AllFamilies_NoParents.csv")
dim(pheno)
#687  33

n2 <- subset(pheno, !is.na(X15NT1))

hist(n2$X15NT1, breaks=30)
hist(n2$X15NT1R, breaks=30)


fit <- get_BLUP(data = pheno, model = X15NT1 ~ (1 | Genotype) + (1 | Rep) + (1 | Longitude) 
         + (1 | Latitude) + (1 | TFamily), which.factor = "Genotype",
         outfile = "data/blup_noparents_famcov_15NT1.csv")
get_H2(fit, numerator="Genotype", denominator=data.frame(f=c("Genotype", "Residual"),
                              df=c(1, 1)))







fit2 <- get_BLUP(data = pheno, model = X15NT2 ~ (1 | Genotype) + (1 | Rep) + (1 | Longitude) 
                + (1 | Latitude) + (1 | TFamily), which.factor = "Genotype",
                outfile = "data/blup_noparents_famcov_15NT2.csv")

get_H2(fit2, numerator="Genotype",
       denominator=data.frame(f=c("Genotype", "Residual"),
                              df=c(1, 1)))


fit3 <- get_BLUP(data = pheno, model = X15NT3 ~ (1 | Genotype) + (1 | Rep) + (1 | Longitude) 
                 + (1 | Latitude) + (1 | TFamily), which.factor = "Genotype",
                 outfile = "data/blup_noparents_famcov_15NT3.csv")

get_H2(fit3, numerator="Genotype",
       denominator=data.frame(f=c("Genotype", "Residual"),
                              df=c(1, 1)))

fit4 <- get_BLUP(data = pheno, model = DT ~ (1 | Genotype) + (1 | Rep) + (1 | Longitude) 
                + (1 | Latitude) + (1 | TFamily), which.factor = "Genotype",
                outfile = "data/blup_noparents_famcov_DT.csv")

get_H2(fit4, numerator="Genotype",
       denominator=data.frame(f=c("Genotype", "Residual"),
                              df=c(1, 1)))

fit5 <- get_BLUP(data = pheno, model = DPS ~ (1 | Genotype) + (1 | Rep) + (1 | Longitude) 
                 + (1 | Latitude) + (1 | TFamily), which.factor = "Genotype",
                 outfile = "data/blup_noparents_famcov_DPS.csv")

get_H2(fit5, numerator="Genotype",
       denominator=data.frame(f=c("Genotype", "Residual"),
                              df=c(1, 1)))

fit6 <- get_BLUP(data = pheno, model = AR ~ (1 | Genotype) + (1 | Rep) + (1 | Longitude) 
                 + (1 | Latitude) + (1 | TFamily), which.factor = "Genotype",
                 outfile = "data/blup_noparents_famcov_AR.csv")

get_H2(fit6, numerator="Genotype",
       denominator=data.frame(f=c("Genotype", "Residual"),
                              df=c(1, 1)))

fit7 <- get_BLUP(data = pheno, model = SC ~ (1 | Genotype) + (1 | Rep) + (1 | Longitude) 
                 + (1 | Latitude) + (1 | TFamily), which.factor = "Genotype",
                 outfile = "data/blup_noparents_famcov_SC.csv")

get_H2(fit7, numerator="Genotype",
       denominator=data.frame(f=c("Genotype", "Residual"),
                              df=c(1, 1)))

fit8 <- get_BLUP(data = pheno, model = PDM ~ (1 | Genotype) + (1 | Rep) + (1 | Longitude) 
                 + (1 | Latitude) + (1 | TFamily), which.factor = "Genotype",
                 outfile = "data/blup_noparents_famcov_PDM.csv")

get_H2(fit8, numerator="Genotype",
       denominator=data.frame(f=c("Genotype", "Residual"),
                              df=c(1, 1)))

fit9 <- get_BLUP(data = pheno, model = PTN ~ (1 | Genotype) + (1 | Rep) + (1 | Longitude) 
                 + (1 | Latitude) + (1 | TFamily), which.factor = "Genotype",
                 outfile = "data/blup_noparents_famcov_PTN.csv")

get_H2(fit9, numerator="Genotype",
       denominator=data.frame(f=c("Genotype", "Residual"),
                              df=c(1, 1)))

fit10 <- get_BLUP(data = pheno, model = PTN_percent ~ (1 | Genotype) + (1 | Rep) + (1 | Longitude) 
                 + (1 | Latitude) + (1 | TFamily), which.factor = "Genotype",
                 outfile = "data/blup_noparents_famcov_PTN_percent.csv")

get_H2(fit10, numerator="Genotype",
       denominator=data.frame(f=c("Genotype", "Residual"),
                              df=c(1, 1)))

fit11 <- get_BLUP(data = pheno, model = GDM ~ (1 | Genotype) + (1 | Rep) + (1 | Longitude) 
                 + (1 | Latitude) + (1 | TFamily), which.factor = "Genotype",
                 outfile = "data/blup_noparents_famcov_GDM.csv")

get_H2(fit11, numerator="Genotype",
       denominator=data.frame(f=c("Genotype", "Residual"),
                              df=c(1, 1)))

fit12 <- get_BLUP(data = pheno, model = GTN ~ (1 | Genotype) + (1 | Rep) + (1 | Longitude) 
                 + (1 | Latitude) + (1 | TFamily), which.factor = "Genotype",
                 outfile = "data/blup_noparents_famcov_GTN.csv")

get_H2(fit12, numerator="Genotype",
       denominator=data.frame(f=c("Genotype", "Residual"),
                              df=c(1, 1)))

fit13 <- get_BLUP(data = pheno, model = GTN_percent ~ (1 | Genotype) + (1 | Rep) + (1 | Longitude) 
                 + (1 | Latitude) + (1 | TFamily), which.factor = "Genotype",
                 outfile = "data/blup_noparents_famcov_GTN_percent.csv")

get_H2(fit13, numerator="Genotype",
       denominator=data.frame(f=c("Genotype", "Residual"),
                              df=c(1, 1)))
