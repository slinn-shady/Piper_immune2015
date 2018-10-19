#Description: Tropical immune data 2015
#Code author: Heather Slinn
#email: heather.slinn@gmail.com

#version details for HS
#vs with GLM code dates back to July 2017
#vs with Ecuador SEMs dates to Aug 27, 2017
#If I want to return to the analyses I did investigating Tara's data look at version Sep 27, 2017
#If I want GLM analyses and Z-scores go to sep 29,2017
#for old splitting databases coce- Oct10
#Nov 16 for bootstrapping
#Nov 28 for data manipulation
#logistic regression figures in parasitism code
#dec 15, 2017 for testing assumptions of linear moodels

## Table of Contents ##
## 1.0 Upload csv file and load libraries ##
## 2.0 Descriptive stats ##
## 3.0 SEM ##
## 4.0 Figures ##

#####################################################################################
## 1.0 Upload csv file and load libraries ##

eois <- read.table("Dropbox/UNR/Research_Design/Immune data/immune data files/eois_Dec15.txt",stringsAsFactors = FALSE,strip.white = TRUE,na.strings = c("NA",""),check.names=FALSE)
str(eois)

quadrus <- read.table("Dropbox/UNR/Research_Design/Immune data/immune data files/quadrus_Dec15.txt",stringsAsFactors = FALSE,strip.white = TRUE,na.strings = c("NA",""),check.names=FALSE)

ecuador <- read.table("Dropbox/UNR/Research_Design/Immune data/immune data files/ecuador_Dec15.txt",stringsAsFactors = FALSE,strip.white = TRUE,na.strings = c("NA",""),check.names=FALSE)

post_cr <- read.csv("Dropbox/OUTMIXEDCR.csv",stringsAsFactors = FALSE,strip.white = TRUE,na.strings = c("NA",""),check.names=FALSE)
  
post_ec <- read.csv("Dropbox/OUTMIXEDEC.csv",stringsAsFactors = FALSE,strip.white = TRUE,na.strings = c("NA",""),check.names=FALSE)

#load libraries
library(visreg) #data visualization
library(piecewiseSEM) #SEM - I prefer the syntax here and it allows you to use lavaan
library(lavaan) #SEM
library(psych) #describe command

#######################################################################################
## Descriptive stats

# mean parasitism and immune response for Eois in CR
describe(eois) # parasit = 0.12, rate.min = 0.03
# mean parasitism for Eois in EC
describe(ecuador) # parasit = 0.04, rate.min = 0.02
# mean parasitism and immune for Quadrus in CR
describe(quadrus) # parasit = 0.34, rate.min = 0.02

lm <- lm(host.parasitism ~ phytchemdiv, data=quadrus)
summary(lm)
visreg(lm)

km2 <- lm(rate.min ~ phytchemdiv, data=quadrus)
summary(km2)
visreg(km2)

km3 <- lm(host.parasitism ~ rate.min, data=quadrus)
summary(km3)
visreg(km3)

#try removing the 4 species that aren't well represented
# six black two pink on lancifolium
ecuadornew <- ecuador[-c(1, 24, 27:28), ]
# two black spots on lancifolium

#pink spots funk on lancifolium

#eight black blur on baezanum


#####################################################################################
## 3.0  Structural equation models ##
######################################################################################
## Structural equation models----
#Hypothesis 1: Original model
sem.rate.minchem <- lm(rate.min ~ phytchemdiv, eois)
sem.paraschem <- lm(host.parasitism ~ dietbreadth, eois)
sem.dietchem <- lm(dietbreadth ~ phytchemdiv, eois)

SEM.list <- list(sem.rate.minchem, sem.paraschem, 
                 sem.dietchem)
#we decided to stick to global estimation methods (as opposed to local estimation)

#lavaan- Global estimation
(lavaan.model = sem.lavaan(SEM.list, eois, estimator="MLMVS"))
varTable(lavaan.model) #diet breadth variance is much larger than the rest
lavaan::summary(lavaan.model, standardized = T, rsq = T)
lavaan::resid(lavaan.model)
lavaan::modindices(lavaan.model)

######################################################################################
#Ecuador - we are not including Ecuador SEMs in paper - Aug 28, 2017
ECsem.rate.minchem <- lm(rate.min ~ phytchemdiv, ecuadornew)
ECsem.paraschem <- lm(host.parasitism ~ dietbreadth, ecuadornew)
ECsem.dietchem <- lm(dietbreadth ~ phytchemdiv, ecuadornew)

ECSEM.list <- list(ECsem.rate.minchem, ECsem.paraschem, 
                   ECsem.dietchem)
(EClavaan.model = sem.lavaan(ECSEM.list, ecuadornew, estimator="MLMVS")) 
lavaan::summary(EClavaan.model, standardized = T, rsq = T)
######################################################################################
#partial correlation using the "car" package
model <- lm(host.parasitism ~ phytchemdiv + rate.min + dietbreadth, data = eois)
avPlots(model, grid=F, frame=F)

model2 <- lm(rate.min ~ phytchemdiv + host.parasitism + dietbreadth, data = eois)
avPlots(model2, grid=F, frame=F)

######################################################################################
#Ecuador
ECmodel <- lm(host.parasitism ~ phytchemdiv + rate.min + dietbreadth, data = ecuadornew)
avPlots(ECmodel, grid=F, frame=F)

ECmodel2 <- lm(rate.min ~ phytchemdiv + host.parasitism + dietbreadth, data = ecuadornew)
avPlots(ECmodel2, grid=F, frame=F)
######################################################################################
#Hypothesis 2: Phytochemical diversity regulation hypothesis
sem.rate.minchem2 <- lm(rate.min ~ phytchemdiv, eois)
sem.paraschem2 <- lm(host.parasitism ~ phytchemdiv, eois)
sem.dietchem2 <- lm(dietbreadth ~ phytchemdiv, eois)
sem.rate.min.paras <- lm(host.parasitism ~ rate.min, eois)
sem.diet.parasitized <- lm(host.parasitism ~ dietbreadth, eois)

SEM.list2 <- list(sem.rate.minchem2, sem.paraschem2, 
                  sem.dietchem2, sem.rate.min.paras, sem.diet.parasitized)

(lavaan.model2 = sem.lavaan(SEM.list2, eois, estimator="MLMVS"))

varTable(lavaan.model2) #diet breadth variance is much larger than the rest
lavaan::summary(lavaan.model2, standardized=T, rsq=T, fit.measures=T)
lavaan::resid(lavaan.model2, type="standardized")
lavaan::modindices(lavaan.model2)
######################################################################################
#hypothesis with updated ecuador data
ECsem.rate.minchem2 <- lm(rate.min ~ phytchemdiv, ecuadornew)
ECsem.paraschem2 <- lm(host.parasitism ~ phytchemdiv, ecuadornew)
ECsem.dietchem2 <- lm(dietbreadth ~ phytchemdiv, ecuadornew)
ECsem.rate.min.paras <- lm(host.parasitism ~ rate.min, ecuadornew)
ECsem.diet.parasitized <- lm(host.parasitism ~ dietbreadth, ecuadornew)

ECSEM.list2 <- list(ECsem.rate.minchem2, ECsem.paraschem2, 
                    ECsem.dietchem2, ECsem.rate.min.paras, ECsem.diet.parasitized)

(EClavaan.model2 = sem.lavaan(ECSEM.list2, ecuadornew, estimator="MLMVS")) 
lavaan::summary(EClavaan.model2, standardized = T) 

hist(filterEC2$phytchemdiv)
######################################################################################
#Hypothesis 3: Simple phytochemical mediation model
sem.rate.minchem3 <- lm(rate.min ~ phytchemdiv, eois)
sem.rate.min.paras2 <- lm(host.parasitism ~ rate.min, eois)

SEM.list3 <- list(sem.rate.minchem3, sem.rate.min.paras2)

(lavaan.model3 = sem.lavaan(SEM.list3, eois, estimator="MLMVS"))
varTable(lavaan.model3) #diet breadth variance is much larger than the rest
lavaan::summary(lavaan.model3, standardized=T, rsq=T, fit.measures=T) 
lavaan::resid(lavaan.model3)
lavaan::modindices(lavaan.model3)
######################################################################################
#Ecuador
ECsem.rate.minchem3 <- lm(rate.min ~ phytchemdiv, ecuadornew)
ECsem.rate.min.paras2 <- lm(host.parasitism ~ rate.min, ecuadornew)

ECSEM.list3 <- list(ECsem.rate.minchem3, ECsem.rate.min.paras2)

(EClavaan.model3 = sem.lavaan(ECSEM.list3, ecuadornew, estimator="MLMVS")) 
######################################################################################
#Hypothesis 4: Diet breadth mediation hypothesis

sem.rate.mindiet <- lm(rate.min ~ dietbreadth, eois)
sem.rate.min.paras2 <- lm(host.parasitism ~ rate.min, eois)

SEM.list4 <- list(sem.rate.mindiet, sem.rate.min.paras2)

(lavaan.model4 = sem.lavaan(SEM.list4, eois, estimator="MLMVS"))
varTable(lavaan.model4) #diet breadth variance is much larger than the rest
lavaan::summary(lavaan.model4, standardized = T, fit.measures=T) 
lavaan::resid(lavaan.model4)
lavaan::modindices(lavaan.model4)

######################################################################################
##hypothesis for Ecuador data

ECsem.rate.mindiet <- lm(rate.min ~ dietbreadth, ecuadornew)
ECsem.rate.min.paras2 <- lm(host.parasitism ~ rate.min, ecuadornew)

ECSEM.list4 <- list(ECsem.rate.mindiet, ECsem.rate.min.paras2)

(EClavaan.model4 = sem.lavaan(ECSEM.list4, ecuadornew, estimator="MLMVS")) 
varTable(EClavaan.model4) 
lavaan::summary(EClavaan.model4, standardized = T) 
lavaan::resid(EClavaan.model4)
lavaan::modindices(EClavaan.model4)
######################################################################################

#Hypothesis 5: Combination model (a)
sem.rate.mindiet <- lm(rate.min ~ dietbreadth, eois)
sem.rate.minchem3 <- lm(rate.min ~ phytchemdiv, eois)
sem.rate.min.paras2 <- lm(host.parasitism ~ rate.min, eois)

SEM.list5 <- list(sem.rate.mindiet, sem.rate.minchem3, sem.rate.min.paras2)

(lavaan.model5 = sem.lavaan(SEM.list5, eois, estimator="MLMVS"))
varTable(lavaan.model5) #diet breadth variance is much larger than the rest
lavaan::summary(lavaan.model5) 
lavaan::resid(lavaan.model5)
lavaan::modindices(lavaan.model5)

######################################################################################
#Ecuador
ECsem.rate.mindiet <- lm(rate.min ~ dietbreadth, ecuadornew)
ECsem.rate.minchem3 <- lm(rate.min ~ phytchemdiv, ecuadornew)
ECsem.rate.min.paras2 <- lm(host.parasitism ~ rate.min, ecuadornew)

ECSEM.list5 <- list(ECsem.rate.mindiet, ECsem.rate.minchem3, ECsem.rate.min.paras2)

(EClavaan.model5 = sem.lavaan(ECSEM.list5, ecuadornew, estimator="MLMVS")) 

######################################################################################
#Hypothesis 6: Combination model (b)
sem.rate.mindiet <- lm(rate.min ~ dietbreadth, eois)
sem.rate.minchem3 <- lm(rate.min ~ phytchemdiv, eois)
sem.rate.min.paras2 <- lm(host.parasitism ~ rate.min, eois)
sem.dietchem2 <- lm(dietbreadth ~ phytchemdiv, eois)

SEM.list6 <- list(sem.rate.mindiet, sem.rate.minchem3, sem.rate.min.paras2, sem.dietchem2)

(lavaan.model6 = sem.lavaan(SEM.list6, eois, estimator="MLMVS"))
varTable(lavaan.model6) #diet breadth variance is much larger than the rest
lavaan::summary(lavaan.model6, standardized=T) 
lavaan::resid(lavaan.model6)
lavaan::modindices(lavaan.model6)

######################################################################################
#Ecuador
ECsem.rate.mindiet <- lm(rate.min ~ dietbreadth, ecuadornew)
ECsem.rate.minchem3 <- lm(rate.min ~ phytchemdiv, ecuadornew)
ECsem.rate.min.paras2 <- lm(host.parasitism ~ rate.min, ecuadornew)
ECsem.dietchem2 <- lm(dietbreadth ~ phytchemdiv, ecuadornew)

ECSEM.list6 <- list(ECsem.rate.mindiet, ECsem.rate.minchem3, ECsem.rate.min.paras2, ECsem.dietchem2)

(EClavaan.model6 = sem.lavaan(ECSEM.list6, ecuadornew, estimator="MLMVS")) 
######################################################################################

#Hypothesis 7: Combination model (c)
sem.rate.minchem2 <- lm(rate.min ~ phytchemdiv, eois)
sem.paraschem2 <- lm(host.parasitism ~ phytchemdiv, eois)
sem.dietchem2 <- lm(dietbreadth ~ phytchemdiv, eois)
sem.rate.min.paras <- lm(host.parasitism ~ rate.min, eois)
sem.diet.rate.min <- lm(rate.min ~ dietbreadth, eois)

SEM.list7 <- list(sem.rate.minchem2, sem.paraschem2, 
                  sem.dietchem2, sem.rate.min.paras, sem.diet.rate.min)

(lavaan.model7 = sem.lavaan(SEM.list7, eois, estimator="MLMVS"))

varTable(lavaan.model7) #diet breadth variance is much larger than the rest
lavaan::summary(lavaan.model7, standardized=T, rsq=T)
lavaan::resid(lavaan.model7, type="standardized")
lavaan::modindices(lavaan.model7)
######################################################################################
#Ecudaor
ECsem.rate.minchem2 <- lm(rate.min ~ phytchemdiv, ecuadornew)
ECsem.paraschem2 <- lm(host.parasitism ~ phytchemdiv, ecuadornew)
ECsem.dietchem2 <- lm(dietbreadth ~ phytchemdiv, ecuadornew)
ECsem.rate.min.paras <- lm(host.parasitism ~ rate.min, ecuadornew)
ECsem.diet.rate.min <- lm(rate.min ~ dietbreadth, ecuadornew)

ECSEM.list7 <- list(ECsem.rate.minchem2, ECsem.paraschem2, 
                    ECsem.dietchem2, ECsem.rate.min.paras, ECsem.diet.rate.min)

(EClavaan.model7 = sem.lavaan(ECSEM.list7, ecuadornew, estimator="MLMVS")) 
lavaan::summary(EClavaan.model7, standardized=T, rsq=T)

######################################################################################
##Quadrus
## Simple phytochemical mediation model
#Can't include dietbreadth here because only 1 species of Quadrus
sem.immchemQU <- lm(rate.min ~ phytchemdiv, quadrus)
sem.paraschemQU <- lm(host.parasitism ~ rate.min, quadrus)

SEM.listQU <- list(sem.immchemQU, sem.paraschemQU)

#global estimation
(lavaan.modelQU = sem.lavaan(SEM.listQU, quadrus, estimator="MLMVS")) 
varTable(lavaan.modelQU) 
lavaan::summary(lavaan.modelQU, standardized= T, rsq=T) 
lavaan::resid(lavaan.modelQU)
lavaan::modindices(lavaan.modelQU)

semPaths(lavaan.modelQU, what="std")
######################################################################################
## 4.0 Figures ##
######################################################################################
#quadplot for Eois in CR
par(mfrow=c(2,2))
plot(eois$dietbreadth, eois$rate.min, xlab = "Diet Breadth (number of host species)", ylab = "Total PO (abs/min)", cex.lab=1.4, cex.axis=1, pch=21, cex=2, col="black", bg="#ffa696", lwd=1) 
pfit <- lm(rate.min ~ dietbreadth, data = eois)
summary(pfit)

plot(eois$rate.min, eois$host.parasitism, xlab = "Total PO (abs/min)", ylab = "Herbivore Parasitism", cex.lab=1.4, cex.axis=1, pch=21, cex=2, col="black", bg="#16a725", lwd=1) 
pfit2 <- lm(host.parasitism ~ rate.min, data = eois)
abline(pfit2)
summary(pfit2)

plot<- plot(eois$phytchemdiv, eois$host.parasitism, xlab = "Phytochemical Diversity (NMR bin diversity)", ylab = "Herbivore Parasitism", cex.lab=1.4, cex.axis=1, pch=21, cex=2, col="black", bg="#6a87ff", lwd=1) 
pfit3 <- lm(host.parasitism ~ phytchemdiv, data = eois)
abline(pfit3)
summary(pfit3)

plot(eois$phytchemdiv, eois$rate.min, xlab = "Phytochemical Diversity (NMR bin diversity)", ylab = "Total PO (abs/min)", cex.lab=1.4, cex.axis=1, pch=21, cex=2, col="black", bg="#a0007d", lwd=1)
lm <- lm(rate.min ~ phytchemdiv, data=eois)
summary(lm)
abline(lm)


############################################################################################################
#quadplot for Eois in EC
par(mfrow=c(2,2))

plot(ecuador$dietbreadth, ecuador$rate.min, xlab = "Diet Breadth (number of host species)", ylab = "Total PO (abs/min)", cex.lab=1.4, cex.axis=1, pch=21, cex=2, col="black", bg="#ffa696", lwd=1) 
pfit <- lm(rate.min ~ dietbreadth, data = ecuador)
summary(pfit)

plot(ecuador$rate.min, ecuador$host.parasitism, xlab = "Total PO (abs/min)", ylab = "Herbivore Parasitism", cex.lab=1.4, cex.axis=1, pch=21, cex=2, col="black", bg="#16a725", lwd=1) 
pfit2 <- lm(host.parasitism ~ rate.min, data = ecuador)
summary(pfit2)

plot<- plot(ecuador$phytchemdiv, ecuador$host.parasitism, xlab = "Phytochemical Diversity (NMR bin diversity)", ylab = "Herbivore Parasitism", cex.lab=1.4, cex.axis=1, pch=21, cex=2, col="black", bg="#6a87ff", lwd=1) 
pfit3 <- lm(host.parasitism ~ phytchemdiv, data = ecuador)
abline(pfit3)
summary(pfit3)

plot(ecuador$phytchemdiv, ecuador$rate.min, xlab = "Phytochemical Diversity (NMR bin diversity)", ylab = "Total PO (abs/min)", cex.lab=1.4, cex.axis=1, pch=21, cex=2, col="black", bg="#a0007d", lwd=1)
lm <- lm(rate.min ~ phytchemdiv, data=ecuador)
summary(lm)
abline(lm)

dev.off()

#to use covariates, they are supposed to be uncorrelated with treatment variable, which here is phytchemidiv and diet breadth and immune response
#correlated with phytchemdiv and immune response but not diet breadth
plot(eois$cov, eois$phytchemdiv)
lm <- lm(phytchemdiv ~cov, data=eois)
abline(lm)
summary(lm)

plot(eois$cov, eois$dietbreadth)
lm <- lm(dietbreadth ~cov, data=eois)
abline(lm)
summary(lm)

##Visualizing phytchem div by total PO slope by new species in EC
pec <- ggplot(ecuador, aes(x=phytchemdiv, y=rate.min))
pec <- qplot(x = phytchemdiv, y = rate.min, color = caterpillar.sp, data = ecuador, geom = "point", xlab = "Phytochemical Diversity", ylab = "Total PO")
pec 
##Visualizing phytchem div by total PO slope by quadrus species
pq <- ggplot(quadrus, aes(x=phytchemdiv, y=rate.min))
pq <- qplot(x = phytchemdiv, y = rate.min, color = caterpillar.sp, data = quadrus, geom = "point", xlab = "Phytochemical Diversity", ylab = "Total PO")
pq 

#############################################################################################################
#posterior distributions
site <- c("Costa Rica")
site2 <- c("Ecuador")
B1 <- post_cr$B1
EC1 <- post_ec$B1
costarica<- data.frame(site, B1, row.names = NULL, check.rows = F, check.names = T)
costarica

ecuador <- data.frame(site2, EC1, row.names = NULL, check.rows = F, check.names = T)
ecuador
colnames(ecuador)<-c("site", "B1")
new <- rbind(costarica, ecuador)

# create value labels 
site.f <- factor(new$site, levels= c("Costa Rica", "Ecuador"),
                labels = c("Costa Rica", "Ecuador")) 
site.f


library(ggplot2)
library(lattice)

densityplot(~ B1, group = site, data = new, auto.key = TRUE)
ggplot(new) + geom_density(aes(x=B1, fill=site), alpha = 0.2) + ylab("Posterior Probability Density") + theme_classic() + theme(legend.position= c(0.8, 0.8)) + scale_fill_discrete(name="Site") + theme(legend.title = element_text(size=20, face="bold")) + theme(legend.text = element_text(size=20)) + theme(axis.title = element_text(size=20)) + theme(axis.text = element_text(size=16)) + scale_fill_manual(values=c("darkblue", "skyblue2"))

