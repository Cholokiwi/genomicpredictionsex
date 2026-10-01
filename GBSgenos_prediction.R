###################################################
####### Breed prediction in KGD       #############
# For reference-based GBS data in sockeye salmon  #
# v1.0  Sep2026                                   #            
# author: Roy Costilla                            #
###################################################

args = commandArgs(trailingOnly=TRUE)
genofile = args[1]
gform = args[2]
myout = args[3]

# libraries
library(RColorBrewer)
require(lattice)
require(dplyr)
require(ggplot2)

#genofile="samples.2025.vcf.gz"; gform="VCF";myout="samples"  # input data
source("~/KGD/GBS-Chip-Gmatrix.R")
source("~/KGD/GBS-PopGen.R")

readGBS()
summary(sampdepth)
summary(callrate)
summary(snpdepth)

# Validation data is already QCed
sampdepth.thresh <- 0 #0.25
snpdepth.thresh  <- 0 #2
cex.pointsize <- 1.2
functions.only <- TRUE
maf.thresh <- 0 #0.01

GBSsummary()
# objects with info per sample
summary(sampdepth)
summary(callrate)
summary(snpdepth)

# Reading Sex SNPs identified using 2024 cohort
sexsnps=read.table("../../sex_snps_sockeye_2024.txt")[,1]
str(sexsnps)
sexsnps
sum(sexsnps %in% SNP_Names)

# Saving genos for sex snps
mydata=mymeta[mymeta$seqID %in% seqID,]
str(mydata)
mydata$sex_num=ifelse(mydata$sex=="Female",1,0)
with(mydata,table(sex,sex_num))
sexsnps_index=which(SNP_Names %in% sexsnps)
genossex=genon[, sexsnps_index]
genossex = data.frame(cbind(seqID=seqID,female=mydata$sex_num, genossex))
head(genossex)
str(genossex)

#########    within-sample breed prediction
# Gender prediction using KGD function genderpred
uM <- which(mydata$sex[match(seqID,mydata$seqID)]=="Male")
uF <- which(mydata$sex[match(seqID,mydata$seqID)]=="Female")
afF <- calcp(uF)
afM <- calcp(uM)
str(cbind(afF,afM))
head(cbind(afF,afM))
plot(afF,afM)
cbind(afF[sexsnps_index],afM[sexsnps_index])
plot(afF[sexsnps_index],afM[sexsnps_index])

genderaf <- as.matrix(cbind(afF,afM)); colnames(genderaf) <- c("Female","Male")

###### breed prediction method
genderpred <- as.data.frame(breedpredict(snpsubset=sexsnps_index,breedaf=genderaf))
head(genderpred)

# logistic transformation for linear probabilities
pF=genderpred[,"pFemale"]
pM=genderpred[,"pMale"]
pF=exp(pF)/(1+(exp(pF)))
pM=exp(pM)/(1+(exp(pM)))
genderpred$pFemale = pF
genderpred$pMale = pM

# Classification
genderpred=genderpred %>% mutate( SexP=ifelse(pFemale>mean(pF),"Female","Male") )

table(genderpred$SexP)
summary(genderpred[,"pFemale"])

alldata=data.frame(seqID,sex=mydata$sex[match(seqID,mydata$seqID)],genderpred,location=mydata$location[match(seqID,mydata$seqID)])
str(alldata)
alldata$sex=factor(alldata$sex, levels = c("Male","Female"))
alldata$SexP=factor(alldata$SexP, levels = c("Male","Female"))
with(alldata,table(sex,SexP)) 
write.csv(alldata,"sex_prediction_2025.csv",quote=FALSE) 
