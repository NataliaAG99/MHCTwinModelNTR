#########################
#
# MHC heritability study
# 
# Data cleaning, preparation and multivariate twin analysis
#
# 28/07/2025; R version 4.3.2
#
#########################

# Preparing the environment
## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
rm(list=ls(all=TRUE)) #Delete every previous object in the environment


#### Install and load packages ####
## First specify the packages of interest
packages = c("OpenMx", "psych", "umx", "haven", "readxl", "dplyr", "ggplot2",
             "lme4", "grid", "reshape2", "here", "moments", "foreign", "gee", "tidyr",
             "ISLR", "polycor", "tidyverse", "broom", "mosaic")

## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)


# Set the directory
#
#
# setwd("C:/Users/vgs055/OneDrive - Vrije Universiteit Amsterdam/2_Projecten/2025_Empathy heritability/2_Data analysis/2.2_Input")
#
# setwd("C:/Users/naa840/OneDrive - Vrije Universiteit Amsterdam/Nueva carpeta/MHC/03. Data/")

setwd("/data/natalia/MHC/")

#########################
# Load the data
dades <- read.spss("PHE_20250701_5055_LL.sav", to.data.frame=T, use.value.labels = F) #use.value.labels=F
describe(dades)
summary(dades)

#Birthorder for 6618 i
# Twinzyg for 6517 i
table(dades$role) #role 1 for twins(multiples)
#role =1 --> 6619 i
table(dades$Extension) #Ext 1 and 2 for twins
# ext =1 --> 3338 i; ext = 2 --> 3191 i; total = 6529i
table(dades$birthorder)
# birthorder = 1 --> 3354; birthorder = 2 --> 3222

twin_df <- dades[!is.na(dades$twzyg),] #6517i

sum(is.na(twin_df$Extension)) # 0 NA in Extension
sum(is.na(twin_df$birthorder)) # 35 NA

table(twin_df$birthorder, twin_df$Extension, useNA = "ifany") #They match so we'll use extension (less NA)

## Create the indices for flo and empathy

# Flourising is the sum up of the items only if they have 1 or none missing values
flo_vars <- c("flo_a14", "flo_b14", "flo_c14", "flo_d14", "flo_e14", "flo_f14", "flo_g14", "flo_h14")
twin_df$flo14 <- ifelse(
  rowSums(is.na(twin_df[flo_vars])) < 2,
  rowSums(twin_df[flo_vars], na.rm = TRUE),
  NA
)
sum(is.na(twin_df$flo14)) # 891 NA in flo14


# MHC is the average of the items
mhc_vars <- c("mhca14", "mhcb14", "mhcc14", "mhcd14", "mhce14", "mhcf14", "mhcg14", "mhch14", "mhci14", "mhcj14", "mhck14", "mhcl14", "mhcm14", "mhcn14")
# # Compute average only when two or less items are missing
twin_df$mhc14 <- ifelse(
  rowSums(is.na(twin_df[mhc_vars])) <= 2,
  rowMeans(twin_df[mhc_vars], na.rm = TRUE),
  NA
)

sum(is.na(twin_df$mhc14)) # 1224 NA in mhc14
subset(twin_df, mhc14 == 0) 


# Descriptives of the wb variables
describe(twin_df[,c("levenc14", "sat14", "hap14", "flo14", "mhc14")])
# Histogram of the wb variables
# cvd1
layo=matrix(1:6,2,3)
layout(layo)
# cvd2
hist(twin_df$levenc14, main="Life satisfaction", xlab="Life satisfaction", col="lightblue")
hist(twin_df$sat14, main="Satisfaction", xlab="Satisfaction", col="lightblue")
hist(twin_df$hap14, main="Happiness", xlab="Happiness", col="lightblue")
hist(twin_df$flo14, main="Flourishing", xlab="Flourishing", col="lightblue")
hist(twin_df$mhc14, main="MHC", xlab="MHC", col="lightblue")
#
# cvd1
hist(twin_df$age14, xlab="age14", col="red")
# cvd2
#
# Scale the variables
twin_df$qol_sc <- scale(twin_df$levenc14)
twin_df$sat_sc <- scale(twin_df$sat14)
twin_df$hap_sc <- scale(twin_df$hap14)
twin_df$flo_sc <- scale(twin_df$flo14)
twin_df$mhc_sc <- scale(twin_df$mhc14)
#
#Re-code Sex
twin_df$sex <- twin_df$sex - 1
#
#Create Age scaled (divided by 100)
twin_df$Age_sc <- (twin_df$age14)/100
#
# Select the relevant columns for the twin data
data_twins <- subset(twin_df, select=-c(role, multiple_type, birthorder,
                                        famzyg, mult_ext, N_mult_fam, mult_id_fam, 
                                        in_AS_14, version14, flo_a14, flo_b14, 
                                        flo_c14, flo_d14, flo_e14, flo_f14, flo_g14, 
                                        flo_h14, levenc14, mhca14, mhcb14, mhcc14, mhcd14, 
                                        mhce14, mhcf14, mhcg14, mhch14, mhci14, mhcj14, 
                                        mhck14, mhcl14, mhcm14, mhcn14, sat14, hap14, flo14, mhc14))

# # Remove the original twin_df to save memory
rm(twin_df) 
rm(dades)
#
#
sum(is.na(data_twins$twzyg))
#
#Remove the multiples
data_twins <- subset(data_twins, !Extension %in% c(3, 4, 5))
#
#Keep the individuals with all data available (complete cases)
data_twins <- data_twins[complete.cases(data_twins[,c('qol_sc', 'sat_sc', 'hap_sc', 'flo_sc', 'mhc_sc')]),]
#
#Keep individuals in the other dataset
fisnumbersMHC <- readRDS(file = "fisnumbersMHC.rds")
data_twins <- data_twins[data_twins$FISNumber %in% fisnumbersMHC, ]
rm(fisnumbersMHC)
#cvd1
#
rph1 = summary(lm(qol_sc~Age_sc + sex, data=data_twins)) # data_twins[,c('qol_sc', 'sat_sc', 'hap_sc', 'flo_sc', 'emp_sc')]))
rph2 = summary(lm(sat_sc~Age_sc + sex, data=data_twins))
rph3 = summary(lm(hap_sc~Age_sc + sex, data=data_twins))
rph4 = summary(lm(flo_sc~Age_sc + sex, data=data_twins))
rph5 = summary(lm(mhc_sc~Age_sc + sex, data=data_twins))
R2s = c(rph1$r.square,rph2$r.square,rph3$r.square,rph4$r.square,rph5$r.square) 
parsb0 = c(
  rph1$coefficients[1,1],
  rph2$coefficients[1,1],
  rph3$coefficients[1,1],
  rph4$coefficients[1,1],
  rph5$coefficients[1,1])
parsbage = c( 
  rph1$coefficients[2,1],
  rph2$coefficients[2,1],
  rph3$coefficients[2,1],
  rph4$coefficients[2,1],
  rph5$coefficients[2,1])
parsbsex = c( 
  rph1$coefficients[3,1],
  rph2$coefficients[3,1],
  rph3$coefficients[3,1],
  rph4$coefficients[3,1],
  rph5$coefficients[3,1])
# cvd2
#
#
length(unique(data_twins$FamilyNumber)) #4195 families, 5227 i

table(as.data.frame(table(data_twins$FamilyNumber))[,2]) #Number of twins per family ## 2 --> 1032 complete pairs

table(data_twins$Extension)

# Renaming the columns for clarity
colnames(data_twins) <- c("FISNumber", "FamilyNumber", "Extension", "Sex", "Zyg", "Age", "QoL", "Sat", "Hap", "Flo", "MHC", "Age_sc")

# how many individuals of each zygosity?
zyg_frequencies <- table(data_twins[,"Zyg"])
zygosities <- c("MZM", "DZM", "MZF", "DZF", "DOSmf", "DOSfm")
rownames(zyg_frequencies) <- zygosities
zyg_frequencies

# MZM   DZM   MZF   DZF DOSmf DOSfm 
# 698   391  2124   921   558   535

describe(data_twins[, c("Age", "Sex", "QoL", "Sat", "Hap", "Flo", "MHC")]) 
#n=6 with a missing sex left

subset(data_twins, is.na(data_twins$Sex)) #We can assign sex to these 6 individuals based on zygosity
data_twins$Sex[is.na(data_twins$Sex) & data_twins$Zyg %in% c(1, 2)] <- 0
data_twins$Sex[is.na(data_twins$Sex) & data_twins$Zyg %in% c(3, 4)] <- 1
data_twins$Sex[is.na(data_twins$Sex) & data_twins$Zyg %in% c(5, 6)] <- NA
#
#Cross-variable correlations
cor.test(data_twins$QoL, data_twins$Sat) #0.7572758 
cor.test(data_twins$QoL, data_twins$Hap) #0.7411566 
cor.test(data_twins$QoL, data_twins$Flo) #0.6902003 
cor.test(data_twins$QoL, data_twins$MHC) #0.5222003 
#
cor.test(data_twins$Sat, data_twins$Hap) #0.7548571 
cor.test(data_twins$Sat, data_twins$Flo) #0.7111875 
cor.test(data_twins$Sat, data_twins$MHC) #0.5024835 
#
cor.test(data_twins$Hap, data_twins$Flo) #0.7151848 
cor.test(data_twins$Hap, data_twins$MHC) #0.5565203 
#
cor.test(data_twins$Flo, data_twins$MHC) #0.6375301 
#
#
##################
#Wide format
#
# Reshape to wide data format 
total_variables <- dim(data_twins)[2] # better to be flexible in case new variables are added (averages, differences, ...)
vars2lineup <- colnames(data_twins)[c(1,3:4,6:total_variables)] # the variables required once for each twin
data_wide <- umx_long2wide(data_twins, famID = "FamilyNumber", twinID = "Extension", 
                           zygosity = "Zyg", vars2keep = vars2lineup)



# Subset by zygosity 
DZ_levels <- c(2,4,5,6)
dataMZ <- subset(data_wide, !(Zyg %in% DZ_levels))
dataDZ <- subset(data_wide, Zyg %in% DZ_levels)


# First check the correlations #2DZ>MZ->ACE
(QoL_mz=cor.test(dataMZ$QoL_T1,dataMZ$QoL_T2)) #0.3758025  
(QoL_dz=cor.test(dataDZ$QoL_T1,dataDZ$QoL_T2)) #0.2651695 

(Sat_mz=cor.test(dataMZ$Sat_T1,dataMZ$Sat_T2)) #0.4042446 
(Sat_dz=cor.test(dataDZ$Sat_T1,dataDZ$Sat_T2)) #0.2127834 

(Hap_mz=cor.test(dataMZ$Hap_T1,dataMZ$Hap_T2)) #0.407134 
(Hap_dz=cor.test(dataDZ$Hap_T1,dataDZ$Hap_T2)) #0.1500732 

(Flo_mz=cor.test(dataMZ$Flo_T1,dataMZ$Flo_T2)) #0.3818246 
(Flo_dz=cor.test(dataDZ$Flo_T1,dataDZ$Flo_T2)) #0.1507913 

(Emp_mz=cor.test(dataMZ$MHC_T1,dataMZ$MHC_T2)) #0.411768 
(Emp_dz=cor.test(dataDZ$MHC_T1,dataDZ$MHC_T2)) #0.1347492 

#
#Cross-twin cross-trait correlations are estimated in the saturated model
#

#Twin models

source("miFunctions.R")
mxOption(NULL, "Default optimizer", "SLSQP")
#
# ----------------------------------------------------------------------------------------------------------------------
# PREPARE DATA

mzData   <- subset(data_wide, (Zyg==3 | Zyg==1))
dzData   <- subset(data_wide, (Zyg==4 | Zyg==2 | Zyg==5 | Zyg==6))

# Select Variables for Analysis
vars      <- c('QoL','Sat', 'Hap', 'Flo', 'MHC')              # list of variables names
nv        <- 5                         # number of variables
ntv       <- nv*2                      # number of total variables
selVars   <- paste(vars,c(rep(1,nv),rep(2,nv)),sep="_T")
defVars   <- c('Age_sc_T1','Sex_T1','Age_sc_T2','Sex_T2')
useVars   <- c(selVars,defVars)
#
# There were no NAs for the covariates (Age and Sex), so as the model doesn't handle NAs, if there is an NA in any of those, it is because it is an uncomplete pair. Therefore, we write 9999
mzData$Age_sc_T1[is.na(mzData$Age_sc_T1)] <- mzData$Age_sc_T2[is.na(mzData$Age_sc_T1)]
mzData$Age_sc_T2[is.na(mzData$Age_sc_T2)] <- mzData$Age_sc_T1[is.na(mzData$Age_sc_T2)]
mzData$Sex_T1[is.na(mzData$Sex_T1)] <- mzData$Sex_T2[is.na(mzData$Sex_T1)]
mzData$Sex_T2[is.na(mzData$Sex_T2)] <- mzData$Sex_T1[is.na(mzData$Sex_T2)]

dzData$Age_sc_T1[is.na(dzData$Age_sc_T1)] <- dzData$Age_sc_T2[is.na(dzData$Age_sc_T1)]
dzData$Age_sc_T2[is.na(dzData$Age_sc_T2)] <- dzData$Age_sc_T1[is.na(dzData$Age_sc_T2)]
dzData$Sex_T1[is.na(dzData$Sex_T1)] <- dzData$Sex_T2[is.na(dzData$Sex_T1)]
dzData$Sex_T2[is.na(dzData$Sex_T2)] <- dzData$Sex_T1[is.na(dzData$Sex_T2)]

# Select Data for Analysis
mzData   <- subset(mzData, (Zyg==3 | Zyg==1), c(useVars))
dzData   <- subset(dzData, (Zyg==4 | Zyg==2 | Zyg==5 | Zyg==6), c(useVars))

# Generate Descriptive Statistics
colMeans(mzData,na.rm=TRUE)
colMeans(dzData,na.rm=TRUE)
cov(mzData,use="complete")
cov(dzData,use="complete")
#
# cvd1 see variances
diag(cov(mzData,use="complete"))
diag(cov(dzData,use="complete"))
# cvd2
#
# Set Starting Values 
svMe      <- c(0)                   # start value for means
svVa      <- c(1)                  # start value for variances
lbVa      <- .0001                     # lower bound for variances
#
# cvd1
svb0 = parsb0
svbsex=parsbsex
svbage=parsbage 
# cvd2
#
#
# Create Labels
labMeMZ   <- labVars("meanMZ",selVars)
labMeDZ   <- labVars("meanDZ",selVars)
#
labMeZ    <- labVars("meanZ",selVars)
labCvMZ   <- labLower("covMZ",ntv)
labCvDZ   <- labLower("covDZ",ntv)
labCvZ    <- labLower("covZ",ntv)
labVaMZ   <- labDiag("covMZ",ntv)
labVaDZ   <- labDiag("covDZ",ntv)
labVaZ    <- labDiag("covZ",ntv)
#
# ------------------------------------------------------------------------------
# PREPARE MODEL
#                                                              Age_sc_T1 Sex_T1 Age_sc_T2 Sex_T2
# cvd1 # parb0 # cvd2
# Saturated Model
defAge    <- mxMatrix( type="Full", nrow=1, ncol=2, free=FALSE, 
                       # cvd1
                       labels=c("data.Age_sc_T1","data.Age_sc_T2"), name="Age" )
# cvd2
# pathB1     <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=.01, 
#               labels=c("bAge_QoL_T1","bAge_Sat_T1", "bAge_Hap_T1","bAge_Flo_T1", "bAge_Emp_T1", 
#                       "bAge_QoL_T2","bAge_Sat_T2", "bAge_Hap_T2","bAge_Flo_T2", "bAge_Emp_T2"), name="b1" )
# cvd1
pathB1     <- mxMatrix( type="Full", nrow=1, ncol=nv, free=TRUE, # TRUE,
                        values=rep(0.05,nv), # c(parsbage),
                        labels=c("bAge_QoL","bAge_Sat", "bAge_Hap","bAge_Flo", "bAge_MHC"), name="b1" )
# cvd2
# cvd1
defSex      <- mxMatrix( type="Full", nrow=1, ncol=2, free=FALSE,
                         labels=c('data.Sex_T1','data.Sex_T2'), name="Sex" )
# cvd2
# pathB2     <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=.2,
#                labels=c("bSex_QoL_T1","bSex_Sat_T1", "bSex_Hap_T1","bSex_Flo_T1", "bSex_Emp_T1",
#                        "bSex_QoL_T2","bSex_Sat_T2", "bSex_Hap_T2","bSex_Flo_T2", "bSex_Emp_T2"), name="b2" ) # Sex effects
#
pathB2     <- mxMatrix( type="Full", nrow=1, ncol=nv, free=TRUE,
                        values=rep(0.05,nv), # c(parsbsex),
                        labels=c("bSex_QoL","bSex_Sat","bSex_Hap","bSex_Flo", "bSex_MHC"), name="b2" ) # Sex effects
#
#
# cvd1 these are parameters b0 : intercepts
# Create Algebra for expected Mean Matrices  
#
# meanMZ    <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=labMeMZ, name="meanMZ" )
# meanDZ    <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=labMeDZ, name="meanDZ" )
# cvd1
b0MZ    <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE,
                     values=c(parsb0,parsb0),
                     labels=c("b0_QoL","b0_Sat", "b0_Hap","b0_Flo", "b0_MHC",
                              "b0_QoL","b0_Sat", "b0_Hap","b0_Flo", "b0_MHC"), name="b0MZ" )
b0DZ    <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE,
                     values=c(parsb0,parsb0),
                     labels=c("b0_QoL","b0_Sat", "b0_Hap","b0_Flo", "b0_MHC",
                              "b0_QoL","b0_Sat", "b0_Hap","b0_Flo", "b0_MHC"), name="b0DZ" )
#
# cvd2
# The order of the data object and the beta object matters
# e.g. for 2 traits the outcome MUST be Twin1-trait1 Twin1-trait2 Twin2-trait2 Twin2-trait2
# Check the matrix algebra
#
expMeanMZ_ <- mxAlgebra( expression = b0MZ + b2%x%Sex + b1%x%Age, name="expMeanMZ")
expMeanDZ_ <- mxAlgebra( expression = b0DZ + b2%x%Sex + b1%x%Age, name="expMeanDZ")
#
#bits   <- list(defSex, defAge, pathB1, pathB2 )
#
## Create Algebra for expected Variance/Covariance Matrices
#covMZ     <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=valDiag(svVa,ntv), lbound=valDiag(lbVa,ntv), labels=labCvMZ, name="covMZ" )
#covDZ     <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=valDiag(svVa,ntv), lbound=valDiag(lbVa,ntv), labels=labCvDZ, name="covDZ" )
#
covMZ_     <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=valDiag(svVa,ntv), labels=labCvMZ, name="covMZ" )
covDZ_     <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=valDiag(svVa,ntv), labels=labCvDZ, name="covDZ" )
#
corMZ  <- mxAlgebra(expression = cov2cor(covMZ),name="corMZ")
corDZ  <- mxAlgebra(expression = cov2cor(covDZ), name="corDZ")
#
# Create Data Objects for Multiple Groups
dataMZ    <- mxData( observed=mzData, type="raw" )
dataDZ    <- mxData( observed=dzData, type="raw" )
#
# Create Expectation Objects for Multiple Groups
# cvd1
expMZ     <- mxExpectationNormal( covariance="covMZ", means="expMeanMZ", dimnames=selVars )
expDZ     <- mxExpectationNormal( covariance="covDZ", means="expMeanDZ", dimnames=selVars )
# cvd2
funML     <- mxFitFunctionML()
#
# Create Model Objects for Multiple Groups
##bits   <- list(defSex, defAge, pathB1, pathB2 )
#
modelMZ   <- mxModel( defSex, defAge, pathB1, pathB2, b0MZ, expMeanMZ_, covMZ_, corMZ, dataMZ, expMZ, funML, name="MZ" )
modelDZ   <- mxModel( defSex, defAge, pathB1, pathB2, b0DZ, expMeanDZ_, covDZ_, corDZ, dataDZ, expDZ, funML, name="DZ" )

multi     <- mxFitFunctionMultigroup( c("MZ","DZ") )
#
# Create Confidence Interval Objects
#
ciCor <- mxCI( c('MZ.corMZ','DZ.corDZ') )
ciCov     <- mxCI( c('MZ.covMZ','DZ.covDZ') )
ciMean    <- mxCI( c('MZ.meanMZ','DZ.meanDZ') )
#
#
# Build Saturated Model with Confidence Intervals
#
modelSAT  <- mxModel("fiveSATc", modelMZ, modelDZ, multi, ciCov, ciCor, ciMean )
#modelSAT  <- mxModel(name="fiveSATc", inclusions, modelMZ, modelDZ, multi) #, ciCov, ciCor, ciMean )
#
# modelSAT  <- mxModel(modelMZ, modelDZ, multi) #, ciCov, ciCor, ciMean )
#
#
# ----------------------------------------------------------------------------------------------------------------------
# RUN MODEL
#mxOption(NULL, "Default optimizer", "NPSOL")
#
# Run Saturated Model
#
fitSAT    <- mxRun(modelSAT)
#fitSAT    <- mxTryHard( modelSAT,10 )
#
#
(sumSAT    <- summary( fitSAT ))
#
#
(corMZ  <- round(cov2cor(fitSAT$MZ$covMZ$values),3))
(corDZ  <- round(cov2cor(fitSAT$DZ$covDZ$values),3))

#--------------------------------------------------------------------------
# Test Significance of Covariate
modelCOV <- mxModel( fitSAT, name="Age" )
modelCOV <- omxSetParameters( modelCOV, label=c("bAge_QoL", "bAge_Sat", "bAge_Hap", "bAge_Flo", "bAge_MHC"), free=FALSE, values=0 )
fitCOV <- mxRun( modelCOV )

summary(fitCOV)
mxCompare(fitSAT,fitCOV) #Significant effect, cannot use the simpler model


modelCOV1 <- mxModel( fitSAT, name="Sex" )
modelCOV1 <- omxSetParameters( modelCOV1, label=c("bSex_QoL", "bSex_Sat", "bSex_Hap", "bSex_Flo", "bSex_MHC"), free=FALSE, values=0 )
fitCOV1 <- mxRun( modelCOV1 ,intervals=F)

mxCompare(fitSAT,fitCOV1) #Significant effect, cannot use the simpler model


summary(fitCOV1)
(corMZ  <- round(cov2cor(fitCOV1$MZ$covMZ$values),3))
(corDZ  <- round(cov2cor(fitCOV1$DZ$covDZ$values),3))

#We should run submodels to see if age or sex do NOT have an effect on any of the wb variables independently

#Age
modelAgeQoL <- mxModel( fitSAT, name="Age" )
modelAgeQoL <- omxSetParameters( modelAgeQoL, label=c("bAge_QoL"), free=FALSE, values=0 )
fitAgeQoL <- mxRun( modelAgeQoL )

summary(fitAgeQoL)
mxCompare(fitSAT,fitAgeQoL) #Significant effect, cannot use the simpler model

modelAgeSat <- mxModel( fitSAT, name="Age" )
modelAgeSat <- omxSetParameters( modelAgeSat, label=c("bAge_Sat"), free=FALSE, values=0 )
fitAgeSat <- mxRun( modelAgeSat )

summary(fitAgeSat)
mxCompare(fitSAT,fitAgeSat) #Significant effect, cannot use the simpler model

modelAgeHap <- mxModel( fitSAT, name="Age" )
modelAgeHap <- omxSetParameters( modelAgeHap, label=c("bAge_Hap"), free=FALSE, values=0 )
fitAgeHap <- mxRun( modelAgeHap )

summary(fitAgeHap)
mxCompare(fitSAT,fitAgeHap) #Significant effect, cannot use the simpler model

modelAgeFlo <- mxModel( fitSAT, name="Age" )
modelAgeFlo <- omxSetParameters( modelAgeFlo, label=c("bAge_Flo"), free=FALSE, values=0 )
fitAgeFlo <- mxRun( modelAgeFlo )

summary(fitAgeFlo)
mxCompare(fitSAT,fitAgeFlo) #Significant effect, cannot use the simpler model

modelAgeMHC <- mxModel( fitSAT, name="Age" )
modelAgeMHC <- omxSetParameters( modelAgeMHC, label=c("bAge_MHC"), free=FALSE, values=0 )
fitAgeMHC <- mxRun( modelAgeMHC )

summary(fitAgeMHC)
mxCompare(fitSAT,fitAgeMHC) #No significant differences (0.067), can use the simpler model

#Sex
modelSexQoL <- mxModel( fitSAT, name="Sex" )
modelSexQoL <- omxSetParameters( modelSexQoL, label=c("bSex_QoL"), free=FALSE, values=0 )
fitSexQoL <- mxRun( modelSexQoL )

summary(fitSexQoL)
mxCompare(fitSAT,fitSexQoL) #Significant effect (0.003), cannot use the simpler model

modelSexSat <- mxModel( fitSAT, name="Sex" )
modelSexSat <- omxSetParameters( modelSexSat, label=c("bSex_Sat"), free=FALSE, values=0 )
fitSexSat <- mxRun( modelSexSat )

summary(fitSexSat)
mxCompare(fitSAT,fitSexSat) #No significant differences(0.2179), can use the simpler model

modelSexHap <- mxModel( fitSAT, name="Sex" )
modelSexHap <- omxSetParameters( modelSexHap, label=c("bSex_Hap"), free=FALSE, values=0 )
fitSexHap <- mxRun( modelSexHap )

summary(fitSexHap)
mxCompare(fitSAT,fitSexHap) #No significant differences(0.563), can use the simpler model

modelSexFlo <- mxModel( fitSAT, name="Sex" )
modelSexFlo <- omxSetParameters( modelSexFlo, label=c("bSex_Flo"), free=FALSE, values=0 )
fitSexFlo <- mxRun( modelSexFlo )

summary(fitSexFlo)
mxCompare(fitSAT,fitSexFlo) #No significant differences (0.603), can use the simpler model

modelSexMHC <- mxModel( fitSAT, name="Sex" )
modelSexMHC <- omxSetParameters( modelSexMHC, label=c("bSex_MHC"), free=FALSE, values=0 )
fitSexMHC <- mxRun( modelSexMHC )

summary(fitSexMHC)
mxCompare(fitSAT,fitSexMHC) #Significant effect, cannot use the simpler model

# ----------------------------------------------------------------------------------------------------------------------
# RUN SUBMODELS

# Constrain expected Means and Variances to be equal across Twin Order
modelEMVO <- mxModel( fitSAT, name="twoEMVOc" )
for (i in 1:nv) {
  modelEMVO <- omxSetParameters( modelEMVO, label=c(labVaMZ[nv+i],labVaMZ[i]), free=TRUE, values=1, newlabels=labVaMZ[i]) 
  modelEMVO <- omxSetParameters( modelEMVO, label=c(labVaDZ[nv+i],labVaDZ[i]), free=TRUE, values=1, newlabels=labVaDZ[i]) }
fitEMVO   <- mxRun( modelEMVO, intervals=F )
fitGofs(fitEMVO); fitEsts(fitEMVO)
mxCompare(fitSAT,fitEMVO) 

# Constrain covariances tw1vstw1 == tw2vstw2
modelcovCons1 <- mxModel( fitEMVO, name="covCons1" )


modelcovCons1 <- omxSetParameters( modelcovCons1, label=c("covMZ21","covMZ76"), free=TRUE, values=fitSAT$MZ$covMZ$values[fitSAT$MZ$covMZ$labels=="covMZ21"][1], newlabels="covMZ21") 
modelcovCons1 <- omxSetParameters( modelcovCons1, label=c("covMZ31","covMZ86"), free=TRUE, values=fitSAT$MZ$covMZ$values[fitSAT$MZ$covMZ$labels=="covMZ31"][1], newlabels="covMZ31") 
modelcovCons1 <- omxSetParameters( modelcovCons1, label=c("covMZ41","covMZ96"), free=TRUE, values=fitSAT$MZ$covMZ$values[fitSAT$MZ$covMZ$labels=="covMZ41"][1], newlabels="covMZ41") 
modelcovCons1 <- omxSetParameters( modelcovCons1, label=c("covMZ51","covMZ106"), free=TRUE, values=fitSAT$MZ$covMZ$values[fitSAT$MZ$covMZ$labels=="covMZ51"][1], newlabels="covMZ51") 
modelcovCons1 <- omxSetParameters( modelcovCons1, label=c("covMZ32","covMZ87"), free=TRUE, values=fitSAT$MZ$covMZ$values[fitSAT$MZ$covMZ$labels=="covMZ32"][1], newlabels="covMZ32") 
modelcovCons1 <- omxSetParameters( modelcovCons1, label=c("covMZ42","covMZ97"), free=TRUE, values=fitSAT$MZ$covMZ$values[fitSAT$MZ$covMZ$labels=="covMZ42"][1], newlabels="covMZ42") 
modelcovCons1 <- omxSetParameters( modelcovCons1, label=c("covMZ52","covMZ107"), free=TRUE, values=fitSAT$MZ$covMZ$values[fitSAT$MZ$covMZ$labels=="covMZ52"][1], newlabels="covMZ52") 
modelcovCons1 <- omxSetParameters( modelcovCons1, label=c("covMZ43","covMZ98"), free=TRUE, values=fitSAT$MZ$covMZ$values[fitSAT$MZ$covMZ$labels=="covMZ43"][1], newlabels="covMZ43") 
modelcovCons1 <- omxSetParameters( modelcovCons1, label=c("covMZ53","covMZ108"), free=TRUE, values=fitSAT$MZ$covMZ$values[fitSAT$MZ$covMZ$labels=="covMZ53"][1], newlabels="covMZ53") 
modelcovCons1 <- omxSetParameters( modelcovCons1, label=c("covMZ54","covMZ109"), free=TRUE, values=fitSAT$MZ$covMZ$values[fitSAT$MZ$covMZ$labels=="covMZ54"][1], newlabels="covMZ54") 

modelcovCons1 <- omxSetParameters( modelcovCons1, label=c("covDZ21","covDZ76"), free=TRUE, values=fitSAT$DZ$covDZ$values[fitSAT$DZ$covDZ$labels=="covDZ21"][1], newlabels="covDZ21") 
modelcovCons1 <- omxSetParameters( modelcovCons1, label=c("covDZ31","covDZ86"), free=TRUE, values=fitSAT$DZ$covDZ$values[fitSAT$DZ$covDZ$labels=="covDZ31"][1], newlabels="covDZ31") 
modelcovCons1 <- omxSetParameters( modelcovCons1, label=c("covDZ41","covDZ96"), free=TRUE, values=fitSAT$DZ$covDZ$values[fitSAT$DZ$covDZ$labels=="covDZ41"][1], newlabels="covDZ41") 
modelcovCons1 <- omxSetParameters( modelcovCons1, label=c("covDZ51","covDZ106"), free=TRUE, values=fitSAT$DZ$covDZ$values[fitSAT$DZ$covDZ$labels=="covDZ51"][1], newlabels="covDZ51") 
modelcovCons1 <- omxSetParameters( modelcovCons1, label=c("covDZ32","covDZ87"), free=TRUE, values=fitSAT$DZ$covDZ$values[fitSAT$DZ$covDZ$labels=="covDZ32"][1], newlabels="covDZ32") 
modelcovCons1 <- omxSetParameters( modelcovCons1, label=c("covDZ42","covDZ97"), free=TRUE, values=fitSAT$DZ$covDZ$values[fitSAT$DZ$covDZ$labels=="covDZ42"][1], newlabels="covDZ42") 
modelcovCons1 <- omxSetParameters( modelcovCons1, label=c("covDZ52","covDZ107"), free=TRUE, values=fitSAT$DZ$covDZ$values[fitSAT$DZ$covDZ$labels=="covDZ52"][1], newlabels="covDZ52") 
modelcovCons1 <- omxSetParameters( modelcovCons1, label=c("covDZ43","covDZ98"), free=TRUE, values=fitSAT$DZ$covDZ$values[fitSAT$DZ$covDZ$labels=="covDZ43"][1], newlabels="covDZ43") 
modelcovCons1 <- omxSetParameters( modelcovCons1, label=c("covDZ53","covDZ108"), free=TRUE, values=fitSAT$DZ$covDZ$values[fitSAT$DZ$covDZ$labels=="covDZ53"][1], newlabels="covDZ53") 
modelcovCons1 <- omxSetParameters( modelcovCons1, label=c("covDZ54","covDZ109"), free=TRUE, values=fitSAT$DZ$covDZ$values[fitSAT$DZ$covDZ$labels=="covDZ54"][1], newlabels="covDZ54") 

fitcovCons1   <- mxRun( modelcovCons1, intervals=F )
mxCompare(fitEMVO,fitcovCons1) 

# Constrain covariances tw1vstw2 == tw2vstw1
modelSATfinal <- mxModel( fitcovCons1, name="SATfinal" )


modelSATfinal <- omxSetParameters( modelSATfinal, label=c("covMZ62","covMZ71"), free=TRUE, values=fitSAT$MZ$covMZ$values[fitSAT$MZ$covMZ$labels=="covMZ62"][1], newlabels="covMZ62") 
modelSATfinal <- omxSetParameters( modelSATfinal, label=c("covMZ63","covMZ81"), free=TRUE, values=fitSAT$MZ$covMZ$values[fitSAT$MZ$covMZ$labels=="covMZ63"][1], newlabels="covMZ63") 
modelSATfinal <- omxSetParameters( modelSATfinal, label=c("covMZ64","covMZ91"), free=TRUE, values=fitSAT$MZ$covMZ$values[fitSAT$MZ$covMZ$labels=="covMZ64"][1], newlabels="covMZ64") 
modelSATfinal <- omxSetParameters( modelSATfinal, label=c("covMZ65","covMZ101"), free=TRUE, values=fitSAT$MZ$covMZ$values[fitSAT$MZ$covMZ$labels=="covMZ65"][1], newlabels="covMZ65") 
modelSATfinal <- omxSetParameters( modelSATfinal, label=c("covMZ73","covMZ82"), free=TRUE, values=fitSAT$MZ$covMZ$values[fitSAT$MZ$covMZ$labels=="covMZ73"][1], newlabels="covMZ73") 
modelSATfinal <- omxSetParameters( modelSATfinal, label=c("covMZ74","covMZ92"), free=TRUE, values=fitSAT$MZ$covMZ$values[fitSAT$MZ$covMZ$labels=="covMZ74"][1], newlabels="covMZ74") 
modelSATfinal <- omxSetParameters( modelSATfinal, label=c("covMZ75","covMZ102"), free=TRUE, values=fitSAT$MZ$covMZ$values[fitSAT$MZ$covMZ$labels=="covMZ75"][1], newlabels="covMZ75") 
modelSATfinal <- omxSetParameters( modelSATfinal, label=c("covMZ84","covMZ93"), free=TRUE, values=fitSAT$MZ$covMZ$values[fitSAT$MZ$covMZ$labels=="covMZ84"][1], newlabels="covMZ84") 
modelSATfinal <- omxSetParameters( modelSATfinal, label=c("covMZ85","covMZ103"), free=TRUE, values=fitSAT$MZ$covMZ$values[fitSAT$MZ$covMZ$labels=="covMZ85"][1], newlabels="covMZ85") 
modelSATfinal <- omxSetParameters( modelSATfinal, label=c("covMZ95","covMZ104"), free=TRUE, values=fitSAT$MZ$covMZ$values[fitSAT$MZ$covMZ$labels=="covMZ95"][1], newlabels="covMZ95") 

modelSATfinal <- omxSetParameters( modelSATfinal, label=c("covDZ62","covDZ71"), free=TRUE, values=fitSAT$DZ$covDZ$values[fitSAT$DZ$covDZ$labels=="covDZ62"][1], newlabels="covDZ62") 
modelSATfinal <- omxSetParameters( modelSATfinal, label=c("covDZ63","covDZ81"), free=TRUE, values=fitSAT$DZ$covDZ$values[fitSAT$DZ$covDZ$labels=="covDZ63"][1], newlabels="covDZ63") 
modelSATfinal <- omxSetParameters( modelSATfinal, label=c("covDZ64","covDZ91"), free=TRUE, values=fitSAT$DZ$covDZ$values[fitSAT$DZ$covDZ$labels=="covDZ64"][1], newlabels="covDZ64") 
modelSATfinal <- omxSetParameters( modelSATfinal, label=c("covDZ65","covDZ101"), free=TRUE, values=fitSAT$DZ$covDZ$values[fitSAT$DZ$covDZ$labels=="covDZ65"][1], newlabels="covDZ65") 
modelSATfinal <- omxSetParameters( modelSATfinal, label=c("covDZ73","covDZ82"), free=TRUE, values=fitSAT$DZ$covDZ$values[fitSAT$DZ$covDZ$labels=="covDZ73"][1], newlabels="covDZ73") 
modelSATfinal <- omxSetParameters( modelSATfinal, label=c("covDZ74","covDZ92"), free=TRUE, values=fitSAT$DZ$covDZ$values[fitSAT$DZ$covDZ$labels=="covDZ74"][1], newlabels="covDZ74") 
modelSATfinal <- omxSetParameters( modelSATfinal, label=c("covDZ75","covDZ102"), free=TRUE, values=fitSAT$DZ$covDZ$values[fitSAT$DZ$covDZ$labels=="covDZ75"][1], newlabels="covDZ75") 
modelSATfinal <- omxSetParameters( modelSATfinal, label=c("covDZ84","covDZ93"), free=TRUE, values=fitSAT$DZ$covDZ$values[fitSAT$DZ$covDZ$labels=="covDZ84"][1], newlabels="covDZ84") 
modelSATfinal <- omxSetParameters( modelSATfinal, label=c("covDZ85","covDZ103"), free=TRUE, values=fitSAT$DZ$covDZ$values[fitSAT$DZ$covDZ$labels=="covDZ85"][1], newlabels="covDZ85") 
modelSATfinal <- omxSetParameters( modelSATfinal, label=c("covDZ95","covDZ104"), free=TRUE, values=fitSAT$DZ$covDZ$values[fitSAT$DZ$covDZ$labels=="covDZ95"][1], newlabels="covDZ95") 

fitSATfinal   <- mxRun( modelSATfinal, intervals=F )
mxCompare(fitcovCons1,fitSATfinal) 

(sumSATfinal    <- summary( fitSATfinal ))
#
#
(corMZfinal  <- round((fitSATfinal$MZ$corMZ$result),2))
(corDZfinal  <- round((fitSATfinal$DZ$corDZ$result),2))

# Constrain expected Means and Variances to be equal across Twin Order and Zygosity
modelEMVZ <- mxModel( fitSATfinal, name="twoEMVZc" )
for (i in 1:nv) {
  modelEMVZ <- omxSetParameters( modelEMVZ, label=c(labVaMZ[i],labVaDZ[i]), free=TRUE, values=svVa, newlabels=labVaZ[i] ) }
fitEMVZ   <- mxRun( modelEMVZ, intervals=F )
fitGofs(fitEMVZ); fitEsts(fitEMVZ)
mxCompare(fitSATfinal,fitEMVZ)

# Constrain covariances across zygosity
modelSAToptimal <- mxModel( fitEMVZ, name="SAToptimal" )


modelSAToptimal <- omxSetParameters( modelSAToptimal, label=c("covMZ21","covDZ21"), free=TRUE, values=fitSAT$MZ$covMZ$values[fitSAT$MZ$covMZ$labels=="covMZ21"][1], newlabels="cov21") 
modelSAToptimal <- omxSetParameters( modelSAToptimal, label=c("covMZ31","covDZ31"), free=TRUE, values=fitSAT$MZ$covMZ$values[fitSAT$MZ$covMZ$labels=="covMZ31"][1], newlabels="cov31") 
modelSAToptimal <- omxSetParameters( modelSAToptimal, label=c("covMZ41","covDZ41"), free=TRUE, values=fitSAT$MZ$covMZ$values[fitSAT$MZ$covMZ$labels=="covMZ41"][1], newlabels="cov41") 
modelSAToptimal <- omxSetParameters( modelSAToptimal, label=c("covMZ51","covDZ51"), free=TRUE, values=fitSAT$MZ$covMZ$values[fitSAT$MZ$covMZ$labels=="covMZ51"][1], newlabels="cov51") 
modelSAToptimal <- omxSetParameters( modelSAToptimal, label=c("covMZ32","covDZ32"), free=TRUE, values=fitSAT$MZ$covMZ$values[fitSAT$MZ$covMZ$labels=="covMZ32"][1], newlabels="cov32") 
modelSAToptimal <- omxSetParameters( modelSAToptimal, label=c("covMZ42","covDZ42"), free=TRUE, values=fitSAT$MZ$covMZ$values[fitSAT$MZ$covMZ$labels=="covMZ42"][1], newlabels="cov42") 
modelSAToptimal <- omxSetParameters( modelSAToptimal, label=c("covMZ52","covDZ52"), free=TRUE, values=fitSAT$MZ$covMZ$values[fitSAT$MZ$covMZ$labels=="covMZ52"][1], newlabels="cov52") 
modelSAToptimal <- omxSetParameters( modelSAToptimal, label=c("covMZ43","covDZ43"), free=TRUE, values=fitSAT$MZ$covMZ$values[fitSAT$MZ$covMZ$labels=="covMZ43"][1], newlabels="cov43") 
modelSAToptimal <- omxSetParameters( modelSAToptimal, label=c("covMZ53","covDZ53"), free=TRUE, values=fitSAT$MZ$covMZ$values[fitSAT$MZ$covMZ$labels=="covMZ53"][1], newlabels="cov53") 
modelSAToptimal <- omxSetParameters( modelSAToptimal, label=c("covMZ54","covDZ54"), free=TRUE, values=fitSAT$MZ$covMZ$values[fitSAT$MZ$covMZ$labels=="covMZ54"][1], newlabels="cov54") 

fitSATOptimal   <- mxRun( modelSAToptimal, intervals=F )

mxCompare(fitEMVZ, fitSATOptimal)

(sumSatOptimal    <- summary( fitSATOptimal ))
#
#
(corMZSAToptimal  <- round((fitSATOptimal$MZ$corMZ$result),2))
(corDZSAToptimal  <- round((fitSATOptimal$DZ$corDZ$result),2))

mxCompare(fitSAT, fitSATOptimal)

save.image("/data/natalia/MHC/MHC_env_Sat_5phe_20260120.RData")
saveRDS(fitSATOptimal, file = "/data/natalia/MHC/Sat_5v_model.rds")

#Pvalue correction for multiple testing (wrong)
p_values_momentary <- c(2.2561969e-11,1.278667e-08,5.1009589e-10,1.9979223e-07,
                        6.249026e-08,7.6066412e-05,0.067421258,0.0035088909,0.21792488,
                        0.56260283, 0.60345615, 1.1434728e-05, 0.42224269, 0.1354094, 
                        0.4140016, 0.04777746, 0.5476674, 0.1302172)
adjusted_p_values_momentary <- p.adjust(p_values_momentary, method = "BH")
print(round(adjusted_p_values_momentary,2))
#With corrected pvalues
# 0.00 0.00 0.00 0.00 0.00 0.00 0.11 0.01 0.31 0.60 0.60 0.00 0.51 0.21 0.51 0.09 0.60
#Still no significant differences:
#AgeMHC, SexSat, SexHap, SexFlo, EqualVariancesTwinOrder, covConstrain1, covConstrain2, EqualVarianceZygosity, EqualCovarianceZygosity


#
dosimcheck=FALSE
# 
if (dosimcheck) { 
  #
  library(MASS)
  nmz=201
  ndz=199
  #
  RA=matrix(.5,5,5)
  RE=matrix(.2,5,5)
  diag(RA) = diag(RE) = 1
  h=diag(sqrt(c(.3,.4,.5,.6,.7)))
  e=diag(sqrt(c(.4,.4,.4,.4,.4)))
  SA=h%*%RA%*%h
  SE=e%*%RE%*%e
  Sph=SA+SE
  Smz = rbind( cbind(Sph,SA), cbind(SA, Sph) ) 
  Sdz = rbind( cbind(Sph,SA/2), cbind(SA/2, Sph) )
  #
  datmz=mvrnorm(nmz, rep(0,10), Sigma=Smz)
  datdz=mvrnorm(ndz, rep(0,10), Sigma=Sdz)
  sexmz=sample(c(1,0), nmz, replace=T) 
  sexdz=sample(c(1,0), ndz, replace=T) 
  agemz = round(rnorm(nmz, 20, 2),2)
  agedz = round(rnorm(ndz, 20, 2),2)
  datmz = cbind(datmz, sexmz, agemz, sexmz, agemz)
  datdz = cbind(datdz, sexdz, agedz, sexdz, agedz)
  colnames(datmz) = colnames(datdz) =c(
    'A1','B1','C1','D1','E1',
    'A2','B2','C2','D2','E2','sex1','age1','sex2','age2')
  datmz = as.data.frame(datmz)
  datdz = as.data.frame(datdz)
  selVars = 
    c('A1','B1','C1','D1','E1',
      'A2','B2','C2','D2','E2')
  #
  #
  nv=5
  ntv=nv*2
  defAge    <- mxMatrix( type="Full", nrow=1, ncol=2, free=FALSE, 
                         labels=c("data.age1","data.age1"), name="Age" )
  #
  pathB1     <- mxMatrix( type="Full", nrow=1, ncol=nv, free=TRUE, # TRUE, 
                          values=rep(0.05,nv), #  
                          labels=c("bAge1","bAge2", "bAge3","bAge4", "bAge5"), name="b1" )
  defSex      <- mxMatrix( type="Full", nrow=1, ncol=2, free=FALSE,
                           labels=c('data.sex1','data.sex1'), name="Sex" )
  pathB2     <- mxMatrix( type="Full", nrow=1, ncol=nv, free=TRUE, 
                          values=rep(0.05,nv), #  
                          labels=c("bSex1","bSex2","bSex3","bSex4", "bSex5"), name="b2" ) # Sex effects
  #
  # cvd1
  b0MZ    <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, 
                       values=c(0), 
                       labels=c("b01","b02", "b03","b04", "b05", 
                                "b01","b02", "b03","b04", "b05"), name="b0MZ" )
  b0DZ    <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, 
                       values=c(0),
                       labels=c("b01","b02", "b03","b04", "b05", 
                                "b01","b02", "b03","b04", "b05"), name="b0DZ" )
  #
  #
  expMeanMZ_ <- mxAlgebra( expression = b0MZ + b2%x%Sex + b1%x%Age, name="expMeanMZ")
  expMeanDZ_ <- mxAlgebra( expression = b0DZ + b2%x%Sex + b1%x%Age, name="expMeanDZ")
  #
  #
  covMZ_     <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=round(Smz,2), labels=labCvMZ, name="covMZ" )
  covDZ_     <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=round(Sdz,2), labels=labCvDZ, name="covDZ" )
  #
  corMZ  <- mxAlgebra(expression = cov2cor(covMZ),name="corMZ")
  corDZ  <- mxAlgebra(expression = cov2cor(covDZ), name="corDZ")
  #
  # Create Data Objects for Multiple Groups
  dataMZ    <- mxData( observed=datmz, type="raw" )
  dataDZ    <- mxData( observed=datdz, type="raw" )
  #
  expMZ     <- mxExpectationNormal( covariance="covMZ", means="expMeanMZ", dimnames=selVars )
  expDZ     <- mxExpectationNormal( covariance="covDZ", means="expMeanDZ", dimnames=selVars )
  # cvd2
  funML     <- mxFitFunctionML()
  #
  # Create Model Objects for Multiple Groups
  ##bits   <- list(defSex, defAge, pathB1, pathB2 )
  #
  modelMZ   <- mxModel( defSex, defAge, pathB1, pathB2, b0MZ, expMeanMZ_, covMZ_, corMZ, dataMZ, expMZ, funML, name="MZ" )
  modelDZ   <- mxModel( defSex, defAge, pathB1, pathB2, b0DZ, expMeanDZ_, covDZ_, corDZ, dataDZ, expDZ, funML, name="DZ" )
  multi     <- mxFitFunctionMultigroup( c("MZ","DZ") )
  #
  #
  # Build Saturated Model with Confidence Intervals
  #
  #modelSAT1  <- mxModel('testx',defSex, defAge, dataMZ, dataDZ, modelMZ, modelDZ, multi,name="fiveSATc") #, ciCov, ciCor, ciMean )
  #
  modelSAT1  <- mxModel('testx',modelMZ, modelDZ, multi,name="fiveSATc") #, ciCov, ciCor, ciMean )
  
  #
  fitSAT1    <- mxRun(modelSAT1)
  # end simcheck
} # dosimcheck
