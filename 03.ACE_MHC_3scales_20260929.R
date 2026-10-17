#########################
#
# MHC heritability study
# 
# Data cleaning, preparation and multivariate ACE twin analysis
#
# 01/08/2025; R version 4.3.2
#
#########################
# Preparing the environment
## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
rm(list=ls(all=TRUE)) #Delete every previous object in the environment


#### Install and load packages ####
## First specify the packages of interest
packages = c("OpenMx", "psych", "umx", "haven", "readxl", "dplyr", "ggplot2",
             "lme4", "grid", "reshape2", "here", "moments", "foreign", "gee", "tidyr",
             "ISLR", "polycor", "tidyverse", "broom", "mosaic", "ggcorrplot")

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
# Define subscale item groups
mhc_e_vars <- c("mhca14", "mhcb14", "mhcc14")                 # Emotional well-being
mhc_s_vars <- c("mhcd14", "mhce14", "mhcf14", "mhcg14", "mhch14") # Social well-being
mhc_p_vars <- c("mhci14", "mhcj14", "mhck14", "mhcl14", "mhcm14", "mhcn14") # Psychological well-being

# Compute subscale means with missingness rule (≤ 2 missing)
twin_df$mhc_e_14 <- ifelse(
  rowSums(is.na(twin_df[mhc_e_vars])) <= 1,
  rowMeans(twin_df[mhc_e_vars], na.rm = TRUE),
  NA
)

twin_df$mhc_s_14 <- ifelse(
  rowSums(is.na(twin_df[mhc_s_vars])) <= 1,
  rowMeans(twin_df[mhc_s_vars], na.rm = TRUE),
  NA
)

twin_df$mhc_p_14 <- ifelse(
  rowSums(is.na(twin_df[mhc_p_vars])) <= 1,
  rowMeans(twin_df[mhc_p_vars], na.rm = TRUE),
  NA
)

# Optional: check how many NAs in each
colSums(is.na(twin_df[, c("mhc_e_14", "mhc_s_14", "mhc_p_14")]))
# mhc_e_14 mhc_s_14 mhc_p_14 
# 1184     1214     1222
# sum(is.na(twin_df$mhc14)) # 1224 NA in mhc14
subset(twin_df, mhc_p_14 == 0) 


# Descriptives of the wb variables
describe(twin_df[,c("levenc14", "sat14", "hap14", "flo14", "mhc_e_14", "mhc_s_14", "mhc_p_14")])
# Histogram of the wb variables
# cvd1
layo=matrix(1:8,2,4)
layout(layo)
# cvd2
hist(twin_df$levenc14, main="QoL", xlab="Quality of Life", col="lightblue")
hist(twin_df$sat14, main="Satisfaction", xlab="Satisfaction", col="lightblue")
hist(twin_df$hap14, main="Happiness", xlab="Happiness", col="lightblue")
hist(twin_df$flo14, main="Flourishing", xlab="Flourishing", col="lightblue")
hist(twin_df$mhc_e_14, main="Emotional MHC", xlab="MHC", col="lightblue")
hist(twin_df$mhc_s_14, main="Social MHC", xlab="MHC", col="lightblue")
hist(twin_df$mhc_p_14, main="Psychological MHC", xlab="MHC", col="lightblue")

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
twin_df$mhc_e_sc <- scale(twin_df$mhc_e_14)
twin_df$mhc_s_sc <- scale(twin_df$mhc_s_14)
twin_df$mhc_p_sc <- scale(twin_df$mhc_p_14)

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
                                        mhck14, mhcl14, mhcm14, mhcn14, sat14, hap14, 
                                        flo14, mhc_e_14, mhc_s_14, mhc_p_14))

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
data_twins <- data_twins[complete.cases(data_twins[,c('qol_sc', 'sat_sc', 'hap_sc', 'flo_sc', 'mhc_e_sc', 'mhc_s_sc', 'mhc_p_sc')]),]
#
#cvd1
#
rph1 = summary(lm(qol_sc~Age_sc + sex, data=data_twins)) # data_twins[,c('qol_sc', 'sat_sc', 'hap_sc', 'flo_sc', 'emp_sc')]))
rph2 = summary(lm(sat_sc~Age_sc + sex, data=data_twins))
rph3 = summary(lm(hap_sc~Age_sc + sex, data=data_twins))
rph4 = summary(lm(flo_sc~Age_sc + sex, data=data_twins))
rph5 = summary(lm(mhc_e_sc~Age_sc + sex, data=data_twins))
rph6 = summary(lm(mhc_s_sc~Age_sc + sex, data=data_twins))
rph7 = summary(lm(mhc_p_sc~Age_sc + sex, data=data_twins))

R2s = c(rph1$r.square,rph2$r.square,rph3$r.square,rph4$r.square,rph5$r.square,rph6$r.square,rph7$r.square) 
parsb0 = c(
  rph1$coefficients[1,1],
  rph2$coefficients[1,1],
  rph3$coefficients[1,1],
  rph4$coefficients[1,1],
  rph5$coefficients[1,1],
  rph6$coefficients[1,1],
  rph7$coefficients[1,1])
parsbage = c( 
  rph1$coefficients[2,1],
  rph2$coefficients[2,1],
  rph3$coefficients[2,1],
  rph4$coefficients[2,1],
  rph5$coefficients[2,1],
  rph6$coefficients[2,1],
  rph7$coefficients[2,1])
parsbsex = c( 
  rph1$coefficients[3,1],
  rph2$coefficients[3,1],
  rph3$coefficients[3,1],
  rph4$coefficients[3,1],
  rph5$coefficients[3,1],
  rph6$coefficients[3,1],
  rph7$coefficients[3,1])
# cvd2
#
#
length(unique(data_twins$FamilyNumber)) #4188 families, 5212 i

table(as.data.frame(table(data_twins$FamilyNumber))[,2]) #Number of twins per family ## 2 --> 1024 complete pairs

table(data_twins$Extension)

# Renaming the columns for clarity
colnames(data_twins) <- c("FISNumber", "FamilyNumber", "Extension", "Sex", "Zyg", "Age", "QoL", "Sat", "Hap", "Flo", "EMHC", "SMHC", "PMHC", "Age_sc")

# how many individuals of each zygosity?
zyg_frequencies <- table(data_twins[,"Zyg"])
zygosities <- c("MZM", "DZM", "MZF", "DZF", "DOSmf", "DOSfm")
rownames(zyg_frequencies) <- zygosities
zyg_frequencies

# MZM   DZM   MZF   DZF DOSmf DOSfm 
# 695   391  2115   918   558   535

describe(data_twins[, c("Age", "Sex", "QoL", "Sat", "Hap", "Flo", "EMHC", "SMHC", "PMHC")]) 
#n=6 with a missing sex left

subset(data_twins, is.na(data_twins$Sex)) #We can assign sex to these 6 individuals based on zygosity
data_twins$Sex[is.na(data_twins$Sex) & data_twins$Zyg %in% c(1, 2)] <- 0
data_twins$Sex[is.na(data_twins$Sex) & data_twins$Zyg %in% c(3, 4)] <- 1
data_twins$Sex[is.na(data_twins$Sex) & data_twins$Zyg %in% c(5, 6)] <- NA
#
#Cross-variable correlations
cor.test(data_twins$QoL, data_twins$Sat) #0.757131 
cor.test(data_twins$QoL, data_twins$Hap) #0.7412574 
cor.test(data_twins$QoL, data_twins$Flo) #0.6906861 
cor.test(data_twins$QoL, data_twins$EMHC) #0.5913329
cor.test(data_twins$QoL, data_twins$SMHC) #0.3732284
cor.test(data_twins$QoL, data_twins$PMHC) #0.4867735

#
cor.test(data_twins$Sat, data_twins$Hap) #0.7550102 
cor.test(data_twins$Sat, data_twins$Flo) #0.7112383 
cor.test(data_twins$Sat, data_twins$EMHC) #0.561041
cor.test(data_twins$Sat, data_twins$SMHC) #0.358656
cor.test(data_twins$Sat, data_twins$PMHC) #0.4729856
#
cor.test(data_twins$Hap, data_twins$Flo) #0.715039 
cor.test(data_twins$Hap, data_twins$EMHC) #0.6283197
cor.test(data_twins$Hap, data_twins$SMHC) #0.3923974
cor.test(data_twins$Hap, data_twins$PMHC) #0.5236602
#
cor.test(data_twins$Flo, data_twins$EMHC) #0.5945883
cor.test(data_twins$Flo, data_twins$SMHC) #0.498382
cor.test(data_twins$Flo, data_twins$PMHC) #0.6167895
#
#
# MHC_comparison <- readRDS("MHC_comparison.rds")
# all(data_twins$FISNumber == MHC_comparison$FISNumber)
# data_twins$MHC <- MHC_comparison$MHC
# 
# cor.test(data_twins$MHC, data_twins$EMHC) #0.7896462 
# cor.test(data_twins$MHC, data_twins$SMHC) #0.89674469 
# cor.test(data_twins$MHC, data_twins$PMHC) #0.93072262 
# #If this part is ran, remove MHC column before continuing

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
(QoL_mz=cor.test(dataMZ$QoL_T1,dataMZ$QoL_T2)) #0.3774608  
(QoL_dz=cor.test(dataDZ$QoL_T1,dataDZ$QoL_T2)) #0.2636428 

(Sat_mz=cor.test(dataMZ$Sat_T1,dataMZ$Sat_T2)) #0.4005449 
(Sat_dz=cor.test(dataDZ$Sat_T1,dataDZ$Sat_T2)) #0.2130364 

(Hap_mz=cor.test(dataMZ$Hap_T1,dataMZ$Hap_T2)) #0.4085048 
(Hap_dz=cor.test(dataDZ$Hap_T1,dataDZ$Hap_T2)) #0.1483032 

(Flo_mz=cor.test(dataMZ$Flo_T1,dataMZ$Flo_T2)) #0.3842038 
(Flo_dz=cor.test(dataDZ$Flo_T1,dataDZ$Flo_T2)) #0.1471827 

(EMHC_mz=cor.test(dataMZ$EMHC_T1,dataMZ$EMHC_T2)) #0.3106738 
(EMHC_dz=cor.test(dataDZ$EMHC_T1,dataDZ$EMHC_T2)) #0.0583326 

(SMHC_mz=cor.test(dataMZ$SMHC_T1,dataMZ$SMHC_T2)) #0.39898 
(SMHC_dz=cor.test(dataDZ$SMHC_T1,dataDZ$SMHC_T2)) #0.1445766 

(PMHC_mz=cor.test(dataMZ$PMHC_T1,dataMZ$PMHC_T2)) #0.3734 
(PMHC_dz=cor.test(dataDZ$PMHC_T1,dataDZ$PMHC_T2)) #0.1221665 

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
vars      <- c('QoL','Sat', 'Hap', 'Flo', 'EMHC', 'SMHC', 'PMHC')              # list of variables names
nv        <- 7                         # number of variables
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
svPa      <- .3                        # start value for path coefficient
svPe      <- .5                        # start value for path coefficient for e

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
labMeZ    <- labVars("mean",selVars)
labCvMZ   <- labLower("covMZ",ntv)
labCvDZ   <- labLower("covDZ",ntv)
labCvZ    <- labLower("covZ",ntv)
labVaMZ   <- labDiag("covMZ",ntv)
labVaDZ   <- labDiag("covDZ",ntv)
labVaZ    <- labDiag("covZ",ntv)
#
# ----------------------------------------------------------------------------------------------------------------------
# PREPARE MODEL
defAge    <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=FALSE,
                       # cvd1
                       labels=c("data.Age_sc_T1","data.Age_sc_T2"), name="Age" )
# #Comparison CFI
# defAge <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=FALSE,
#                     labels=rep(c("data.Age_sc_T1","data.Age_sc_T2"), each=nv),
#                     name="Age" )
# cvd2
#pathB1     <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=.01, 
#              labels=c("bAge_QoL_T1","bAge_Sat_T1", "bAge_Hap_T1","bAge_Flo_T1", "bAge_Emp_T1", 
#                      "bAge_QoL_T2","bAge_Sat_T2", "bAge_Hap_T2","bAge_Flo_T2", "bAge_Emp_T2"), name="b1" )
# cvd1
pathB1     <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=c(T,T,T,T,T,T,T), # TRUE, 
                        values=rep(0,nv), # c(parsbage), 
                        labels=c("bAge_QoL","bAge_Sat", "bAge_Hap","bAge_Flo", "bAge_EMHC", "bAge_SMHC", "bAge_PMHC"), name="b1" )
# cvd2
# cvd1
defSex      <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=FALSE,
                         labels=c('data.Sex_T1','data.Sex_T2'), name="Sex" )
# #Comparison CFI:
# defSex <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=FALSE,
#                     labels=rep(c("data.Sex_T1","data.Sex_T2"), each=nv),
#                     name="Sex" )
# cvd2
#pathB2     <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=.2, 
#                labels=c("bSex_QoL_T1","bSex_Sat_T1", "bSex_Hap_T1","bSex_Flo_T1", "bSex_Emp_T1", 
#                        "bSex_QoL_T2","bSex_Sat_T2", "bSex_Hap_T2","bSex_Flo_T2", "bSex_Emp_T2"), name="b2" ) # Sex effects
#
pathB2     <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=c(T,T,F,F,F,T,T), 
                        values=c(0.05, 0.5, 0, 0, 0, 0.05,0.05), # c(parsbsex), 
                        labels=c("bSex_QoL","bSex_Sat","bSex_Hap","bSex_Flo", "bSex_EMHC", "bSex_SMHC", "bSexPMHC"), name="b2" ) # Sex effects

# Create Algebra for expected Mean Matrices
meanG     <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=labVars("mean",vars), name="meanG" )

expMeanG <- mxAlgebra( meanG + Sex*b2 + Age*b1, name="expMeanG")

inclusions   <- list (meanG, expMeanG, defAge, defSex, pathB1, pathB2 )

# Create Matrices for Variance Components
covA      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPa,nv), label=labLower("VA",nv), name="VA" ) 
covC      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPa,nv), label=labLower("VC",nv), name="VC" )
covE      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPa,nv), label=labLower("VE",nv), name="VE" )

# Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
covP      <- mxAlgebra( expression= VA+VC+VE, name="V" )
covMZ     <- mxAlgebra( expression= VA+VC, name="cMZ" )
covDZ     <- mxAlgebra( expression= 0.5%x%VA+ VC, name="cDZ" )
expCovMZ  <- mxAlgebra( expression= rbind( cbind(V, cMZ), cbind(t(cMZ), V)), name="expCovMZ" )
expCovDZ  <- mxAlgebra( expression= rbind( cbind(V, cDZ), cbind(t(cDZ), V)), name="expCovDZ" )

# Create Data Objects for Multiple Groups
dataMZ    <- mxData( observed=mzData, type="raw" )
dataDZ    <- mxData( observed=dzData, type="raw" )

# Create Expectation Objects for Multiple Groups
expMZ     <- mxExpectationNormal( covariance="expCovMZ", means="expMeanG", dimnames=selVars )
expDZ     <- mxExpectationNormal( covariance="expCovDZ", means="expMeanG", dimnames=selVars )
funML     <- mxFitFunctionML()

# Create Model Objects for Multiple Groups
pars      <- list( covA, covC, covE, covP )
modelMZ   <- mxModel( inclusions, pars, covMZ, expCovMZ, dataMZ, expMZ, funML, name="MZ" )
modelDZ   <- mxModel( inclusions, pars, covDZ, expCovDZ, dataDZ, expDZ, funML, name="DZ" )
multi     <- mxFitFunctionMultigroup( c("MZ","DZ") )

# Create Algebra for Standardization
matI      <- mxMatrix( type="Iden", nrow=nv, ncol=nv, name="I" )
invSD     <- mxAlgebra( expression=solve(sqrt(I*V)), name="iSD" )

# Calculate genetic and environmental correlations
corA      <- mxAlgebra( expression=solve(sqrt(I*VA))%&%VA, name ="rA" ) #cov2cor()
corC      <- mxAlgebra( expression=solve(sqrt(I*VC))%&%VC, name ="rC" )
corE      <- mxAlgebra( expression=solve(sqrt(I*VE))%&%VE, name ="rE" )
corP      <- mxAlgebra( expression=solve(sqrt(I*V))%&%V, name ="rP") 

# Create Algebra for Unstandardized and Standardized Variance Components
rowUV     <- rep('UV',nv)
colUV     <- rep(c('VA','VC','VE','SA','SC','SE'),each=nv)
estUV     <- mxAlgebra( expression=cbind(VA,VC,VE,VA/V,VC/V,VE/V), name="UV", dimnames=list(rowUV,colUV) )

rowSV     <- rep('SV',nv)
colSV     <- rep(c('SA','SC','SE'),each=nv)
estSV     <- mxAlgebra( expression=cbind(VA/V,VC/V,VE/V), name="SV", dimnames=list(rowSV,colSV) )

rowcorr     <- rep('corr',nv)
colcorr     <- rep(c('rA','rC','rE','rP'),each=nv)
estcorr     <- mxAlgebra( expression=cbind(rA,rC,rE,rP), name="corr", dimnames=list(rowcorr,colcorr) )

# # Create Confidence Interval Objects
# odd       <- seq(1+3*nv,2*3*nv,nv)
# even      <- seq(2+3*nv,2*3*nv,nv)
# ciACE     <- mxCI( c("US[1,odd]","US[2,odd]","US[2,even]") )

## Create Confidence Interval Objects
ciUV      <- mxCI( c("UV") )
ciSV      <- mxCI( c("SV") )
cicorr      <- mxCI( c("corr") )

# Build Model with Confidence Intervals
calc      <- list( matI, invSD, corA, corC, corE, corP, estUV, estSV, ciUV, ciSV,estcorr,cicorr )
modelACE  <- mxModel( "fiveACEvc", pars, modelMZ, modelDZ, multi, calc )

# ----------------------------------------------------------------------------------------------------------------------
# RUN MODEL

# Run ACE Model
fitACE    <- mxRun( modelACE, intervals=F )
(sumACE    <- summary( fitACE ))

# Compare with Saturated Model
#mxCompare( fitSATOptimal, fitACE )
# base comparison ep minus2LL    df      AIC    diffLL diffdf  p
# 1 SAToptimal       <NA> 60 56405.92 26075 56525.92        NA     NA NA
# 2 SAToptimal  fiveACEvc 56 56216.42 26079 56328.42 -189.4999      4  1

fitACE$algebras$UV
fitACE$algebras$SV

# Print Goodness-of-fit Statistics & Parameter Estimates
# fitGofs(fitACE)
# fitEstCis(fitACE)

# ----------------------------------------------------------------------------------------------------------------------
# RUN SUBMODELS

# Run AE model
modelAE   <- mxModel( fitACE, name="twoAEvc" )
modelAE   <- omxSetParameters( modelAE, labels=labLower("VC",nv), free=FALSE, values=0 )
fitAE     <- mxRun( modelAE, intervals=T )
fitGofs(fitAE); fitEstCis(fitAE)

mxCompare( fitACE, fitAE )

summary(fitAE)

# library(openxlsx)
# write.xlsx(as.data.frame(summary(fitAE)$parameters), "SummaryAE_7v.xlsx")

fitAE$algebras$UV
fitAE$algebras$SV

fitAE$algebras$rA
fitAE$algebras$rE
fitAE$algebras$rP

#If variance MZ is lower than DZ (massive) it implies GxE interaction
save.image("/data/natalia/MHC/MHC_env_ACEAE_20251103.RData")
saveRDS(fitAE, file = "/data/natalia/MHC/AE_7v_model.rds")

ACor_df <- as.data.frame(fitAE$algebras$rA@result)
rownames(ACor_df) <- c("QoL", "LS", "SH", "FL", "EMHC", "SMHC", "PMHC")  
colnames(ACor_df) <- c("QoL", "LS", "SH", "FL", "EMHC", "SMHC", "PMHC")  

PCor_df <- as.data.frame(fitAE$algebras$rP@result)
rownames(PCor_df) <- c("QoL", "LS", "SH", "FL", "EMHC", "SMHC", "PMHC")  
colnames(PCor_df) <- c("QoL", "LS", "SH", "FL", "EMHC", "SMHC", "PMHC")  

ECor_df <- as.data.frame(fitAE$algebras$rE@result)
rownames(ECor_df) <- c("QoL", "LS", "SH", "FL", "EMHC", "SMHC", "PMHC")  
colnames(ECor_df) <- c("QoL", "LS", "SH", "FL", "EMHC", "SMHC", "PMHC")  

# Plot significant correlations
# ggcorrplot(ACor_df, 
#            hc.order = FALSE,  
#            type = "lower",   
#            lab = TRUE,
#            lab_size = 3,
#            colors = c("white", "lightgreen", "darkgreen"))  
# 
# ggcorrplot(PCor_df, 
#            hc.order = FALSE,  
#            type = "upper",   
#            lab = TRUE,
#            lab_size = 3,
#            colors = c("white", "lightblue", "darkblue"))

# Plot upper triangle with one color scheme
library(corrplot)
library(ggplot2)

layo=matrix(1,1)
layout(layo)

corrplot(as.matrix(PCor_df), method = "circle", 
         type = "upper", 
         #         col = colorRampPalette(c("#b3e0db", "#5dc1b9", "#325c58"))(200), 
         col = colorRampPalette(c("#275317", "#275317", "#d9f2d0"))(200), 
         diag = FALSE, 
         tl.col = "black",       # labels in black
         tl.pos = "lt",
         tl.srt = 45,            # rotate labels 45 degrees (adjust angle)
         tl.offset = 0.5,        # move labels a bit out (adjust as needed)
         cl.pos = "r",           # color legend bottom here
         cl.offset = 0.15,       # move legend horizontally (adjust)
         cl.align.text = "l",    # align legend text to left
         addCoef.col = "black",  # add correlation values in black
         number.cex = 1)       # size of correlation numbers

# save the image #To guess the inches, divide the px shown in the url of the external viewer by 96
#width 1142 / 96 ≈ 11.9 inches; height =  975 / 96 ≈ 10.2 inches
#ggsave("PCorrelationMatrix.tiff",last_plot(), path = "C:/Users/naa840/OneDrive - Vrije Universiteit Amsterdam/Nueva carpeta/MHC/06.Plots/", width = 12, height = 11, dpi=300)
#Saved as tiff with dimensions 670*568
# Add lower triangle with different colors and legend on right with adjusted position
# corrplot(as.matrix(ACor_df), method = "circle", 
#          type = "lower", 
#          #  col = colorRampPalette(c("#a6dac2", "#3eb489", "#265643"))(200),
#          col = colorRampPalette(c("#084F6A", "#084F6A", "#CAEEFB"))(200), 
#          diag = FALSE, 
#          tl.pos = "l",           # labels here
#          tl.col = "black",       # labels in black
#          tl.srt = 45,            # rotate labels 45 degrees (adjust angle)
#          tl.offset = 0.5,        # move labels a bit out (adjust as needed)
#          cl.pos = "r",           # color legend on right
#          cl.offset = 0.15,       # move legend horizontally (adjust)
#          cl.align.text = "l",    # align legend text to left
#          #         add = TRUE,
#          addCoef.col = "black",
#          number.cex = 1)
pdf("ACorrelationMatrix_PDF.pdf", width = 9, height = 8)
corrplot(as.matrix(ACor_df), method = "circle", 
         type = "upper", 
         #         col = colorRampPalette(c("#b3e0db", "#5dc1b9", "#325c58"))(200), 
         col = colorRampPalette(c("#084F6A", "#084F6A", "#CAEEFB"))(200), 
         diag = FALSE, 
         tl.col = "black",       # labels in black
         tl.pos = "lt",
         tl.srt = 45,            # rotate labels 45 degrees (adjust angle)
         tl.offset = 0.5,        # move labels a bit out (adjust as needed)
         cl.pos = "r",           # color legend bottom here
         cl.offset = 0.15,       # move legend horizontally (adjust)
         cl.align.text = "l",    # align legend text to left
         addCoef.col = "black",  # add correlation values in black
         number.cex = 1)       # size of correlation numbers

dev.off()
#width 860 / 96 ≈ 8.96 inches; height =  751 / 96 ≈ 7.83 inches
ggsave("ACorrelationMatrix.tiff",last_plot(), path = "/data/natalia/MHC/", width = 9, height = 8, dpi=300)

pdf("ECorrelationMatrix_PDF.pdf", width = 9, height = 8)
corrplot(as.matrix(ECor_df), method = "circle", 
         type = "lower", 
         #  col = colorRampPalette(c("#a6dac2", "#3eb489", "#265643"))(200),
         col = colorRampPalette(c("#80350E", "#80350E", "#FBE3D6"))(200), 
         diag = FALSE, 
         tl.pos = "l",           # labels here
         tl.col = "black",       # labels in black
         tl.srt = 45,            # rotate labels 45 degrees (adjust angle)
         tl.offset = 0.5,        # move labels a bit out (adjust as needed)
         cl.pos = "r",           # color legend on right
         cl.offset = 0.15,       # move legend horizontally (adjust)
         cl.align.text = "l",    # align legend text to left
         #         add = TRUE,
         addCoef.col = "black",
         number.cex = 1)
dev.off()

#width 860 / 96 ≈ 8.96 inches; height =  751 / 96 ≈ 7.83 inches
ggsave("ECorrelationMatrix.tiff",last_plot(), path = "/data/natalia/MHC/", width = 9, height = 8, dpi=300)

# corrplot(as.matrix(ECor_df), method = "circle", 
#          type = "upper", 
#          #         col = colorRampPalette(c("#b3e0db", "#5dc1b9", "#325c58"))(200), 
#          col = colorRampPalette(c("#80350E", "#80350E", "#FBE3D6"))(200), 
#          diag = FALSE, 
#          tl.col = "black",       # labels in black
#          tl.pos = "lt",
#          tl.srt = 45,            # rotate labels 45 degrees (adjust angle)
#          tl.offset = 0.5,        # move labels a bit out (adjust as needed)
#          cl.pos = "r",           # color legend bottom here
#          cl.offset = 0.15,       # move legend horizontally (adjust)
#          cl.align.text = "l",    # align legend text to left
#          addCoef.col = "black",  # add correlation values in black
#          number.cex = 1)       # size of correlation numbers


#Graph Lianne
library(epiDisplay)
library(gridExtra)
library(ggplot2)
# Build the dataframe
#To calculate the value, it is the Pcor multiplied by the standardized variance (either A or E). To calculate the CI, it is the correlation multiplied by the upper or lower bound of the standardized covariance. For the CI of A, as inthe graph they have to appear on top of E, you sum up the value of BivE.
#fitAE$algebras$rP
#fitAE$algebras$SV
#fitAE$output$confidenceIntervals
# df <- data.frame(
#   Correlation = rep(c("Biv A", "Biv E"), times = 4),
#   Measure = rep(c(rep("QoL", 2), 
#                   rep("SWLS", 2), 
#                   rep("SHS", 2), 
#                   rep("Flourishing", 2)), each = 1),
#   Value = c(0.53*0.52, 0.53*0.48, 0.50*0.52, 0.50*0.48,0.56*0.51, 0.56*0.49, 0.64*0.49, 0.64*0.51),
#   L_ci = c(0.53*0.48+0.53*0.43, 0.53*0.40, 0.50*0.48+0.50*0.43, 0.50*0.40,0.56*0.49+0.56*0.43, 0.56*0.41, 0.64*0.51+0.64*0.42, 0.64*0.44),
#   U_ci = c(0.53*0.48+0.53*0.60, 0.53*0.57, 0.50*0.48+0.50*0.60, 0.50*0.57,0.56*0.49+0.56*0.59, 0.56*0.57, 0.64*0.51+0.64*0.56, 0.64*0.58),
#   Cor = c(NA, 0.53, NA, 0.50, NA, 0.56, NA, 0.64)
# )
# df$Correlation <- factor(df$Correlation, levels = c('Biv A', 'Biv E'))
# df$Measure <- factor(df$Measure, levels = c('QoL','SWLS', 'SHS', 'Flourishing'))
# 
# # tiff("Figure 2 Plot_gencor.tiff",  width = 10, height = 6, units = 'in', res = 300)
# ggplot(data = df, aes(x = factor(Measure), y = Value, fill = factor(Correlation))) + 
#   geom_bar(stat = "identity") + 
#   xlab(" ") + ylab("Correlation") + 
#   scale_fill_manual(name = " ", values = c('#CAEEFB', '#FBE3D6')) +
#   theme_classic(base_size = 12) + 
#   geom_text(aes(label = round(Cor, 2)), 
#             #            vjust = c(9,8,5,5,-1,-1,-1,-1, 9,10,6,4,-2,-2,-.5,-1), 
#             position = "stack") + 
#   ylim(0, 0.75) + 
#   theme(legend.position = "bottom", axis.text.x = element_text(angle = 0)) +
#   geom_errorbar(aes(ymax = L_ci, ymin = U_ci),
#                 position = position_dodge(0.9), width = 0.10)
#rP
rP_E <- fitAE$rP$result[1:4,5]
rP_S <- fitAE$rP$result[1:4,6]
rP_P <- fitAE$rP$result[1:4,7]

# Standardized A and E components
SVA <- fitAE$VA$values / (fitAE$VA$values+fitAE$VE$values) 
#SVA <- SVA[2:6,1]
SVA_E <- SVA[1:4,5]
SVE_E <- 1 - SVA_E
SVA_S <- SVA[1:4,6]
SVE_S <- 1 - SVA_S
SVA_P <- SVA[1:4,7]
SVE_P <- 1 - SVA_P

#CIs 
CIs_E <- as.data.frame(fitAE$output$confidenceIntervals[grep("^twoAEvc.SV\\[5", rownames(fitAE$output$confidenceIntervals)), ])
CIs_S <- as.data.frame(fitAE$output$confidenceIntervals[grep("^twoAEvc.SV\\[6", rownames(fitAE$output$confidenceIntervals)), ])
CIs_P <- as.data.frame(fitAE$output$confidenceIntervals[grep("^twoAEvc.SV\\[7", rownames(fitAE$output$confidenceIntervals)), ])

# CIs for A
CI_lbA_E <- CIs_E$lbound[1:4]
CI_ubA_E <- CIs_E$ubound[1:4]

CI_lbA_S <- CIs_S$lbound[1:4]
CI_ubA_S <- CIs_S$ubound[1:4]

CI_lbA_P <- CIs_P$lbound[1:4]
CI_ubA_P <- CIs_P$ubound[1:4]

# CIs for E
CI_lbE_E <- CIs_E$lbound[15:18]
CI_ubE_E <- CIs_E$ubound[15:18]

CI_lbE_S <- CIs_S$lbound[15:18]
CI_ubE_S <- CIs_S$ubound[15:18]

CI_lbE_P <- CIs_P$lbound[15:18]
CI_ubE_P <- CIs_P$ubound[15:18]

# Build the dataframe
# To calculate the value, it is the Pcor multiplied by the standardized variance (either A or E). 
# To calculate the CI, it is the correlation multiplied by the upper or lower bound of the standardized covariance. 
# For the CI of A, as inthe graph they have to appear on top of E, you sum up the value of BivE.
df_E <- data.frame(
  Correlation = rep(c("Biv A", "Biv E"), times = 4),
  Measure = rep(c(rep("QoL", 2),
                  rep("LS", 2),
                  rep("SH", 2),
                  rep("FL", 2)), each = 1),
  Value = c(rP_E[1]*SVA_E[1], rP_E[1]*SVE_E[1], 
            rP_E[2]*SVA_E[2], rP_E[2]*SVE_E[2],
            rP_E[3]*SVA_E[3], rP_E[3]*SVE_E[3], 
            rP_E[4]*SVA_E[4], rP_E[4]*SVE_E[4]),
  # (cor * SVE + cor*LLA, cor * LLE )
  L_ci = c(rP_E[1]*SVE_E[1]+rP_E[1]*CI_lbA_E[1], rP_E[1]*CI_lbE_E[1], 
           rP_E[2]*SVE_E[2]+rP_E[2]*CI_lbA_E[2], rP_E[2]*CI_lbE_E[2],
           rP_E[3]*SVE_E[3]+rP_E[3]*CI_lbA_E[3], rP_E[3]*CI_lbE_E[3], 
           rP_E[4]*SVE_E[4]+rP_E[4]*CI_lbA_E[4], rP_E[4]*CI_lbE_E[4]),
  # (cor * SVE + cor*ULA, cor * ULE )
  U_ci = c(rP_E[1]*SVE_E[1]+rP_E[1]*CI_ubA_E[1], rP_E[1]*CI_ubE_E[1], 
           rP_E[2]*SVE_E[2]+rP_E[2]*CI_ubA_E[2], rP_E[2]*CI_ubE_E[2],
           rP_E[3]*SVE_E[3]+rP_E[3]*CI_ubA_E[3], rP_E[3]*CI_ubE_E[3], 
           rP_E[4]*SVE_E[4]+rP_E[4]*CI_ubA_E[4], rP_E[4]*CI_ubE_E[4]),
  Cor = c(NA, rP_E[1], NA, rP_E[2], NA, rP_E[3], NA, rP_E[4])
)

df_E$Correlation <- factor(df_E$Correlation, levels = c('Biv A', 'Biv E'))
df_E$Measure <- factor(df_E$Measure, levels = c('QoL','LS', 'SH', 'FL'))

df_S <- data.frame(
  Correlation = rep(c("Biv A", "Biv E"), times = 4),
  Measure = rep(c(rep("QoL", 2),
                  rep("LS", 2),
                  rep("SH", 2),
                  rep("FL", 2)), each = 1),
  Value = c(rP_S[1]*SVA_S[1], rP_S[1]*SVE_S[1], 
            rP_S[2]*SVA_S[2], rP_S[2]*SVE_S[2],
            rP_S[3]*SVA_S[3], rP_S[3]*SVE_S[3], 
            rP_S[4]*SVA_S[4], rP_S[4]*SVE_S[4]),
  # (cor * SVE + cor*LLA, cor * LLE )
  L_ci = c(rP_S[1]*SVE_S[1]+rP_S[1]*CI_lbA_S[1], rP_S[1]*CI_lbE_S[1], 
           rP_S[2]*SVE_S[2]+rP_S[2]*CI_lbA_S[2], rP_S[2]*CI_lbE_S[2],
           rP_S[3]*SVE_S[3]+rP_S[3]*CI_lbA_S[3], rP_S[3]*CI_lbE_S[3], 
           rP_S[4]*SVE_S[4]+rP_S[4]*CI_lbA_S[4], rP_S[4]*CI_lbE_S[4]),
  # (cor * SVE + cor*ULA, cor * ULE )
  U_ci = c(rP_S[1]*SVE_S[1]+rP_S[1]*CI_ubA_S[1], rP_S[1]*CI_ubE_S[1], 
           rP_S[2]*SVE_S[2]+rP_S[2]*CI_ubA_S[2], rP_S[2]*CI_ubE_S[2],
           rP_S[3]*SVE_S[3]+rP_S[3]*CI_ubA_S[3], rP_S[3]*CI_ubE_S[3], 
           rP_S[4]*SVE_S[4]+rP_S[4]*CI_ubA_S[4], rP_S[4]*CI_ubE_S[4]),
  Cor = c(NA, rP_S[1], NA, rP_S[2], NA, rP_S[3], NA, rP_S[4])
)

df_S$Correlation <- factor(df_S$Correlation, levels = c('Biv A', 'Biv E'))
df_S$Measure <- factor(df_S$Measure, levels = c('QoL','LS', 'SH', 'FL'))

df_P <- data.frame(
  Correlation = rep(c("Biv A", "Biv E"), times = 4),
  Measure = rep(c(rep("QoL", 2),
                  rep("LS", 2),
                  rep("SH", 2),
                  rep("FL", 2)), each = 1),
  Value = c(rP_P[1]*SVA_P[1], rP_P[1]*SVE_P[1], 
            rP_P[2]*SVA_P[2], rP_P[2]*SVE_P[2],
            rP_P[3]*SVA_P[3], rP_P[3]*SVE_P[3], 
            rP_P[4]*SVA_P[4], rP_P[4]*SVE_P[4]),
  # (cor * SVE + cor*LLA, cor * LLE )
  L_ci = c(rP_P[1]*SVE_P[1]+rP_P[1]*CI_lbA_P[1], rP_P[1]*CI_lbE_P[1], 
           rP_P[2]*SVE_P[2]+rP_P[2]*CI_lbA_P[2], rP_P[2]*CI_lbE_P[2],
           rP_P[3]*SVE_P[3]+rP_P[3]*CI_lbA_P[3], rP_P[3]*CI_lbE_P[3], 
           rP_P[4]*SVE_P[4]+rP_P[4]*CI_lbA_P[4], rP_P[4]*CI_lbE_P[4]),
  # (cor * SVE + cor*ULA, cor * ULE )
  U_ci = c(rP_P[1]*SVE_P[1]+rP_P[1]*CI_ubA_P[1], rP_P[1]*CI_ubE_P[1], 
           rP_P[2]*SVE_P[2]+rP_P[2]*CI_ubA_P[2], rP_P[2]*CI_ubE_P[2],
           rP_P[3]*SVE_P[3]+rP_P[3]*CI_ubA_P[3], rP_P[3]*CI_ubE_P[3], 
           rP_P[4]*SVE_P[4]+rP_P[4]*CI_ubA_P[4], rP_P[4]*CI_ubE_P[4]),
  Cor = c(NA, rP_P[1], NA, rP_P[2], NA, rP_P[3], NA, rP_P[4])
)

df_P$Correlation <- factor(df_P$Correlation, levels = c('Biv A', 'Biv E'))
df_P$Measure <- factor(df_P$Measure, levels = c('QoL','LS', 'SH', 'FL'))

# tiff("Figure 2 Plot_gencor.tiff",  width = 10, height = 6, units = 'in', res = 300)

p1 <- ggplot(data = df_E, aes(x = factor(Measure), y = Value, fill = factor(Correlation))) +
  geom_bar(stat = "identity") +
  xlab(" ") + ylab("Phenotypic correlation") +
  ggtitle("EMHC") +
  scale_fill_manual(name = " ", values = c('#CAEEFB', '#FBE3D6')) +
  theme_classic(base_size = 12) +
  geom_text(
    data = subset(df_E, !is.na(Cor)),
    aes(label = sub("^0", "", sprintf("%.2f", Cor))),
    position = position_stack(vjust = 1.15)
  ) +
  ylim(0, 0.75) +
  theme(legend.position = "", axis.text.x = element_text(angle = 0),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
  geom_errorbar(aes(ymax = L_ci, ymin = U_ci),
                position = position_dodge(0.9), width = 0.10)

p2 <- ggplot(data = df_S, aes(x = factor(Measure), y = Value, fill = factor(Correlation))) +
  geom_bar(stat = "identity") +
  xlab(" ") + ylab("Phenotypic correlation") +
  ggtitle("SMHC") +
  scale_fill_manual(name = " ", values = c('#CAEEFB', '#FBE3D6')) +
  theme_classic(base_size = 12) +
  geom_text(
    data = subset(df_S, !is.na(Cor)),
    aes(label = sub("^0", "", sprintf("%.2f", Cor))),
    position = position_stack(vjust = 1.15)
  ) +
  ylim(0, 0.75) +
  theme(legend.position = "", axis.text.x = element_text(angle = 0),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
  geom_errorbar(aes(ymax = L_ci, ymin = U_ci),
                position = position_dodge(0.9), width = 0.10)

p3 <- ggplot(data = df_P, aes(x = factor(Measure), y = Value, fill = factor(Correlation))) +
  geom_bar(stat = "identity") +
  xlab(" ") + ylab("Phenotypic correlation") +
  ggtitle("PMHC") +
  scale_fill_manual(name = " ", values = c('#CAEEFB', '#FBE3D6')) +
  theme_classic(base_size = 12) +
  geom_text(
    data = subset(df_P, !is.na(Cor)),
    aes(label = sub("^0", "", sprintf("%.2f", Cor))),
    position = position_stack(vjust = 1.15)
  ) +
  ylim(0, 0.75) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 0),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
  geom_errorbar(aes(ymax = L_ci, ymin = U_ci),
                position = position_dodge(0.9), width = 0.10)

library(patchwork)
p1 / p2 / p3
#width 654 / 96 ≈ 6.81 inches; height =  900 / 96 ≈ 9.38 inches
ggsave("BivariateCorrelation.tiff",last_plot(), path = "/data/natalia/MHC/", width = 7, height = 10, dpi=300)

# modelAE2   <- mxModel( fitAE, name="twoAEvc2" )
# modelAE2 <- omxSetParameters( modelAE2, labels=c("VA32", "VA42", "VA43"), free=FALSE, values=0 ) #doesnot work
# fitAE2     <- mxRun( modelAE2, intervals=F )
# 
# rAc1 <- mxConstraint(abs(rA[1,2]) == 1)
# rAc2 <- mxConstraint(abs(rA[1,3]) == 1)
# rAc3 <- mxConstraint(abs(rA[1,4]) == 1)
# rAc4 <- mxConstraint(abs(rA[2,3]) == 1)
# rAc5 <- mxConstraint(abs(rA[2,4]) == 1)
# rAc6 <- mxConstraint(abs(rA[3,4]) == 1)
# modelAE2 <- mxModel(modelAE, name="rG1", rAc1, rAc2, rAc3, rAc4, rAc5, rAc6)
# fitAE2 <- mxRun( modelAE2, intervals=F)
# mxCompare(fitAE,fitAE2)
# 
# fitAE2$algebras$rA
# fitAE2$algebras$UV
# 
# cov2cor(fitAE2$algebras$UV$result[1:5,1:5])



# Run CE model
modelCE   <- mxModel( fitACE, name="twoCEvc" )
modelCE   <- omxSetParameters( modelCE, labels=labLower("VA",nv), free=FALSE, values=0 )
modelCE   <- omxSetParameters( modelCE, labels=labLower("VC",nv), free=TRUE, values=.6 )
fitCE     <- mxRun( modelCE, intervals=F )
fitGofs(fitCE); fitEstCis(fitCE)




# Run E model
modelE    <- mxModel( fitAE, name="twoEvc" )
modelE    <- omxSetParameters( modelE, labels=labLower("VA",nv), free=FALSE, values=0 )
fitE      <- mxRun( modelE, intervals=F )
fitGofs(fitE); fitEstCis(fitE)

# Print Comparative Fit Statistics
mxCompare( fitACE, nested <- list(fitAE, fitCE, fitE) )
#round(rbind(fitACE$US$result,fitAE$US$result,fitCE$US$result,fitE$US$result),4)

# Add genetic and environmental correlations
modelAE_r <- mxModel(fitAE, name="AE_r" )
# Set up rA confidence intervals
ciVarA <- mxCI( c("rA[6,5]")  )
modelAE_rA <-mxModel(modelAE_r,ciVarA)
#fitAE_rA <- omxRunCI(modelAE_rA, optimizer = "NPSOL")
fitAE_rA <- omxRunCI(modelAE_rA, optimizer = "CSOLNP")

options(max.print = 99999)
summary(fitAE_rA)