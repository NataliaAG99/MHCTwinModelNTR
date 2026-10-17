#########################
#
# MHC heritability study
# 
# Revision: Indep model creation and comparison for RMSEA
#
# 23/09/2026; R version 4.5.1
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

# Flourishing is the sum up of the items only if they have 1 or none missing values
flo_vars <- c("flo_a14", "flo_b14", "flo_c14", "flo_d14", "flo_e14", "flo_f14", "flo_g14", "flo_h14")
twin_df$flo14 <- ifelse(
  rowSums(is.na(twin_df[flo_vars])) < 2,
  rowSums(twin_df[flo_vars], na.rm = TRUE),
  NA
)
sum(is.na(twin_df$flo14)) # 891 NA in flo14


# MHC is the average of the items
# mhc_vars <- c("mhca14", "mhcb14", "mhcc14", "mhcd14", "mhce14", "mhcf14", "mhcg14", "mhch14", "mhci14", "mhcj14", "mhck14", "mhcl14", "mhcm14", "mhcn14")
# # # Compute average only when two or less items are missing
# twin_df$mhc14 <- ifelse(
#   rowSums(is.na(twin_df[mhc_vars])) <= 2,
#   rowMeans(twin_df[mhc_vars], na.rm = TRUE),
#   NA
# )
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

# fisnumbersMHC <- data_twins$FISNumber
# saveRDS(fisnumbersMHC, file = "/data/natalia/MHC/fisnumbersMHC.rds")
# rm(fisnumbersMHC)

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
                        labels=c("bAge_QoL","bAge_Sat", "bAge_Hap","bAge_Flo", "bAge_EMHC", "bAge_SMHC", "bAge_PMHC"), name="b1" )
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
                        labels=c("bSex_QoL","bSex_Sat","bSex_Hap","bSex_Flo", "bSex_EMHC", "bSex_SMHC", "bSex_PMHC"), name="b2" ) # Sex effects
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
                     labels=c("b0_QoL","b0_Sat", "b0_Hap","b0_Flo", "b0_EMHC", "b0_SMHC", "b0_PMHC",
                              "b0_QoL","b0_Sat", "b0_Hap","b0_Flo", "b0_EMHC", "b0_SMHC", "b0_PMHC"), name="b0MZ" )
b0DZ    <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE,
                     values=c(parsb0,parsb0),
                     labels=c("b0_QoL","b0_Sat", "b0_Hap","b0_Flo", "b0_EMHC", "b0_SMHC", "b0_PMHC",
                              "b0_QoL","b0_Sat", "b0_Hap","b0_Flo", "b0_EMHC", "b0_SMHC", "b0_PMHC"), name="b0DZ" )
#
# cvd2
# The order of the data object and the beta object matters
# e.g. for 2 traits the outcome MUST be Twin1-trait1 Twin1-trait2 Twin2-trait2 Twin2-trait2
# Check the matrix algebra
#
# expMeanMZ_ <- mxAlgebra( expression = b0MZ + b2%x%Sex + b1%x%Age, name="expMeanMZ")
# expMeanDZ_ <- mxAlgebra( expression = b0DZ + b2%x%Sex + b1%x%Age, name="expMeanDZ")
expMeanMZ_ <- mxAlgebra( expression = b0MZ + Sex%x%b2 + Age%x%b1, name="expMeanMZ")
expMeanDZ_ <- mxAlgebra( expression = b0DZ + Sex%x%b2 + Age%x%b1, name="expMeanDZ")

#
#bits   <- list(defSex, defAge, pathB1, pathB2 )
#
## Create Algebra for expected Variance/Covariance Matrices
#covMZ     <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=valDiag(svVa,ntv), lbound=valDiag(lbVa,ntv), labels=labCvMZ, name="covMZ" )
#covDZ     <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=valDiag(svVa,ntv), lbound=valDiag(lbVa,ntv), labels=labCvDZ, name="covDZ" )
#
# covMZ_     <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=valDiag(svVa,ntv), labels=labCvMZ, name="covMZ" )
# covDZ_     <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=valDiag(svVa,ntv), labels=labCvDZ, name="covDZ" )
#
covMZ_ <- mxMatrix(
  type   = "Symm",
  nrow   = ntv,
  ncol   = ntv,
  free   = diag(ntv) == 1,
  values = diag(svVa, ntv),
  name   = "covMZ"
)

covDZ_ <- mxMatrix(
  type   = "Symm",
  nrow   = ntv,
  ncol   = ntv,
  free   = diag(ntv) == 1,
  values = diag(svVa, ntv),
  name   = "covDZ"
)
#
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
# ciMean    <- mxCI( c('MZ.meanMZ','DZ.meanDZ') )
#
#
# Build Saturated Model with Confidence Intervals
#
modelSAT  <- mxModel("fiveSATc", modelMZ, modelDZ, multi, ciCov, ciCor) #ciMean )

# Check number of parameters
check_params <- omxGetParameters(modelSAT)
length(names(check_params))

modelSAT$MZ$covMZ
modelSAT$DZ$covDZ


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
fitIND    <- mxRun(modelSAT)
#fitSAT    <- mxTryHard( modelSAT,10 )
#
#
(sumIND    <- summary( fitIND ))

saveRDS(fitIND, file = "/data/natalia/MHC/Indep_model_revision_comparison.rds")


#Comparison

INDEP <- readRDS("/data/natalia/MHC/revision/Indep_model_revision_comparison.rds")
SAT <- readRDS("/data/natalia/MHC/revision/Sat_model_revision_comparison.rds")
MOD1 <- readRDS("/data/natalia/MHC/revision/AE_7v_model_CFI.rds")
MOD2 <- readRDS("/data/natalia/MHC/revision/AE_cP_model_CFI.rds")
MOD3 <- readRDS("/data/natalia/MHC/revision/AE_iP_model_CFI.rds")


CFI_RMSEA <- function(indep,sat,mod) {
  sum_ind <- summary(indep)
  sum_mod <- summary(mod)
  sum_sat <- summary(sat)
  chi <- sum_mod$Minus2LogLikelihood - sum_sat$Minus2LogLikelihood
  chiDoF <- sum_mod$degreesOfFreedom - sum_sat$degreesOfFreedom
  numObs <- sum_mod$numObs
  
  RMSEA <- sqrt( (chi / (chiDoF) - 1) / numObs)
  
  num <- (sum_ind$Minus2LogLikelihood - sum_ind$degreesOfFreedom - sum_mod$Minus2LogLikelihood + sum_mod$degreesOfFreedom)
  denom <- (sum_ind$Minus2LogLikelihood - sum_ind$degreesOfFreedom - sum_sat$Minus2LogLikelihood + sum_sat$degreesOfFreedom)
  CFI <- num / denom
  
  return(list(CFI=CFI,RMSEA=RMSEA))
}

CFI_RMSEA(indep = INDEP, sat = SAT, mod = MOD1)
CFI_RMSEA(indep = INDEP, sat = SAT, mod = MOD2)
CFI_RMSEA(indep = INDEP, sat = SAT, mod = MOD3)
