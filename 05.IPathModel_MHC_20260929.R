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
#CFI Comparison
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
#CFI Comparison
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

svPa <- .4 # start value for path coefficient
svPe <- .5 # start value for path coefficient for e

# INDEPENDENT PATHWAY 1 FACTOR

# Fit Independent Pathway ACE Model
# ------------------------------------------------------------------------------  

nf        <- 1 # number of factors for each variance component

# ------------------------------------------------------------------------------
# PREPARE MODEL

# Factor loading matrix (5 variables x 1 factor)
svfl      <- rep(.2,nv)
frfl      <- TRUE


# Matrices Al, Dl, and El for variance components for latent phenotype(s)
covAl     <- mxMatrix( type="Full", nrow=nv, ncol=nf, free=TRUE, values=rep(0.2,nv), label=labDiag("Al",nv), name="Al" ) 
covCl     <- mxMatrix( type="Full", nrow=nv, ncol=nf, free=TRUE, values=rep(0,nv), label=labDiag("Cl",nv), name="Cl" )
covEl     <- mxMatrix( type="Full", nrow=nv, ncol=nf, free=TRUE, values=rep(0.4,nv), label=labDiag("El",nv), name="El" )

# Matrices As, Ds, and Es for variance components for specific factors
covAs     <- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=0.1, labels=labDiag("As",nv), name="As" )
covCs     <- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=0.1, labels=labDiag("Cs",nv), name="Cs" )
covEs     <- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=0.2, labels=labDiag("Es",nv), name="Es" )

# Matrices A, D, and E compute variance components
covA      <- mxAlgebra( expression=( Al %*% t(Al)  + As), name="A" )
covC      <- mxAlgebra( expression=( Cl %*% t(Cl)  + Cs), name="C" )
covE      <- mxAlgebra( expression=( El %*% t(El)  + Es), name="E" )

# Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
#0.25xD for ADE model
covP      <- mxAlgebra( expression= A+C+E, name="V" )
covMZ     <- mxAlgebra( expression= A+C, name="cMZ" )
#covDZ     <- mxAlgebra( expression= 0.5%x%A+ 0.25%x%D, name="cDZ" ) # for dominance
covDZ     <- mxAlgebra( expression= 0.5%x%A+ C, name="cDZ" )
expCovMZ  <- mxAlgebra( expression= rbind( cbind(V, cMZ), 
                                           cbind(t(cMZ), V)), name="expCovMZ" )
expCovDZ  <- mxAlgebra( expression= rbind( cbind(V, cDZ), cbind(t(cDZ), V)), name="expCovDZ" )

# Create Data Objects for Multiple Groups
dataMZ <- mxData( observed=mzData, type="raw" )
dataDZ <- mxData( observed=dzData, type="raw" )

# Create Expectation Objects for Multiple Groups
expMZ     <- mxExpectationNormal( covariance="expCovMZ", means="expMeanG", dimnames=selVars  )
expDZ     <- mxExpectationNormal( covariance="expCovDZ", means="expMeanG", dimnames=selVars  )
funML     <- mxFitFunctionML()

# Create Model Objects for Multiple Groups

pars      <- list(covA, covC, covE, covAl, covCl, covEl,
                  covP, covAs, covCs, covEs)
modelMZ   <- mxModel( name="MZ", inclusions, pars, expMeanG, covMZ, covDZ, expCovMZ, dataMZ, expMZ, funML )
modelDZ   <- mxModel( name="DZ", inclusions, pars, expMeanG, covDZ, expCovDZ, dataDZ, expDZ, funML )
multi     <- mxFitFunctionMultigroup( c("MZ","DZ") )

# Build & Run Model 
modelIP   <- mxModel( "mulIP", pars, modelMZ, modelDZ, multi )

check_params <- omxGetParameters(modelIP)
length(names(check_params))

fitIP <- mxRun(modelIP, intervals = F)
#fitIP <- mxTryHard(modelIP, extraTries = 10,bestInitsOutput=T, intervals=FALSE)

summary(fitIP)

fitIP$output$matrices

# Run AE model (fix both general C and specific C to zero)
modelIP_AE <- mxModel( fitIP, name="IP_AE" )
modelIP_AE <- omxSetParameters( modelIP_AE, labels=c(fitIP$Cl$labels,
                                                     diag(fitIP$Cs$labels)), free=FALSE, values=0 )

# Check number of parameters
check_params <- omxGetParameters(modelIP_AE)
length(names(check_params))

# Run model
tick <- proc.time()[3]
fitIP_AE <- mxRun(modelIP_AE, intervals = F)
#fitIP_AE <- mxTryHard( modelIP_AE, extraTries = 10,bestInitsOutput=T, intervals=FALSE)
tock-tick

mxCompare(fitIP, fitIP_AE)

summary(fitIP_AE)

fitIP_AE$output$matrices

save.image("/data/natalia/MHC/IPath_env_20251124.RData")
saveRDS(fitIP_AE, file = "/data/natalia/MHC/AE_iP_model.rds")

