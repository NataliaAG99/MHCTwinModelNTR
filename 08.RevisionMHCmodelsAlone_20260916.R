#########################
#
# MHC heritability study
# 
# Revision
#
# 10/07/2026; R version 4.5.1
#
#########################
# Preparing the environment
## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
rm(list=ls(all=TRUE)) #Delete every previous object in the environment


#### Install and load packages ####
## First specify the packages of interest
packages = c("OpenMx", "psych", "umx", "haven", "readxl", "dplyr", "ggplot2",
             "lme4", "grid", "reshape2", "here", "moments", "foreign", "gee", "tidyr",
             "ISLR", "polycor", "tidyverse", "broom", "mosaic", "ggcorrplot", "corrplot")

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


#Cor matrix among the MHC items
# Subset your data (replace `df` with your actual dataframe name)
mhc_items <- twin_df[, c("mhca14", "mhcb14", "mhcc14", "mhcd14", "mhce14", "mhcf14", "mhcg14", "mhch14", "mhci14", "mhcj14", "mhck14", "mhcl14", "mhcm14", "mhcn14")]
MHC_describe <- describe(twin_df[, c("mhca14", "mhcb14", "mhcc14", "mhcd14", "mhce14", "mhcf14", "mhcg14", "mhch14", "mhci14", "mhcj14", "mhck14", "mhcl14", "mhcm14", "mhcn14")])
library(openxlsx)
write.xlsx(as.data.frame(MHC_describe), "/revision/SummaryMHC_items.xlsx")
rm(MHC_describe)

# Compute correlation matrix (use pairwise complete obs to handle missing data)
PCor_items <- cor(mhc_items, use = "pairwise.complete.obs")

# Plot
layo <- matrix(1, 1)
layout(layo)
pdf("MHCItemsCorrelationMatrix_PDF.pdf", width = 9, height = 8)
corrplot(as.matrix(PCor_items), method = "circle",
         type = "upper",
         col = colorRampPalette(c("#5D0201", "#F28805", "#FEEB2D"))(200),
         diag = FALSE,
         tl.col = "black",
         tl.pos = "lt",
         tl.srt = 45,
         tl.offset = 0.5,
         cl.pos = "r",
         cl.offset = 0.15,
         cl.align.text = "l",
         addCoef.col = "black",
         number.cex = 0.7)   # smaller than 1 since you now have 14 items, not 3
dev.off()

# MHC is the average of the items
# Define subscale item groups
mhc_e_vars <- c("mhca14", "mhcb14", "mhcc14")                 # Emotional well-being
mhc_s_vars <- c("mhcd14", "mhce14", "mhcf14", "mhcg14", "mhch14") # Social well-being
mhc_p_vars <- c("mhci14", "mhcj14", "mhck14", "mhcl14", "mhcm14", "mhcn14") # Psychological well-being

#Correlations within and between scales
all_vars <- colnames(PCor_items)

subscales <- list(
  Emotional     = mhc_e_vars,
  Social        = mhc_s_vars,
  Psychological = mhc_p_vars
)

subscale_cor_summary <- do.call(rbind, lapply(names(subscales), function(nm) {
  inside  <- subscales[[nm]]
  outside <- setdiff(all_vars, inside)
  
  # Within: unique pairs among the subscale's own items
  sub_in <- PCor_items[inside, inside]
  within_r <- mean(sub_in[lower.tri(sub_in)])
  
  # Between: every correlation between a subscale item and a non-subscale item
  outside_r <- mean(PCor_items[inside, outside])
  
  data.frame(
    Subscale = nm,
    Mean_r_within = within_r,
    N_within = choose(length(inside), 2),
    Mean_r_inside_vs_outside = outside_r,
    N_inside_vs_outside = length(inside) * length(outside)
  )
}))

subscale_cor_summary[, grepl("Mean_r", names(subscale_cor_summary))] <-
  round(subscale_cor_summary[, grepl("Mean_r", names(subscale_cor_summary))], 3)

print(subscale_cor_summary)


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

# # Scale the variables
# twin_df$qol_sc <- scale(twin_df$levenc14)
# twin_df$sat_sc <- scale(twin_df$sat14)
# twin_df$hap_sc <- scale(twin_df$hap14)
# twin_df$flo_sc <- scale(twin_df$flo14)
# twin_df$mhc_e_sc <- scale(twin_df$mhc_e_14)
# twin_df$mhc_s_sc <- scale(twin_df$mhc_s_14)
# twin_df$mhc_p_sc <- scale(twin_df$mhc_p_14)

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
                                        flo14))

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
data_twins <- data_twins[complete.cases(data_twins[,c('mhc_e_14', 'mhc_s_14', 'mhc_p_14')]),]
#
#cvd1
#
rph1 = summary(lm(mhc_e_14~Age_sc + sex, data=data_twins))
rph2 = summary(lm(mhc_s_14~Age_sc + sex, data=data_twins))
rph3 = summary(lm(mhc_p_14~Age_sc + sex, data=data_twins))

R2s = c(rph1$r.square,rph2$r.square,rph3$r.square) 
parsb0 = c(
  rph1$coefficients[1,1],
  rph2$coefficients[1,1],
  rph3$coefficients[1,1])
parsbage = c( 
  rph1$coefficients[2,1],
  rph2$coefficients[2,1],
  rph3$coefficients[2,1])
parsbsex = c( 
  rph1$coefficients[3,1],
  rph2$coefficients[3,1],
  rph3$coefficients[3,1])
# cvd2
#
#
length(unique(data_twins$FamilyNumber)) #4207 families, 5248 i

table(as.data.frame(table(data_twins$FamilyNumber))[,2]) #Number of twins per family ## 2 --> 1041 complete pairs

table(data_twins$Extension)

# Renaming the columns for clarity
colnames(data_twins) <- c("FISNumber", "FamilyNumber", "Extension", "Sex", "Zyg", "Age", "EMHC", "SMHC", "PMHC", "Age_sc")

# how many individuals of each zygosity?
zyg_frequencies <- table(data_twins[,"Zyg"])
zygosities <- c("MZM", "DZM", "MZF", "DZF", "DOSmf", "DOSfm")
rownames(zyg_frequencies) <- zygosities
zyg_frequencies

# MZM   DZM   MZF   DZF DOSmf DOSfm 
# 700   394  2127   925   562   540

describe(data_twins[, c("Age", "Sex", "EMHC", "SMHC", "PMHC")]) 
#n=6 with a missing sex left

subset(data_twins, is.na(data_twins$Sex)) #We can assign sex to these 6 individuals based on zygosity
data_twins$Sex[is.na(data_twins$Sex) & data_twins$Zyg %in% c(1, 2)] <- 0
data_twins$Sex[is.na(data_twins$Sex) & data_twins$Zyg %in% c(3, 4)] <- 1
data_twins$Sex[is.na(data_twins$Sex) & data_twins$Zyg %in% c(5, 6)] <- NA
#
#Cross-variable correlations
cor.test(data_twins$EMHC, data_twins$SMHC) #0.5735568 
cor.test(data_twins$EMHC, data_twins$PMHC) #0.6738831 
cor.test(data_twins$SMHC, data_twins$PMHC) #0.7205525 
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
(EMHC_mz=cor.test(dataMZ$EMHC_T1,dataMZ$EMHC_T2)) #0.3095693  
(EMHC_dz=cor.test(dataDZ$EMHC_T1,dataDZ$EMHC_T2)) #0.06001571  

(SMHC_mz=cor.test(dataMZ$SMHC_T1,dataMZ$SMHC_T2)) #0.394007  
(SMHC_dz=cor.test(dataDZ$SMHC_T1,dataDZ$SMHC_T2)) #0.1452804  

(PMHC_mz=cor.test(dataMZ$PMHC_T1,dataMZ$PMHC_T2)) #0.368907  
(PMHC_dz=cor.test(dataDZ$PMHC_T1,dataDZ$PMHC_T2)) #0.1166189  

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
vars      <- c('EMHC', 'SMHC', 'PMHC')              # list of variables names
nv        <- 3                         # number of variables
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
svMe      <- c(3)                   # start value for means
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
# cvd2
#pathB1     <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=.01, 
#              labels=c("bAge_QoL_T1","bAge_Sat_T1", "bAge_Hap_T1","bAge_Flo_T1", "bAge_Emp_T1", 
#                      "bAge_QoL_T2","bAge_Sat_T2", "bAge_Hap_T2","bAge_Flo_T2", "bAge_Emp_T2"), name="b1" )
# cvd1
pathB1     <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=c(T,T,T), # TRUE, 
                        values=rep(0,nv), # c(parsbage), 
                        labels=c("bAge_EMHC", "bAge_SMHC", "bAge_PMHC"), name="b1" )
# cvd2
# cvd1
defSex      <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=FALSE,
                         labels=c('data.Sex_T1','data.Sex_T2'), name="Sex" )
# cvd2
#pathB2     <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=.2, 
#                labels=c("bSex_QoL_T1","bSex_Sat_T1", "bSex_Hap_T1","bSex_Flo_T1", "bSex_Emp_T1", 
#                        "bSex_QoL_T2","bSex_Sat_T2", "bSex_Hap_T2","bSex_Flo_T2", "bSex_Emp_T2"), name="b2" ) # Sex effects
#
pathB2     <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=c(F,T,T), 
                        values=c(0, 0.05,0.05), # c(parsbsex), 
                        labels=c("bSex_EMHC", "bSex_SMHC", "bSexPMHC"), name="b2" ) # Sex effects

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
#fitGofs(fitAE); fitEstCis(fitAE)

mxCompare( fitACE, fitAE )

summary(fitAE)

library(openxlsx)
write.xlsx(as.data.frame(summary(fitAE)$parameters), "/data/natalia/MHC/revision/SummaryAE_MHConly.xlsx")
write.xlsx(as.data.frame(summary(fitACE)$parameters), "/data/natalia/MHC/revision/SummaryACE_MHConly.xlsx")
write.xlsx(as.data.frame(mxCompare( fitACE, fitAE )), "/data/natalia/MHC/revision/ComparisonAE_MHConly.xlsx")

fitAE$algebras$UV
fitAE$algebras$SV

fitAE$algebras$rA
fitAE$algebras$rE
fitAE$algebras$rP

#If variance MZ is lower than DZ (massive) it implies GxE interaction
save.image("/data/natalia/MHC/revision/MHCrev_env_ACEAE_20260916.RData")
saveRDS(fitAE, file = "/data/natalia/MHC/revision/AE_MHConly_model.rds")

ACor_df <- as.data.frame(fitAE$algebras$rA@result)
rownames(ACor_df) <- c("EMHC", "SMHC", "PMHC")  
colnames(ACor_df) <- c("EMHC", "SMHC", "PMHC")  

PCor_df <- as.data.frame(fitAE$algebras$rP@result)
rownames(PCor_df) <- c("EMHC", "SMHC", "PMHC")  
colnames(PCor_df) <- c("EMHC", "SMHC", "PMHC")  

ECor_df <- as.data.frame(fitAE$algebras$rE@result)
rownames(ECor_df) <- c("EMHC", "SMHC", "PMHC")  
colnames(ECor_df) <- c("EMHC", "SMHC", "PMHC")  

#######################################################################################
# Independent pathway
#######################################################################################

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

mxCompare(fitACE, fitIP)

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
fitIP_AE <- mxRun(modelIP_AE, intervals = T)
#fitIP_AE <- mxTryHard( modelIP_AE, extraTries = 10,bestInitsOutput=T, intervals=FALSE)

mxCompare(fitIP, fitIP_AE)
#mxCompare(fitIP_AE, fitCP_AE)
write.xlsx(as.data.frame(mxCompare( fitIP, fitIP_AE )), "/data/natalia/MHC/revision/ComparisonIP_MHConly.xlsx")

summary(fitIP_AE)

write.xlsx(as.data.frame(summary(fitIP)$parameters), "/data/natalia/MHC/revision/SummaryIP_MHConly.xlsx")
write.xlsx(as.data.frame(summary(fitIP_AE)$parameters), "/data/natalia/MHC/revision/SummaryAEIP_MHConly.xlsx")

fitIP_AE$output$matrices

mxCompare(fitIP_AE, fitCP_AE) # ambiguous because p value says CP, AIC says IP
mxCompare(fitIP_AE, fitAE) # cannot be determined
mxCompare(fitAE, fitCP_AE) # ambiguous because p value says CP, AIC says IP

write.xlsx(as.data.frame(mxCompare( fitIP, fitACE )), "/data/natalia/MHC/revision/ComparisonIP_ACE_MHConly.xlsx")
write.xlsx(as.data.frame(mxCompare( fitAE, fitIP_AE )), "/data/natalia/MHC/revision/ComparisonIPAE_AE_MHConly.xlsx")

save.image("/data/natalia/MHC/revision/MHCrev_env_IPath_20260916.RData")
saveRDS(fitIP_AE, file = "/data/natalia/MHC/revision/AE_iP_MHConly_model.rds")

########################################################################################

# Fit Common Pathway ACE Model

svPa <- .4 # start value for path coefficient
svPe <- .5 # start value for path coefficient for e

# Number of latent variables (common MHC)
nl        <- 1

# Create Labels
labMe     <- paste("mean",vars,sep="_")
# ------------------------------------------------------------------------------
# PREPARE MODEL

# ACE Model

# Create Matrices for Variance Components latent variables
covA <- mxMatrix( type="Symm", nrow=nl, ncol=nl, free=TRUE, values=valDiag(svPa,nl), label=labLower("VA",nl), name="VA" )
covC <- mxMatrix( type="Symm", nrow=nl, ncol=nl, free=TRUE, values=valDiag(svPa,nl), label=labLower("VC",nl), name="VC" )
covE <- mxMatrix( type="Symm", nrow=nl, ncol=nl, free=TRUE, values=valDiag(svPa,nl), label=labLower("VE",nl), name="VE" )

# Matrix and Algebra for constraint on variance of latent phenotype
covarLP <- mxAlgebra( expression= VA + VC + VE, name="CovarLP" )
varLP <- mxAlgebra( expression= diag2vec(CovarLP), name="VarLP" )
unit <- mxMatrix( type="Unit", nrow=nl, ncol=1, name="Unit")
varLP1 <- mxConstraint( expression=VarLP == Unit, name="varLP1")

# Create Matrices for Variance Components for specific factors
covAs <- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPa,nv), label=labDiag("VAs",nv), name="VAs" )
covCs <- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPa,nv), label=labDiag("VCs",nv), name="VCs" )
covEs <- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPa,nv), label=labDiag("VEs",nv), name="VEs" )


svfl      <- matrix(c(1,0.5,0.5),3,1)

frfl      <- c(TRUE, TRUE, TRUE)
# Matrix for factor loadings on latent phenotype
pathFl    <- mxMatrix( type="Full", nrow=nv, ncol=nl, free=frfl, 
                       values=svfl, labels=labFull("fl",nv,nl), name="fl" )

# Matrices A, C, and E compute variance components of domain indicators
covAd <- mxAlgebra( expression = fl %&% VA + VAs, name="A" )
covCd <- mxAlgebra( expression = fl %&% VC + VCs, name="C" )
covEd <- mxAlgebra( expression = fl %&% VE + VEs, name="E" )

# Matrices for standardized common variance components
covSAc      <- mxAlgebra( expression=VA /(VA+VC+VE), name="SAc" )
covSCc      <- mxAlgebra( expression=VC /(VA+VC+VE), name="SCc" )
covSEc      <- mxAlgebra( expression=VE /(VA+VC+VE), name="SEc" )

# Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
covP      <- mxAlgebra( expression= A+C+E, name="V" )
covMZ     <- mxAlgebra( expression= A+C, name="cMZ" )
#covDZ     <- mxAlgebra( expression= 0.5%x%A+ 0.25%x%D, name="cDZ" ) # for dominance
covDZ     <- mxAlgebra( expression= 0.5%x%A+ C, name="cDZ" )
expCovMZ <- mxAlgebra( expression= rbind( cbind(V, cMZ), cbind(t(cMZ), V)), name="expCovMZ" )
expCovDZ <- mxAlgebra( expression= rbind( cbind(V, cDZ), cbind(t(cDZ), V)), name="expCovDZ" )


# Create Data Objects for Multiple Groups
dataMZ <- mxData( observed=mzData, type="raw" )
dataDZ <- mxData( observed=dzData, type="raw" )

# Create Expectation Objects for Multiple Groups
expMZ     <- mxExpectationNormal( covariance="expCovMZ", means="expMeanG", dimnames=selVars)
expDZ     <- mxExpectationNormal( covariance="expCovDZ", means="expMeanG", dimnames=selVars)
funML     <- mxFitFunctionML()

# Create Model Objects for Multiple Groups
pars      <- list(pathFl, covA, covC, covE, covAd, covCd, covEd, 
                  covP, covAs, covCs, covEs, covSAc, covSCc, covSEc, covarLP, varLP, unit)
modelMZ   <- mxModel( name="MZ", inclusions, pars, covMZ, expCovMZ, dataMZ, expMZ, funML )
modelDZ   <- mxModel( name="DZ", inclusions, pars, covDZ, expCovDZ, dataDZ, expDZ, funML )
multi     <- mxFitFunctionMultigroup( c("MZ","DZ") )

# Build & Run Model 
modelCP   <- mxModel( "mulCP_SWLc", pars, modelMZ, modelDZ, multi, varLP1 )

check_params <- omxGetParameters(modelCP)
length(names(check_params))
# always good to check whether you estimate what you expect


fitCP <- mxRun(modelCP, intervals = F)
#fitCP <- mxTryHard(modelCP, extraTries = 10,bestInitsOutput=T, intervals=FALSE)


summary(fitCP)

# Run AE model (fix both general C and specific C to zero)
modelCP_AE <- mxModel( fitCP, name="CP_AE" )
modelCP_AE <- omxSetParameters( modelCP_AE, labels=c(unique(c(fitCP$VC$labels)) ,
                                                     diag(fitCP$VCs$labels)), free=FALSE, values=0 )

# Check number of parameters
check_params <- omxGetParameters(modelCP_AE)
length(names(check_params))

0.579198436-0.135412936+0.556214500
0.45035101381 + 0.58442318579

# Set options
#modelCP_AE <- mxOption(model= modelCP_AE, key="Number of Threads", value= 32 )

# Run model
fitCP_AE <- mxRun(modelCP_AE, intervals = T)
#fitCP_AE <- mxTryHard( modelCP_AE, extraTries = 10,bestInitsOutput=T, intervals=FALSE)


summary(fitCP_AE)

write.xlsx(as.data.frame(summary(fitCP)$parameters), "/data/natalia/MHC/revision/SummaryCP_MHConly.xlsx")
write.xlsx(as.data.frame(summary(fitCP_AE)$parameters), "/data/natalia/MHC/revision/SummaryAECP_MHConly.xlsx")

mxCompare(fitCP, fitCP_AE)
write.xlsx(as.data.frame(mxCompare( fitCP, fitCP_AE )), "/data/natalia/MHC/revision/ComparisonCP_MHConly.xlsx")

# LATENT LEVEL (MAIN RESULTS)
Tot_var <- fitCP_AE$VA$values+fitCP_AE$VE$values
fitCP_AE$VA$values/Tot_var
fitCP_AE$VE$values/Tot_var
cov2cor(fitCP_AE$VA$values)
cov2cor(fitCP_AE$VE$values)

# Calculate standardized first loading
Lat_var <- Tot_var[1,1] 
sqrt(Lat_var) / sqrt(Lat_var + diag(Vs)[1])

write.xlsx(as.data.frame(mxCompare( fitIP, fitCP )), "/data/natalia/MHC/revision/ComparisonIP_CP_MHConly.xlsx")
write.xlsx(as.data.frame(mxCompare( fitAE, fitCP_AE )), "/data/natalia/MHC/revision/ComparisonAE_CPAE_MHConly.xlsx")
write.xlsx(as.data.frame(mxCompare( fitACE, fitCP )), "/data/natalia/MHC/revision/ComparisonACE_CP_MHConly.xlsx")
write.xlsx(as.data.frame(mxCompare( fitIP_AE, fitCP_AE )), "/data/natalia/MHC/revision/ComparisonIPAE_CPAE_MHConly.xlsx")

save.image("/data/natalia/MHC/revision/MHCrev_env_CPath_20260916.RData")
saveRDS(fitCP_AE, file = "/data/natalia/MHC/revision/AE_cP_MHConly_model.rds")
