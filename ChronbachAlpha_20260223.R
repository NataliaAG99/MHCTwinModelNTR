# Preparing the environment
## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
rm(list=ls(all=TRUE)) #Delete every previous object in the environment


#### Install and load packages ####
## First specify the packages of interest

packages = c("foreign", "lavaan")

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
setwd("/data/natalia/MHC/")

#########################
# Load the data
dades <- read.spss("PHE_20260130_5055_LL.sav", to.data.frame=T, use.value.labels = F) #use.value.labels=F
#describe(dades)
summary(dades)
colnames(dades)

twin_df <- dades[!is.na(dades$twzyg),] #6517i

#Remove the multiples
twin_df <- subset(twin_df, !Extension %in% c(3, 4, 5))

# Select the relevant columns for the twin data
data_twins <- subset(twin_df, select=c(FISNumber, FamilyNumber, welba14,
                                       welbb14, welbc14, welbd14, welbe14, 
                                       welbf14, welbg14, welbh14, welbi14, 
                                       flo_a14, flo_b14, flo_c14, flo_d14, 
                                       flo_e14, flo_f14, flo_g14, flo_h14, 
                                       levenc14, mhca14, mhcb14, mhcc14, mhcd14, 
                                       mhce14, mhcf14, mhcg14, mhch14, mhci14, 
                                       mhcj14, mhck14, mhcl14, mhcm14, mhcn14))

# # Remove the original twin_df to save memory
rm(twin_df) 
rm(dades)
#
#
sum(is.na(data_twins$twzyg))
#
#Keep individuals in the other dataset
fisnumbersMHC <- readRDS(file = "fisnumbersMHC.rds")
data_twins <- data_twins[data_twins$FISNumber %in% fisnumbersMHC, ]
rm(fisnumbersMHC)


#Make all variables in a way that the higher the score, the more positive the response
#Variables/questions g and i
#Survey 8
#welbgX is SHS2 and welbiX is SHS4
VarToInv <- c('welbg14', 'welbi14')
# table(dades$welbg14, useNA="ifany")
#Define function to invert the score
inv <- function(x){ ((x - max(x, na.rm = TRUE)) * -1) + min(x, na.rm = TRUE) }

# Loop through each variable in VarToInv and apply inv() function
for (var in VarToInv) {
  data_twins[[paste0(var, "_inv")]] <- inv(data_twins[[var]])
}


#RELIABILITY
names(data_twins)

library(psych)

#Happiness
psych:::alpha(data_twins[,c("welbf14", "welbg14_inv", "welbh14", "welbi14_inv")])

#SWLS
psych:::alpha(data_twins[,c("welba14", "welbb14", "welbc14", "welbd14", "welbe14")])

#Flourishing
psych:::alpha(data_twins[,c("flo_a14", "flo_b14", "flo_c14", "flo_d14", "flo_e14", "flo_f14", "flo_g14", "flo_h14")])

#MHC
psych:::alpha(data_twins[,c("mhca14", "mhcb14", "mhcc14", "mhcd14", "mhce14", "mhcf14", "mhcg14", "mhch14", "mhci14", "mhcj14", "mhck14", "mhcl14", "mhcm14", "mhcn14")])

#EMHC
psych:::alpha(data_twins[,c("mhca14", "mhcb14", "mhcc14")])

#SMHC
psych:::alpha(data_twins[,c("mhcd14", "mhce14", "mhcf14", "mhcg14", "mhch14")])

#MHC
psych:::alpha(data_twins[,c("mhci14", "mhcj14", "mhck14", "mhcl14", "mhcm14", "mhcn14")])
