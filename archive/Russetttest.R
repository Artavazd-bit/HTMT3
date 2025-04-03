library(lavaan)
library(semTools)
library(cSEM)
library(stringr)
library(doParallel)
library(foreach)
library(dplyr)
library(boot)

source("C:/Forschung/HTMT3/2025_19_02_models_3_4.R")
source("C:/Forschung/HTMT3/2024_01_08_functions.R")
source("C:/Forschung/HTMT3/2024_01_10_gradient_analytically_of_htmt.R")
Russett <- as.data.frame(readxl::read_excel("C:/Forschung/Bootstrap-versus-Delta/01_Data/Russett.xlsx"))

model_Russett = ' # Specify composite models
              AgrIneq =~ gini + farm +rent
              IndDev  =~ gnpr + labo
              PolInst =~ inst + ecks + deat + stab + dict 
              
              # Specify the relation among 
              # the emergent variables
              PolInst ~ AgrIneq + IndDev
              '


#### Beispiel Russet mit meinem HTMT ###########

out=sem(model_Russett,data = Russett1)
summary(out)


Russett1[,"gnpr"]<- Russett[,"gnpr"]*-1

cor(Russett)
cor(Russett1)

vc_r <- calculate_cov_cov(data = Russett)

gradient_htmt_1 <- calc_grad_htmt_ana(data = Russett1, model = model_Russett, latent1 = "AgrIneq", latent2 = "IndDev", scale = FALSE)

Gradient_htmt <- as.matrix(gradient_htmt_1$output$gradient)

se_htmt_1 = sqrt(t(Gradient_htmt) %*% vc_r %*% Gradient_htmt / n)

z_value_htmt_1 = (gradient_htmt_1$HTMT - 1)/se_htmt_1
