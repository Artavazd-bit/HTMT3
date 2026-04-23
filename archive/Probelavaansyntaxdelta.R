library(lavaan)
library(doParallel)
library(foreach)
library(dplyr)

source("setup.R")

data <- lavaan::simulateData(model = simModels$model[1],
                             sample.nobs = 10000, # Number of observations.
                             seed = 1234123, # Set random seed.
                             empirical = FALSE, # Logical. If TRUE, the implied moments (Mu and Sigma) specify the empirical not population mean and covariance matrix.
                             return.type = "data.frame"
)

write.table(x = data,file = choose.files(),sep = "\t",row.names = F,col.names = F)
colnames(data)=paste0('x',1:6)

model_est<- '
              #  latent variables
                xi_1 =~ x1 + x2 + x3
                xi_2 =~ x4 + x5 + x6 
                
                xi_1 ~~ xi_2
              ' 


derivhtmt(data = data, model = model_est, latent1 = "xi_1", latent2 = "xi_2", scale = FALSE, htmt2 = FALSE)

delta <- deltamethod(data = data, model = model_est, alpha = c(0.05, 0.10), latent1 = "xi_1", latent2 = "xi_2", scale = FALSE, htmt2 = FALSE)

htmtmodel = ' 
x1~~r12*x2 + r13*x3 + r14*x4 + r15*x5 + r16*x6
x2~~r23*x3 + r24*x4 + r25*x5 + r26*x6
x3~~r34*x4 + r35*x5 + r36*x6
x4~~r45*x5 + r46*x6
x5~~r56*x6

rxx:=1/3*(r12+r13+r23)
ryy:=1/3*(r45+r46+r56)
rxy:=1/9*(r14+r15+r16+r24+r25+r26+r34+r35+r36)
HTMT:=rxy/(sqrt(ryy*rxx))

'

out=sem(model = htmtmodel,data = data, sample.cov.rescale = TRUE)
summary(out)

fit_cfa <- lavaan::cfa(model  = model_est, 
                       data = data)

lavdata <- fit_cfa@Data
lavoptions <- lavInspect(fit_cfa, "options")
lavoptions$gamma.unbiased <- FALSE
lavoptions$gamma.n.minus.one <- FALSE
lavoptions$correlation <- FALSE
lavoptions$conditional.x <- FALSE
lavoptions$fixed.x <- FALSE
lavoptions$meanstructure <- FALSE

test <- lavaan::lav_samplestats_from_data(lavdata = lavdata, lavoptions = lavoptions, NACOV = TRUE)

out1 <- sem(NACOV = test@NACOV, model = model_est)
summary(out1)



qnorm(0.95)
