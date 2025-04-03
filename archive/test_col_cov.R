popmodel <- 'eta1 =~ 1*x1 + 1.5*x2 + 2*x3
              x1 ~~ 3 * x1 
              x2 ~~ 6.75 * x2
              x3 ~~ 12 * x3
              
              eta1~~1*eta1
              
              eta2 =~ 1*z2 + 1.5*z3 + 2*z1
              z2 ~~ 3 * z2 
              z3 ~~ 6.75 * z3
              z1 ~~ 12 * z1
              
              eta2~~1*eta2
              
              eta1~~0.5*eta2

'



# parallel
popmodel <- 'eta1 =~ 1.5*x1 + 1.5*x2 + 1.5*x3
              x1 ~~ 2 * x1 
              x2 ~~ 2 * x2
              x3 ~~ 2 * x3
              
              eta1~~1*eta1
              
              eta2 =~ 2*z2 + 2*z3 + 2*z1
              z2 ~~ 1 * z2 
              z3 ~~ 1 * z3
              z1 ~~ 1 * z1
              
              eta2~~1*eta2
              
              eta1~~0.5*eta2

'



# tau
popmodel <- 'eta1 =~ 1.5*x1 + 1.5*x2 + 1.5*x3
              x1 ~~ 2.3 * x1 
              x2 ~~ 2.7 * x2
              x3 ~~ 2 * x3
              
              eta1~~1*eta1
              
              eta2 =~ 2*z2 + 2*z3 + 2*z1
              z2 ~~ 1.8 * z2 
              z3 ~~ 0.1 * z3
              z1 ~~ 2 * z1
              
              eta2~~1*eta2
              
              eta1~~0.5*eta2

'


mod <- 'eta1 =~ NA*x1 +x2 +x3
            eta1~~1*eta1          
              eta2 =~ NA*z2 + z3 + z1
            eta2~~1*eta2  
'

library(lavaan)
datapop=simulateData(model = popmodel,empirical = T)
out=sem(model = modelest,data = datapop)
summary(out,standardized=T)

cor(datapop)
semTools::htmt(out)

modelhtmt <- 'eta1 =~ x1 +x2 +x3
            eta1~~1*eta1          
              eta2 =~ z2 + z3 + z1
            eta2~~1*eta2  
'

calc_htmt(datapop,model = modelhtmt,latent1='eta1',latent2='eta2',htmt2 = T,scale = T)
