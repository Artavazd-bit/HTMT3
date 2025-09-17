library(dplyr)
library(ggplot2)
library(patchwork)
df <- read.csv2("simresults/confidenceintervals/confruns/HTMTsimresults100200400.csv")
df2 <- read.csv2("simresults/confidenceintervals/confruns/HTMTsimresults50800.csv")
df3 <- read.csv2("simresults/confidenceintervals/confruns/HTMTsimresults251600.csv")

dfall <- rbind(df, df2, df3)

rm(df, df2, df3)

missing <- dfall[!is.na(dfall$missing) & dfall$alpha == 0.05, ]

summary(missing)

missingfilt <- missing[,c("seed", "method", "missing", "sim_runs", "datatype", "correlation", "n")]


ggplot(missingfilt[missingfilt$datatype == "normal" & missingfilt$method == "boot" & missingfilt$missing != 0, ], aes(x = missing)) + 
  geom_bar() + 
  facet_grid(cols = vars(correlation), rows  = vars(n)) +  
  theme(legend.position = "bottom") 

missingfilt$missingint[missingfilt$missing == 0] <- "None"
missingfilt$missingint[1 <= missingfilt$missing& missingfilt$missing <10] <- "[1,10)"
missingfilt$missingint[10 <= missingfilt$missing& missingfilt$missing <20] <- "[10,20)"
missingfilt$missingint[20 <= missingfilt$missing& missingfilt$missing <30] <- "[20,30)"
missingfilt$missingint[30 <= missingfilt$missing& missingfilt$missing <40] <- "[30,40)"
missingfilt$missingint[40 <= missingfilt$missing& missingfilt$missing <50] <- "[40,50)"
missingfilt$missingint[50 <= missingfilt$missing& missingfilt$missing <60] <- "[50,60)"
missingfilt$missingint[60 <= missingfilt$missing& missingfilt$missing <70] <- "[60,70)"
missingfilt$missingint[70 <= missingfilt$missing& missingfilt$missing <80] <- "[70,80)"
missingfilt$missingint[80 <= missingfilt$missing& missingfilt$missing <90] <- "[80,90)"
missingfilt$missingint[90 <= missingfilt$missing& missingfilt$missing <100] <- "[90,100)"
missingfilt$missingint[100 <= missingfilt$missing& missingfilt$missing <110] <- "[100,110)"
missingfilt$missingint[110 <= missingfilt$missing& missingfilt$missing <120] <- "[110,120)"
missingfilt$missingint[120 <= missingfilt$missing& missingfilt$missing <130] <- "[120,130)"
missingfilt$missingint[130 <= missingfilt$missing& missingfilt$missing <140] <- "[130,140)"
missingfilt$missingint[140 <= missingfilt$missing& missingfilt$missing <150] <- "[140,150)"
missingfilt$missingint[150 <= missingfilt$missing& missingfilt$missing <160] <- "[150,160)"
missingfilt$missingint[160 <= missingfilt$missing& missingfilt$missing <170] <- "[160,170)"
missingfilt$missingint[170 <= missingfilt$missing& missingfilt$missing <180] <- "[170,180)"

missingfilt$missingint <- ordered(missingfilt$missingint, levels=c("None", "[1,10)", "[10,20)", "[20,30)",  "[30,40)", "[40,50)", "[50,60)"
                                                                  , "[60,70)", "[70,80)", "[80,90)", "[90,100)", "[100,110)", "[110,120)", 
                                                                  "[120,130)", "[130,140)", "[140,150)", "[150,160)", "[160,170)", "[170,180)"))

missingfilttest <- missingfilt %>% 
  filter(method == "boot") %>%
  group_by(correlation, missingint, datatype, method) %>%
  select(n) %>%
  table()


ggplot(missingfilt[missingfilt$datatype == "normal" & missingfilt$method == "boot" & missingfilt$missing != 0, ], aes(x = missingint)) + 
  geom_bar() + 
  facet_grid(cols = vars(correlation), rows  = vars(n)) +  
  theme(legend.position = "bottom") 


missingfiltwithou0 <- missingfilt[missingfilt$missing != 0, ]
ggplot(missingfilt[missingfilt$datatype == "normal" & missingfilt$method == "boot" , ], aes(x = missing)) + 
  stat_bin(breaks = c(1, seq(from = 10, to = 180, by =10)), aes(y = after_stat(density*width))) + 
  facet_grid(cols = vars(correlation), rows  = vars(n)) +  
  theme(legend.position = "bottom") + 
  scale_x_continuous(breaks = c(2, seq(from = 10, to = 180, by =10)))

ggplot(missingfilt[missingfilt$datatype == "normal" & missingfilt$method == "boot" & missingfilt$missing != 0, ], aes(x = missing)) + 
  geom_histogram(aes(y= pr)) + 
  facet_grid(cols = vars(correlation), rows  = vars(n)) +  
  theme(legend.position = "bottom") + 
  scale_x_continuous(breaks = c(1, seq(from = 10, to = 30, by =10)))


test <- missingfilt %>% 
  group_by(correlation, n, datatype, method) %>%
  table()


head(missingfilt)

################################################################################

test <- ggplot(missingfilt[missingfilt$datatype == "normal" & missingfilt$n == 25 & missingfilt$method == "boot" & missingfilt$missing != 0, ], aes(x = missing)) + 
  geom_histogram(aes(y=after_stat(count/10000)), breaks = c(1, seq(from = 10, to = 180, by =10))) + 
  facet_grid(cols = vars(correlation), rows  = vars(n)) +  
  theme(legend.position = "bottom") + 
  scale_y_continuous(name = "Proportion") + 
  scale_x_continuous(name = "Inadmissible solutions")

ggplot(missingfilt[missingfilt$datatype == "normal" & missingfilt$method == "bcaboot" &missingfilt$missing != 0, ], aes(x = missing)) + 
  geom_histogram(aes(y=after_stat(count/10000))) + 
  facet_grid(cols = vars(correlation), rows  = vars(n)) +  
  theme(legend.position = "bottom") + 
  scale_y_continuous(name = "Proportion") + 
  scale_x_continuous(name = "Inadmissible solutions")



table <- missingfilt[missingfilt$datatype == "normal" & missingfilt$method == "boot" & missingfilt$missing != 0, ]

missingfiltbca[missingfiltbca$missing == max(missingfiltbca$missing), ]

missingfilt[missingfilt$missing == max(missingfilt$missing), ]

psych::describe(missingfilt)

psych::describe(missingfilt$missing)

psych::describe(missingfiltbca$missing)


test <- missingfilt[missingfilt$missing != 0,]

psych::describe(test)

nrow(missingfiltdist) 



