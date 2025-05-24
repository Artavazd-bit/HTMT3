library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(patchwork)

res <- read.csv2("simresultsconfidence.csv")

missingcases <- res %>% filter(is.na(HTMT))
head(missingcases, n = 10)
nrow(missingcases)
table(missingcases$correlation)
table(missingcases$n)
table(missingcases$datatype)

reswona <- res %>% filter(!is.na(HTMT))

head(reswona)

reswona$upperwithin <- reswona$correlation < reswona$upperbound
reswona$lowerwithin <- reswona$correlation > reswona$lowerbound



resag <- reswona %>% 
  group_by(correlation, n, datatype, alpha, method) %>%
  summarize(upperwithin = mean(upperwithin)*100,
            lowerwithin = mean(lowerwithin)*100
  )

head(resag)

resag$method2[resag$method == "boot"] <- "CI_boot" 
resag$method2[resag$method == "delta"] <- "CI_normal" 


y_breaks <- seq(70, 100, by = 5)

p_upper_n <- ggplot(resag[resag$alpha == 0.05 & resag$datatype == "normal",], aes(x = as.factor(n), y = upperwithin, group = method)) +
  geom_line(aes(linetype = method2)) +
  geom_point(aes(shape = method2)) +
  facet_grid(cols = vars(correlation), labeller = label_parsed) + 
  geom_hline(yintercept = 97.5) + 
  scale_y_continuous(limits = c(75, 100), name = "Prop. with pop. correlation smaller upper limit (%)", breaks = y_breaks) +
  theme_minimal() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) 

p_lower_n <- ggplot(resag[resag$alpha == 0.05 & resag$datatype == "normal",], aes(x = as.factor(n), y = lowerwithin, group = method)) +
  geom_line(aes(linetype = method2)) +
  geom_point(aes(shape = method2)) +
  facet_grid(cols = vars(correlation), labeller = label_parsed, ) + 
  geom_hline(yintercept = 97.5) +
  scale_y_reverse(name = "Prop. with pop. correlation higher lower limit (%)", limits = c(100, 75),  breaks = y_breaks) +
  theme_minimal() +
  theme(legend.position = "bottom", strip.text.x = element_blank()) +
  labs(x = "Sample size") + 
  scale_linetype_discrete(name = "Type of CI:") +
  scale_shape_discrete(name = "Type of CI:")

p_upper_n / p_lower_n

##############################################################################
p_upper_nn <- ggplot(resag[resag$alpha == 0.05 & resag$datatype == "nonnormal",], aes(x = as.factor(n), y = upperwithin, group = method)) +
  geom_line(aes(linetype = method2)) +
  geom_point(aes(shape = method2)) +
  facet_grid(cols = vars(correlation), labeller = label_parsed) + 
  geom_hline(yintercept = 97.5) + 
  scale_y_continuous(limits = c(75, 100), name = "Prop. with pop. correlation smaller upper limit (%)", breaks = y_breaks) +
  theme_minimal() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) 

p_lower_nn <- ggplot(resag[resag$alpha == 0.05 & resag$datatype == "nonnormal",], aes(x = as.factor(n), y = lowerwithin, group = method)) +
  geom_line(aes(linetype = method2)) +
  geom_point(aes(shape = method2)) +
  facet_grid(cols = vars(correlation), labeller = label_parsed, ) + 
  geom_hline(yintercept = 97.5) +
  scale_y_reverse(name = "Prop. with pop. correlation higher lower limit (%)", limits = c(100, 75),  breaks = y_breaks) +
  theme_minimal() +
  theme(legend.position = "bottom", strip.text.x = element_blank()) +
  labs(x = "Sample size") + 
  scale_linetype_discrete(name = "Type of CI:") +
  scale_shape_discrete(name = "Type of CI:")

p_upper_nn / p_lower_nn
##############################################################################
nobs <- 50
method <- "delta"
correl <- 0.8
valpha <- 0.05
vdata <- "nonnormal"

df <- reswona[reswona$correlation == correl & reswona$method == method & reswona$alpha == valpha & reswona$n == nobs & reswona$datatype == vdata,]

ggplot(df, aes(x = X)) +
  geom_line(aes(y=upperbound, colour = "red")) + 
  geom_line(aes(y=lowerbound, colour = "green")) 

ggplot(df, aes(x = X)) + 
  geom_ribbon(aes(ymin = lowerbound, ymax = upperbound), fill = "lightblue", alpha = 0.5) +
  geom_line(aes(y = lowerbound), color = "red") +
  geom_line(aes(y = upperbound), color = "green")

ggplot(df, aes(x = X)) + 
  geom_smooth(aes(y = lowerbound), color = "red", se = FALSE) +
  geom_smooth(aes(y = upperbound), color = "green", se = FALSE)

ggplot(df, aes(x = X)) +
  stat_density_2d(aes(y = lowerbound, fill = after_stat(density)), geom = "raster", contour = FALSE) +
  stat_density_2d(aes(y = upperbound, fill = after_stat(density)), geom = "raster", contour = FALSE) +
  scale_fill_viridis_c()
###############################################################################


ggplot(reswona[reswona$correlation == correl & reswona$method == method & reswona$alpha == valpha & reswona$n == nobs & reswona$datatype == vdata,], aes(x = lowerbound)) + 
  geom_density(alpha = 0.9)

###############################################################################
nobs <- 200
method <- "delta"
correl <- 0.7
valpha <- 0.05
vdata <- "normal"

df <- reswona[reswona$correlation == correl & reswona$method == method & reswona$alpha == valpha & reswona$n == nobs & reswona$datatype == vdata,]

plot(density(df$HTMT))
psych::describe(df$HTMT)

psych::describe(df$HTMT)


dfss <- sample(df, size = 1000)
ggplot(df, aes(sim_runs, HTMT))+ 
  geom_line() +
   geom_ribbon(data=df,aes(ymin=lowerbound,ymax=upperbound),alpha=0.3)

ggplot(df, aes(sim_runs, HTMT))+ 
  geom_point() +
  geom_ribbon(data=df,aes(ymin=lowerbound,ymax=upperbound),alpha=0.3)
  
ggplot(df, aes(sim_runs, HTMT)) + 
  geom_point() +
  geom_smooth(method = "lm", formula = df$HTMT ~ 1) 
  

library(zoo)
# Calculate rolling means with window size
window_size <- 50
df$lower_roll_mean <- rollmean(df$lowerbound, k = window_size, fill = NA, align = "center")
df$upper_roll_mean <- rollmean(df$upperbound, k = window_size, fill = NA, align = "center")
df$lower_roll_min <- rollmin(df$lowerbound, k = window_size, fill = NA, align = "center")
df$upper_roll_max <- rollmax(df$upperbound, k = window_size, fill = NA, align = "center")

# Plot rolling statistics
ggplot(df, aes(x = X)) +
  geom_line(aes(y = lower_roll_mean), color = "red", size = 1) +
  geom_line(aes(y = upper_roll_mean), color = "green", size = 1) +
  geom_ribbon(aes(ymin = lower_roll_min, ymax = upper_roll_max), fill = "gray", alpha = 0.3)
 
################################################################################

