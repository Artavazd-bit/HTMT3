library(dplyr)
library(ggplot2)
library(patchwork)
df <- read.csv2("simresults/confidenceintervals/confruns/HTMTsimresults100200400.csv")
df2 <- read.csv2("simresults/confidenceintervals/confruns/HTMTsimresults50800.csv")
df3 <- read.csv2("simresults/confidenceintervals/confruns/HTMTsimresults251600.csv")

dfall <- rbind(df, df2, df3)

dfall$upperwithin <- dfall$correlation < dfall$upperbound
dfall$lowerwithin <- dfall$correlation > dfall$lowerbound



resag <- dfall %>% 
  group_by(correlation, n, datatype, alpha, method) %>%
  summarize(upperwithin = mean(upperwithin)*100,
            lowerwithin = mean(lowerwithin)*100
  )

resag$method2[resag$method == "boot"] <- "CI_boot" 
resag$method2[resag$method == "delta"] <- "CI_normal" 
resag$method2[resag$method == "bcaboot"] <- "CI_bcaboot"

y_breaks <- seq(85, 100, by = 5)

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
