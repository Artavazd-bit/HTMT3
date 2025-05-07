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

p_upper <- ggplot(resag[resag$alpha == 0.05 & resag$datatype == "normal",], aes(x = as.factor(n), y = upperwithin, group = method)) +
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

p_lower <- ggplot(resag[resag$alpha == 0.05 & resag$datatype == "normal",], aes(x = as.factor(n), y = lowerwithin, group = method)) +
  geom_line(aes(linetype = method2)) +
  geom_point(aes(shape = method2)) +
  facet_grid(cols = vars(correlation), labeller = label_parsed, ) + 
  geom_hline(yintercept = 97.5) +
  scale_y_reverse(name = "Prop. with pop. correlation lower lower limit (%)", limits = c(100, 75),  breaks = y_breaks) +
  theme_minimal() +
  theme(legend.position = "none", strip.text.x = element_blank())

p_upper / p_lower
