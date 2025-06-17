library(dplyr)
library(ggplot2)
library(patchwork)
df <- read.csv2("simresults/confidenceintervals/confruns/HTMTsimresults100200400.csv")
df2 <- read.csv2("simresults/confidenceintervals/confruns/HTMTsimresults50800.csv")
df3 <- read.csv2("simresults/confidenceintervals/confruns/HTMTsimresults251600.csv")

dfall <- rbind(df, df2, df3)

rm(df, df2, df3)

dfall$upperwithin <- dfall$correlation < dfall$upperbound
dfall$lowerwithin <- dfall$correlation > dfall$lowerbound

resag <- dfall %>% 
  group_by(correlation, n, datatype, alpha, method) %>%
  summarize(upperwithin = mean(upperwithin)*100,
            lowerwithin = mean(lowerwithin)*100,
            covagoneag = mean(coverageone)*100,
            covagcorrag = mean(coveragecorr)*100
  )

resag$correlation <- format(resag$correlation, nsmall = 2)

resag$method2[resag$method == "boot"] <- "percentile" 
resag$method2[resag$method == "delta"] <- "symmetric" 
resag$method2[resag$method == "bcaboot"] <- "BCa"

################################################################################
## Coverage of corr
################################################################################
alpha = c(0.05)
datatype = "nonnormal"
lowertick <- 80

y_breaks <- c(seq(lowertick, 100, by = 5), (1-alpha/2)*100)

decimal_labeller <- function(x) {
  sprintf("%.2f", as.numeric(x))
}

p_upper_n <- ggplot(resag[resag$alpha == alpha & resag$datatype == datatype,], aes(x = as.factor(n), y = upperwithin, group = method)) +
  geom_line(aes(linetype = method2)) +
  geom_point(aes(shape = method2)) +
  facet_grid(cols = vars(correlation), labeller = labeller(col_var = decimal_labeller)) + 
  geom_hline(yintercept = (1-alpha/2)*100) + 
  scale_y_continuous(name = "pop. correlation value below upper limit (%)", breaks = y_breaks, 
                     limits = c(lowertick, 100) )  +
  theme_minimal() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) 

p_lower_n <- ggplot(resag[resag$alpha == alpha & resag$datatype == datatype,], aes(x = as.factor(n), y = lowerwithin, group = method)) +
  geom_line(aes(linetype = method2)) +
  geom_point(aes(shape = method2)) +
  facet_grid(cols = vars(correlation), labeller = label_parsed, ) + 
  geom_hline(yintercept = (1-alpha/2)*100) +
  scale_y_reverse(name = "pop. correlation value above lower limit (%)", breaks = y_breaks, 
                  limits = c(100, lowertick)) +
  theme_minimal() +
  theme(legend.position = "bottom", strip.text.x = element_blank()) +
  labs(x = "sample size") + 
  scale_linetype_discrete(name = "Type of CI:") +
  scale_shape_discrete(name = "Type of CI:")

popcov005nonnormal <- p_upper_n / p_lower_n
ggsave("popcoverage_005_nonnormal.png", plot = popcov005nonnormal, width = 12.375, height = 9.15625)
###############################################################################alpha = c(0.05)
alpha = c(0.05)
datatype = "normal"
lowertick <- 85

y_breaks <- c(seq(lowertick, 100, by = 5), (1-alpha/2)*100)

decimal_labeller <- function(x) {
  sprintf("%.2f", as.numeric(x))
}

p_upper_n <- ggplot(resag[resag$alpha == alpha & resag$datatype == datatype,], aes(x = as.factor(n), y = upperwithin, group = method)) +
  geom_line(aes(linetype = method2)) +
  geom_point(aes(shape = method2)) +
  facet_grid(cols = vars(correlation), labeller = labeller(col_var = decimal_labeller)) + 
  geom_hline(yintercept = (1-alpha/2)*100) + 
  scale_y_continuous(name = "pop. correlation value below upper limit (%)", breaks = y_breaks, 
                     limits = c(lowertick, 100) )  +
  theme_minimal() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) 

p_lower_n <- ggplot(resag[resag$alpha == alpha & resag$datatype == datatype,], aes(x = as.factor(n), y = lowerwithin, group = method)) +
  geom_line(aes(linetype = method2)) +
  geom_point(aes(shape = method2)) +
  facet_grid(cols = vars(correlation), labeller = label_parsed, ) + 
  geom_hline(yintercept = (1-alpha/2)*100) +
  scale_y_reverse(name = "pop. correlation value above lower limit (%)", breaks = y_breaks, 
                  limits = c(100, lowertick)) +
  theme_minimal() +
  theme(legend.position = "bottom", strip.text.x = element_blank()) +
  labs(x = "sample size") + 
  scale_linetype_discrete(name = "Type of CI:") +
  scale_shape_discrete(name = "Type of CI:")

popcov005normal <- p_upper_n / p_lower_n
ggsave("popcoverage_005_normal.png", plot = popcov005normal, width = 12.375, height = 9.15625)
################################################################################
alpha = c(0.10)
datatype = "nonnormal"
lowertick <- 75

y_breaks <- c(seq(lowertick, 100, by = 5), (1-alpha/2)*100)

decimal_labeller <- function(x) {
  sprintf("%.2f", as.numeric(x))
}

p_upper_n <- ggplot(resag[resag$alpha == alpha & resag$datatype == datatype,], aes(x = as.factor(n), y = upperwithin, group = method)) +
  geom_line(aes(linetype = method2)) +
  geom_point(aes(shape = method2)) +
  facet_grid(cols = vars(correlation), labeller = labeller(col_var = decimal_labeller)) + 
  geom_hline(yintercept = (1-alpha/2)*100) + 
  scale_y_continuous(name = "pop. correlation value below upper limit (%)", breaks = y_breaks, 
                     limits = c(lowertick, 100) )  +
  theme_minimal() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) 

p_lower_n <- ggplot(resag[resag$alpha == alpha & resag$datatype == datatype,], aes(x = as.factor(n), y = lowerwithin, group = method)) +
  geom_line(aes(linetype = method2)) +
  geom_point(aes(shape = method2)) +
  facet_grid(cols = vars(correlation), labeller = label_parsed, ) + 
  geom_hline(yintercept = (1-alpha/2)*100) +
  scale_y_reverse(name = "pop. correlation value above lower limit (%)", breaks = y_breaks, 
                  limits = c(100, lowertick)) +
  theme_minimal() +
  theme(legend.position = "bottom", strip.text.x = element_blank()) +
  labs(x = "sample size") + 
  scale_linetype_discrete(name = "Type of CI:") +
  scale_shape_discrete(name = "Type of CI:")

popcov010nonnormal <- p_upper_n / p_lower_n
ggsave("popcoverage_010_nonnormal.png", plot = popcov010nonnormal, width = 12.375, height = 9.15625)
################################################################################
alpha = c(0.10)
datatype = "normal"
lowertick <- 80

y_breaks <- c(seq(lowertick, 100, by = 5), (1-alpha/2)*100)

decimal_labeller <- function(x) {
  sprintf("%.2f", as.numeric(x))
}

p_upper_n <- ggplot(resag[resag$alpha == alpha & resag$datatype == datatype,], aes(x = as.factor(n), y = upperwithin, group = method)) +
  geom_line(aes(linetype = method2)) +
  geom_point(aes(shape = method2)) +
  facet_grid(cols = vars(correlation), labeller = labeller(col_var = decimal_labeller)) + 
  geom_hline(yintercept = (1-alpha/2)*100) + 
  scale_y_continuous(name = "pop. correlation value below upper limit (%)", breaks = y_breaks, 
                     limits = c(lowertick, 100) )  +
  theme_minimal() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) 

p_lower_n <- ggplot(resag[resag$alpha == alpha & resag$datatype == datatype,], aes(x = as.factor(n), y = lowerwithin, group = method)) +
  geom_line(aes(linetype = method2)) +
  geom_point(aes(shape = method2)) +
  facet_grid(cols = vars(correlation), labeller = label_parsed, ) + 
  geom_hline(yintercept = (1-alpha/2)*100) +
  scale_y_reverse(name = "pop. correlation value above lower limit (%)", breaks = y_breaks, 
                  limits = c(100, lowertick)) +
  theme_minimal() +
  theme(legend.position = "bottom", strip.text.x = element_blank()) +
  labs(x = "sample size") + 
  scale_linetype_discrete(name = "Type of CI:") +
  scale_shape_discrete(name = "Type of CI:")

popcov010normal <- p_upper_n / p_lower_n
ggsave("popcoverage_010_normal.png", plot = popcov010normal, width = 12.375, height = 9.15625)
################################################################################
alpha = c(0.01)
datatype = "nonnormal"
lowertick <- 85

y_breaks <- c(seq(lowertick, 100, by = 5), (1-alpha/2)*100)

decimal_labeller <- function(x) {
  sprintf("%.2f", as.numeric(x))
}

p_upper_n <- ggplot(resag[resag$alpha == alpha & resag$datatype == datatype,], aes(x = as.factor(n), y = upperwithin, group = method)) +
  geom_line(aes(linetype = method2)) +
  geom_point(aes(shape = method2)) +
  facet_grid(cols = vars(correlation), labeller = labeller(col_var = decimal_labeller)) + 
  geom_hline(yintercept = (1-alpha/2)*100) + 
  scale_y_continuous(name = "pop. correlation value below upper limit (%)", breaks = y_breaks, 
                     limits = c(lowertick, 100) )  +
  theme_minimal() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) 

p_lower_n <- ggplot(resag[resag$alpha == alpha & resag$datatype == datatype,], aes(x = as.factor(n), y = lowerwithin, group = method)) +
  geom_line(aes(linetype = method2)) +
  geom_point(aes(shape = method2)) +
  facet_grid(cols = vars(correlation), labeller = label_parsed, ) + 
  geom_hline(yintercept = (1-alpha/2)*100) +
  scale_y_reverse(name = "pop. correlation value above lower limit (%)", breaks = y_breaks, 
                  limits = c(100, lowertick)) +
  theme_minimal() +
  theme(legend.position = "bottom", strip.text.x = element_blank()) +
  labs(x = "sample size") + 
  scale_linetype_discrete(name = "Type of CI:") +
  scale_shape_discrete(name = "Type of CI:")

popcov001nonnormal <- p_upper_n / p_lower_n
ggsave("popcoverage_001_nonnormal.png", plot = popcov001nonnormal, width = 12.375, height = 9.15625)
################################################################################
alpha = c(0.01)
datatype = "normal"
lowertick <- 90

y_breaks <- c(seq(lowertick, 100, by = 5), (1-alpha/2)*100)

decimal_labeller <- function(x) {
  sprintf("%.2f", as.numeric(x))
}

p_upper_n <- ggplot(resag[resag$alpha == alpha & resag$datatype == datatype,], aes(x = as.factor(n), y = upperwithin, group = method)) +
  geom_line(aes(linetype = method2)) +
  geom_point(aes(shape = method2)) +
  facet_grid(cols = vars(correlation), labeller = labeller(col_var = decimal_labeller)) + 
  geom_hline(yintercept = (1-alpha/2)*100) + 
  scale_y_continuous(name = "pop. correlation value below upper limit (%)", breaks = y_breaks, 
                     limits = c(lowertick, 100) )  +
  theme_minimal() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) 

p_lower_n <- ggplot(resag[resag$alpha == alpha & resag$datatype == datatype,], aes(x = as.factor(n), y = lowerwithin, group = method)) +
  geom_line(aes(linetype = method2)) +
  geom_point(aes(shape = method2)) +
  facet_grid(cols = vars(correlation), labeller = label_parsed, ) + 
  geom_hline(yintercept = (1-alpha/2)*100) +
  scale_y_reverse(name = "pop. correlation value above lower limit (%)", breaks = y_breaks, 
                  limits = c(100, lowertick)) +
  theme_minimal() +
  theme(legend.position = "bottom", strip.text.x = element_blank()) +
  labs(x = "sample size") + 
  scale_linetype_discrete(name = "Type of CI:") +
  scale_shape_discrete(name = "Type of CI:")

popcov001normal <- p_upper_n / p_lower_n
ggsave("popcoverage_001_normal.png", plot = popcov001normal, width = 12.375, height = 9.15625)
################################################################################
## Coverage of one 
################################################################################
datatype = "nonnormal"
resag$covagone2 <- 100 - resag$covagoneag
resag$hline <- 80
resag$hline[resag$correlation == "1.00"] <- (resag$alpha[resag$correlation == "1.00"]) * 100

y_breaks2 <- seq(0, 100, by = 10)

ggplot(resag[resag$datatype == datatype,], aes(x = as.factor(n), y = covagone2, group = method)) + 
  geom_line(aes(linetype = method2)) + 
  facet_grid(cols = vars(correlation), labeller = labeller(col_var = decimal_labeller), rows  = vars(alpha)) + 
  geom_point(aes(shape = method2)) + 
  theme(legend.position = "bottom") +
  labs(x = "sample size") + 
  scale_linetype_discrete(name = "Type of CI:") +
  scale_shape_discrete(name = "Type of CI:") + 
  geom_hline(data = resag, aes(yintercept = hline)) + 
  scale_y_continuous(breaks = y_breaks2, name = "rejection rate")

################################################################################
datatype = "normal"
resag$covagone2 <- 100 - resag$covagoneag
resag$hline <- 80
resag$hline[resag$correlation == "1.00"] <- (resag$alpha[resag$correlation == "1.00"]) * 100

y_breaks2 <- seq(0, 100, by = 10)

ggplot(resag[resag$datatype == datatype,], aes(x = as.factor(n), y = covagone2, group = method)) + 
  geom_line(aes(linetype = method2)) + 
  facet_grid(cols = vars(correlation), labeller = labeller(col_var = decimal_labeller), rows  = vars(alpha)) + 
  geom_point(aes(shape = method2)) + 
  theme(legend.position = "bottom") +
  labs(x = "sample size") + 
  scale_linetype_discrete(name = "Type of CI:") +
  scale_shape_discrete(name = "Type of CI:") + 
  geom_hline(data = resag, aes(yintercept = hline)) + 
  scale_y_continuous(breaks = y_breaks2, name = "rejection rate")
################################################################################

cmpdf <- dfall %>% filter(alpha == 0.05, datatype == "nonnormal")

cmpdf$tme <- as.numeric(cmpdf$time)
cmpdfag <- cmpdf %>% group_by(method, n) %>% summarize(time  = mean(tme))


compcmpdf <- tidyr::pivot_wider(data = cmpdfag, names_from = method, values_from = time)







