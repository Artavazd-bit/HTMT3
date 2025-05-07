library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)

res <- read.csv2("simresultsconfidence.csv")

head(res)

missingcases <- res %>% filter(is.na(HTMT))
nrow(missingcases)
table(missingcases$correlation)
table(missingcases$n)

missingcases005 <- missingcases %>% filter(alpha == 0.05)

reswona <- res %>% filter(!is.na(HTMT))

head(reswona)

resag <- reswona %>% 
  group_by(correlation, n, datatype, alpha, method) %>%
  summarize(covagcorr = mean(coveragecorr), 
            covagone = mean(coverageone)
  )

head(resag)

resag$intercept <- 1 - resag$alpha
resag$correlation <- paste("phi ==", resag$correlation)
resag$alpha <- paste("gamma ==", resag$alpha)

resag$intercept <- 1 - resag$alpha

##################### non normal coverage corr ##################################
ggplot(resag[resag$datatype == "nonnormal",], aes(n, covagcorr)) + 
  geom_line(aes(linetype = method), show.legend = TRUE) + 
  geom_point() + facet_grid(cols = vars(alpha), rows = vars(correlation), labeller = label_parsed, scales = "free_y")  + 
  geom_hline(data = resag, aes(yintercept = intercept)) + 
  theme(legend.position="left", legend.title=element_blank()) + ylab("Coverage Rate") + xlab("Sample Size (n)") + 
  scale_x_continuous(breaks=c(12, 100, 200, 400, 800))

##################### non normal coverage one ##################################
ggplot(resag[resag$datatype == "nonnormal",], aes(n, covagone)) + 
  geom_line(aes(linetype = method), show.legend = TRUE) + 
  geom_point() + facet_grid(cols = vars(alpha), rows = vars(correlation), labeller = label_parsed, scales = "free_y")  + 
  #geom_hline(data = resag, aes(yintercept = intercept)) + 
  theme(legend.position="left", legend.title=element_blank()) + ylab("Coverage Rate") + xlab("Sample Size (n)") + 
  scale_x_continuous(breaks=c(12, 100, 200, 400, 800))

##################### normal coverage corr #####################################
ggplot(resag[resag$datatype == "normal",], aes(n, covagcorr)) + 
  geom_line(aes(linetype = method), show.legend = TRUE) + 
  geom_point() + facet_grid(cols = vars(alpha), rows = vars(correlation), labeller = label_parsed, scales = "free_y")  + 
  geom_hline(data = resag, aes(yintercept = intercept)) + 
  theme(legend.position="left", legend.title=element_blank()) + ylab("Coverage Rate") + xlab("Sample Size (n)") + 
  scale_x_continuous(breaks=c(12, 100, 200, 400, 800))

##################### normal coverage one #####################################
ggplot(resag[resag$datatype == "normal",], aes(n, covagone)) + 
  geom_line(aes(linetype = method), show.legend = TRUE) + 
  geom_point() + facet_grid(cols = vars(alpha), rows = vars(correlation), labeller = label_parsed, scales = "free_y")  + 
  #geom_hline(data = resag, aes(yintercept = intercept)) + 
  theme(legend.position="left", legend.title=element_blank()) + ylab("Coverage Rate") + xlab("Sample Size (n)") + 
scale_x_continuous(breaks=c(12, 100, 200, 400, 800))

################### Cmp time #################################################

cmpdf <- reswona %>% filter(alpha == 0.05, datatype == "nonnormal")

cmpdf$tme <- as.numeric(cmpdf$time)
cmpdfag <- cmpdf %>% group_by(method, n) %>% summarize(time  = mean(tme))


compcmpdf <- pivot_wider(data = cmpdfag, names_from = method, values_from = time)

compcmpdf$savingsinpercent <- (1-(compcmpdf$delta / compcmpdf$boot))*100
