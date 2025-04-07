normalsim <- read.csv2("normalsim.csv")
colnames(normalsim)


test <- normalsim[, c("correlation", "n" , "htmt", "sehtmt")]

test$nt <- (test$htmt - 0) / test$sehtmt

test$dec <- test$nt < qnorm(p = 0.95, mean = 0, sd = 1)


test2 <- test %>% group_by(correlation, n) %>% summarize(t = mean(dec))
