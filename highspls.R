library(tidyverse)

highspls <- read.csv2("simresults/testingwith1/Normalsimaghighspls.csv")
highnonnormalspls <- read.csv2("simresults/testingwith1/nonnormalsimaghighspls.csv") 

simdatahighspls <- rbind(highspls, highnonnormalspls)

simdatahighsplslong <- simdatahighspls %>%
  pivot_longer(
    cols = starts_with("rr"),    # alle Spalten, die mit "RR_" beginnen
    names_to = "Method",       # neue Spalte für die Nummern
    values_to = "rr",             # neue Spalte für die Werte
    names_prefix = "rr"          # entfernt "RR_" aus den Werten in RR_nummer
  )

simdatahighsplslong$correlation <- paste("phi ==", simdatahighsplslong$correlation)


ggplot(simdatahighsplslong[simdatahighsplslong$correlation == "phi == 1",], aes(n, rr)) + 
  geom_line(aes(color = as.factor(Method))) + 
  geom_point() + facet_grid(cols = vars(data), rows = vars(correlation), labeller = label_parsed)  + 
  geom_hline(yintercept = c(0.01, 0.05, 0.10), color = "red") + 
  theme(legend.position="left", legend.title=element_blank()) + ylab("Rejectionrate")
