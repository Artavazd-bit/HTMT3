library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)

normalsim <- read.csv2("simresults/testingwith1/normalsim.csv")
normalsimag <- read.csv2("simresults/testingwith1/normalsimag.csv")

nonnormalsim <- read.csv2("simresults/testingwith1/nonnormalsim.csv")
nonormalsimag <- read.csv2("simresults/testingwith1/nonnormalsimag.csv")


normalsim$data <- "normal"
normalsimag$data <- "normal"
nonnormalsim$data <- "nonnormal"
nonormalsimag$data <- "nonnormal"

testcommit


simdata <- rbind(normalsim, nonnormalsim)
simdataag <- rbind(normalsimag, nonormalsimag)

simdataag_long <- simdataag %>%
  pivot_longer(
    cols = starts_with("rr"),    # alle Spalten, die mit "RR_" beginnen
    names_to = "Method",       # neue Spalte für die Nummern
    values_to = "rr",             # neue Spalte für die Werte
    names_prefix = "rr"          # entfernt "RR_" aus den Werten in RR_nummer
  )

simdataag_long$correlation <- paste("phi ==", simdataag_long$correlation)


ggplot(simdataag_long[simdataag_long$correlation != "phi == 1",], aes(n, rr)) + 
  geom_line(aes(linetype = as.factor(Method))) + 
  geom_point() + facet_grid(cols = vars(data), rows = vars(correlation), labeller = label_parsed)  + 
  geom_hline(yintercept = c( 0.8), color = "red") + 
  theme(legend.position="left", legend.title=element_blank()) + ylab("Rejectionrate")


ggplot(simdataag_long[simdataag_long$correlation == "phi == 1",], aes(n, rr)) + 
  geom_line(aes(color = as.factor(Method))) + 
  geom_point() + facet_grid(cols = vars(data), rows = vars(correlation), labeller = label_parsed)  + 
  geom_hline(yintercept = c(0.01, 0.05, 0.10), color = "red") + 
  theme(legend.position="left", legend.title=element_blank()) + ylab("Rejectionrate")


