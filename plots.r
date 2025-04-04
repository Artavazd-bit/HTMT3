library(ggplot2)
library(dplyr)
library(tidyr)

normalsim <- read.csv2("normalsim.csv")
normalsimag <- read.csv2("normalsimag.csv")

nonnormalsim <- read.csv2("nonnormalsim.csv")
nonormalsimag <- read.csv2("nonnormalsimag.csv")

######################### NORMAL ################################################
normalsimaglong <- normalsimag %>%
  pivot_longer(
    cols = starts_with("rr"),    # alle Spalten, die mit "RR_" beginnen
    names_to = "Method",       # neue Spalte für die Nummern
    values_to = "rr",             # neue Spalte für die Werte
    names_prefix = "rr"          # entfernt "RR_" aus den Werten in RR_nummer
  )


##################### fuer 10 prozent NORMAL ####################################
plot001powernormal <- ggplot(normalsimaglong[normalsimaglong$correlation != 1 & normalsimaglong$Method %in% c("delta001", "boot001"),], aes(n, rr)) + 
  geom_line(aes(linetype = as.factor(Method))) + 
  geom_point() + facet_grid(rows = vars(correlation), labeller = label_parsed)  + 
  geom_hline(yintercept = c( 0.8), color = "red") + 
  theme(legend.position="left", legend.title=element_blank())


plot001errornormal <- ggplot(normalsimaglong[normalsimaglong$correlation == 1 & normalsimaglong$Method %in% c("delta001", "boot001"),], aes(n, rr)) + 
  geom_line(aes(linetype = as.factor(Method))) + 
  geom_point() + facet_grid(rows = vars(correlation), labeller = label_parsed)  + 
  geom_hline(yintercept = c( 0.01), color = "red") + 
  theme(legend.position="left", legend.title=element_blank())

##################### fuer 5 prozent NORMAL ####################################
plot005powernormal <- ggplot(normalsimaglong[normalsimaglong$correlation != 1 & normalsimaglong$Method %in% c("delta005", "boot005"),], aes(n, rr)) + 
  geom_line(aes(linetype = as.factor(Method))) + 
  geom_point() + facet_grid(rows = vars(correlation), labeller = label_parsed)  + 
  geom_hline(yintercept = c( 0.8), color = "red") + 
  theme(legend.position="left", legend.title=element_blank())


plot005errornormal <- ggplot(normalsimaglong[normalsimaglong$correlation == 1 & normalsimaglong$Method %in% c("delta005", "boot005"),], aes(n, rr)) + 
  geom_line(aes(linetype = as.factor(Method))) + 
  geom_point() + facet_grid(rows = vars(correlation), labeller = label_parsed)  + 
  geom_hline(yintercept = c( 0.05), color = "red") + 
  theme(legend.position="left", legend.title=element_blank())


##################### fuer 10 prozent NORMAL ####################################
plot010powernormal <- ggplot(normalsimaglong[normalsimaglong$correlation != 1 & normalsimaglong$Method %in% c("delta010", "boot010"),], aes(n, rr)) + 
  geom_line(aes(linetype = as.factor(Method))) + 
  geom_point() + facet_grid(rows = vars(correlation), labeller = label_parsed)  + 
  geom_hline(yintercept = c( 0.8), color = "red") + 
  theme(legend.position="left", legend.title=element_blank())


plot010errornormal <- ggplot(normalsimaglong[normalsimaglong$correlation == 1 & normalsimaglong$Method %in% c("delta010", "boot010"),], aes(n, rr)) + 
  geom_line(aes(linetype = as.factor(Method))) + 
  geom_point() + facet_grid(rows = vars(correlation), labeller = label_parsed)  + 
  geom_hline(yintercept = c( 0.10), color = "red") + 
  theme(legend.position="left", legend.title=element_blank())

####################### NON NORMAL ############################################

##################### fuer 1 prozent NONNORMAL #################################
nonnormalsimaglong <- nonormalsimag %>%
  pivot_longer(
    cols = starts_with("rr"),    # alle Spalten, die mit "RR_" beginnen
    names_to = "Method",       # neue Spalte für die Nummern
    values_to = "rr",             # neue Spalte für die Werte
    names_prefix = "rr"          # entfernt "RR_" aus den Werten in RR_nummer
  )

plot001powernonnormal <- ggplot(nonnormalsimaglong[nonnormalsimaglong$correlation != 1 & nonnormalsimaglong$Method %in% c("delta001", "boot001"),], aes(n, rr)) + 
  geom_line(aes(linetype = as.factor(Method))) + 
  geom_point() + facet_grid(rows = vars(correlation), labeller = label_parsed)  + 
  geom_hline(yintercept = c( 0.8), color = "red") + 
  theme(legend.position="left", legend.title=element_blank())


plot001errornonnormal <- ggplot(nonnormalsimaglong[nonnormalsimaglong$correlation == 1 & nonnormalsimaglong$Method %in% c("delta001", "boot001"),], aes(n, rr)) + 
  geom_line(aes(linetype = as.factor(Method))) + 
  geom_point() + facet_grid(rows = vars(correlation), labeller = label_parsed)  + 
  geom_hline(yintercept = c( 0.01), color = "red") + 
  theme(legend.position="left", legend.title=element_blank())


##################### fuer 5 prozent NONNORMAL #################################

plot005powernonnormal <- ggplot(nonnormalsimaglong[nonnormalsimaglong$correlation != 1 & nonnormalsimaglong$Method %in% c("delta005", "boot005"),], aes(n, rr)) + 
  geom_line(aes(linetype = as.factor(Method))) + 
  geom_point() + facet_grid(rows = vars(correlation), labeller = label_parsed)  + 
  geom_hline(yintercept = c( 0.8), color = "red") + 
  theme(legend.position="left", legend.title=element_blank())


plot005errornonnormal <- ggplot(nonnormalsimaglong[nonnormalsimaglong$correlation == 1 & nonnormalsimaglong$Method %in% c("delta005", "boot005"),], aes(n, rr)) + 
  geom_line(aes(linetype = as.factor(Method))) + 
  geom_point() + facet_grid(rows = vars(correlation), labeller = label_parsed)  + 
  geom_hline(yintercept = c( 0.05), color = "red") + 
  theme(legend.position="left", legend.title=element_blank())



##################### fuer 10 prozent NONNORMAL ################################# 
plot010powernonnormal <- ggplot(nonnormalsimaglong[nonnormalsimaglong$correlation != 1 & nonnormalsimaglong$Method %in% c("delta010", "boot010"),], aes(n, rr)) + 
  geom_line(aes(linetype = as.factor(Method))) + 
  geom_point() + facet_grid(rows = vars(correlation), labeller = label_parsed)  + 
  geom_hline(yintercept = c( 0.8), color = "red") + 
  theme(legend.position="left", legend.title=element_blank())


plot010errornonnormal <- ggplot(nonnormalsimaglong[nonnormalsimaglong$correlation == 1 & nonnormalsimaglong$Method %in% c("delta010", "boot010"),], aes(n, rr)) + 
  geom_line(aes(linetype = as.factor(Method))) + 
  geom_point() + facet_grid(rows = vars(correlation), labeller = label_parsed)  + 
  geom_hline(yintercept = c( 0.10), color = "red") + 
  theme(legend.position="left", legend.title=element_blank())

