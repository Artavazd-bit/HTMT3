library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)

normalsim <- read.csv2("Normalsim0704.csv")
normalsimag <- read.csv2("Normalsimag0704.csv")

nonnormalsim <- read.csv2("nonnormal0704.csv")
nonnormalsimag <- read.csv2("nonnormalsimag0704.csv")

normalsim$data <- "normal"
normalsimag$data <- "normal"

nonnormalsim$data <- "nonnormal"
nonnormalsimag$data <- "nonnormal"


simdata <- rbind(normalsim, nonnormalsim)


simdataag <- rbind(normalsimag, nonnormalsimag)

simdataag_long <- simdataag %>%
  pivot_longer(
    cols = starts_with("rr"),    # alle Spalten, die mit "RR_" beginnen
    names_to = "Method",       # neue Spalte für die Nummern
    values_to = "rr",             # neue Spalte für die Werte
    names_prefix = "rr"          # entfernt "RR_" aus den Werten in RR_nummer
  )

simdataag_long$correlation <- paste("phi ==", simdataag_long$correlation)

simdataag_long$alpha <- substr(simdataag_long$Method, nchar(simdataag_long$Method)-2,  nchar(simdataag_long$Method))

simdataag_long$alpha <- as.numeric(paste0(substr(simdataag_long$alpha, 1, 1), "." , substr(simdataag_long$alpha, 2, 3)))

simdataag_long$method <- substr(simdataag_long$Method, 1,  nchar(simdataag_long$Method)-3) 


p <- ggplot(simdataag_long[simdataag_long$correlation ==  "phi == 1",], aes(n, rr)) + 
  geom_line(aes(color = as.factor(method))) + 
  geom_point() + facet_grid(cols = vars(data), rows = vars(correlation), labeller = label_parsed)  + 
  geom_hline(yintercept = vars(alpha), color = "red") + 
  theme(legend.position="left", legend.title=element_blank()) + ylab("Rejectionrate")

p + facet_grid(cols= vars(), rows = vars(alpha))



range_act <- range(range(simdataag_long$rr[simdataag_long$correlation == "phi == 1",]))

type1er <- ggplot(simdataag_long[simdataag_long$correlation == "phi == 1", ], aes(n, rr)) + 
        geom_line(aes(linetype = as.factor(method)), show.legend = FALSE) + 
        geom_point() + 
        facet_grid(cols = vars(data), rows = vars(alpha), scales = "free_y") + 
        geom_hline(data = simdataag_long, aes(yintercept = alpha)) + 
        theme(legend.position="left", legend.title=element_blank()) + ylab("Rejectionrate")


powernormal <- ggplot(simdataag_long[simdataag_long$correlation != "phi == 1" & simdataag_long$n < 400 & simdataag_long$data == "normal",], aes(n, rr)) + 
  geom_line(aes(linetype = method), show.legend = FALSE) + 
  geom_point() + facet_grid(cols = vars(alpha), rows = vars(correlation), labeller = label_parsed, scales = "free_y")  + 
  geom_hline(yintercept = c(0.8)) + 
  theme(legend.position="left", legend.title=element_blank()) + ylab("Rejectionrate")



powernonnormal <- ggplot(simdataag_long[simdataag_long$correlation != "phi == 1" & simdataag_long$n < 400 & simdataag_long$data == "nonnormal",], aes(n, rr)) + 
  geom_line(aes(linetype = method), show.legend = FALSE) + 
  geom_point() + facet_grid(cols = vars(alpha), rows = vars(correlation), labeller = label_parsed, scales = "free_y")  + 
  geom_hline(yintercept = c(0.8)) + 
  theme(legend.position="left", legend.title=element_blank()) + ylab("Rejectionrate")


################################################################################
simdataag_long$intercept <- simdataag_long$alpha
simdataag_long$intercept[simdataag_long$correlation != "phi == 1"] <- 0.8

simdataag_long$cor <- factor(simdataag_long$correlation, levels = c("phi == 1", "phi == 0.9", "phi == 0.8", "phi == 0.7")) 

powerandtypeonenonnormal <- ggplot(simdataag_long[simdataag_long$data == "nonnormal" & simdataag_long$n < 1600 ,], aes(n, rr)) + 
  geom_line(aes(linetype = method), show.legend = FALSE) + 
  geom_point() + facet_grid(cols = vars(alpha), rows = vars(cor), labeller = label_parsed, scales = "free_y")  + 
  geom_hline(data = simdataag_long, aes(yintercept = intercept)) + 
  theme(legend.position="left", legend.title=element_blank()) + ylab("Rejection Rate") + xlab("Sample Size (n)")



simdatajj <- simdata %>% 
  group_by(n) %>%
  summarize(
            tdelta = mean(as.numeric(tdelta)), 
            tboot = mean(as.numeric(tboot))
  )