library(ggplot2)

sim_overview_2 <- read.csv2("2025_03_14_results")

sim_overview_delta <- sim_overview_2[, c("loading1" , "loading2", "correlation", "n" , 
                                        "Rejection_rate_htmt_09", "Rejection_rate_htmt_095",
                                        "Rejection_rate_htmt_099", "comp_time_delta")]

colnames(sim_overview_delta) <- c("loading1" , "loading2", "correlation", "n" , 
                                  "RR_09", "RR_095",
                                  "RR_099", "comp_time")


sim_overview_delta$Method <- "Delta"

###############################################################################

sim_overview_boot <- sim_overview_2[, c("loading1" , "loading2", "correlation", "n" , 
                                         "Rejection_boot_09", "Rejection_boot_095",
                                         "Rejection_boot_099", "comp_time_boot")]

colnames(sim_overview_boot) <- c("loading1" , "loading2", "correlation", "n" , 
                                  "RR_09", "RR_095",
                                  "RR_099", "comp_time")


sim_overview_boot$Method <- "Boot"

sim_overview_com <-  rbind(sim_overview_delta, sim_overview_boot)

sim_overview_com$loading_com <- paste0("lambda[1]==", sim_overview_com$loading1, "~lambda[2]==" , sim_overview_com$loading2)
sim_overview_com$corr_new <- paste("phi ==", sim_overview_com$correlation)

ggplot(sim_overview_com[sim_overview_com$correlation != 1,], aes(n, RR_09)) + 
  geom_line(aes(linetype = as.factor(Method))) + 
  geom_point() + facet_grid(cols = vars(loading_com), rows = vars(corr_new), labeller = label_parsed)  + 
  geom_hline(yintercept = c( 0.8), color = "red") + 
  theme(legend.position="left", legend.title=element_blank())


ggplot(sim_overview_com[sim_overview_com$correlation == 1,], aes(n, RR_09)) + 
  geom_line(aes(linetype = as.factor(Method))) + 
  geom_point() + facet_grid(cols = vars(loading_com), rows = vars(corr_new), labeller = label_parsed)  + 
  geom_hline(yintercept = c(0.1), color = "red") + 
  theme(legend.position="left", legend.title=element_blank())
######################################################################################################

ggplot(sim_overview_com[sim_overview_com$correlation != 1,], aes(n, RR_099)) + 
  geom_line(aes(linetype = as.factor(Method))) + geom_point() + 
  facet_grid(cols = vars(loading_com), rows = vars(corr_new), labeller = label_parsed)  + 
  geom_hline(yintercept = c(0.8), color = "red") + 
  theme(legend.position="left", legend.title=element_blank())

ggplot(sim_overview_com[sim_overview_com$correlation == 1,], aes(n, RR_099)) + 
  geom_line(aes(linetype = as.factor(Method))) + geom_point() + 
  facet_grid(cols = vars(loading_com), rows = vars(corr_new), labeller = label_parsed)  + 
  geom_hline(yintercept = c(0.01), color = "red") + 
  theme(legend.position="left", legend.title=element_blank()) 
####################################################################################################

ggplot(sim_overview_com[sim_overview_com$correlation != 1,], aes(n, RR_095)) + 
  geom_line(aes(linetype = as.factor(Method))) + geom_point() + 
  facet_grid(cols = vars(loading_com), rows = vars(corr_new), labeller = label_parsed)  + 
  geom_hline(yintercept = c(0.8), color = "red") + 
  theme(legend.position="left", legend.title=element_blank()) 


ggplot(sim_overview_com[sim_overview_com$correlation == 1,], aes(n, RR_095)) + 
  geom_line(aes(linetype = as.factor(Method))) + geom_point() + 
  facet_grid(cols = vars(loading_com), rows = vars(corr_new), labeller = label_parsed)  + 
  geom_hline(yintercept = c(0.05), color = "red") + 
  theme(legend.position="left", legend.title=element_blank())   



