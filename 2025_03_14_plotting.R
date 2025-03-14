library(ggplot2)

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

sim_overview_com$loading_com <- paste( "lmabda_1:", sim_overview_com$loading1, "and", "lambda_2:", sim_overview_com$loading2)

ggplot(sim_overview_com, aes(n, RR_09)) + geom_line(aes(linetype = as.factor(sim_overview_com$Method))) + geom_point() + facet_grid(cols = vars(loading_com), rows = vars(correlation))  + geom_hline(yintercept = c(0.1, 0.8), color = "red")+ theme(legend.position="none")  


ggplot(sim_overview_com, aes(n, RR_099)) + geom_line(aes(linetype = as.factor(sim_overview_com$Method))) + geom_point() + facet_grid(cols = vars(loading_com), rows = vars(correlation))  + geom_hline(yintercept = c(0.01, 0.8), color = "red")+ theme(legend.position="none")  


ggplot(sim_overview_com, aes(n, RR_095)) + geom_line(aes(linetype = as.factor(sim_overview_com$Method))) + geom_point() + facet_grid(cols = vars(loading_com), rows = vars(correlation))  + geom_hline(yintercept = c(0.05, 0.8), color = "red")+ theme(legend.position="none")  
