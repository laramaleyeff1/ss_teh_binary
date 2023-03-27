#
# Want to compare sigma_4^2 assuming no het to sigma_4_sq^2 assuming het
# 
library(dplyr)
source("~/Dropbox/dissertation/paper2/rcode/helper_boots.R")
B = 10

beta_1 = log(0.14/(1-0.14))
beta_2 = log(0.14/(1-0.14)) - log(0.1/(1-0.1))
beta_3 = 0
beta_4 = 0.05*beta_2
K = 40
icc_yx = c(0.01,0.1,0.2)
p_x = 0.5
bar_W = 0.5
results = expand.grid(cv_n = c(0.25,0.5,0.75),
                      n = c(50,100),
                      icc_x = c(0,0.1,0.2,0.3,0.5,0.75,1),
                      sigma_gamma_sq = (pi^2/3*icc_yx)/(1-icc_yx))
results$binary_fixed = apply(results,
                               1,
                               function(x) {
                                 return(get_var_general(beta_1,
                                                        beta_2,
                                                        beta_3,
                                                        beta_4,
                                                        n = x[["n"]],
                                                        p_x = p_x,
                                                        icc_x = x[["icc_x"]],
                                                        bar_W,
                                                        sigma_gamma_sq = x[["sigma_gamma_sq"]],
                                                        sigma_sq_logit))
                               })


results$binary_vary = apply(results,
                            1,
                            function(x) {
                              return(get_monte_carlo_var(beta_1,
                                                         beta_2,
                                                         beta_3,
                                                         beta_4,
                                                         n = x[["n"]],
                                                         cv_n = x[["cv_n"]],
                                                         icc_x = x[["icc_x"]],
                                                         bar_W,
                                                         sigma_gamma_sq = x[["sigma_gamma_sq"]],
                                                         moments_x = p_x,
                                                         generate_binx,
                                                         sigma_sq_logit,
                                                         K,
                                                         B))
                            })



results$B = B
results$ratio = results$binary_vary/results$binary_fixed
#results_1000_vary = results
#save(results_1000_vary, file="results_1000_vary.Rda")
