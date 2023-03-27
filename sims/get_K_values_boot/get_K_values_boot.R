source("../../helper_boots.R")
#-----------------------------------------------------#
#              Model params                           #
#-----------------------------------------------------#

#-----------------------------------------------------
# Loading/defining all important functions
#-----------------------------------------------------
start_time = Sys.time()

library(foreach)
library(doParallel)
library(lme4)
library(ICCbin)
library(dplyr)

set.seed(123, kind = "L'Ecuyer-CMRG")



B = 10
beta_1 <- 0.1
beta_2 <- -0.1
beta_3 <- 0.2
bar_W <- 0.5
mu_x <- 0
sigma_x_sq <- 1

get_K_values = expand.grid(n = c(50,100,200),
                           icc_x = c(0.1,1),
                           icc_yx = c(0.01,0.1),
                           beta_4 = c(0.25,0.5,0.75))
get_K_values$sigma_gamma_sq <- (pi^2/3*get_K_values$icc_yx)/(1-get_K_values$icc_yx)

get_K_values$K = sapply(1:nrow(get_K_values), function(x) {
  return(get_K_boot(beta_1,
                    beta_2,
                    beta_3,
                    get_K_values$beta_4[x],
                    get_K_values$n[x],
                    cv_n = 0,
                    get_K_values$icc_x[x],
                    bar_W,
                    get_K_values$sigma_gamma_sq[x],
                    moments_x = c(mu_x, sigma_x_sq),
                    generate_contx,
                    sigma_sq_logit,
                    B))
})

get_K_values
save(get_K_values, file="get_K_values.Rda")
