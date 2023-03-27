source("../../helper_boots.R")
load("../get_K_values_boot/get_K_values.Rda")
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

params <- commandArgs(trailingOnly=TRUE)

# Number of time points observed
n <- as.numeric(params[1])

# Covariate ICC
icc_x <- as.numeric(params[2])

# Outcome ICC
icc_yx <- as.numeric(params[3])

# beta
beta_4_nonzero <- as.numeric(params[4])


# Index
index <- as.numeric(params[6])


#########################################
#   Set seed for one-time generation    #
#   of random effects (normal and non). #
#   Then, set seed to unique simulation #
#   iteration (to ensure variability)   #
#########################################
set.seed(123, kind = "L'Ecuyer-CMRG")
seeds <- runif(125, 1,10000)


# Reset seed (different seed for each job array) so that
# simulations are reproducible
set.seed(seeds[index],kind = "L'Ecuyer-CMRG")


B = 10
beta_1 <- 0.1
beta_2 <- -0.1
beta_3 <- 0.2
beta_4 <- c(0, beta_4_nonzero)
bar_W <- 0.5
mu_x <- 0
sigma_x_sq <- 1
sigma_gamma_sq <- (pi^2/3*icc_yx)/(1-icc_yx)
K <- get_K_values[get_K_values$n == n &
                    get_K_values$icc_x == icc_x &
                    get_K_values$icc_yx == icc_yx &
                    get_K_values$beta_4 == beta_4_nonzero]$K



generateData <- function(beta_4) {
  cluster_eff  = rnorm(K, 0, sqrt(icc_x*sigma_x_sq))
  group.ranef <- rnorm(K,0,sqrt(sigma_gamma_sq))
  group = rep(1:K,each=n)
  W = c(rep(0,floor((1-bar_W)*K*n)),
        rep(1,K*n - floor((1-bar_W)*K*n)))
  X =  mu_x + 
    rnorm(n*K, 0, sqrt(sigma_x_sq*(1-icc_x))) +
    cluster_eff[group]
  
  W_X = W*X
  data <- data.frame(
    group,
    W,
    X,
    W_X)
  data$g <- beta_1 +
    beta_2*data$W +
    beta_3*data$X + 
    beta_4*data$W_X +
    apply(data,
          1,
          function(x) {
            # Random effect for group 
            group.ranef[as.numeric(x["group"])]
          })
  
  data$prob <- expit(data$g)
  data$Y <- rbinom(nrow(data),1,data$prob)
  return(data)
  
}


one_iter <- function(beta_4) {
  data <- generateData(beta_4)
  model <- glmer(Y ~ W*X + (1|group), 
                 data=data,
                 family=binomial())
  sigma_4_sq = get_monte_carlo_var(beta_1,
                                   beta_2,
                                   beta_3,
                                   beta_4,
                                   n,
                                   cv_n = 0,
                                   icc_x,
                                   bar_W,
                                   sigma_gamma_sq,
                                   moments_x = c(mu_x, sigma_x_sq),
                                   generate_contx,
                                   sigma_sq_logit,
                                   K,
                                   B)
  
  
  
  
  returned = cbind( 
    data.frame( n = n,
                K = K,
                mu_x = mu_x,
                sigma_x_sq = sigma_x_sq,
                icc_x = icc_x,
                beta_1 = beta_1,
                beta_2 = beta_2,
                beta_3 = beta_3,
                beta_4 = beta_4,
                bar_W = bar_W,
                sigma_gamma_sq = sigma_gamma_sq,
                icc_yx = icc_yx,
                p_value_model = summary(model)$coefficients[4,4],
                sd_model = summary(model)$coefficients[4,2],
                sd_boot = sqrt(sigma_4_sq/K),
                power_boot = get_power_general(K, 
                                               sigma_4_sq, 
                                               beta_4),
                B = B
                ))
  return(returned) 

}


one_run <- function() {
  ret <- foreach(i = 1:2, .combine=rbind) %do% 
    one_iter(beta_4[i])
  return(ret)
}


registerDoParallel(cores=10)
out <- foreach(iter = 1:100, .combine=rbind) %dopar% one_run()

end_time = Sys.time()
out$time = difftime(end_time, start_time, units = "mins")

# Creating and writing output dataframe
setwd("results")
out.name <- paste0(paste("cont_power", 
                         n,
                         K,
                         icc_yx,
                         icc_x,
                         beta_4_nonzero,
                         index, sep="_"), ".csv")
write.csv(out, out.name, row.names=FALSE)

