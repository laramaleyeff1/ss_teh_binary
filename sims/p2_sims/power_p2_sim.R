source("../../helper.R")
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

params <- commandArgs(trailingOnly=TRUE)

beta_4 <- as.numeric(params[1])

icc_yx <- as.numeric(params[2])

# Index
index <- as.numeric(params[3])


p_x = 0.5
beta_1 <- 0.9
beta_2 <- 0.1
beta_3 <- -0.1
bar_W <- 0.5
sigma_gamma_sq <- (pi^2/3*icc_yx)/(1-icc_yx)

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



generateData <- function(n, p_x, icc_x, K, beta_4) {
  group.ranef <- rnorm(K, 0, sqrt(sigma_gamma_sq))
  group = rep(1:K,each=n)
  W = c(rep(0,floor((1-bar_W)*K*n)),
        rep(1,K*n - floor((1-bar_W)*K*n)))
  X = rcbin(p_x, 0, K, n, 0, icc_x)$y
  if (icc_x == 1) {
    X = rep(c(rep(0,n),rep(1,n)),K/2)
  }
  
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


one_iter <- function(scenario) {
  icc_x = expgrid[scenario,]$icc_x
  n = expgrid[scenario,]$n

  K = get_K(beta_1,
            beta_2,
            beta_3,
            beta_4,
            n,
            p_x,
            icc_x,
            bar_W,
            sigma_gamma_sq,
            sigma_sq_logit)

  data <- generateData(n, p_x, icc_x, K, beta_4)
  model <- glmer(Y ~ W*X + (1|group), 
                 data=data,
                 family=binomial())
  

  data.null <- generateData(n, p_x, icc_x, K, 0)
  model.null <- glmer(Y ~ W*X + (1|group), 
                 data=data.null,
                 family=binomial())

  returned = cbind(
    data.frame( n = n,
                K = K,
                p_x = p_x,
                icc_x = icc_x,
                beta_1 = beta_1,
                beta_2 = beta_2,
                beta_3 = beta_3,
                beta_4 = beta_4,
                bar_W = bar_W,
                sigma_gamma_sq = sigma_gamma_sq,
                icc_yx = icc_yx,
                p_val = summary(model)$coefficients[4,4],
                sd_model = summary(model)$coefficients[4,2],
                p_val_null = summary(model.null)$coefficients[4,4],
                sd_model_null = summary(model.null)$coefficients[4,2],
                sigma_gamma_sq_model = as.numeric(summary(model)$varcor$group[1]),
                est_power = get_power(K, 
                                      beta_1,
                                      beta_2,
                                      beta_3,
                                      beta_4,
                                      n,
                                      p_x,
                                      icc_x,
                                      bar_W,
                                      sigma_gamma_sq,
                                      sigma_sq_logit)))
  return(returned) 
}


expgrid = expand.grid(n = c(30,50,100),
                      icc_x = c(0.1,0.25,0.5,1)
                      )

one_run <- function() {
  ret <- foreach(scenario = 1:nrow(expgrid), .combine=rbind) %do% 
    one_iter(scenario)
  return(ret)
}


registerDoParallel(cores=10)
out <- foreach(iter = 1:10, .combine=rbind) %dopar% one_run()

end_time = Sys.time()
out$time = difftime(end_time, start_time, units = "mins")

# Creating and writing output dataframe
setwd("results")
out.name <- paste0(paste("bin_power", 
                         beta_4,
                         icc_yx,
                         index, sep="_"), ".csv")
write.csv(out, out.name, row.names=FALSE)

