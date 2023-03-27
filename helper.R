#############################################################
#            Sample size requirements for                   #
# an interaction test of a single binary effect modifier    #
#     on either the cluster or individual level             #
#############################################################


expit <- function(x) {
  return(exp(x)/(1+exp(x)))
}

logit <- function(x) {
  return(log(x/(1-x)))
}

## General function that takes as input
## @param K   number of clusters
## @param sigma_4_sq    variance of beta_4
## @param beta_4   effect size of interest
##
## And returns the power to detect such an interaction
## effect
get_power_general = function(K, sigma_4_sq, beta_4, type_1_error = 0.05) {
  return((pnorm(sqrt(K/sigma_4_sq)*abs(beta_4) - qnorm(1-type_1_error/2)) + 
        pnorm(-1*sqrt(K/sigma_4_sq)*abs(beta_4) - qnorm(1-type_1_error/2)))*100 )
}

## Helper function for sigma_{w,x}^2 term in paper 
## i.e. the approximated individual-level variance term
## This function assumes param of interest is OR with a 
## single binary effect modifier
## @param beta_1  value of intercept
## @param beta_2  main effect of treatment
## @param beta_3  main covariate effect (for univariate binary variable)
## @param beta_4  interaction effect - this is the main 
##                effect size of interest!
## @param W       value of treatment (0/1)
## @param X       value of X (0/1)
sigma_sq_logit <- function(beta_1, beta_2, beta_3, beta_4, W, X) {
  mu <- expit(beta_1+beta_2*W+beta_3*X+beta_4*X*W)
  return(1/(mu*(1-mu)))
}

# Helper function for sigma_{w,x}^2 term in paper 
# i.e. the approximated individual-level variance term
# This function assumes a param of interest is RD and a
# single binary effect modifier
## @param beta_1  value of intercept
## @param beta_2  main effect of treatment
## @param beta_3  main covariate effect (for univariate binary variable)
## @param beta_4  interaction effect - this is the main 
##                effect size of interest!
## @param W       value of treatment (0/1)
## @param X       value of X (0/1)
sigma_sq_identity <- function(beta_1, beta_2, beta_3, beta_4, W, X) {
  mu <- beta_1+beta_2*W+beta_3*X+beta_4*X*W
  return(mu*(1-mu))
}

## Approximated variance to be used in the formula of 
## Yang et al. - returns an approximate CONSTANT variance
##
## @param beta_1  value of intercept
## @param beta_2  main effect of treatment
## @param beta_3  main covariate effect (for univariate binary variable)
## @param beta_4  interaction effect - this is the main 
##                effect size of interest!
## @param W       value of treatment (0/1)
## @param X       value of X (0/1)
sigma_sq_linear_approx <- function(beta_1, beta_2, beta_3, beta_4, W, X) {
  X <- c(0,0,1,1)
  W <- c(0,1,0,1)
  mu <- beta_1+beta_2*W+beta_3*X+beta_4*X*W
  return(mean(mu*(1-mu)))
}


##
## Helper function when covariate is on the cluster level
## corresponds to Section 2.4 in text
##
## @param beta_1  value of intercept
## @param beta_2  main effect of treatment
## @param beta_3  main covariate effect (for univariate binary variable)
## @param beta_4  interaction effect - this is the main 
##                effect size of interest!
## @param n   number of individuals per cluster, assumed to be constant
## @param p_x   prevalence of binary variable X
## @param bar_W   proportion of cluster receiving treatment
## @param sigma_gamma_sq    cluster-level heterogeneity term 
## @param sigma_sq    function which computes sigma_{w,x}^2
get_var_icc_1 <- function(beta_1,
                            beta_2,
                            beta_3,
                            beta_4,
                            n,
                            p_x,
                            bar_W,
                            sigma_gamma_sq,
                            sigma_sq) {
  kappa <- function(W,X,n) {
    sigma_sq_inv = 1/sigma_sq(W,X)
    d = -1*sigma_gamma_sq/(1+n*sigma_gamma_sq*sigma_sq_inv)
    return(sigma_sq_inv+sigma_sq_inv^2*n*d)
  }
  num1 = (1-p_x)*kappa(0,0,n) + (p_x)*kappa(0,1,n)
  denom1 = (1-bar_W)*p_x*(1-p_x)*kappa(0,1,n)*kappa(0,0,n)
  num2 = (1-p_x)*kappa(1,0,n) + (p_x)*kappa(1,1,n)
  denom2 = (bar_W)*p_x*(1-p_x)*kappa(1,1,n)*kappa(1,0,n)
  
  returned = (1/n * (num1/denom1 + num2/denom2))
  return(returned)
}

## Returns variance expression for beta_4 based on the formula of 
## Section 2.3 and 2.4 for individual- and cluster-level covaraites,
## respectively.
##
## @param beta_1  value of intercept
## @param beta_2  main effect of treatment
## @param beta_3  main covariate effect (for univariate binary variable)
## @param beta_4  interaction effect - this is the main 
##                effect size of interest!
## @param n   number of individuals per cluster, assumed to be constant
## @param p_x   prevalence of binary variable X
## @param icc_x   covariate ICC for binary variable X
## @param bar_W   proportion of cluster receiving treatment
## @param sigma_gamma_sq    cluster-level heterogeneity term 
## @param sigma_sq_original_fn    function which computes sigma_{w,x}^2
## @param use_exact   (defaults to TRUE) if TRUE, the exact formula is used 
##                    when the covariate is on the cluster-level
get_var_general <- function(beta_1,
                      beta_2,
                      beta_3,
                      beta_4,
                      n,
                      p_x,
                      icc_x,
                      bar_W,
                      sigma_gamma_sq,
                      sigma_sq_original_fn,
                      use_exact = TRUE
                      ) {
  ## Here we adapt the functional form of "sigma_sq_original_fn"
  ## for ease of use
  sigma_sq = function(W, X) {
    return(sigma_sq_original_fn(beta_1, beta_2, beta_3, beta_4, W,X))
  }

  if (icc_x == 1 & use_exact) {
    return(get_var_icc_1(beta_1,
                           beta_2,
                           beta_3,
                           beta_4,
                           n,
                           p_x,
                           bar_W,
                           sigma_gamma_sq,
                           sigma_sq))
  }
  eta_2 = (1/n)*((1+(n-1)*icc_x)*p_x + (n-1)*((1-icc_x)*p_x^2))
  P_1 <- function(s) {
    num1 = sigma_gamma_sq*((1-p_x)/sigma_sq(s,0 )+p_x/sigma_sq(s,1))^2
    denom1 <- 1+sigma_gamma_sq*n*((1-p_x)/sigma_sq(s,0)+p_x/sigma_sq(s,1))
    num2 = sigma_gamma_sq*(1/sigma_sq(s,0)-1/sigma_sq(s,1))^2*(eta_2-p_x^2)
    denom2 = (sigma_gamma_sq*(p_x-1)*n/sigma_sq(s,0)-
                sigma_gamma_sq*p_x*n/sigma_sq(s,1) - 1)^3
    return((-1*num1/denom1+num2/denom2))
  }
  
  P_2 <- function(s) {
    num1 = sigma_gamma_sq*p_x*((1-p_x)/sigma_sq(s,0)+p_x/sigma_sq(s,1))/sigma_sq(s,1)
    denom1 <- 1+sigma_gamma_sq*n*((1-p_x)/sigma_sq(s,0)+p_x/sigma_sq(s,1))
    num2 = sigma_gamma_sq*(1/sigma_sq(s,0)-1/sigma_sq(s,1))*
      (1/sigma_sq(s,1))*((sigma_gamma_sq * n / sigma_sq(s,0)) + 1)*(eta_2-p_x^2)
    denom2 = (sigma_gamma_sq*(p_x-1)*n/sigma_sq(s,0)-
                sigma_gamma_sq*p_x*n/sigma_sq(s,1) - 1)^3
    return((-1*num1/denom1+num2/denom2))
  }
  
  P_3 <- function(s) {
    num1 = sigma_gamma_sq*eta_2*(1/sigma_sq(s,1))^2
    denom1 <- 1+sigma_gamma_sq*n*((1-p_x)/sigma_sq(s,0)+p_x/sigma_sq(s,1))
    
    return((-1*num1/denom1))
  }
  
  a_1 = bar_W*p_x/sigma_sq(1,1) + bar_W*(1-p_x)/sigma_sq(1,0)  + 
    (1-bar_W)*p_x/sigma_sq(0,1)  + 
    (1-bar_W)*(1-p_x)/sigma_sq(0,0)  + (n*(bar_W*P_1(1) + (1-bar_W)*P_1(0)))
  
  a_2 = bar_W*(p_x/sigma_sq(1,1) + (1-p_x)/sigma_sq(1,0)) +
    n*bar_W*P_1(1) 
  b_1 = p_x*((1-bar_W)/sigma_sq(0,1)  + bar_W/sigma_sq(1,1)) +
    (n*(bar_W*P_2(1) + (1-bar_W)*P_2(0)))
  b_2 = bar_W*p_x/sigma_sq(1,1) + n*bar_W*P_2(1)
  
  d_1 = p_x*((1-bar_W)/sigma_sq(0,1)  + bar_W/sigma_sq(1,1)) + 
    (n*(bar_W*P_3(1) + (1-bar_W)*P_3(0)))
  d_2 = bar_W*p_x/sigma_sq(1,1) + n*bar_W*P_3(1)
  
  U = n*matrix(c(a_1,a_2,b_1,b_2,
                 a_2,a_2,b_2,b_2,
                 b_1,b_2,d_1,d_2,
                 b_2,b_2,d_2,d_2),byrow=TRUE,nrow=4)
  
  num1 = a_1-a_2
  denom1 = (a_1-a_2)*(d_1-d_2)-(b_1-b_2)^2
  
  num2 = a_2
  denom2 = a_2*d_2-b_2^2
  returned = (num1/denom1 + num2/denom2) * (1/(n))
  return(returned)
}


## @param beta_1  value of intercept
## @param beta_2  main effect of treatment
## @param beta_3  main covariate effect (for univariate binary variable)
## @param beta_4  interaction effect - this is the main 
##                effect size of interest!
## @param n   number of individuals per cluster, assumed to be constant
## @param p_x   prevalence of binary variable X
## @param icc_x   covariate ICC for binary variable X
## @param bar_W   proportion of cluster receiving treatment
## @param sigma_gamma_sq    cluster-level heterogeneity term 
## @param sigma_sq    function which computes sigma_{w,x}^2
## @param type_1_error     (defaults to 0.05) type 1 error 
get_power = function(K, 
                     beta_1,
                     beta_2,
                     beta_3,
                     beta_4,
                     n,
                     p_x,
                     icc_x,
                     bar_W,
                     sigma_gamma_sq,
                     sigma_sq,
                     type_1_error = 0.05
                     ) {
  sigma_4_sq = get_var_general(beta_1,
                         beta_2,
                         beta_3,
                         beta_4,
                         n,
                         p_x,
                         icc_x,
                         bar_W,
                         sigma_gamma_sq,
                         sigma_sq,
                         TRUE
                         )
  return = get_power_general(K, sigma_4_sq, beta_4, type_1_error)
  return(return) 
}

## Returns required number of clusters for trial based on
## parameters below - rounds up to next integer
##
## @param beta_1  value of intercept
## @param beta_2  main effect of treatment
## @param beta_3  main covariate effect (for univariate binary variable)
## @param beta_4  interaction effect - this is the main 
##                effect size of interest!
## @param n   number of individuals per cluster, assumed to be constant
## @param p_x   prevalence of binary variable X
## @param icc_x   covariate ICC for binary variable X
## @param bar_W   proportion of cluster receiving treatment
## @param sigma_gamma_sq    cluster-level heterogeneity term 
## @param sigma_sq    function which computes sigma_{w,x}^2
## @param power    (defaults to 0.8) desired power
## @param type_1_error     (defaults to 0.05) type 1 error 
get_K = function(beta_1,
                       beta_2,
                       beta_3,
                       beta_4,
                       n,
                       p_x,
                       icc_x,
                       bar_W,
                       sigma_gamma_sq,
                       sigma_sq,
                       power = 0.8,
                       type_1_error = 0.05) {
  sigma_4_sq = get_var_general(beta_1,
                         beta_2,
                         beta_3,
                         beta_4,
                         n,
                         p_x,
                         icc_x,
                         bar_W,
                         sigma_gamma_sq,
                         sigma_sq,
                         TRUE)
  K_raw = sigma_4_sq*(qnorm(1-type_1_error/2)+qnorm(power))^2/beta_4^2
  return(ceiling(K_raw))
}

## @param K   number of clusters
## @param beta_1  value of intercept
## @param beta_2  main effect of treatment
## @param beta_3  main covariate effect (for univariate binary variable)
## @param n   number of individuals per cluster, assumed to be constant
## @param p_x   prevalence of binary variable X
## @param icc_x   covariate ICC for binary variable X
## @param bar_W   proportion of cluster receiving treatment
## @param sigma_gamma_sq    cluster-level heterogeneity term 
## @param sigma_sq    function which computes sigma_{w,x}^2
## @param power    (defaults to 0.8) desired power
## @param type_1_error     (defaults to 0.05) type 1 error 
get_beta4 = function(K, 
                     beta_1,
                     beta_2,
                     beta_3,
                     beta_4,
                     n,
                     p_x,
                     icc_x,
                     bar_W,
                     sigma_gamma_sq,
                     sigma_sq,
                     power = 0.8,
                     type_1_error = 0.05) {
  sigma_4_sq = get_var_general(beta_1,
                         beta_2,
                         beta_3,
                         beta_4,
                         n,
                         p_x,
                         icc_x,
                         bar_W,
                         sigma_gamma_sq,
                         sigma_sq, 
                         TRUE)
  return(sqrt(K/(sigma_4_sq*(qnorm(1-type_1_error/2)+qnorm(power)))^2))
}

n = 1584
K = 26
p_x = 0.4418

beta_1 = log(0.139/(1-0.139))
beta_2 = log(0.139/(1-0.139)) - log(0.104/(1-0.104))
beta_3 = 0
beta_4 = log(1.4)
sigma_gamma_sq = 0.32
#icc_x = 0.007
bar_W = 0.5
get_var_general(beta_1,
                beta_2,
                beta_3,
                beta_4,
                n,
                p_x,
                icc_x = 0.84,
                bar_W,
                sigma_gamma_sq,
                sigma_sq_logit)
