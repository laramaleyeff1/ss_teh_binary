source("/Users/laramaleyeff/Dropbox/dissertation/paper2/rcode/helper.R")

## Function to generate binary covariate X_ki with covariate ICC icc_x
## and mean p_x
##
## @param n_vector    vector of length K containing the number of individuals in 
##                    each cluster
## @param K           number of clusters
## @param icc_x       covariate ICC
## @param moments_x   vector of length 1 containing the first moment of X_ki
generate_binx <- function(n_vector,
                          K,
                          icc_x,
                          moments_x
                          ) {
  p_x = moments_x[1]
  Z_k <- rep(rbinom(K, 1, p_x),n_vector)
  
  Y_ki <-  unlist(lapply(1:K, function(k) {rbinom(n_vector[k],1,p_x)}))
  U_ki <- unlist(lapply(1:K, function(k) {rbinom(n_vector[k],1,sqrt(icc_x))}))
  (1-U_ki)*Y_ki + U_ki*Z_k
  return((1-U_ki)*Y_ki + U_ki*Z_k)
}

## Function to generate continuous covariate X_ki with covariate ICC icc_x,
## mean mu_x, and variance sigma_x_sq
##
## @param n_vector    vector of length K containing the number of individuals in 
##                    each cluster
## @param K           number of clusters
## @param icc_x       covariate ICC
## @param moments_x   vector of length 2 containing the first two moments of X_ki
generate_contx <- function(n_vector,
                           K,
                           icc_x,
                           moments_x) {
  mu_x = moments_x[1]
  sigma_x_sq = moments_x[2]
  cluster_eff = rnorm(K, 0, sqrt(icc_x*sigma_x_sq))
  return(mu_x + rnorm(sum(n_vector),0,sqrt(sigma_x_sq*(1-icc_x))) +
           unlist(lapply(1:K, function(k) {rep(cluster_eff[k],n_vector[k])})))
}

## @param beta_1  value of intercept
## @param beta_2  main effect of treatment
## @param beta_3  main covariate effect (for univariate binary variable)
## @param beta_4  interaction effect - this is the main 
##                effect size of interest!
## @param n   average number of individuals per cluster
## @param cv_n    coefficient of variation for n (0 for constant size)
## @param icc_x   covariate ICC for binary variable X
## @param bar_W   proportion of cluster receiving treatment
## @param sigma_gamma_sq    cluster-level heterogeneity term 
## @param moments_x     vector of length 1 containing P(X_ki = 1) for 
##                      binary X_ki and vector of length 2 containing
##                      mean and variance for continuous X_ki
## @param generate_x_fn   function which generates covariates X_ki
## @param sigma_sq_fn    function which computes sigma_{w,x}^2
## @param K         number of clusters
## @param B         number of Monte Carlo samples
get_monte_carlo_var <- function(beta_1,
                                 beta_2,
                                 beta_3,
                                 beta_4,
                                 n,
                                 cv_n,
                                 icc_x,
                                 bar_W,
                                 sigma_gamma_sq,
                                 moments_x,
                                 generate_x_fn,
                                 sigma_sq_fn,
                                 K,
                                 B) {

  treatment <- data.frame(group = 1:K,
                          W = c(rep(1, round(bar_W*K)),
                                rep(0, K-round(bar_W*K))))
  oneStep <- function() {
    if (cv_n == 0) {
      n_vector <- rep(n,K)
    } else {
      n_vector <- round(rgamma(K, 1/cv_n^2, (cv_n^2*n)^(-1)))
    }
    
    n_vector[n_vector < 2] = 2
    data = data.frame(id = unlist(lapply(1:K, function(k) {1:n_vector[k]})),
                      group=rep(1:K,n_vector))
    data = merge(data, treatment, by = c("group"), all.x = TRUE) 
    data$X = generate_x_fn(n_vector,
                           K,
                           icc_x,
                           moments_x
                           )
    
    data$W_X = data$W * data$X
    data$sigma_sq = sigma_sq_fn(beta_1, 
                             beta_2, 
                             beta_3, 
                             beta_4, 
                             data$W, 
                             data$X)
    XWX_k = lapply(1:K, function(k) {
      ones <- rep(1,n_vector[k])
      data.group = data[data$group == k,]
      X = t(matrix(c(ones,
                     data.group$W, 
                     data.group$X, 
                     data.group$W_X),
                   byrow=TRUE, nrow=4))
      V <- diag(data.group$sigma_sq) + sigma_gamma_sq * ones %*% t(ones)
      
      return((t(X) %*% solve(V) %*% X))
    })
    
    solve(Reduce('+', XWX_k))[4,4]*K
  }
  
  vars = unlist(lapply(1:B, function(i) {
    oneStep()
  }))
  
  return(mean(vars))
}

## @param beta_1  value of intercept
## @param beta_2  main effect of treatment
## @param beta_3  main covariate effect (for univariate binary variable)
## @param beta_4  interaction effect - this is the main 
##                effect size of interest!
## @param n   average number of individuals per cluster
## @param cv_n    coefficient of variation for n (0 for constant size)
## @param icc_x   covariate ICC for binary variable X
## @param bar_W   proportion of cluster receiving treatment
## @param sigma_gamma_sq    cluster-level heterogeneity term 
## @param moments_x     vector of length 1 containing P(X_ki = 1) for 
##                      binary X_ki and vector of length 2 containing
##                      mean and variance for continuous X_ki
## @param generate_x_fn   function which generates covariates X_ki
## @param sigma_sq_fn    function which computes sigma_{w,x}^2
## @param K         number of clusters
## @param B         (defaults to 1000) number of Monte Carlo samples
## @param type_1_error     (defaults to 0.05) type 1 error 
get_power_boot = function(beta_1,
                     beta_2,
                     beta_3,
                     beta_4,
                     n,
                     cv_n,
                     icc_x,
                     bar_W,
                     sigma_gamma_sq,
                     moments_x,
                     generate_x_fn,
                     sigma_sq_fn,
                     K,
                     B = 1000,
                     type_1_error = 0.05
) {
  sigma_4_sq = get_monte_carlo_var(beta_1,
                            beta_2,
                            beta_3,
                            beta_4,
                            n,
                            cv_n,
                            icc_x,
                            bar_W,
                            sigma_gamma_sq,
                            moments_x,
                            generate_x_fn,
                            sigma_sq_fn,
                            K,
                            B)
  get_power_general(K, sigma_4_sq, beta_4, type_1_error)
}

get_K_boot = function(beta_1,
                      beta_2,
                      beta_3,
                      beta_4,
                      n,
                      cv_n,
                      icc_x,
                      bar_W,
                      sigma_gamma_sq,
                      moments_x,
                      generate_x_fn,
                      sigma_sq_fn,
                      B = 1000,
                      power = 0.8,
                      type_1_error = 0.05) {
  opt = optimize(function(K) {
    K = ceiling(K)
    est_power = get_power_boot(beta_1,
                               beta_2,
                               beta_3,
                               beta_4,
                               n,
                               cv_n,
                               icc_x,
                               bar_W,
                               sigma_gamma_sq,
                               moments_x,
                               generate_x_fn,
                               sigma_sq_fn,
                               K,
                               B,
                               type_1_error)
    return(ifelse(round(est_power,0) >= (power*100), abs(est_power - (power*100)), 10000))
  }, interval = c(2,100), maximum = FALSE, tol=2)
  return(opt$minimum)
}
