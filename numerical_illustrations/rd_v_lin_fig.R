library(viridis)
library(doParallel)
library(dplyr)
library(tidyr)
library(ggplot2)
source("/Users/laramaleyeff/Dropbox/dissertation/paper2/rcode/helper_copy.R")

p_x = 0.5
bar_W <- 0.5
beta_1 <- 0.01
beta_2 <- 0.2
beta_3 <- 0
n = 50
K = 20
sigma_gamma_sq <- 0.1

expgrid = expand.grid(icc_x = c(seq(0,1,by=0.01)),
                      icc_yx = seq(0,0.3,by=0.01),
                      beta_4 = c(0.1, 0.2))
expgrid$sigma_gamma_sq = (pi^2/3*expgrid$icc_yx)/(1-expgrid$icc_yx)

expgrid$rd_var <- foreach(scenario = 1:nrow(expgrid), .combine=rbind) %do% 
  get_var_general(beta_1,
            beta_2,
            beta_3,
            expgrid[scenario,]$beta_4,
            n,
            p_x,
            expgrid[scenario,]$icc_x, 
            bar_W,
            expgrid[scenario,]$sigma_gamma_sq,
            sigma_sq_identity
    )

expgrid$linear_var <- foreach(scenario = 1:nrow(expgrid), .combine=rbind) %do% 
  get_var_general(beta_1,
             beta_2,
             beta_3,
             expgrid[scenario,]$beta_4,
             n,
             p_x,
             expgrid[scenario,]$icc_x, 
             bar_W,
             expgrid[scenario,]$sigma_gamma_sq,
             sigma_sq_linear_approx)
expgrid$ratio = expgrid$rd_var/expgrid$linear_var

expgrid %>%
  filter(beta_4 != -0.05)%>%
ggplot(aes(x=icc_yx,y=icc_x,fill=ratio)) +
  geom_tile() +
  scale_fill_viridis(discrete=FALSE) +
  labs(x=expression(paste("Outcome ICC (", rho["y|x"], ")", sep="")), 
       y=expression(paste("Covariate ICC (", rho[x] ,")", sep="")),
       fill="Ratio of variances") +
  facet_wrap(~beta_4)

ggsave(file = "/Users/laramaleyeff/Dropbox/dissertation/paper2/figures/rd_v_lin.eps",
       width=8,
       height=6,
       device=cairo_ps)
# Plot for presentation
expgrid %>%
  filter(beta_4 == 0.1)%>%
ggplot(aes(x=icc_yx,y=icc_x,fill=ratio)) +
  geom_tile() +
  scale_fill_viridis(discrete=FALSE) +
  labs(x=expression(paste("Outcome ICC (", rho["y|x"], ")", sep="")), 
       y=expression(paste("Covariate ICC (", rho[x] ,")", sep="")),
       fill="Ratio of variances") 

ggsave(file = "/Users/laramaleyeff/Dropbox/dissertation/paper2/figures/rd_v_lin_pres.jpg",
       width=4,
       height=5)

