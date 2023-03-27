library(ggplot2)
library(dplyr)
library(cowplot)
library(ggthemes)
library(tidyr)
library(viridis)
source("/Users/laramaleyeff/Dropbox/dissertation/paper2/rcode/helper_copy.R")


p_x = 0.5
beta_1 <- log(0.14/(1-0.14))
beta_2 <- log(0.14/(1-0.14)) - log(0.1/(1-0.1))
beta_3 <- 0.1
beta_4 <- 0.5
bar_W <- 0.5


p_yx = c(0,0.01,0.05,0.1,0.2,0.5)
expgrid = expand.grid(n = c(30,50,100),
                      icc_x = c(0, 0.01,0.05,0.1,0.2,0.5,1),
                      sigma_gamma_sq = (pi^2/3*p_yx)/(1-p_yx)
)

expgrid$var = apply(expgrid,
                    1,
                    function(x) {
                      return(get_var_general(beta_1,
                             beta_2,
                             beta_3,
                             beta_4,
                             n = x[[1]],
                             p_x = p_x,
                             icc_x = x[[2]],
                             bar_W,
                             sigma_gamma_sq = x[[3]],
                             sigma_sq_logit))
                    })


p1 <- expgrid %>%
  filter(icc_x != 1 & icc_x != 0) %>%
  mutate(p_yx = sigma_gamma_sq/(sigma_gamma_sq+pi^2/3)) %>%
  ggplot(aes(x=p_yx,y=var,group=factor(icc_x), col=factor(icc_x))) + geom_line() +
  facet_wrap(~n, labeller = as_labeller(c(
    `30` = "n=30",
    `50` = "n=50",
    `100` = "n=100"
  )))+
   #scale_color_gradient()+
   labs(col=expression(rho[x]), 
       y=expression(sigma[4]^2),
       x=expression(paste("Outcome ICC (", rho['y|x'],")", sep="")))
p2 <- expgrid %>%
  mutate(p_yx = round(sigma_gamma_sq/(sigma_gamma_sq+pi^2/3),2)) %>%
  filter(p_yx < 0.5 & p_yx > 0) %>%
  ggplot(aes(x=icc_x,y=var,group=factor(p_yx),col=factor(p_yx))) + 
  geom_line(linetype="dashed") +
  facet_wrap(~n, labeller = as_labeller(c(
    `30` = "n=30",
    `50` = "n=50",
    `100` = "n=100"
  )))+
  #scale_color_gradient()+
  labs(x=expression(paste("Covariate ICC (",
               rho[x],
               ")",sep="")), 
       y=expression(sigma[4]^2),
       col=expression(rho['y|x']))


plot_grid(p1, p2, labels = c('A', 'B'),nrow=2)
ggsave(file = "/Users/laramaleyeff/Dropbox/dissertation/paper2/figures/icc_sigma_same_sign.eps",
       width=8,
       height=6,
       device=cairo_ps)

beta_4 <- -0.5
expgrid = expand.grid(n = c(30,50,100),
                      icc_x = c(0, 0.01,0.05,0.1,0.2,0.5,1),
                      sigma_gamma_sq = (pi^2/3*p_yx)/(1-p_yx)
)

expgrid$var = apply(expgrid,
                    1,
                    function(x) {
                      return(get_var_general(beta_1,
                                             beta_2,
                                             beta_3,
                                             beta_4,
                                             n = x[[1]],
                                             p_x = p_x,
                                             icc_x = x[[2]],
                                             bar_W,
                                             sigma_gamma_sq = x[[3]],
                                             sigma_sq_logit))
                    })
p1 <- expgrid %>%
  filter(icc_x != 1 & icc_x != 0) %>%
  mutate(p_yx = sigma_gamma_sq/(sigma_gamma_sq+pi^2/3)) %>%
  ggplot(aes(x=p_yx,y=var,group=factor(icc_x), col=factor(icc_x))) + geom_line() +
  facet_wrap(~n, labeller = as_labeller(c(
    `30` = "n=30",
    `50` = "n=50",
    `100` = "n=100"
  )))+
  #scale_color_gradient()+
  labs(col=expression(rho[x]), 
       y=expression(sigma[4]^2),
       x=expression(paste("Outcome ICC (", rho['y|x'],")", sep="")))
p2 <- expgrid %>%
  mutate(p_yx = round(sigma_gamma_sq/(sigma_gamma_sq+pi^2/3),2)) %>%
  filter(p_yx < 0.5 & p_yx > 0) %>%
  ggplot(aes(x=icc_x,y=var,group=factor(p_yx),col=factor(p_yx))) + 
  geom_line(linetype="dashed") +
  facet_wrap(~n, labeller = as_labeller(c(
    `30` = "n=30",
    `50` = "n=50",
    `100` = "n=100"
  )))+
  #scale_color_gradient()+
  labs(x=expression(paste("Covariate ICC (",
                          rho[x],
                          ")",sep="")), 
       y=expression(sigma[4]^2),
       col=expression(rho['y|x']))


plot_grid(p1, p2, labels = c('A', 'B'),nrow=2)
ggsave(file = "/Users/laramaleyeff/Dropbox/dissertation/paper2/figures/icc_sigma_opp_sign.eps",
       width=8,
       height=6,
       device=cairo_ps)