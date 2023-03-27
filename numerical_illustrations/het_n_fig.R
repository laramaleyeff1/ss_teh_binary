library(ggplot2)
library(dplyr)
library(cowplot)
load("/Users/laramaleyeff/Dropbox/dissertation/paper2/rcode/numerical_illustrations/results_1000_vary.Rda")

p1 <- results_1000_vary %>%
  filter(n==50)%>%
  mutate(icc_yx = sigma_gamma_sq/(pi^2/3 + sigma_gamma_sq)) %>%
  ggplot(aes(x=icc_x, y=ratio, group=factor(cv_n),col=factor(cv_n))) +
  geom_point() +
  geom_line()+
  labs(col="CV",x=expression(paste("Covariate ICC (", rho[x] ,")", sep="")),y="Ratio", title=expression(paste(bar(n), "=50",sep="")))+
  facet_wrap(~icc_yx)

p2 <- results_1000_vary %>%
  filter(n==100) %>%
  mutate(icc_yx = sigma_gamma_sq/(pi^2/3 + sigma_gamma_sq)) %>%
  ggplot(aes(x=icc_x, y=ratio, group=factor(cv_n),col=factor(cv_n))) +
  geom_point() +
  geom_line()+
  labs(col="CV",x=expression(paste("Covariate ICC (", rho[x] ,")", sep="")),y="Ratio", title=expression(paste(bar(n), "=100",sep="")))+
  facet_wrap(~icc_yx)

plot_grid(p1, p2, labels=c("A","B"), nrow=2)

ggsave(file = "/Users/laramaleyeff/Dropbox/dissertation/paper2/figures/cv_vary_n.eps",
       width=8,
       height=6,
       device=cairo_ps)

results_1000_vary %>%
  mutate(icc_yx = sigma_gamma_sq/(pi^2/3 + sigma_gamma_sq)) %>%
  filter(n==50, icc_yx == 0.1)%>%
  mutate(icc_yx = sigma_gamma_sq/(pi^2/3 + sigma_gamma_sq)) %>%
  ggplot(aes(x=icc_x, y=ratio, group=factor(cv_n),col=factor(cv_n))) +
  geom_point() +
  geom_line()+
  labs(col="CV",x=expression(paste("Covariate ICC (", rho[x] ,")", sep="")),y="Ratio", title=expression(paste(bar(n), "=50",sep="")))


ggsave(file = "/Users/laramaleyeff/Dropbox/dissertation/paper2/figures/cv_vary_n_pres.jpg",
       width=4,
       height=5)
