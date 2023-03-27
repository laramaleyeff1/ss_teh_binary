library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(xtable)
results <- list.files(path = "/Users/laramaleyeff/Dropbox/dissertation/paper2/rcode/sims/p2_sims/results",     
                            pattern = "", full.names = TRUE) %>% 
  lapply(read_csv) %>%   
  bind_rows  


set.seed(3)
sum_bin = results %>%
  filter(n!=30) %>%
  group_by(n, icc_x, K, beta_4, icc_yx) %>%
  sample_n(1000) %>%
  dplyr::summarize(power.empirical=sum(p.val <= 0.05)*100/n(),
                   power.ss=mean(power.ss),
                   type.1.error = sum(p.val.null <= 0.05)*100/n(),
                   sd.model = mean(sd.model),
                   sd.ss = mean(sd.ss),
                   sigma.gamma.sq.model = mean(sigma.gamma.sq.model),
                   n_sim=n()) 
quantile((sum_bin$power.empirical-sum_bin$power.ss), c(0.25,0.75))
mean((sum_bin$power.empirical-sum_bin$power.ss))

sum_bin$n_sim

View(sum_bin)

df_output = sum_bin %>%
  select(n, icc_x, icc_yx, K, beta_4, power.empirical, power.ss, type.1.error)%>%
  gather("type", "power", K, type.1.error,power.empirical, power.ss) %>%
  unite(temp, beta_4, type) %>%
  spread(temp, power) 
  

mean(abs(sum_bin$power.empirical-sum_bin$power.ss))
IQR(abs(sum_bin$power.empirical-sum_bin$power.ss))
mean(abs(sum_bin$power.empirical-sum_bin$power.ss))
median((sum_bin$power.empirical-sum_bin$power.ss))

nrow(sum_bin[abs(sum_bin$power.empirical-sum_bin$power.ss)<=2,])
nrow(sum_bin[abs(sum_bin$power.empirical-sum_bin$power.ss)>2,])

nrow(sum_bin[])

sum_bin %>%
  filter(icc_yx == 0.01,beta_4 ==0.25) %>%
  ggplot(aes(x=icc_x,y=K)) +
  geom_point()+
  geom_line()+
  facet_wrap(~n)

sum_bin_subset = sum_bin %>%
  filter(beta_4==0.75, n==50) 

df_output = sum_bin_subset %>%
  select(n, icc_x, icc_yx, K, beta_4, power.empirical, power.ss, type.1.error)%>%
  gather("type", "power", K, type.1.error,power.empirical, power.ss) %>%
  unite(temp, beta_4, type) %>%
  spread(temp, power) 

quantile((sum_bin_subset$power.empirical-sum_bin_subset$power.ss), c(0.25,0.75))
mean((sum_bin_subset$power.empirical-sum_bin_subset$power.ss))
print(xtable(df_output[,c(1:4,7,5:6)], digits = c(0,0,2,2, rep(c(0,1,1,1),1))), include.rownames=FALSE)

results %>%
  group_by(bar_W, cond_icc, beta_4, n, icc_x, K) %>%
  dplyr::summarize(power.empirical=sum(p.val <= 0.05)*100/n(),
                   power.ss=mean(power.ss),
                   sd.model = mean(sd.model),
                   sd.model.approx = mean(sd.model.approx),
                   sd.ss = mean(sd.ss),
                   sigma.gamma.sq.model = mean(sigma.gamma.sq.model),
                   n_sim=n()) %>%
  View()


  



old = results %>%
  group_by(n, K, p_x, icc_x,beta_4) %>%
  dplyr::summarize(power.empirical=sum(p.val <= 0.05)/n(),
                   power.logit = mean(power.logit),
                   n_sim=n()) 

old %>% View()

model1 <- glmer(Y ~ W*X  + (1|group), family="binomial", data=data)
fixef(model1) <- c(beta_1,beta_2,beta_3,beta_4)
powerSim(model1, test="W:X") 


results_cont <- list.files(path = "/Users/laramaleyeff/Dropbox/dissertation/paper2/rcode/sims/results_cont/",     
                      pattern = "", full.names = TRUE) %>% 
  lapply(read_csv) %>%   
  bind_rows  

sum_cont = results_cont %>%
  group_by(n,icc_x,K,beta_4,R,sigma_gamma_sq)%>%
  dplyr::summarize(power.empirical=sum(p.val <= 0.05)*100/n(),
                   power.ss=mean(power.ss),
                   sd.model = mean(sd.model),
                   sd.ss = mean(sd.ss),
                   n_sim=n()) 
sum_cont$n_sim
#View(sum_cont)

results_cont_0.1 = sum_cont %>%
  filter(R > 200 & sigma_gamma_sq == 0.1) %>%
  select(n, icc_x, K, beta_4, power.empirical, power.ss, R) %>%
  gather("type", "power", power.empirical, power.ss) %>%
  unite(temp, beta_4, type) %>%
  spread(temp, power) %>%
  select(-`0_power.ss`) %>%
  dplyr::rename( "type1error" = `0_power.empirical`,
          "power.empirical" = `0.5_power.empirical`,
          "power.estimated" = `0.5_power.ss`) 

mean(abs(results_cont_0.1$power.empirical-results_cont_0.1$power.estimated))
IQR(abs(results_cont_0.1$power.empirical-results_cont_0.1$power.estimated))
mean((results_cont_0.1$power.empirical-results_cont_0.1$power.estimated))
median((results_cont_0.1$power.empirical-results_cont_0.1$power.estimated))

quantile((results_cont_0.1$power.empirical-results_cont_0.1$power.estimated), c(0.25,0.75))
# 
# %>%
#   filter(power.estimated > 80) %>%
#   group_by(n, icc_x) %>% 
#   slice(which.min(power.estimated)) %>%
#   View()




View(results_cont_0.1)


print(xtable(results_cont_0.1[,-4], 
             digits = c(0,0,2,0,1,1,1)), include.rownames=FALSE)

  
sum_cont %>%

  gather("type", "power", power.empirical, power.ss) %>%
  ggplot(aes(x=K,y=power,group=interaction(type,beta_4),linetype=type)) +
  geom_line() +
  geom_hline(aes(yintercept = 80), col="red") +
  scale_linetype_manual(name = "",
                        values=c("solid", "dashed"), 
                        labels=c("Empirical power", "Estimated power")) +
  facet_wrap(~icc_x*n, scales = "free_x")+
  labs(x="Number of clusters (K)", y="Power")
