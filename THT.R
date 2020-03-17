library(tidyverse)
# (a)
theta<-rep(NA, 201)
ts<-rep(NA, 201)
theta[1]<-0
omega<-rnorm(201, 0, 1)
nu<-rnorm(201, 0, 1)

for (t in 1:200) {
  theta[t+1]<-0.9*theta[t]+omega[t+1]
  ts[t+1]<-theta[t+1]+nu[t+1]
}

ts<-ts[-1]
theta<-theta[-1]
ts<-ts-mean(ts)
plot(ts, type="l")
plot(theta, type="l")

# (b)
# placeholder
phi_p<-vector(mode="list", length = 15)
resid_p<-vector(mode="list", length = 15)
s_square_p<-vector(mode="list", length = 15)
AIC_p<-rep(NA, 15)
BIC_p<-rep(NA, 15)

n<-length(ts)

for(p in 1:15) {
  y<-ts[(p+1):(n)]
  y_prev<-matrix(NA, ncol = n-p, nrow = p)
  for(i in 1:p) {
    y_prev[i,]<-ts[(p-i+1):(n-i)]
  }
  
  phi_p[[p]]<-solve(y_prev %*% t(y_prev))%*%y_prev%*%y
  resid_p[[p]]<-y-t(y_prev)%*% phi_p[[p]]
  s_square_p[[p]]<-(t(resid_p[[p]]) %*%(resid_p[[p]]))/(n-p)
  AIC_p[[p]]<-2*(p)+(n-p)*log(s_square_p[[p]])
  BIC_p[[p]]<-log(n-p)*p+(n-p)*log(s_square_p[[p]])
}

AICBIC<-data.frame(p=1:15,AIC=AIC_p, BIC=BIC_p)
ggplot(AICBIC)+
  geom_point(aes(x=p, y=AIC), col="red", pch="a", cex=5)+
  geom_point(aes(x=p, y=BIC), col="blue", pch="b", cex=5)+
  ylab("Score")+
  xlab("Model order")+
  theme_classic()

phi_p[[2]]
s_square_p[[2]]

p<-2
# Reference prior
y<-ts[(p+1):(n)]
y_prev<-matrix(NA, ncol = n-p, nrow = p)
for(i in 1:p) {
  y_prev[i,]<-ts[(p-i+1):(n-i)]
}

beta_hat<-solve(y_prev %*% t(y_prev))%*%y_prev%*%y
resid<-y-t(y_prev)%*% beta_hat
s_square<-(t(resid) %*%(resid))/(n-p)

v_posterior<-LaplacesDemon::rinvgamma(500, (n-p)/2, (n-p)*s_square/2)
phi_posterior<-matrix(NA, nrow=500, ncol=p)
for (i in 1:500) {
  phi_posterior[i,]<-MASS::mvrnorm(1, mu = beta_hat, Sigma = v_posterior[i]*solve(y_prev %*% t(y_prev)))
}

# # Conjugate prior
# 
# # Prior
# m_0<-rep(0, p)
# C_0<-diag(p)
# n_0<-1
# d_0<-1/100
# # v_prior<-LaplacesDemon::rinvgamma(500, n_0/2, d_0/2)
# # hist(v_prior)
# 
# # Data
# y<-matrix(ts[(p+1):(n)])
# y_prev<-matrix(NA, ncol = n-p, nrow = p)
# for(i in 1:p) {
#   y_prev[i,]<-ts[(p-i+1):(n-i)]
# }
# 
# # Posterior
# 
# n_star<-length(y)+n_0
# e<-y-t(y_prev)%*%m_0
# Q<-t(y_prev)%*%C_0%*%y_prev+diag(length(y))
# 
# d_star<-t(e)%*%solve(Q)%*%e+d_0
# 
# v_posterior<-LaplacesDemon::rinvgamma(500, n_star/2, d_star/2)
# 
# m<-m_0+C_0%*%y_prev%*%solve(t(y_prev)%*%C_0%*%y_prev+diag(length(y)))%*%(y-t(y_prev)%*%m_0)
# C<-C_0-C_0%*%y_prev%*%solve(t(y_prev)%*%C_0%*%y_prev+diag(length(y)))%*%t(y_prev)%*%C_0
# 
# phi_posterior<-matrix(NA, nrow=500, ncol=p)
# for (i in 1:500) {
#   phi_posterior[i,]<-MASS::mvrnorm(1, mu = m, Sigma = v_posterior[i]*C)
# }


## summary
v_posterior_summary <- v_posterior_df %>% 
  summarize(mean=mean(v_posterior, na.rm=T),
            lower=quantile(v_posterior, prob=0.025, na.rm=T),
            upper=quantile(v_posterior, prob=0.975, na.rm=T))
# v_posterior_summary

v_posterior_df<-data.frame(v_posterior)
ggplot(v_posterior_df, aes(x=v_posterior))+
  geom_histogram()+
  geom_vline(data=v_posterior_summary, aes(xintercept=mean), col="blue")+
  geom_vline(data=v_posterior_summary, aes(xintercept=lower), col="blue", lty=2)+
  geom_vline(data=v_posterior_summary, aes(xintercept=upper), col="blue", lty=2)+
  theme_classic()
# summary(v_posterior)

phi_posterior_summary <- data.frame(phi_posterior) %>% 
  rename(phi_1=X1, phi_2=X2) %>% 
  gather() %>% 
  group_by(key) %>% 
  summarize(mean=mean(value, na.rm=T),
            lower=quantile(value, prob=0.025, na.rm=T),
            upper=quantile(value, prob=0.975, na.rm=T))
phi_posterior_summary

phi_posterior_df<-data.frame(phi_posterior) %>% 
  rename(phi_1=X1, phi_2=X2) %>% 
  gather()

ggplot(phi_posterior_df, aes(x=value))+
  geom_histogram()+
  geom_vline(data=phi_posterior_summary, aes(xintercept=mean), col="blue")+
  geom_vline(data=phi_posterior_summary, aes(xintercept=lower), col="blue", lty=2)+
  geom_vline(data=phi_posterior_summary, aes(xintercept=upper), col="blue", lty=2)+
  facet_wrap(. ~ key)+
  theme_classic()
summary(phi_posterior)

poly_coef<-c(1, -phi_posterior[201,])
poly_coef
recroots<-data.frame(recroot=1/polyroot(poly_coef)) %>% 
  mutate(im=Im(recroot)) %>% 
  mutate(ReCom=if_else(abs(im)>0.0000000001, "Complex", "Real")) %>% 
  # filter(abs(im)>0.0000000001) %>% # remove real roots 
  mutate(modulus=if_else(ReCom=="Complex", Mod(recroot), Re(recroot))) %>% 
  mutate(arg=abs(Arg(recroot))) %>% 
  mutate(wavelength=2*pi/arg) %>% 
  arrange(desc(wavelength)) %>% 
  arrange(ReCom) %>% 
  distinct(round(modulus, digits = 10),.keep_all = T)
recroots
# unlist(recroots %>% filter(ReCom=="Complex") %>% select(modulus))
# unlist(recroots %>% filter(ReCom=="Complex") %>% select(wavelength))
unlist(recroots %>% filter(ReCom=="Real") %>% select(modulus))

# complex_modulus_posterior<-data.frame(modulus_1=rep(NA, 500), modulus_2=rep(NA, 500), modulus_3=rep(NA, 500), modulus_4=rep(NA, 500))
# complex_wavelength_posterior<-data.frame(wavelength_1=rep(NA, 500), wavelength_2=rep(NA, 500), wavelength_3=rep(NA, 500), wavelength_4=rep(NA, 500))
real_modulus_posterior<-data.frame(modulus_1_r=rep(NA, 500), modulus_2_r=rep(NA, 500))

for(i in 1:500) {
  poly_coef<-c(1, -phi_posterior[i,])
  recroots<-data.frame(recroot=1/polyroot(poly_coef)) %>% 
    mutate(im=Im(recroot)) %>% 
    mutate(ReCom=if_else(abs(im)>0.0000000001, "Complex", "Real")) %>% 
    # filter(abs(im)>0.0000000001) %>% # remove real roots 
    mutate(modulus=if_else(ReCom=="Complex", Mod(recroot), Re(recroot))) %>% 
    mutate(arg=abs(Arg(recroot))) %>% 
    mutate(wavelength=2*pi/arg) %>% 
    arrange(desc(wavelength)) %>% 
    arrange(ReCom) %>% 
    distinct(round(modulus, digits = 10),.keep_all = T)
  # complex_modulus<-unlist(recroots %>% filter(ReCom=="Complex") %>% select(modulus))
  # length(complex_modulus)<-4
  # complex_wavelength<-unlist(recroots %>% filter(ReCom=="Complex") %>% select(wavelength))
  # length(complex_wavelength)<-4
  real_modulus<-unlist(recroots %>% filter(ReCom=="Real") %>% select(modulus))
  length(real_modulus)<-2
  
  # complex_modulus_posterior[i,]<-complex_modulus
  # complex_wavelength_posterior[i,]<-complex_wavelength
  real_modulus_posterior[i,]<-real_modulus
}

ggplot(real_modulus_posterior %>% gather(), aes(x=value))+
  geom_histogram()+
  facet_wrap(. ~ key)+
  theme_classic()
summary(real_modulus_posterior)
# 
# ggplot(complex_wavelength_posterior %>% gather(), aes(x=value))+
#   geom_histogram()+
#   facet_wrap(. ~ key)+
#   theme_classic()
# summary(complex_wavelength_posterior)

modulus_posterior_summary <- real_modulus_posterior %>% 
  gather() %>% 
  group_by(key) %>% 
  summarize(mean=mean(value, na.rm=T),
            lower=quantile(value, prob=0.025, na.rm=T),
            upper=quantile(value, prob=0.975, na.rm=T))
modulus_posterior_summary

# wavelength_posterior_summary <- complex_wavelength_posterior %>% 
#   gather() %>% 
#   group_by(key) %>% 
#   summarize(mean=mean(value, na.rm=T),
#             lower=quantile(value, prob=0.025, na.rm=T),
#             upper=quantile(value, prob=0.975, na.rm=T))
# wavelength_posterior_summary

y_pred_posterior<-t(y_prev)%*%t(phi_posterior)
y_pred<-rowMeans(y_pred_posterior)
plot(y, type="l")
lines(y_pred, type="l", col="blue")

y_resid<-y_pred-y

qqnorm(y_resid)
qqline(y_resid, col="red")

# (c)
AIC_arma<-matrix(NA, nrow=5, ncol=5)

for (p in 1:5) {
  for (q in 1:5) {
    AIC_arma[p, q]<-AIC(arima(ts, order = c(p, 0, q), include.mean = F))
  }
}

image(AIC_arma)

AIC_arma_df<-data.frame(AIC_arma) %>% 
  mutate(p=c(1:5)) %>% 
  gather("q", "AIC",-p) %>% 
  mutate(q=substr(q, 2,2)) %>% 
  mutate(q=as.numeric(q))
AIC_arma_df
ggplot(AIC_arma_df, aes(x=p, y=q, fill=AIC))+
  geom_tile()+
  scale_fill_viridis_c(option="magma", direction = -1)+
  theme_minimal()

arma_coef<-coef(arima(ts, order = c(2, 0, 2), include.mean = F))
arma_coef_sd<-sqrt(diag(vcov(arima(ts, order = c(2, 0, 2), include.mean = F))))
# sigma2<-arima(ts, order = c(2, 0, 2), include.mean = F)$sigma2
arma_summary<-data.frame(Coefficient=names(arma_coef), Mean=arma_coef, SD=arma_coef_sd)

ggplot(arma_summary)+
  geom_point(aes(x=Coefficient, y=Mean))+
  geom_errorbar(aes(x=Coefficient, ymin=Mean-1.96*SD,ymax=Mean+1.96*SD), width=0.2)+
  geom_hline(yintercept = 0, lty=2)+
  theme_classic()
# (d)
# ARMA (1, 1)
p<-2
q<-1
psi<-matrix(0, nrow=500, ncol=q)
for (k in 1:500) {
  for (i in 2:p) {
    psi[k,1]<-psi[k,1]+real_modulus_posterior[k, i]
  }
}

head(psi)
hist(psi[,1])

head(real_modulus_posterior)

psi_df<-data.frame(ar_1=modulus_posterior[,1], ma_1=psi[,1]) %>% 
  gather()
psi_df_summary<-psi_df %>%
  group_by(key) %>%
  summarize(mean=mean(value, na.rm=T),
            lower=quantile(value, prob=0.025, na.rm=T),
            upper=quantile(value, prob=0.975, na.rm=T))

psi_df_summary

ggplot(data.frame(psi_df), aes(x=value))+
  geom_histogram()+
  geom_vline(data=psi_df_summary, aes(xintercept=mean), col="blue")+
  geom_vline(data=psi_df_summary, aes(xintercept=lower), col="blue", lty=2)+
  geom_vline(data=psi_df_summary, aes(xintercept=upper), col="blue", lty=2)+
  facet_wrap(. ~ key)+
  theme_classic()


