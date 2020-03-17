#Chap 2 Problem 22

library(tidyverse)
# eeg<-scan("./STAT 223/eeg_F3.txt")
eeg<-scan("./Data/eeg_subset.txt")
str(eeg)
n<-length(eeg)

plot(eeg, type="l")
lines(fitted(loess(eeg~c(1:n))), col="red")

eeg<-eeg-mean(eeg) # adjust mean to 0

plot(eeg, type="l")
lines(fitted(loess(eeg~c(1:n))), col="red")

# (a)
p<-8
# Conjugate prior

# Prior
m_0<-rep(0, p)
C_0<-diag(p)
n_0<-1
d_0<-1/100
v_prior<-LaplacesDemon::rinvgamma(500, n_0/2, d_0/2)
hist(v_prior)

# Data
y<-matrix(eeg[(p+1):(n)])
y_prev<-matrix(NA, ncol = n-p, nrow = p)
for(i in 1:p) {
  y_prev[i,]<-eeg[(p-i+1):(n-i)]
}

# Posterior

n_star<-length(y)+n_0

e<-y-t(y_prev)%*%m_0
Q<-t(y_prev)%*%C_0%*%y_prev+diag(length(y))

d_star<-t(e)%*%solve(Q)%*%e+d_0

v_posterior<-LaplacesDemon::rinvgamma(500, n_star/2, d_star/2)

v_posterior_df<-data.frame(v_posterior)
ggplot(v_posterior_df, aes(x=v_posterior))+
  geom_histogram()+
  theme_classic()

m<-m_0+C_0%*%y_prev%*%solve(t(y_prev)%*%C_0%*%y_prev+diag(length(y)))%*%(y-t(y_prev)%*%m_0)
C<-C_0-C_0%*%y_prev%*%solve(t(y_prev)%*%C_0%*%y_prev+diag(length(y)))%*%t(y_prev)%*%C_0

phi_posterior<-matrix(NA, nrow=500, ncol=p)
for (i in 1:500) {
  phi_posterior[i,]<-MASS::mvrnorm(1, mu = m, Sigma = v_posterior[i]*C)
}

phi_posterior_df<-data.frame(phi_posterior) %>% 
  rename(phi_1=X1, phi_2=X2, phi_3=X3, phi_4=X4, phi_5=X5, phi_6=X6, phi_7=X7, phi_8=X8) %>% 
  gather()

ggplot(phi_posterior_df, aes(x=value))+
  geom_histogram()+
  facet_wrap(. ~ key)+
  theme_classic()
summary(phi_posterior)

# (b)
poly_coef<-c(1, -phi_posterior[201,])
poly_coef
recroots<-data.frame(recroot=1/polyroot(poly_coef)) %>% 
  mutate(im=Im(recroot)) %>% 
  mutate(ReCom=if_else(abs(im)>0.0000000001, "Complex", "Real")) %>% 
  # filter(abs(im)>0.0000000001) %>% # remove real roots 
  mutate(modulus=Mod(recroot)) %>% 
  mutate(arg=abs(Arg(recroot))) %>% 
  mutate(wavelength=2*pi/arg) %>% 
  arrange(desc(wavelength)) %>% 
  arrange(ReCom) %>% 
  distinct(round(modulus, digits = 10),.keep_all = T)
recroots
unlist(recroots %>% filter(ReCom=="Complex") %>% select(modulus))
unlist(recroots %>% filter(ReCom=="Complex") %>% select(wavelength))
unlist(recroots %>% filter(ReCom=="Real") %>% select(modulus))

complex_modulus_posterior<-data.frame(modulus_1=rep(NA, 500), modulus_2=rep(NA, 500), modulus_3=rep(NA, 500), modulus_4=rep(NA, 500))
complex_wavelength_posterior<-data.frame(wavelength_1=rep(NA, 500), wavelength_2=rep(NA, 500), wavelength_3=rep(NA, 500), wavelength_4=rep(NA, 500))
real_modulus_posterior<-data.frame(modulus_1_r=rep(NA, 500), modulus_2_r=rep(NA, 500))

for(i in 1:500) {
  poly_coef<-c(1, -phi_posterior[i,])
  recroots<-data.frame(recroot=1/polyroot(poly_coef)) %>% 
    mutate(im=Im(recroot)) %>% 
    mutate(ReCom=if_else(abs(im)>0.0000000001, "Complex", "Real")) %>% 
    # filter(abs(im)>0.0000000001) %>% # remove real roots 
    mutate(modulus=Mod(recroot)) %>% 
    mutate(arg=abs(Arg(recroot))) %>% 
    mutate(wavelength=2*pi/arg) %>% 
    arrange(desc(wavelength)) %>% 
    arrange(ReCom) %>% 
    distinct(round(modulus, digits = 10),.keep_all = T)
  complex_modulus<-unlist(recroots %>% filter(ReCom=="Complex") %>% select(modulus))
  length(complex_modulus)<-4
  complex_wavelength<-unlist(recroots %>% filter(ReCom=="Complex") %>% select(wavelength))
  length(complex_wavelength)<-4
  real_modulus<-unlist(recroots %>% filter(ReCom=="Real") %>% select(modulus))
  length(real_modulus)<-2
  
  complex_modulus_posterior[i,]<-complex_modulus
  complex_wavelength_posterior[i,]<-complex_wavelength
  real_modulus_posterior[i,]<-real_modulus
}

ggplot(complex_modulus_posterior %>% gather(), aes(x=value))+
  geom_histogram()+
  facet_wrap(. ~ key)+
  theme_classic()
summary(complex_modulus_posterior)

ggplot(complex_wavelength_posterior %>% gather(), aes(x=value))+
  geom_histogram()+
  facet_wrap(. ~ key)+
  theme_classic()
summary(complex_wavelength_posterior)

modulus_posterior_summary <- complex_modulus_posterior %>% 
  gather() %>% 
  group_by(key) %>% 
  summarize(mean=mean(value, na.rm=T),
            lower=quantile(value, prob=0.025, na.rm=T),
            upper=quantile(value, prob=0.975, na.rm=T))
modulus_posterior_summary

wavelength_posterior_summary <- complex_wavelength_posterior %>% 
  gather() %>% 
  group_by(key) %>% 
  summarize(mean=mean(value, na.rm=T),
            lower=quantile(value, prob=0.025, na.rm=T),
            upper=quantile(value, prob=0.975, na.rm=T))
wavelength_posterior_summary

# (c)
modulus_posterior_station<-cbind(complex_modulus_posterior, real_modulus_posterior) %>% 
  mutate(nonstationary=((modulus_1>1)+(modulus_2>1)+(modulus_3>1)+(modulus_4>1)+(modulus_1_r>1)+(modulus_2_r>1))>1)
sum(modulus_posterior_station$nonstationary, na.rm = T)

# (d)

y_pred_posterior<-t(y_prev)%*%t(phi_posterior)
y_pred<-rowMeans(y_pred_posterior)
plot(y, type="l")
lines(y_pred, type="l", col="blue")

y_resid<-y_pred-y

qqnorm(y_resid)
qqline(y_resid, col="red")

# (e)
# ARMA (2, 8)
q<-8
psi<-matrix(0, nrow=500, ncol=q)
for (k in 1:500) {
  for (i in 3:p) {
  psi[k,1]<-psi[k,1]+phi_posterior[k, i]
    for (j in 2:q) {
      psi[k,j]<-psi[k,j]+phi_posterior[k, i]*psi[k,(j-1)]
    }
  }
}

head(psi)
plot(psi[,1])

psi_df<-data.frame(psi) %>% 
  gather() %>% 
  group_by(key) %>% 
  summarize(mean=mean(value, na.rm=T),
            lower=quantile(value, prob=0.025, na.rm=T),
            upper=quantile(value, prob=0.975, na.rm=T)) %>% 
  mutate(index=1:8)

psi_df

ggplot(psi_df)+
  geom_point(aes(x=index, y=mean))+
  geom_errorbar(aes(x=index, ymin = lower, ymax = upper), width=0.2)+
  geom_hline(yintercept = 0)+
  ylab("coefficient")+
  theme_classic()
