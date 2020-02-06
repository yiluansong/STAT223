library(tidyverse)
# eeg<-scan("./STAT 223/eeg_F3.txt")
eeg<-scan("./STAT 223/eeg_subset.txt")
str(eeg)
plot(eeg, type="l")

n<-length(eeg)
AIC_p<-rep(NA, 25)
BIC_p<-rep(NA, 25)
phi_p<-vector(mode="list", length = 25)
s_square_p<-vector(mode="list", length = 25)
for(p in 1:25) {
  y<-eeg[(p+1):(n)]
  y_prev<-matrix(NA, ncol = n-p, nrow = p)
  for(i in 1:p) {
    y_prev[i,]<-iso_detrend[(p-i+1):(n-i)]
  }
  
  phi_p[[p]]<-solve(y_prev %*% t(y_prev))%*%y_prev%*%y
  s_square_p[[p]]<-(t(y-t(y_prev)%*% phi_p[[p]]) %*%(y-t(y_prev)%*% phi_p[[p]]))/(n-p)
  AIC_p[[p]]<- -(n-p)/2*log2(s_square_p[[p]]*(n-p))
  # BIC_p[[p]]<-log10(n)*p+n*log10(s_square_p[[p]])
}
plot(AIC_p)
y<-eeg[9:(length(eeg))]
str(y)


phi<-rep(NA, 8)

y_back1<-eeg[8:(length(eeg)-1)]
y_back2<-eeg[7:(length(eeg)-2)]
y_back3<-eeg[6:(length(eeg)-3)]
y_back4<-eeg[5:(length(eeg)-4)]
y_back5<-eeg[4:(length(eeg)-5)]
y_back6<-eeg[3:(length(eeg)-6)]
y_back7<-eeg[2:(length(eeg)-7)]
y_back8<-eeg[1:(length(eeg)-8)]

y_prev<-rbind(y_back1, y_back2, y_back3, y_back4, y_back5, y_back6, y_back7, y_back8)

dim(y_prev)

phi_hat<- solve(y_prev %*% t(y_prev))%*%y_prev%*%y

R<-t(y-t(y_prev)%*% phi_hat) %*%(y-t(y_prev)%*% phi_hat)

s_square<-as.numeric(R/(length(y)-8))
s_square
v_hat<-R/length(y)
v_hat

# reference prior
library(MASS)
phi_posterior<-mvrnorm(n = 500, mu=phi_hat, Sigma=s_square*solve(y_prev%*%t(y_prev)))

library(LaplacesDemon)
phi_posterior<-rmvt(n=500, mu=phi_hat, S=s_square*solve(y_prev%*%t(y_prev)), df=length(y)-8)
# ?????

library(invgamma)
v_posterior<-invgamma::rinvgamma(n=500, shape=(length(y)-8)/2, rate=(length(y)-8)*s_square/2 )
head(v_posterior)

phi_posterior_df<-data.frame(phi_posterior) %>% 
  rename(phi_1=y_back1, phi_2=y_back2, phi_3=y_back3, phi_4=y_back4, phi_5=y_back5, phi_6=y_back6, phi_7=y_back7, phi_8=y_back8) %>% 
  gather()

ggplot(phi_posterior_df, aes(x=value))+
  geom_histogram()+
  facet_wrap(. ~ key)+
  theme_classic()
summary(phi_posterior)

# roots
poly_coef<-c(1, -phi_posterior[1,])
poly_coef
roots<-data.frame(root=polyroot(poly_coef)) %>% 
  # mutate(re=Re(root)) %>% 
  mutate(modulus=Mod(root)) %>% 
  mutate(arg=abs(Arg(root))) %>% 
  mutate(wavelength=2*pi/arg) %>% 
  arrange(modulus)
roots

modulus_posterior<-data.frame(modulus_1=rep(NA, 500), modulus_2=rep(NA, 500), modulus_3=rep(NA, 500), modulus_4=rep(NA, 500))
wavelength_posterior<-data.frame(wavelength_1=rep(NA, 500), wavelength_2=rep(NA, 500), wavelength_3=rep(NA, 500), wavelength_4=rep(NA, 500))

for(i in 1:500) {
  poly_coef<-c(1, -phi_posterior[i,])
  roots<-data.frame(root=polyroot(poly_coef)) %>% 
    # mutate(re=Re(root)) %>% 
    mutate(modulus=Mod(root)) %>% 
    mutate(arg=abs(Arg(root))) %>% 
    mutate(wavelength=2*pi/arg) %>% 
    arrange(desc(modulus))
  modulus_posterior[i,]<-roots$modulus[c(1, 3, 5, 7)]
  wavelength_posterior[i,]<-roots$wavelength[c(1, 3, 5, 7)]
}

modulus_posterior
modulus_posterior_df<-data.frame(modulus_posterior) %>% 
  gather()

ggplot(modulus_posterior_df, aes(x=value))+
  geom_histogram()+
  facet_wrap(. ~ key)+
  theme_classic()
summary(modulus_posterior)

wavelength_posterior
wavelength_posterior_df<-data.frame(wavelength_posterior) %>% 
  gather()

ggplot(wavelength_posterior_df, aes(x=value))+
  geom_histogram()+
  facet_wrap(. ~ key)+
  theme_classic()
summary(wavelength_posterior)


v_posterior
v_posterior_df<-data.frame(v_posterior) %>% 
  gather()

ggplot(v_posterior_df, aes(x=value))+
  geom_histogram()+
  facet_wrap(. ~ key)+
  theme_classic()
summary(v_posterior)
