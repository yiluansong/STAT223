library(tidyverse)
library(LaplacesDemon)
ts<-read_table("./Data/arplusnoise.txt", col_names =F)$X1

plot(ts, type="l")

mean(ts)

n<-length(ts)
p<-6
# Reference prior
y<-ts[(p+1):(n)]
y_prev<-matrix(NA, ncol = n-p, nrow = p)
for(i in 1:p) {
  y_prev[i,]<-ts[(p-i+1):(n-i)]
}

beta_hat<-solve(y_prev %*% t(y_prev))%*%y_prev%*%y
resid<-y-t(y_prev)%*% beta_hat
s_square<-(t(resid) %*%(resid))/(n-p)

fitted<-t(y_prev)%*% beta_hat
plot(y, type="l")
lines(fitted, col="blue")

v_posterior<-LaplacesDemon::rinvgamma(500, (n-p)/2, (n-p)*s_square/2)
phi_posterior<-matrix(NA, nrow=500, ncol=p)
for (i in 1:500) {
  phi_posterior[i,]<-MASS::mvrnorm(1, mu = beta_hat, Sigma = v_posterior[i]*solve(y_prev %*% t(y_prev)))
}

complex_modulus_posterior<-data.frame(modulus_1=rep(NA, 500), modulus_2=rep(NA, 500))
complex_wavelength_posterior<-data.frame(wavelength_1=rep(NA, 500), wavelength_2=rep(NA, 500))
real_modulus_posterior<-data.frame(reciprocal_root_1=rep(NA, 500), reciprocal_root_2=rep(NA, 500))

for(i in 1:500) {
  poly_coef<-c(1, -phi_posterior[i,])
  recroots<-data.frame(recroot=1/polyroot(poly_coef)) %>% 
    mutate(im=Im(recroot)) %>% 
    mutate(ReCom=if_else(abs(im)>0.0000000001, "Complex", "Real")) %>% 
    mutate(modulus=if_else(ReCom=="Complex", Mod(recroot), Re(recroot))) %>% 
    mutate(arg=abs(Arg(recroot))) %>% 
    mutate(wavelength=2*pi/arg) %>% 
    arrange(desc(modulus)) %>% 
    arrange(ReCom) %>% 
    distinct(round(modulus, digits = 10),.keep_all = T)
  complex_modulus<-unlist(recroots %>% filter(ReCom=="Complex") %>% select(modulus))
  length(complex_modulus)<-2
  complex_wavelength<-unlist(recroots %>% filter(ReCom=="Complex") %>% select(wavelength))
  length(complex_wavelength)<-2
  real_modulus<-unlist(recroots %>% filter(ReCom=="Real") %>% select(modulus))
  length(real_modulus)<-2
  
  complex_modulus_posterior[i,]<-complex_modulus
  complex_wavelength_posterior[i,]<-complex_wavelength
  real_modulus_posterior[i,]<-real_modulus
}

# ggplot(real_modulus_posterior %>% gather(), aes(x=value))+
#   geom_histogram()+
#   facet_wrap(. ~ key)+
#   theme_classic()
# 
# ggplot(complex_modulus_posterior %>% gather(), aes(x=value))+
#   geom_histogram()+
#   facet_wrap(. ~ key)+
#   theme_classic()
# 
# ggplot(complex_wavelength_posterior %>% gather(), aes(x=value))+
#   geom_histogram()+
#   facet_wrap(. ~ key)+
#   theme_classic()

real_modulus_posterior_df<-real_modulus_posterior %>% 
  gather()
real_modulus_posterior_summary <- real_modulus_posterior_df%>% 
  group_by(key) %>% 
  summarize(median=quantile(value, prob=0.5, na.rm=T),
            lower=quantile(value, prob=0.025, na.rm=T),
            upper=quantile(value, prob=0.975, na.rm=T))
ggplot(real_modulus_posterior_df, aes(x=value))+
  geom_histogram()+
  geom_vline(data=real_modulus_posterior_summary, aes(xintercept=median), col="blue")+
  geom_vline(data=real_modulus_posterior_summary, aes(xintercept=lower), col="blue", lty=2)+
  geom_vline(data=real_modulus_posterior_summary, aes(xintercept=upper), col="blue", lty=2)+
  facet_wrap(. ~ key)+
  theme_classic()


complex_modulus_posterior_df<-complex_modulus_posterior %>% 
  gather()
complex_modulus_posterior_summary <- complex_modulus_posterior_df%>% 
  group_by(key) %>% 
  summarize(median=quantile(value, prob=0.5, na.rm=T),
            lower=quantile(value, prob=0.025, na.rm=T),
            upper=quantile(value, prob=0.975, na.rm=T))
ggplot(complex_modulus_posterior_df, aes(x=value))+
  geom_histogram()+
  geom_vline(data=complex_modulus_posterior_summary, aes(xintercept=median), col="blue")+
  geom_vline(data=complex_modulus_posterior_summary, aes(xintercept=lower), col="blue", lty=2)+
  geom_vline(data=complex_modulus_posterior_summary, aes(xintercept=upper), col="blue", lty=2)+
  facet_wrap(. ~ key)+
  theme_classic()


complex_wavelength_posterior_df<-complex_wavelength_posterior %>% 
  gather()
complex_wavelength_posterior_summary <- complex_wavelength_posterior_df%>% 
  group_by(key) %>% 
  summarize(median=quantile(value, prob=0.5, na.rm=T),
            lower=quantile(value, prob=0.025, na.rm=T),
            upper=quantile(value, prob=0.975, na.rm=T))
ggplot(complex_wavelength_posterior_df, aes(x=value))+
  geom_histogram()+
  geom_vline(data=complex_wavelength_posterior_summary, aes(xintercept=median), col="blue")+
  geom_vline(data=complex_wavelength_posterior_summary, aes(xintercept=lower), col="blue", lty=2)+
  geom_vline(data=complex_wavelength_posterior_summary, aes(xintercept=upper), col="blue", lty=2)+
  facet_wrap(. ~ key)+
  theme_classic()


###
# Prado pp 86

y<-ts

period<-12
omega<-2*pi/period

F_har<-matrix(nrow=2, ncol=n)
for(t in 1:n) {
  F_har[,t]<-c(cos(omega*t),sin(omega*t))
}

# seq(4, n-1, length.out = 50)
omega.grid<-seq(0+0.001, pi-0.001, length.out = 500)
period.grid<-2*pi/omega.grid

  
log_likelihood<-rep(NA, 500)
for(i in 1:500) {
  omega<-omega.grid[[i]]
  F_har<-matrix(nrow=2, ncol=n)
  for(t in 1:n) {
    F_har[,t]<-c(cos(omega*t),sin(omega*t))
  }
  beta_hat<-solve(F_har%*%t(F_har))%*%F_har%*%y
  resid<-y-t(F_har)%*% beta_hat
  R<-as.numeric(t(resid) %*%(resid))
  
  log_likelihood[[i]]<- -1/2*log(det(F_har%*%t(F_har)))+(2-n)/2*log(R)
}
  
plot(x=omega.grid, y=log_likelihood, type="l")

plot(x=period.grid, y=log_likelihood, type="l", xlim=c(0, 100))

# exp(log_likelihood+1205)

omega_posterior<-sample(x = omega.grid, size = 500, 
       prob = exp(log_likelihood+1205), replace=T)

hist(omega_posterior)

v_posterior<-rep(NA, 500)
beta_posterior<-data.frame(a=rep(NA, 500), b=rep(NA, 500))
y_pred_posterior<-matrix(nrow=n, ncol=500)

for (i in 1:500) {
  omega<-omega_posterior[[i]]
  F_har<-matrix(nrow=2, ncol=n)
  for(t in 1:n) {
    F_har[,t]<-c(cos(omega*t),sin(omega*t))
  }
  beta_hat<-solve(F_har%*%t(F_har))%*%F_har%*%y
  resid<-y-t(F_har)%*% beta_hat
  s_square<-as.numeric(t(resid) %*%(resid))/(n-2)
  
  v_posterior[[i]]<-LaplacesDemon::rinvgamma(1, (n-2)/2, (n-2)*s_square/2)
  beta_posterior[i,]<-rmvt(n=1, mu=as.numeric(beta_hat), S=s_square*round(solve(F_har %*% t(F_har)),10), df=n-2)
  
  y_pred_posterior[,i]<-t(F_har)%*%matrix(as.numeric(beta_posterior[i,]),ncol=1)+rnorm(n, sd=v_posterior[[i]])
}


hist(v_posterior)

beta_posterior_df<-beta_posterior %>% 
  gather()
beta_posterior_summary <- beta_posterior_df%>% 
  group_by(key) %>% 
  summarize(median=quantile(value, prob=0.5, na.rm=T),
            lower=quantile(value, prob=0.025, na.rm=T),
            upper=quantile(value, prob=0.975, na.rm=T))
ggplot(beta_posterior_df, aes(x=value))+
  geom_histogram()+
  geom_vline(data=beta_posterior_summary, aes(xintercept=median), col="blue")+
  geom_vline(data=beta_posterior_summary, aes(xintercept=lower), col="blue", lty=2)+
  geom_vline(data=beta_posterior_summary, aes(xintercept=upper), col="blue", lty=2)+
  facet_wrap(. ~ key)+
  theme_classic()

y_pred_posterior_df<-data.frame(y_pred_posterior) %>% 
  mutate(index=row_number()) %>% 
  gather(key = "run", value="value",-index) %>% 
  group_by(index) %>% 
  summarize(median=quantile(value, prob=0.5, na.rm=T),
            lower=quantile(value, prob=0.025, na.rm=T),
            upper=quantile(value, prob=0.975, na.rm=T)) %>% 
  mutate(y=y) 
  
head(y_pred_posterior_df)

ggplot(data = y_pred_posterior_df)+
  geom_line(aes(x=index, y=y))+
  geom_line(aes(x=index, y=median), col="blue")+
  geom_line(aes(x=index, y=lower), col="blue", lty=2)+
  geom_line(aes(x=index, y=upper), col="blue", lty=2)+
  theme_classic()

### pp 104
phi_estimate<-data.frame(phi_posterior) %>% 
  gather() %>% 
  group_by(key) %>% 
  summarize(median=median(value)) %>% 
  select(median) %>% 
  unlist() %>% 
  as.numeric()

spectral.grid<-
  var(y)/2/pi*
  Mod(1
      -phi_estimate[1]*complex(length.out = 1, argument = -omega.grid)
      -phi_estimate[2]*complex(length.out = 1, argument = -2*omega.grid)
      -phi_estimate[3]*complex(length.out = 1, argument = -3*omega.grid)
      -phi_estimate[4]*complex(length.out = 1, argument = -4*omega.grid)
      -phi_estimate[5]*complex(length.out = 1, argument = -5*omega.grid)
      -phi_estimate[6]*complex(length.out = 1, argument = -6*omega.grid))^(-2)

spectral.grid2<-
  var(y)/2/pi*
  Mod(1
      -0.95*complex(length.out = 1, argument = -omega.grid)
      +0.9*complex(length.out = 1, argument = -2*omega.grid))^(-2)

# phi_1<- 0.95
# phi_2<- -0.9
# spectral.grid2<-
#   var(y)/2/pi/(1+phi_1^2+2*phi_2+phi_2^2+2*(phi_1*phi_2-phi_1)*cos(omega.grid)-4*phi_2*(cos(omega.grid))^2)

par(mfrow=c(3,1))
plot(x=omega.grid, y=log_likelihood+1205, type="l")
plot(x=omega.grid, y=log(spectral.grid), type="l")
plot(x=omega.grid, y=log(spectral.grid2), type="l")
par(mfrow=c(1,1))

###

N<-length(ts)
x<-matrix(ncol=N+1, nrow = 1000)
phi_1<-rep(NA, 1000)
phi_2<-rep(NA, 1000)
phi_1[[1]]<- 0.5
phi_2[[1]]<- 0.5

m<-vector(mode = "list", length=N+1)
m[[1]]<-matrix(rep(1, 2), ncol=1)

C<-vector(mode = "list", length=N+1)
library(ltsa)
gamma_0<-tacvfARMA(phi = c(phi_1[[1]], -phi_2[[1]]),  maxLag = 1, sigma2 = 1)[1]
gamma_1<-tacvfARMA(phi = c(phi_1[[1]], -phi_2[[1]]),  maxLag = 1, sigma2 = 1)[2]
C[[1]]<-matrix(c(gamma_0, gamma_1, gamma_1, gamma_0), ncol=2)

a<-vector(mode = "list", length=N+1)
R<-vector(mode = "list", length=N+1)
f<-vector(mode = "list", length=N+1)
Q<-vector(mode = "list", length=N+1)
e<-vector(mode = "list", length=N+1)
A<-vector(mode = "list", length=N+1)
B<-vector(mode="list", length=N+1)

a_back<-vector(mode="list", length=N+1)
R_back<-vector(mode="list", length=N+1)


for(j in 1:1000) {
  ### FFBS
  F_dlm<-matrix(c(1,0), ncol=1)
  F_dlm
  
  G_dlm<-matrix(c(phi_1[j], 1, -phi_2[j], 0), ncol = 2)
  G_dlm
  
  W_dlm<-matrix(c(1, 0, 0, 0), ncol = 2)
  v_dlm<-4
  
  ## filtering
  for(i in 1:N) {
    a[[i+1]]<-G_dlm%*%m[[i]]
    R[[i+1]]<-G_dlm%*%C[[i]]%*%t(G_dlm)+W_dlm
    f[[i+1]]<-as.numeric(t(F_dlm)%*%a[[i+1]])
    Q[[i+1]]<-as.numeric(t(F_dlm)%*%R[[i+1]]%*%F_dlm+v_dlm)
    e[[i+1]]<-ts[[i]]-f[[i+1]]
    A[[i+1]]<-R[[i+1]]%*%F_dlm/Q[[i+1]]
    m[[i+1]]<-a[[i+1]]+A[[i+1]]*e[[i+1]]
    C[[i+1]]<-R[[i+1]]-A[[i+1]]%*%t(A[[i+1]])*Q[[i+1]]
    
    B[[i]]<-C[[i]]%*%t(G_dlm)%*%solve(R[[i+1]])
  }
  
  ### smoothing
  
  a_back[[N+1]]<-m[[N+1]]
  R_back[[N+1]]<-C[[N+1]]
  
  for (i in N:1) {
    a_back[[i]]<-m[[i]]-B[[i]]%*%(a[[i+1]]-a_back[[i+1]])
    R_back[[i]]<-C[[i]]-B[[i]]%*%(R[[i+1]]-R_back[[i+1]])%*%t(B[[i]])
  }
  
  m[[1]]<-a_back[[1]]
  C[[1]]<-R_back[[1]]
  
  theta<-matrix(NA, ncol=N+1, nrow = 2)
  for (i in 1:(N+1)) {
    # theta[,i]<-rmvn(mu=as.numeric(a_back[[i]]), Sigma = round(R_back[[i]],10))
    theta[,i]<-as.numeric(a_back[[i]])
  }
  
  x[j,]<-theta[1,]
  ### sample phi
  p<-2
  y<-theta[1,(p+1+1):(N+1)]
  F_ar<-matrix(NA, ncol = N-2, nrow = p)
  F_ar[1,]<- theta[2,(p+1+1):(N+1)]
  F_ar[2,]<- -theta[2,(p-1+1+1):(N-1+1)]
  
  beta_hat<-solve(F_ar %*% t(F_ar))%*%F_ar%*%y
  beta_var<-1*solve(F_ar %*% t(F_ar))
  
  phi_sample<-rmvn(mu=as.numeric(beta_hat), Sigma = round(beta_var,10))
  phi_1[j+1]<-phi_sample[1]
  phi_2[j+1]<-phi_sample[2]
}

plot(phi_1, type = "l")
summary(phi_1[501:1001])
hist(phi_1[501:1001])
acf(phi_1[501:1001])
pacf(phi_1[501:1001])

plot(phi_2, type = "l")
summary(phi_2[501:1001])
hist(phi_2[501:1001])
acf(phi_2[501:1001])
pacf(phi_2[501:1001])

plot(ts, type = "l")
lines(x[500,2:(n+1)], col="blue")

