library(tidyverse)
library(LaplacesDemon)
flu<-read_csv("./Data/fludata.csv", col_names =F)$X2[1:192]


plot(flu, type="l")

F_dlm<-matrix(c(1,0), ncol=1)
for (r in 1:6) {
  F_dlm<-rbind(F_dlm,matrix(c(1,0), ncol = 1))
}
F_dlm<-rbind(F_dlm,matrix(c(1)))

F_dlm

G_dlm<-matrix(c(1, 0, 1, 1), ncol = 2)
p<-192
for (r in c(2,4, 8,16, 32, 64)) {
  G_dlm<-bdiag(G_dlm,matrix(c(cos(2*pi*r/p), -sin(2*pi*r/p), sin(2*pi*r/p), cos(2*pi*r/p)), ncol = 2))
}
G_dlm<-bdiag(G_dlm,matrix(c(-1)))
G_dlm

N<-length(flu)

delta.grid<-seq(0.99,1, length.out = 50)
delta.log.lik<-rep(NA, 50)

for (j in 1:50) {
  delta<-delta.grid[[j]]
  
  m<-vector(mode = "list", length=N+1)
  m[[1]]<-matrix(rep(0, 15), ncol=1)
  
  C_star<-vector(mode = "list", length=N+1)
  C_star[[1]]<-diag(rep(1, 15))
  
  n<-vector(mode = "list", length=N+1)
  n[[1]]<-1
  
  S<-vector(mode = "list", length=N+1)
  S[[1]]<-10
  
  a<-vector(mode = "list", length=N+1)
  R_star<-vector(mode = "list", length=N+1)
  f<-vector(mode = "list", length=N+1)
  Q_star<-vector(mode = "list", length=N+1)
  e<-vector(mode = "list", length=N+1)
  A<-vector(mode = "list", length=N+1)
  
  R<-vector(mode = "list", length=N+1)
  Q<-vector(mode = "list", length=N+1)
  C<-vector(mode = "list", length=N+1)
  
  likelihood<-vector(mode = "list", length=N+1)
  
  for(i in 1:N) {
    a[[i+1]]<-G_dlm%*%m[[i]]
    R_star[[i+1]]<-G_dlm%*%C_star[[i]]%*%t(G_dlm)/delta
    f[[i+1]]<-as.numeric(t(F_dlm)%*%a[[i+1]])
    Q_star[[i+1]]<-as.numeric(1+t(F_dlm)%*%R_star[[i+1]]%*%F_dlm)
    e[[i+1]]<-flu[[i]]-f[[i+1]]
    A[[i+1]]<-R_star[[i+1]]%*%F_dlm/Q_star[[i+1]]
    m[[i+1]]<-a[[i+1]]+A[[i+1]]*e[[i+1]]
    C_star[[i+1]]<-R_star[[i+1]]-A[[i+1]]%*%t(A[[i+1]])*Q_star[[i+1]]
    
    n[[i+1]]<-n[[i]]+1
    S[[i+1]]<-(n[[i]]*S[[i]]+e[[i+1]]^2/Q_star[[i+1]])/n[[i+1]]
    
    R[[i+1]]<-S[[i]]*R_star[[i+1]]
    Q[[i+1]]<-S[[i]]*Q_star[[i+1]]
    C[[i+1]]<-S[[i+1]]*C_star[[i+1]]
    
    likelihood[[i]]<-dst(flu[[i]], mu=f[[i+1]], sigma=Q[[i+1]], nu=n[[i]], log=T)
  }
  delta.log.lik[[j]]<-sum(unlist(cbind(likelihood)))
}

plot(unlist(cbind(likelihood)), type="l")
plot(unlist(cbind(delta.log.lik)), type="l")


#####
delta<-1

m<-vector(mode = "list", length=N+1)
m[[1]]<-matrix(rep(0, 15), ncol=1)

C_star<-vector(mode = "list", length=N+1)
C_star[[1]]<-diag(rep(1, 15))

n<-vector(mode = "list", length=N+1)
n[[1]]<-1

S<-vector(mode = "list", length=N+1)
S[[1]]<-100

a<-vector(mode = "list", length=N+1)
R_star<-vector(mode = "list", length=N+1)
f<-vector(mode = "list", length=N+1)
Q_star<-vector(mode = "list", length=N+1)
e<-vector(mode = "list", length=N+1)
A<-vector(mode = "list", length=N+1)

R<-vector(mode = "list", length=N+1)

Q<-vector(mode = "list", length=N+1)

C<-vector(mode = "list", length=N+1)
C[[1]]<-S[[1]]*C_star[[1]]
# likelihood<-vector(mode = "list", length=N+1)

for(i in 1:N) {
  a[[i+1]]<-G_dlm%*%m[[i]]
  R_star[[i+1]]<-G_dlm%*%C_star[[i]]%*%t(G_dlm)/delta
  f[[i+1]]<-as.numeric(t(F_dlm)%*%a[[i+1]])
  Q_star[[i+1]]<-as.numeric(1+t(F_dlm)%*%R_star[[i+1]]%*%F_dlm)
  e[[i+1]]<-flu[[i]]-f[[i+1]]
  A[[i+1]]<-R_star[[i+1]]%*%F_dlm/Q_star[[i+1]]
  m[[i+1]]<-a[[i+1]]+A[[i+1]]*e[[i+1]]
  C_star[[i+1]]<-R_star[[i+1]]-A[[i+1]]%*%t(A[[i+1]])*Q_star[[i+1]]
  
  n[[i+1]]<-n[[i]]+1
  S[[i+1]]<-(n[[i]]*S[[i]]+e[[i+1]]^2/Q_star[[i+1]])/n[[i+1]]
  
  R[[i+1]]<-S[[i]]*R_star[[i+1]]
  Q[[i+1]]<-S[[i]]*Q_star[[i+1]]
  C[[i+1]]<-S[[i+1]]*C_star[[i+1]]
  # 
  # likelihood[[i]]<-dst(flu[[i]], mu=f[[i+1]], sigma=Q[[i+1]], nu=n[[i]], log=T)
}


plot(flu, type="l")
lines(unlist(cbind(f)), col="blue")

#####

theta_summary<-vector(mode = "list", length=N)
for (i in 1:N) {
  theta_sample<-rmvt(n=100, mu=as.numeric(m[[i+1]]), S=round(C[[i+1]],8), df=n[[i+1]]) # sometimes not symmetric? computational issue?

theta_sample<-data.frame(theta_sample)
colnames(theta_sample)<-paste("theta_",seq(1:15), sep="")

theta_summary[[i]]<-theta_sample %>%
  gather(key="param") %>% 
  group_by(param) %>% 
  summarize(median=quantile(value, probs = 0.5),
            lower=quantile(value, probs = 0.025),
            upper=quantile(value, probs = 0.975)) %>% 
  mutate(index=i) %>% 
  ungroup()
}


theta_summary<-bind_rows(theta_summary)

# trend
ggplot(theta_summary %>% filter(param=="theta_1") %>% mutate(flu=flu))+
  geom_line(aes(x=index, y=flu), col="grey")+
  geom_line(aes(x=index, y=median))+
  geom_line(aes(x=index, y=lower), lty=2)+
  geom_line(aes(x=index, y=upper), lty=2)+
  theme_classic()

# harmonic 1
ggplot(theta_summary %>% filter(param=="theta_3") %>% mutate(flu=flu))+
  geom_line(aes(x=index, y=flu), col="grey")+
  geom_line(aes(x=index, y=median))+
  geom_line(aes(x=index, y=lower), lty=2)+
  geom_line(aes(x=index, y=upper), lty=2)+
  theme_classic()

# harmonic 2
ggplot(theta_summary %>% filter(param=="theta_5") %>% mutate(flu=flu))+
  geom_line(aes(x=index, y=flu), col="grey")+
  geom_line(aes(x=index, y=median))+
  geom_line(aes(x=index, y=lower), lty=2)+
  geom_line(aes(x=index, y=upper), lty=2)+
  theme_classic()

# harmonic 3
ggplot(theta_summary %>% filter(param=="theta_7") %>% mutate(flu=flu))+
  geom_line(aes(x=index, y=flu), col="grey")+
  geom_line(aes(x=index, y=median))+
  geom_line(aes(x=index, y=lower), lty=2)+
  geom_line(aes(x=index, y=upper), lty=2)+
  theme_classic()

# harmonic 4
ggplot(theta_summary %>% filter(param=="theta_9") %>% mutate(flu=flu))+
  geom_line(aes(x=index, y=flu), col="grey")+
  geom_line(aes(x=index, y=median))+
  geom_line(aes(x=index, y=lower), lty=2)+
  geom_line(aes(x=index, y=upper), lty=2)+
  theme_classic()

# harmonic 5
ggplot(theta_summary %>% filter(param=="theta_11") %>% mutate(flu=flu))+
  geom_line(aes(x=index, y=flu), col="grey")+
  geom_line(aes(x=index, y=median))+
  geom_line(aes(x=index, y=lower), lty=2)+
  geom_line(aes(x=index, y=upper), lty=2)+
  theme_classic()

# harmonic 6
ggplot(theta_summary %>% filter(param=="theta_13") %>% mutate(flu=flu))+
  geom_line(aes(x=index, y=flu), col="grey")+
  geom_line(aes(x=index, y=median))+
  geom_line(aes(x=index, y=lower), lty=2)+
  geom_line(aes(x=index, y=upper), lty=2)+
  theme_classic()

# nyquist
ggplot(theta_summary %>% filter(param=="theta_15") %>% mutate(flu=flu))+
  geom_line(aes(x=index, y=flu), col="grey")+
  geom_line(aes(x=index, y=median))+
  geom_line(aes(x=index, y=lower), lty=2)+
  geom_line(aes(x=index, y=upper), lty=2)+
  theme_classic()

composite<-
  theta_summary %>% filter(param=="theta_1") %>% select(median)+
  theta_summary %>% filter(param=="theta_3") %>% select(median)+
  theta_summary %>% filter(param=="theta_5") %>% select(median)+
  theta_summary %>% filter(param=="theta_7") %>% select(median)+
  theta_summary %>% filter(param=="theta_9") %>% select(median)+
  theta_summary %>% filter(param=="theta_11") %>% select(median)+
  theta_summary %>% filter(param=="theta_13") %>% select(median)+
  theta_summary %>% filter(param=="theta_15") %>% select(median)

plot(flu, type="l")
lines(unlist(cbind(f)), col="blue")
lines(composite$median, col="grey")


# one-step ahead

##### smoothing pp130
a_back<-vector(mode="list", length=N+1)
a_back[[N+1]]<-m[[N+1]] #Check!

R_back_star<-vector(mode="list", length=N+1)
R_back_star[[N+1]]<-C[[N+1]]
R_back<-vector(mode="list", length=N+1)
R_back[[N+1]]<-C[[N+1]]

f_back<-vector(mode="list", length=N+1)
f_back[[N+1]]<-f[[N+1]]

Q_back_star<-vector(mode="list", length=N+1)
Q_back_star[[N+1]]<-Q[[N+1]]
Q_back<-vector(mode="list", length=N+1)
Q_back[[N+1]]<-Q[[N+1]]

for (i in N:1) {
  a_back[[i]]<-(1-delta)*m[[i]]+delta*solve(G_dlm)%*%a_back[[i+1]]
  R_back_star[[i]]<-(1-delta)*C[[i]]+delta^2*solve(G_dlm)%*%R_back_star[[i+1]]%*%solve(t(G_dlm))
  R_back[[i]]<-R_back_star[[i]]*S[[N+1]]/S[[i]]
  
  f_back[[i]]<-t(F_dlm)%*%a_back[[i]]
  Q_back_star[[i]]<-t(F_dlm)%*%R_back_star[[i]]%*%F_dlm
  Q_back[[i]]<-Q_back_star[[i]]*S[[N+1]]/S[[i]] # pp 115
}

plot(flu, type="l")
lines(unlist(cbind(f_back))[2:(N+1)], col="blue")

#####

theta_back_summary<-vector(mode = "list", length=N)
for (i in 1:N) {
  theta_back_sample<-rmvt(n=100, mu=as.numeric(a_back[[i+1]]), S=round(R_back[[i+1]],8), df=n[[N+1]]) # sometimes not symmetric? computational issue?
  
  theta_back_sample<-data.frame(theta_back_sample)
  colnames(theta_back_sample)<-paste("theta_",seq(1:15), sep="")
  
  theta_back_summary[[i]]<-theta_back_sample %>%
    gather(key="param") %>% 
    group_by(param) %>% 
    summarize(median=quantile(value, probs = 0.5),
              lower=quantile(value, probs = 0.025),
              upper=quantile(value, probs = 0.975)) %>% 
    mutate(index=i) %>% 
    ungroup()
}


theta_back_summary<-bind_rows(theta_back_summary)

ggplot(theta_back_summary %>% filter(param=="theta_1") %>% mutate(flu=flu))+
  geom_line(aes(x=index, y=flu), col="grey")+
  geom_line(aes(x=index, y=median))+
  geom_line(aes(x=index, y=lower), lty=2)+
  geom_line(aes(x=index, y=upper), lty=2)+
  theme_classic()

ggplot(theta_back_summary %>% filter(param=="theta_3") %>% mutate(flu=flu))+
  geom_line(aes(x=index, y=flu), col="grey")+
  geom_line(aes(x=index, y=median))+
  geom_line(aes(x=index, y=lower), lty=2)+
  geom_line(aes(x=index, y=upper), lty=2)+
  theme_classic()

ggplot(theta_back_summary %>% filter(param=="theta_5") %>% mutate(flu=flu))+
  geom_line(aes(x=index, y=flu), col="grey")+
  geom_line(aes(x=index, y=median))+
  geom_line(aes(x=index, y=lower), lty=2)+
  geom_line(aes(x=index, y=upper), lty=2)+
  theme_classic()

ggplot(theta_back_summary %>% filter(param=="theta_7") %>% mutate(flu=flu))+
  geom_line(aes(x=index, y=flu), col="grey")+
  geom_line(aes(x=index, y=median))+
  geom_line(aes(x=index, y=lower), lty=2)+
  geom_line(aes(x=index, y=upper), lty=2)+
  theme_classic()

ggplot(theta_back_summary %>% filter(param=="theta_9") %>% mutate(flu=flu))+
  geom_line(aes(x=index, y=flu), col="grey")+
  geom_line(aes(x=index, y=median))+
  geom_line(aes(x=index, y=lower), lty=2)+
  geom_line(aes(x=index, y=upper), lty=2)+
  theme_classic()

ggplot(theta_back_summary %>% filter(param=="theta_11") %>% mutate(flu=flu))+
  geom_line(aes(x=index, y=flu), col="grey")+
  geom_line(aes(x=index, y=median))+
  geom_line(aes(x=index, y=lower), lty=2)+
  geom_line(aes(x=index, y=upper), lty=2)+
  theme_classic()

ggplot(theta_back_summary %>% filter(param=="theta_13") %>% mutate(flu=flu))+
  geom_line(aes(x=index, y=flu), col="grey")+
  geom_line(aes(x=index, y=median))+
  geom_line(aes(x=index, y=lower), lty=2)+
  geom_line(aes(x=index, y=upper), lty=2)+
  theme_classic()

ggplot(theta_back_summary %>% filter(param=="theta_15") %>% mutate(flu=flu))+
  geom_line(aes(x=index, y=flu), col="grey")+
  geom_line(aes(x=index, y=median))+
  geom_line(aes(x=index, y=lower), lty=2)+
  geom_line(aes(x=index, y=upper), lty=2)+
  theme_classic()

##### residual analysis

plot(x=flu, y=unlist(cbind(f)))
abline(a=0, b=1, col="red")
cor(x=flu, y=unlist(cbind(f)), method = "pearson")
filter_error<-flu-unlist(cbind(f))
filter_MAD<-mean(abs(filter_error))
filter_MAD
filter_MSE<-mean((filter_error)^2)
filter_MSE
plot(filter_error, type = "l")
acf(filter_error)
pacf(filter_error)

plot(x=flu, y=unlist(cbind(f_back))[2:(N+1)])
abline(a=0, b=1, col="red")
cor(x=flu, y=unlist(cbind(f_back))[2:(N+1)], method = "pearson")
smooth_error<-flu-unlist(cbind(f_back))[2:(N+1)]
smooth_MAD<-mean(abs(smooth_error))
smooth_MAD
smooth_MSE<-mean((smooth_error)^2)
smooth_MSE
plot(smooth_error, type = "l")
acf(smooth_error)
pacf(smooth_error)

##### forecast
a_forecast<-G_dlm%*%m[[N+1]]
R_forecast<-(G_dlm%*%C[[N+1]]%*%t(G_dlm))/delta
f_forecast<-t(F_dlm)%*%a_forecast
Q_forecast<-t(F_dlm)%*%R_forecast%*%F_dlm+S[[N+1]]

forecast_sample<-rst(n=1000, mu=f_forecast, sigma=Q_forecast, nu=n[[N+1]])
hist(forecast_sample)

# 12 steps
forecast_sample_summary<-vector(mode="list", length = 12)
# forecast_mean<-rep(NA, 12)
for (j in 1:12) {
  a_forecast<-vector(mode="list", length=12+1)
  a_forecast[[1]]<-m[[N-12+j+1]]
  R_forecast<-vector(mode="list", length=12+1)
  R_forecast[[1]]<-C[[N-12+j+1]]
  f_forecast<-vector(mode="list", length=12+1)
  Q_forecast<-vector(mode="list", length=12+1)
  
  for (i in 2:13) {
    a_forecast[[i]]<-G_dlm%*%a_forecast[[i-1]]
    R_forecast[[i]]<-(G_dlm%*%R_forecast[[i-1]]%*%t(G_dlm))/delta
    f_forecast[[i]]<-t(F_dlm)%*%a_forecast[[i]]
    Q_forecast[[i]]<-t(F_dlm)%*%R_forecast[[i]]%*%F_dlm+S[[N-12+j+1]]
  }
  # forecast_mean[[j]]<-f_forecast[[13]]
  forecast_sample<-rst(n=1000, mu=f_forecast[[13]], sigma=Q_forecast[[13]], nu=n[[N-12+j+1]])
  forecast_sample_summary[[j]]<-data.frame(forecast_sample) %>% 
    select(value=forecast_sample) %>% 
    summarize(median=quantile(value, probs = 0.5),
                                     lower=quantile(value, probs = 0.025),
                                     upper=quantile(value, probs = 0.975)) %>% 
    mutate(index=j+N)
}

forecast_sample_summary<-bind_rows(forecast_sample_summary) %>% 
  bind_rows(data.frame(index=c(1:N), flu=flu))
forecast_sample_summary

forecast_sample_summary
ggplot(forecast_sample_summary)+
  geom_line(aes(x=index, y=median))+
  geom_line(aes(x=index, y=lower), lty=2)+
  geom_line(aes(x=index, y=upper), lty=2)+
  geom_line(aes(x=index, y=flu))+
  theme_classic()

### Pradp pp 134
h<-diag(c(10,2,rep(0,13)))%*%a[[N+1]]
H<-diag(c(10,2,rep(0,13)))*diag(R[[N+1]])
a_forecast_int<-a[[N+1]]+h
R_forecast_int<-R[[N+1]]+H
f_forecast_int<-t(F_dlm)%*%a_forecast_int
Q_forecast_int<-t(F_dlm)%*%R_forecast_int%*%F_dlm+S[[N+1]]
forecast_int_sample<-rst(n=1000, mu=f_forecast_int, sigma=Q_forecast_int, nu=n[[N+1]])
hist(forecast_int_sample)
hist(forecast_sample)

forecast_df<-data.frame(no_intervention=forecast_sample, intervention=forecast_int_sample) %>% 
  gather(key=type)
ggplot(forecast_df)+
  geom_histogram(aes(value, fill=type), alpha=0.5, position = "identity")+
  theme_classic()
