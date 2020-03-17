#Chap 2 Problem 22 (This is a draft. Final version in R markdown)

library(tidyverse)

iso<-read.table("./STAT 223/oxygen_isotope.txt", skip=6, header=F)
str(iso)
iso<-iso$V2
n<-length(iso)

plot(iso, type="l")
lines(fitted(loess(iso~c(1:n))), col="red")


loess(iso~c(1:n))
iso_detrend<-residuals(loess(iso~c(1:n)))
plot(iso_detrend, type="l")

# a

# placeholder
phi_p<-vector(mode="list", length = 15)
resid_p<-vector(mode="list", length = 15)
s_square_p<-vector(mode="list", length = 15)
AIC_p<-rep(NA, 15)
BIC_p<-rep(NA, 15)

for(p in 1:15) {
  y<-iso_detrend[(p+1):(n)]
  y_prev<-matrix(NA, ncol = n-p, nrow = p)
  for(i in 1:p) {
    y_prev[i,]<-iso_detrend[(p-i+1):(n-i)]
  }
  
  phi_p[[p]]<-solve(y_prev %*% t(y_prev))%*%y_prev%*%y
  resid_p[[p]]<-y-t(y_prev)%*% phi_p[[p]]
  s_square_p[[p]]<-(t(resid_p[[p]]) %*%(resid_p[[p]]))/(n-p)
  AIC_p[[p]]<-2*(p)+(n-p)*log(s_square_p[[p]])
  BIC_p[[p]]<-log(n-p)*p+(n-p)*log(s_square_p[[p]])
}
# lm(y~y_prev[1,])
# AIC(arima(iso_detrend, c(1,0,0), include.mean = F, method = "ML"))
AICBIC<-data.frame(p=1:15,AIC=AIC_p, BIC=BIC_p)
head(AICBIC)

ggplot(AICBIC)+
  geom_point(aes(x=p, y=AIC), col="red", pch="a", cex=5)+
  geom_point(aes(x=p, y=BIC), col="blue", pch="b", cex=5)+
  ylab("score")+
  theme_classic()

# b
phi_p[[3]]

# resid_p[[3]]<-residuals(arima(iso_detrend, c(3,0,0), method="ML"))
plot(resid_p[[3]], type="l")

acf(resid_p[[3]], main="ACF")
pacf(resid_p[[3]],main="PACF")
qqnorm(resid_p[[3]])
qqline(resid_p[[3]], col="red")

# c
poly_coef<-c(1, -phi_p[[3]])
invroot<-data.frame(invroot=1/polyroot(poly_coef)) %>% 
  mutate(modulus=Mod(invroot)) %>% 
  mutate(argument=abs(Arg(invroot))) %>% 
  mutate(wavelength=2*pi/argument) %>% 
  arrange(desc(wavelength)) 
invroot

# d
phi<-c(phi_p[[3]])
v<-c(s_square_p[[3]])
library(ltsa)
gamma_0<-tacvfARMA(phi = phi,  maxLag = 0, sigma2 = v)
gamma_0

plot(tacvfARMA(phi = phi,  maxLag = 16, sigma2 = v))

autocov<-rep(NA, 16)
autocov[1]<-cov(y, y)
for (i in 1:15) {
  autocov[i+1]<-cov(y, y_prev[i,])
}
plot(autocov)

GAMMA_t<-matrix(NA, ncol=n, nrow=n)
for (j in 1:n) {
  GAMMA_t[j,c(j:n)]<-tacvfARMA(phi = phi,  maxLag = 1+n-2, sigma2 =v)[1:(n-j+1)]
  if(j>1) {
    GAMMA_t[j, c(1:(j-1))]<-rev(tacvfARMA(phi = phi,  maxLag = 1+n-2, sigma2 = v)[2:j])
  }
}
GAMMA_t

gamma_t<-matrix(NA, ncol=100, nrow=n)
for (h in 1:100) {
  gamma_t[,h]<-matrix(tacvfARMA(phi = phi,  maxLag = n+h-1, sigma2 = v)[(h+1):(n+h)])
}
str(gamma_t)

pacf_t<-matrix(NA, ncol=100, nrow=n)
for (h in 1:100) {
  pacf_t[,h]<-solve(a=GAMMA_t, b=gamma_t[,h])
}
str(pacf_t)

plot(pacf_t[,1], type="l")
round(pacf_t[1:10,1],5)
round(phi,5)

y_pred<-rep(NA, 100)
for (h in 1:100) {
    y_pred[[h]]<-t(matrix(rev(iso_detrend)))%*%matrix(pacf_t[,h])
}
plot(y_pred, type="l", ylab="Predicted")

gamma_0<-tacvfARMA(phi = phi,  maxLag = 0, sigma2 = v)
MSE_pred<-rep(NA, 100)
for (h in 1:100) {
  MSE_pred[[h]]<-gamma_0-t(gamma_t[,h])%*%solve(GAMMA_t)%*%gamma_t[,h]
}
plot(MSE_pred, type="l", ylab="MSE")

# y_pred
# plot(y_pred, type="l")
# lines(y_pred-1.96*MSE_pred, type="l", lty=2)
# lines(y_pred+1.96*MSE_pred, type="l", lty=2)

pred_df<-data.frame(time=n+c(1:100),predicted=y_pred, lower=y_pred-1.96*MSE_pred, upper=y_pred+1.96*MSE_pred) 

obs_df<-data.frame(time=1:n, observed=iso_detrend)

iso_trend<-coef(lm(iso~c(1:n)))[1]+coef(lm(iso~c(1:n)))[2]*c(1:(n+100))

combine_df<-full_join(pred_df, obs_df, by="time") %>% 
  arrange(time) %>% 
  mutate(predicted_trend=predicted+iso_trend,
         lower_trend=lower+iso_trend,
         upper_trend=upper+iso_trend,
         observed_trend=observed+iso_trend)

ggplot(combine_df)+
  geom_line(aes(x=time, y=observed))+
  geom_line(aes(x=time, y=predicted), col="blue")+
  geom_line(aes(x=time, y=lower), col="blue", lty=2)+
  geom_line(aes(x=time, y=upper), col="blue", lty=2)+
  ylab("")+
  theme_classic()

ggplot(combine_df)+
  geom_line(aes(x=time, y=observed_trend))+
  geom_line(aes(x=time, y=predicted_trend), col="blue")+
  geom_line(aes(x=time, y=lower_trend), col="blue", lty=2)+
  geom_line(aes(x=time, y=upper_trend), col="blue", lty=2)+
  ylab("")+
  theme_classic()

##### Method 2
# placeholders
y_forecast<-matrix(NA, nrow=n+100+1, ncol=n+100+1)
y_new<-rep(NA, n+100+1)
MSE_forecast<-matrix(NA, nrow=n+100+1, ncol=n+100+1)
b<-matrix(NA, nrow=n+100+1, ncol=n+100+1)

# initialize
y_forecast[1, 2]<-0
MSE_forecast[1, 2]<-tacvfARMA(phi = phi,  maxLag = 0, sigma2 = v)
y_new[2:(n+1)]<-iso_detrend


# calculate y and MSE for h=1
h<-1
t<-1
while (t<=n) {
  # update b
  for (j in 0:(t+h-1)) {
    gamma<-tacvfARMA(phi = phi,  maxLag = t+h-1-j, sigma2 = v)[(t+h-1-j+1)]
    b[(t+h-1+1),(t+h-1-j+1) ]<-gamma
    if (j>0) {
      for (l in 0:(j-1)) {
        b[(t+h-1+1),(t+h-1-j+1) ]<-b[(t+h-1+1),(t+h-1-j+1) ]-b[(j+1), (j-l+1)]*b[(t+h-1+1), (t+h-1-l+1)]*MSE_forecast[(l+1), (l+1+1)]
      }
    }
    b[(t+h-1+1),(t+h-1-j+1) ]<-b[(t+h-1+1),(t+h-1-j+1) ]/MSE_forecast[(j+1), (j+1+1)] 
  }
  
  # update y_forecast
  y_forecast[(t+1),(t+h+1)]<-0
  for (j in h:(t+h-1)){
    y_forecast[(t+1),(t+h+1)]<-y_forecast[(t+1),(t+h+1)]+b[(t+h-1+1), (j+1)]*(y_new[t+h-j+1]-y_forecast[(t+h-j-1+1),(t+h-j+1)])
  }
  
  # update MSE_forecast
  MSE_forecast[(t+1),(t+h+1)]<-tacvfARMA(phi = phi,  maxLag = 0, sigma2 = v)
  for (j in h:(t+h-1)) {
    MSE_forecast[(t+1),(t+h+1)]<-MSE_forecast[(t+1),(t+h+1)]-(b[(t+h-1+1),(j+1)])^2*MSE_forecast[(t+h-j-1+1), (t+h-j+1)]
  }
  t<-t+1
}

# set up parallel computing
library(doParallel)
library(foreach)

cores<-detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)


# calculate y and MSE for all h
inno_pred<-
  foreach (h=1:100, .combine=rbind) %dopar% {
    library(ltsa)
    t<-1
    while (t<=n) {
      # update b
      for (j in 0:(t+h-1)) {
        gamma<-tacvfARMA(phi = phi,  maxLag = t+h-1-j, sigma2 = v)[(t+h-1-j+1)]
        b[(t+h-1+1),(t+h-1-j+1) ]<-gamma
        if (j>0) {
          for (l in 0:(j-1)) {
            b[(t+h-1+1),(t+h-1-j+1) ]<-b[(t+h-1+1),(t+h-1-j+1) ]-b[(j+1), (j-l+1)]*b[(t+h-1+1), (t+h-1-l+1)]*MSE_forecast[(l+1), (l+1+1)]
          }
        }
        b[(t+h-1+1),(t+h-1-j+1) ]<-b[(t+h-1+1),(t+h-1-j+1) ]/MSE_forecast[(j+1), (j+1+1)] 
      }
      
      # update y_forecast
      y_forecast[(t+1),(t+h+1)]<-0
      for (j in h:(t+h-1)){
        y_forecast[(t+1),(t+h+1)]<-y_forecast[(t+1),(t+h+1)]+b[(t+h-1+1), (j+1)]*(y_new[t+h-j+1]-y_forecast[(t+h-j-1+1),(t+h-j+1)])
      }
      
      # update MSE_forecast
      MSE_forecast[(t+1),(t+h+1)]<-tacvfARMA(phi = phi,  maxLag = 0, sigma2 = v)
      for (j in h:(t+h-1)) {
        MSE_forecast[(t+1),(t+h+1)]<-MSE_forecast[(t+1),(t+h+1)]-(b[(t+h-1+1),(j+1)])^2*MSE_forecast[(t+h-j-1+1), (t+h-j+1)]
      }
      t<-t+1
    }
    c(h, y_forecast[n+1, n+h+1],MSE_forecast[n+1, n+h+1])
  }
# stop parallel computing
stopCluster(cl)

# save output
colnames(inno_pred)<-c("h", "y_pred", "MSE_pred")
write_csv(data.frame(inno_pred), "./STAT 223/inno_pred.csv")

inno_pred<-read_csv("./STAT 223/inno_pred.csv")

compare_df<-data.frame(y_pred_Method_1=y_pred, MSE_pred_Method_1=MSE_pred, y_pred_Method_2=inno_pred$y_pred, MSE_pred_Method_2=inno_pred$MSE_pred)
compare_df[1:10,]
