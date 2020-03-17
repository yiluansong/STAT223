# (a)
ts<-arima.sim(n = 400, list(ar = c(0.9), ma=c(0.6)), sd=1)
plot(ts, type="l")

# (b)
y<-ts

epsilon<-rep(NA, length(y))
z<-matrix(NA, ncol=2, nrow = length(y))

beta<-c(1,0)

count<-0

while (count<=1000) {
  num<-0
  denom<-0
  
  for (t in 1:length(y)) {
    if(t==1) {
      epsilon[t]<-y[t]
    } else {
      epsilon[t]<-y[t]-beta[1]*y[t-1]-beta[2]*epsilon[t-1]
    }
    
    if(t==1) {
      z[t,1]<-0
      z[t,2]<-0
    } else {
      z[t,1]<-y[t-1]-beta[2]*z[t-1,1]
      z[t,2]<-epsilon[t-1]-beta[2]*z[t-1,2]
    }
    
    if(t>=2) {num<-num+z[t,]*epsilon[t]} else {num<-num}
    if(t>=2) {denom<-denom+t(matrix(z[t,]))%*%matrix(z[t,])} else {denom<-denom}
  }
  beta<-beta+c(num)/c(denom)
  count<-count+1
}

beta

# (c)

