#Chap 2 Problem 21

r1<-0.9
w1<-5
phi1<-2*r1*cos(w1)
phi2<--r1^2

r2<-0.75
w2<-1.35
phi3<-2*r2*cos(w2)
phi4<--r2^2

ts<-arima.sim(n = 500, list(ar = c(phi1, phi2, phi3, phi4)), sd=1)

plot(ts, type="l")

acf(ts)

acf(ts, type = "partial")

# 