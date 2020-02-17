soi<-read.table("soi.txt")[,1]
class(soi)
plot(soi, type="l")
mean(soi)

t<-length(soi)
p<-length(soi)/2
m<-floor(p/2)

y<-soi

f<-matrix(NA, nrow= 2*m, ncol=t)
for (i in 1:t) {
  for (k in 1:m) {
  f[(2*k-1):(2*k), i]<-c(cos(2*pi*i*k/p), sin(2*pi*i*k/p))
}
}

f<-t(f)

beta_hat<-coef(lm(soi~f))[-1]
residual<-resid(lm(soi~f))




#####
r1<-0.9
r2<- -0.95


r1<-complex(modulus = 0.95, argument = 2*pi/8)
r2<-complex(modulus = 0.95, argument = -2*pi/8)

r1<-complex(modulus = 0.5, argument = 2*pi/8)
r2<-complex(modulus = 0.5, argument = -2*pi/8)

phi_1<- r1+r2
phi_2<- -r1*r2

v<-1

f<-function(w) {
  v/2/pi/(Mod(1-phi_1*exp(-1i*w)-phi_2*exp(-2*1i*w)))^2
}

plot(f(1:100), type="l")


