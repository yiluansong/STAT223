# W&H Chapter 2
mu[(0+1)]<-25
nu<-rnorm(101, 0, 1)
Y<-rep(NA, 101)
W<-rnorm(101, 0, 0.05)
W<-rnorm(101, 0, 0.5)

for (t in 1:100) {
  mu[t+1]<-mu[t+1-1]+W[t+1]
  Y[t+1]<-mu[t+1]+nu[t+1]
}

plot(Y, type="l")

# W&H Chapter 2
f<-matrix(c(8,13,10), ncol=3)
g<-matrix(c(4,-1,2,3,9,3,1,5,5),ncol=3, byrow=T)
f%*%g

det(matrix(c(1,1,1,8,13,10,81,159,105),ncol=3, byrow=T))

t <-(matrix(c(1,0,1,1,1,0.5,1,2,0.25),ncol=3, byrow=T))
t_1 <-(matrix(c(1,1,1,1,2,2.5,1,3,4.25),ncol=3, byrow=T))
solve(t)%*%t_1

