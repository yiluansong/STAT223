# W&H Chapter 8 Problem 5

y<-c(6, 13, 36, 46, 99, 99, 89, 82, 64, 43, 19, 9)

X<-matrix(nrow=12, ncol=12)
for (t in 1:12){
  X[t,]<-c(1, 
           cos(2*pi*1*t/12),sin(2*pi*1*t/12),
           cos(2*pi*2*t/12),sin(2*pi*2*t/12),
           cos(2*pi*3*t/12),sin(2*pi*3*t/12),
           cos(2*pi*4*t/12),sin(2*pi*4*t/12),
           cos(2*pi*5*t/12),sin(2*pi*5*t/12),
           (-1)^t)
}
X

summary(lm(y~X-1))


plot(y, type="l", ylim=c(-100,100))
abline(a=coef(lm(y~X-1))[[1]], b=0)
lines(x=c(1:12), y=coef(lm(y~X-1))[[2]]*cos(2*pi*1*c(1:12)/12)+coef(lm(y~X-1))[[3]]*sin(2*pi*1*c(1:12)/12))
lines(x=c(1:12), y=coef(lm(y~X-1))[[4]]*cos(2*pi*1*c(1:12)/12)+coef(lm(y~X-1))[[5]]*sin(2*pi*1*c(1:12)/12))
lines(x=c(1:12), y=coef(lm(y~X-1))[[6]]*cos(2*pi*1*c(1:12)/12)+coef(lm(y~X-1))[[7]]*sin(2*pi*1*c(1:12)/12))
lines(x=c(1:12), y=coef(lm(y~X-1))[[8]]*cos(2*pi*1*c(1:12)/12)+coef(lm(y~X-1))[[9]]*sin(2*pi*1*c(1:12)/12))
lines(x=c(1:12), y=coef(lm(y~X-1))[[10]]*cos(2*pi*1*c(1:12)/12)+coef(lm(y~X-1))[[11]]*sin(2*pi*1*c(1:12)/12))
lines(x=c(1:12), y=coef(lm(y~X-1))[[12]]*(-1)^c(1:12))

harmonics_df<-data.frame(A=rep(NA, 6), gamma=rep(NA, 6))
for (r in 1:5) {
  harmonics_df$A[r]<-sqrt(coef(lm(y~X-1))[[2*(r-1)+2]]^2+coef(lm(y~X-1))[[2*(r-1)+3]]^2)
  harmonics_df$gamma[r]<-atan(-coef(lm(y~X-1))[[2*(r-1)+3]]/coef(lm(y~X-1))[[2*(r-1)+2]])
}
harmonics_df$A[6]<-sqrt(coef(lm(y~X-1))[[2*(6-1)+2]]^2+0^2)
harmonics_df$gamma[6]<-atan(-0/coef(lm(y~X-1))[[2*(6-1)+2]])
harmonics_df

harmonics_df$percentage<-harmonics_df$A^2/sum(harmonics_df$A^2)*100
harmonics_df

### using eqn on pp 247
a_0<-1/12*sum(y)
a_1<-2/12*sum(y*cos(2*pi/12*1*c(1:12)))
b_1<-2/12*sum(y*sin(2*pi/12*1*c(1:12)))
a_6<-1/12*sum(y*(-1)^c(1:12))
