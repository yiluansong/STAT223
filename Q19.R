roots<-polyroot(c(1, -0.9, 0.9))
Mod(roots)

a<-1/roots

A<-rbind(a, a^2)
A

b<-rbind(0.9, (0.9^2-0.9))
b 

psi_coef<-data.frame(coefficient=solve(A, b)) %>% 
  mutate(modulus=Mod(coefficient)) %>% 
  mutate(argument=Arg(coefficient))
psi_coef

psi<-rep(NA, 50)
for(j in 1:50){
  psi[[j]]<-Re(a[[1]]^j*psi_coef$coefficient[[1]]+a[[2]]^j*psi_coef$coefficient[[2]])
}
psi
plot(psi, type="l")
