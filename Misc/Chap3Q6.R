# Chapter 3 Problem 6

# (a)
alpha_1<-0.9
alpha_2<--0.95
phi_1<- alpha_1+alpha_2
phi_2<- -alpha_1*alpha_2
poly_coef<-c(1, -phi_1, -phi_2)
1/polyroot(poly_coef)

spectral<-function(omega) {
  f<-1/2/pi/(1+phi_1^2+phi_2^2+2*((-phi_1+phi_1*phi_2)*cos(omega)+(-phi_2)*cos(2*omega)))
  return (f)
}

spectral_df<-data.frame(x=seq(0, pi,length.out = 200))
spectral_df$y<-spectral(spectral_df$x)

ggplot(spectral_df)+
  geom_line(aes(x=x, y=y))+
  theme_classic()

# (b)
alpha_1<-complex(modulus = 0.95, argument = 2*pi/8)
alpha_2<-complex(modulus = 0.95, argument = -2*pi/8)
phi_1<- alpha_1+alpha_2
phi_2<- -alpha_1*alpha_2
poly_coef<-c(1, -phi_1, -phi_2)
Mod(1/polyroot(poly_coef))
Arg(1/polyroot(poly_coef))

spectral<-function(omega) {
  f<-1/2/pi/(1+phi_1^2+phi_2^2+2*((-phi_1+phi_1*phi_2)*cos(omega)+(-phi_2)*cos(2*omega)))
  return (f)
}

spectral_df<-data.frame(x=seq(0, pi,length.out = 200))
spectral_df$y<-Re(spectral(spectral_df$x))

ggplot(spectral_df)+
  geom_line(aes(x=x, y=y))+
  theme_classic()

# (c)
alpha_1<-complex(modulus = 0.5, argument = 2*pi/8)
alpha_2<-complex(modulus = 0.5, argument = -2*pi/8)
phi_1<- alpha_1+alpha_2
phi_2<- -alpha_1*alpha_2
poly_coef<-c(1, -phi_1, -phi_2)
Mod(1/polyroot(poly_coef))
Arg(1/polyroot(poly_coef))

spectral<-function(omega) {
  f<-1/2/pi/(1+phi_1^2+phi_2^2+2*((-phi_1+phi_1*phi_2)*cos(omega)+(-phi_2)*cos(2*omega)))
  return (f)
}

spectral_df<-data.frame(x=seq(0, pi,length.out = 200))
spectral_df$y<-Re(spectral(spectral_df$x))

ggplot(spectral_df)+
  geom_line(aes(x=x, y=y))+
  theme_classic()
