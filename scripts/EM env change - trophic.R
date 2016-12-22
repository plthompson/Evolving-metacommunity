#packages####
library(dplyr)
library(ggplot2)
library(igraph)

#functions####
Env_perform<-function(env,z,zmax=NA,sig_p){
  Tmat<-matrix(env,length(env),species) 
  wT<-exp(-((Tmat-rep(z,each=patches))/2*rep(sig_p,each=patches))^2)
  #wT2<-1-((Tmat-z)/(z-z+zmax))^2
  #wT[Tmat>=z]<-wT2[Tmat>=z]
  wT[wT<0]<-0
  wT<-wT-1
  return(wT)
}

#simulation parameters####
species<-40
patches<-50

nplants<-species*0.5
nherb<-species*0.3
npred<-species*0.2

#Environmental change####
burn_in<-1000
change<-7000
change_mag<-7

ChangeV<-c(rep(0,burn_in),seq(0,change_mag-1,length=change),rep(change_mag-1,burn_in))
Tmax<-length(ChangeV)

#Species interactions####
weight=1/80*3

plantV<-1:nplants
herbV<-(nplants+1):(nplants+nherb)
predV<-(species-npred+1):(species)
trophicV<-factor(c(rep("plant",nplants),rep("herbivore",nherb),rep("predator",npred)),levels=c("plant","herbivore","predator"),ordered = T)

b11=-0.05#-0.1
b12=-0.3
b21=0.1
b23=-.1
b32=.08
bdiag1=-.2
bdiag2=-.15

#tritrophic BB Matrix####
B11=b11*matrix(runif(nplants*nplants),nplants,nplants)
B12=b12*matrix(runif(nplants*nherb),nplants,nherb)
B13=matrix(0,nplants,npred)
B21=b21*matrix(runif(nherb*nplants),nherb,nplants)
B22=matrix(0,nherb,nherb)
B23=b23*matrix(runif(nherb*npred),nherb,npred)
B31=matrix(0,npred,nplants)
B32=b32*matrix(runif(npred*nherb),npred,nherb)
B33=matrix(0,npred,npred)
BB=rbind(cbind(B11 ,B12, B13),cbind(B21,B22, B23),cbind(B31, B32, B33))
diag(BB)<-bdiag1
diag(BB[(nplants+nherb+1):species,(nplants+nherb+1):species])<-bdiag2
BB=weight*BB

BB

C3<-c(rep(0.2,nplants),rep(0,species-nplants))

#Environment
Temp_I<-c(1:(1+patches/2),(patches/2):2)

#Environmental optima
z<-runif(n = species,min = min(Temp_I),max=max(Temp_I))
z<-c(seq(min(Temp_I),max(Temp_I),length=nplants),seq(min(Temp_I),max(Temp_I),length=nherb),seq(min(Temp_I),max(Temp_I),length=npred))

sig_p<-0.15

#dispersal####
edges<-rep(1:patches,each=2)
edges<-c(edges[-1],edges[1])
graph<-make_graph(edges, directed = FALSE)
dist_mat<-distances(graph)

dd<-0.5#rnorm(n = species,mean=0.5,sd=0.1)
disp_matrix<-exp(-dd*dist_mat)
diag(disp_matrix)<-0
disp_matrix<-disp_matrix/rowSums(disp_matrix)

ddV<-rnorm(species,mean = dd,sd = 0.5*0.1)

disp_array<-array(NA,dim=c(patches,patches,species))
for(s in 1:species){
  disp_matrix<-exp(-ddV[s]*dist_mat)
  diag(disp_matrix)<-0
  disp_matrix<-disp_matrix/rowSums(disp_matrix)
  disp_array[,,s]<-disp_matrix
}
disp_array_initial<-disp_array
#matplot(disp_array[,25,],type='l')

N<-matrix(c(rep(5,nplants),rep(3,nherb),rep(1,npred)),species,patches)
Nsave<-array(NA,dim=c(species,patches,Tmax))


disp<-0.0001
for(l in 1:Tmax){
  Temp<-Temp_I+ChangeV[l]
  A<-Env_perform(env = Temp,z=z,sig_p = sig_p)
  Nt<-N*exp(C3+BB%*%N+t(A))
  Nt<-Nt+Nt%*%disp_matrix*disp-Nt*disp
  Nt[Nt<10^-3]<-0
  N<-Nt
  Nsave[,,l]<-N
}

trophic_abundance<-data.frame(N=c(N),Patch=rep(1:patches,each=species),T_level=trophicV) %>% 
  group_by(T_level,Patch) %>% 
  summarise(Abundance=sum(N),S=sum(N>0))

ggplot(trophic_abundance,aes(x=Patch,y=Abundance,group=T_level,color=T_level))+
  geom_line()

ggplot(trophic_abundance,aes(x=Patch,y=S,group=T_level,color=T_level))+
  geom_line()

species_abundance<-data.frame(N=c(N),Patch=rep(1:patches,each=species),T_level=trophicV,Species=1:species)

ggplot(species_abundance,aes(x=Patch,y=N,group=Species,color=T_level))+
  geom_line()

filled.contour(Nsave[,,burn_in])
filled.contour(N)

data.frame(N=c(N),Patch=rep(1:patches,each=species),T_level=trophicV) %>% 
  filter(N>0) %>% 
  group_by(T_level) %>% 
  summarise(N=mean(N))

# plot(colSums(Nsave[,,burn_in]>0))
# plot(colSums(N>0))
# plot(rowSums(N>0),col=trophicV,pch=19)
# 
# plot(rowSums(Nsave[,,burn_in]),col=trophicV,pch=19)
# plot(rowSums(N),col=trophicV,pch=19)
# 
# 
# matplot(t(Nsave[,25,]),type='l',col=trophicV)
# matplot(t(Nsave[,2,]),type='l',col=trophicV)
# matplot(t(Nsave[,15,]),type='l',col=trophicV)
