#Evolving metacommunity simulation
#Code by Patrick Thompson
#2016

V<-0.005 #1 #additive genetic variation in thermal optimum
d<-0.001

species<-15
patches<-50

Tmax<-10000

#specify the environment
env_min<-0
env_max_i<-0.8
T<-c(seq(env_min,env_max_i,length=(1+patches/2)),rev(seq(env_min,env_max_i,length=(1+patches/2))[-c(1,(1+patches/2))]))

#initial species optima means
#z<-matrix(rep(runif(species,min = env_min,max=env_max_i),each=patches),ncol=species) #temperature optima
z<-matrix(rep(seq(from = env_min,to = env_max_i,length=species),each=patches),ncol=species) #temperature optima


N<-matrix(10,patches,species) #population size
r_max<-rep(0.2,species) #max intrinsic rate of increase
w<-0.1 #width of selection on thermal performance trait

#interactions
weight<-1/(species/2)

a<-matrix(-.15*runif(species*species),species,species)*weight#competitive matrix
diag(a)<--0.02

sig_p<-5 #rise in performance as temp increases
zmax<-0.2 #distance from temp opt to zero growth

T_perform<-function(T,z,zmax,sig_p){
  Tmat<-matrix(T,length(T),species) 
  wT<-exp(-((Tmat-z)/2*sig_p)^2)
  #wT2<-1-((Tmat-z)/(z-z+zmax))^2
  #wT[Tmat>=z]<-wT2[Tmat>=z]
  wT[wT<0]<-0
  return(wT)
}

disp_matrix<-matrix(0,patches,patches)
for(j in 1:patches){
  if(j<patches){
  disp_matrix[j,j+1]<-0.5
  disp_matrix[j+1,j]<-0.5
  } else {
    disp_matrix[j,1]<-0.5
    disp_matrix[1,j]<-0.5
  }
}

z_store<-array(NA,dim = c(patches,species,Tmax))
N_store<-array(NA,dim = c(patches,species,Tmax))
for(i in 1:Tmax){
  
  g<-exp(
    T_perform(T,z,zmax,sig_p)*rep(r_max,each=patches)+ #temperature dependent growth
      N%*%a) #competition
  
  Nt1<-g*N
  
  #change in trait z
  z_change<-z_up<-(exp(T_perform(T,z+0.01,z+zmax+0.01,sig_p)*rep(r_max,each=patches)+N%*%a)-
    exp(T_perform(T,z,z+zmax,sig_p)*rep(r_max,each=patches)+N%*%a))
  
  z_down<-(exp(T_perform(T,z-0.01,z+zmax-0.01,sig_p)*rep(r_max,each=patches)+N%*%a)-
    exp(T_perform(T,z,z+zmax,sig_p)*rep(r_max,each=patches)+N%*%a))
  
  z_up_down<-(z_up>0)*1
  z_up_down[z_down>0]<-((z_down>0)*-1)[z_down>0]
  
  z_change[z_down>z_up]<-z_down[z_down>z_up]
  z_change[z_change<0]<-0
  
  zt<-z+(Nt1*V*z_change*z_up_down)
  
  gene_matrix<-d*disp_matrix
  diag(gene_matrix)<-1-(d*1/(patches-1))*(patches-1)
  
  Nt<-Nt1-Nt1*d+d*disp_matrix%*%Nt1
  
  zt1<-(gene_matrix%*%(zt*Nt1))/Nt
  zt1[is.na(zt1)]<-zt[is.na(zt1)]
  zt<-zt1
  
  z_store[,,i]<-z
  
  Nt[Nt<10^-2]<-0
  N_store[,,i]<-N
  N<-Nt
  z<-zt
}

library(tidyr)
library(ggplot2)
library(ggExtra)
Trait.df<-gather(data.frame(time=1:Tmax,t(apply(z_store,3,colMeans))),key = Species,value = Trait,X1:X15)

Trait.df$max<-c(t(apply(z_store,3,function(x){apply(x,2,max)})))
Trait.df$min<-c(t(apply(z_store,3,function(x){apply(x,2,min)})))
Trait.df$Species<-factor(Trait.df$Species,levels = paste("X",1:species,sep=""),ordered = T)

ggplot(Trait.df,aes(x=time,y=Trait,color=Species,fill=Species))+
  geom_line()+
  geom_ribbon(aes(ymin=min,ymax=max),alpha=0.2,color=NA)+
  theme_bw()+
  removeGrid()
ggsave(filename = "./figures/Trait evolution.pdf",width = 8,height = 6)

N.df<-data.frame(patch=1:patches,N)
names(N.df)<-c("patch",1:species)
N.df<-gather(N.df,key=Species,value=N,-1)
N.df$Species<-factor(N.df$Species,levels = 1:species,ordered = T)

ggplot(N.df,aes(x=Species,y=patch,z=N,color=N, fill=N))+
  geom_tile()
ggsave(filename = "./figures/Abundance spatial.pdf",width = 8,height = 6)

matplot(t(apply(N_store,3,rowSums)),type='l', lty=1, ylim=c(0,25))
plot(T,rowSums(N))

