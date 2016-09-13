#Evolving metacommunity simulation
#Code by Patrick Thompson
#2016
species<-15
patches<-50

Tmax<-20000

#specify the environment
env_min<-0
env_max_i<-0.8
T<-c(seq(env_min,env_max_i,length=(1+patches/2)),rev(seq(env_min,env_max_i,length=(1+patches/2))[-c(1,(1+patches/2))]))


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

V_all<-c(0,0.0001,0.001,0.01,0.1,1) #additive genetic variation in thermal optimum
dispV<-c(0,0.0001,0.001,0.01,0.1,1)

for(V in V_all){
  for(d in dispV){
    
    #initial species optima means
    z<-matrix(rep(seq(from = env_min,to = env_max_i,length=species),each=patches),ncol=species) #temperature optima
    z_initial<-z
    
    N<-matrix(10,patches,species) #population size
    
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
    library(RColorBrewer)
    library(viridis)
    Trait.df<-gather(data.frame(time=1:Tmax,t(apply(z_store,3,colMeans))),key = Species,value = Trait,X1:X15)
    
    Trait.df$max<-c(t(apply(z_store,3,function(x){apply(x,2,max)})))
    Trait.df$min<-c(t(apply(z_store,3,function(x){apply(x,2,min)})))
    Trait.df$Species<-factor(Trait.df$Species,levels = paste("X",1:species,sep=""),ordered = T)
    
    ggplot(Trait.df,aes(x=time,y=Trait,color=Species,fill=Species))+
      geom_line()+
      geom_ribbon(aes(ymin=min,ymax=max),alpha=0.2,color=NA)+
      theme_bw()+
      removeGrid()
    ggsave(filename = paste("./figures/Trait evolution, V = ",V,", d = ",d,".pdf",sep=""),width = 8,height = 6)
    
    N.df<-data.frame(patch=1:patches,N)
    names(N.df)<-c("patch",1:species)
    N.df<-gather(N.df,key=Species,value=N,-1)
    N.df$Species<-factor(N.df$Species,levels = 1:species,ordered = T)
    
    z.df<-data.frame(patch=1:patches,z)
    names(z.df)<-c("patch",1:species)
    z.df<-gather(z.df,key=Species,value=z,-1)
    
    N.df$Env_trait<-z.df$z
    N.df$z_initial<-c(z_initial)
    N.df$N[N.df$N==0]<-NA
    
    
    matplot(t(apply(N_store,3,rowSums)),type='l', lty=1, ylim=c(0,25))
    matplot(t(apply(N_store,3,colMeans)),type='l', lty=1, ylim=c(0,10))
    
    plot(z_store[1,,1],colSums(N), pch=19, ylab="Regional abundance",xlab="Initial trait")
    
    plot(T,rowSums(N), pch=19)
    
    ggplot(N.df,aes(x=Species,y=Env_trait,size=N,fill=Env_trait))+
      geom_point(pch=21)+
      #scale_fill_gradientn(colors=brewer.pal(9,name = "BuGn"))+
      scale_fill_viridis(option = "D")+
      geom_point(aes(x=Species,y=z_initial),col="red",shape=4,size=2)+
      removeGrid()+
      theme_bw()
    ggsave(filename = paste("./figures/Trait distributions, V = ",V,", d = ",d,".pdf",sep=""),width = 8,height = 6)
    
    
    ggplot(N.df,aes(x=Species,y=patch,size=N,fill=Env_trait))+
      geom_point(pch=21)+
      #scale_fill_gradientn(colors=rev(brewer.pal(9,name = "YlGnBu")))+
      scale_fill_viridis()+
      removeGrid()+
      theme_bw()+
      ylim(25,50)
    ggsave(filename = paste("./figures/Abundance and traits spatial, V = ",V,", d = ",d,".pdf",sep=""),width = 8,height = 6)
    
    z_present<-z
    z_present[N==0]<-NA
    
    mean_z_change<-rep(NA,species)
    for(j in 1:species){
      if(sum(N[,j])>0){
        max_patches<-which(N[,j]==max(N[,j]))
        mean_z_change[j]<-mean(z[max_patches,j]-z_initial[1,j])
      }
    }
  
  hold<-apply(z_present,2,sd,na.rm=T)
  
  response.data1<-data.frame(Dispersal=d,
                             Adapt_potential=V,
                             Value=c(mean(rowSums(N>0)),
                                     sum(colSums(N)>0),
                                     mean(rowSums(N)),
                                     mean(colSums(N>0)[colSums(N>0)>0]),
                                     mean(hold[!is.na(hold)]),
                                     mean(mean_z_change[!is.na(mean_z_change)])
                                     ),
                             Response=c("Local richness","Regional richness","Local biomass","Range size","Optima sd","Optima change"))
  if(d==dispV[1] & V==V_all[1]){
    response.df<-response.data1
  } else {response.df<-rbind(response.df,response.data1)}
}
}

ggplot(response.df,aes(x=Dispersal,y=Value,group=Adapt_potential, color=Adapt_potential))+
  scale_color_viridis(trans="log",breaks=V_all)+
  geom_point()+
  geom_line()+
  facet_wrap(~Response,scales = "free")+
  scale_x_log10()+
  theme_bw()+
  removeGrid()
