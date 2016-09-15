#Evolving metacommunity simulation - changing environment
#Code by Patrick Thompson
#2016

#packages####
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggExtra)
library(RColorBrewer)
library(viridis)

#simulation model####
reps<-10

species<-20
patches<-50

Burn_in<-50000
Change_time<-20000
Tmax<-Burn_in+Change_time

#specify the environment
env_min<-0
env_max_i<-0.8
env_max<-1
Temp_initial<-c(seq(env_min,env_max_i,length=(1+patches/2)),rev(seq(env_min,env_max_i,length=(1+patches/2))[-c(1,(1+patches/2))]))
Temp_changeV<-c(rep(0,Burn_in),seq(0,env_max-env_max_i,length=Change_time))

r_max<-rep(0.2,species) #max intrinsic rate of increase
w<-0.1 #width of selection on thermal performance trait

#interactions
weight<-1/(species/2)

sig_p<-5 #rise in performance as temp increases
zmax<-0.2 #distance from temp opt to zero growth

T_perform<-function(Temp,z,zmax,sig_p){
  Tmat<-matrix(Temp,length(Temp),species) 
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

V_all<-c(0.00001,0.0001,0.001,0.01,0.1,1,10) #additive genetic variation in thermal optimum
dispV<-c(0.00001,0.0001,0.001,0.01,0.1,1)

for(r in 1:reps){
  a<-matrix(-.15*runif(species*species),species,species)*weight#competitive matrix
  diag(a)<--0.02
  
  for(V in V_all){
    for(d in dispV){
      print(paste("rep =",r,", V = ",V,", d = ",d,sep=""))
      
      #initial species optima means
      z<-matrix(rep(seq(from = env_min,to = env_max_i,length=species),each=patches),ncol=species) #temperature optima
      z_initial<-z
      
      N<-matrix(10,patches,species) #population size
      
      z_store<-array(NA,dim = c(patches,species,Tmax))
      N_store<-array(NA,dim = c(patches,species,Tmax))
      for(i in 1:Tmax){
        
        Temp<-Temp_initial+Temp_changeV[i]
        
        g<-exp(
          T_perform(Temp,z,zmax,sig_p)*rep(r_max,each=patches)+ #temperature dependent growth
            N%*%a) #competition
        
        Nt1<-g*N
        
        #change in trait z
        z_change<-z_up<-(exp(T_perform(Temp,z+0.01,z+zmax+0.01,sig_p)*rep(r_max,each=patches)+N%*%a)-
                           exp(T_perform(Temp,z,z+zmax,sig_p)*rep(r_max,each=patches)+N%*%a))
        
        z_down<-(exp(T_perform(Temp,z-0.01,z+zmax-0.01,sig_p)*rep(r_max,each=patches)+N%*%a)-
                   exp(T_perform(Temp,z,z+zmax,sig_p)*rep(r_max,each=patches)+N%*%a))
        
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
      
      
      # Trait.df<-gather(data.frame(time=1:Tmax,t(apply(z_store,3,colMeans))),key = Species,value = Trait,X1:X15)
      # 
      # Trait.df$max<-c(t(apply(z_store,3,function(x){apply(x,2,max)})))
      # Trait.df$min<-c(t(apply(z_store,3,function(x){apply(x,2,min)})))
      # Trait.df$Species<-factor(Trait.df$Species,levels = paste("X",1:species,sep=""),ordered = T)
      # 
      # ggplot(Trait.df,aes(x=time,y=Trait,color=Species,fill=Species))+
      #   geom_line()+
      #   geom_ribbon(aes(ymin=min,ymax=max),alpha=0.2,color=NA)+
      #   theme_bw()+
      #   removeGrid()
      # ggsave(filename = paste("./figures/Trait evolution, V = ",V,", d = ",d,".pdf",sep=""),width = 8,height = 6)
      
      # N.df<-data.frame(patch=1:patches,N)
      # names(N.df)<-c("patch",1:species)
      # N.df<-gather(N.df,key=Species,value=N,-1)
      # N.df$Species<-factor(N.df$Species,levels = 1:species,ordered = T)
      # 
      # z.df<-data.frame(patch=1:patches,z)
      # names(z.df)<-c("patch",1:species)
      # z.df<-gather(z.df,key=Species,value=z,-1)
      # 
      # N.df$Env_trait<-z.df$z
      # N.df$z_initial<-c(z_initial)
      # N.df$N[N.df$N==0]<-NA
      
      
      #matplot(t(apply(N_store,3,rowSums)),type='l', lty=1, ylim=c(0,25))
      #matplot(t(apply(N_store,3,colMeans)),type='l', lty=1, ylim=c(0,10))
      
      #plot(z_store[1,,1],colSums(N), pch=19, ylab="Regional abundance",xlab="Initial trait")
      
      #plot(Temp,rowSums(N), pch=19)
      
      # ggplot(N.df,aes(x=Species,y=Env_trait,size=N,fill=Env_trait))+
      #   geom_point(pch=21)+
      #   #scale_fill_gradientn(colors=brewer.pal(9,name = "BuGn"))+
      #   scale_fill_viridis(option = "D")+
      #   geom_point(aes(x=Species,y=z_initial),col="red",shape=4,size=2)+
      #   removeGrid()+
      #   theme_bw()
      # ggsave(filename = paste("./figures/Trait distributions, V = ",V,", d = ",d,".pdf",sep=""),width = 8,height = 6)
      # 
      # 
      # ggplot(N.df,aes(x=Species,y=patch,size=N,fill=Env_trait))+
      #   geom_point(pch=21)+
      #   #scale_fill_gradientn(colors=rev(brewer.pal(9,name = "YlGnBu")))+
      #   scale_fill_viridis()+
      #   removeGrid()+
      #   theme_bw()+
      #   ylim(25,50)
      # ggsave(filename = paste("./figures/Abundance and traits spatial, V = ",V,", d = ",d,".pdf",sep=""),width = 8,height = 6)
      # 
      
      z_present<-z
      z_present[N==0]<-NA
      z_present[is.infinite(z_present)]<-NA
      
      z_prechange<-z_store[,,Burn_in]
      z_prechange[N_store[,,Burn_in]==0]<-NA
      z_prechange[is.infinite(z_prechange)]<-NA
      
      mean_z_change<-rep(NA,species)
      for(j in 1:species){
        if(sum(N[,j])>0){
          mean_z_change[j]<-abs(weighted.mean(z_present[,j],N[,j])-weighted.mean(z_prechange[,j],N_store[,j,Burn_in]))
        }
      }

      
      
      hold<-apply(z_present,2,sd,na.rm=TRUE)
      
      response.data1<-data.frame(Dispersal=d,
                                 Adapt_potential=V,
                                 Rep=r,
                                 Value=c(mean(rowSums(N>0))/mean(rowSums(N_store[,,Burn_in]>0)),
                                         sum(colSums(N)>0)/sum(colSums(N_store[,,Burn_in])>0),
                                         mean(rowSums(N))/mean(rowSums(N_store[,,Burn_in])),
                                         mean(colSums(N>0)[colSums(N>0)>0])/mean(colSums(N_store[,,Burn_in]>0)[colSums(N_store[,,Burn_in]>0)>0]),
                                         mean(apply(z_present,2,sd,na.rm=TRUE),na.rm=TRUE)/mean(apply(z_prechange,2,sd,na.rm=TRUE),na.rm=TRUE),
                                         mean(mean_z_change, na.rm=TRUE)
                                 ),
                                 Response=c("Local richness","Regional richness","Local biomass","Range size","Optima sd","Optima change"))
      if(d==dispV[1] & V==V_all[1] & r == 1){
        response.df<-response.data1
      } else {response.df<-rbind(response.df,response.data1)}
    }
  }
}

response.df$Response<-factor(response.df$Response,levels=c("Local richness","Regional richness","Local biomass","Range size","Optima sd","Optima change"),ordered = TRUE)

response_means<-response.df %>% 
  group_by(Response,Dispersal,Adapt_potential) %>% 
  summarise(Mean=mean(Value, na.rm=T),Lower=quantile(Value,probs = 0.25),Upper = quantile(Value,probs = 0.75))

ggplot(response_means,aes(x=Dispersal,y=Mean,group=Adapt_potential, color=Adapt_potential))+
  scale_color_viridis(trans="log",breaks=V_all)+
  geom_point()+
  geom_line()+
  facet_wrap(~Response,scales = "free")+
  scale_x_log10()+
  theme_bw()+
  removeGrid()

ggplot(filter(response_means, Response!="Optima change" & Response!="Optima sd"),aes(x=Dispersal,y=Adapt_potential,fill=Mean))+
  scale_fill_gradient2(low = brewer.pal(5,name = "RdBu")[5],mid = brewer.pal(5,name = "RdBu")[3],high = brewer.pal(5,name = "RdBu")[1],midpoint = 1)+
  geom_tile()+
  facet_wrap(~Response,scales = "free")+
  scale_y_log10()+
  scale_x_log10()+
  theme_bw()+
  removeGrid()

ggplot(filter(response_means, Response=="Optima change"),aes(x=Dispersal,y=Adapt_potential,fill=Mean))+
  scale_fill_viridis()+
  geom_tile()+
  facet_wrap(~Response,scales = "free")+
  scale_y_log10()+
  scale_x_log10()+
  theme_bw()+
  removeGrid()

ggplot(filter(response_means, Response=="Optima sd"),aes(x=Dispersal,y=Adapt_potential,fill=Mean))+
  scale_fill_gradient2(low = brewer.pal(5,name = "RdBu")[5],mid = brewer.pal(5,name = "RdBu")[3],high = brewer.pal(5,name = "RdBu")[1],midpoint = 1)+
  geom_tile()+
  facet_wrap(~Response,scales = "free")+
  scale_y_log10()+
  scale_x_log10()+
  theme_bw()+
  removeGrid()


save(response.df,file = "./workspace/Evolving MC - change.RData")
