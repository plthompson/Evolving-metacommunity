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
library(igraph)
library(betalink)
library(NetIndices)

source("./functions/EM env change functions.R")

#simulation model####
reps<-1

species<-20
patches<-50

Burn_in<-50000
Change_time<-5000
Burn_out<-10000
Tmax<-Burn_in+Change_time+Burn_out

#specify the environment
env_min<-0
env_max_i<-0.8
env_max<-1.2
Temp_initial<-c(seq(env_min,env_max_i,length=(1+patches/2)),rev(seq(env_min,env_max_i,length=(1+patches/2))[-c(1,(1+patches/2))]))
Temp_changeV<-c(rep(0,Burn_in),seq(0,env_max-env_max_i,length=Change_time),rep(env_max-env_max_i,Burn_out))

r_max<-rep(0.2,species) #max intrinsic rate of increase
w<-0.1 #width of selection on thermal performance trait

#interactions
weight<-1/(species/2)

sig_p<-5 #rise in performance as temp increases
zmax<-0.2 #distance from temp opt to zero growth

edges<-rep(1:patches,each=2)
edges<-c(edges[-1],edges[1])
graph<-make_graph(edges, directed = FALSE)
dist_mat<-distances(graph)


V_all<-c(0.00001,0.0001,0.001,0.005,0.01,0.1,1) #additive genetic variation in thermal optimum
dispV<-c(0.00001,0.0001,0.001,0.01,0.1,0.5)

for(r in 1:reps){
  a<-matrix(-.15*runif(species*species),species,species)*weight#competitive matrix
  diag(a)<--0.02
  
  dd<-0.5#rnorm(n = species,mean=0.5,sd=0.1)
  disp_matrix<-exp(-dd*dist_mat)
  diag(disp_matrix)<-0
  disp_matrix<-disp_matrix/rowSums(disp_matrix)
  # disp_array<-array(NA,dim=c(patches,patches,species))
  # for(s in 1:species){
  #   disp_matrix<-exp(-dd[s]*dist_mat)
  #   diag(disp_matrix)<-0
  #   disp_matrix<-disp_matrix/rowSums(disp_matrix)
  #   disp_array[,,s]<-disp_matrix
  # }
  
  for(V in V_all){
    V_species<-rnorm(n = species,mean = V,sd=V*0.1)
    for(d in dispV){
      print(paste("rep = ",r,", V = ",V,", d = ",d,sep=""))
      
      dV<-rnorm(n = species,mean = d,sd=d*0.1)
      
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
        
        zt<-z+(rep(V_species,each=patches)*z_change*z_up_down)#+rnorm(species*patches,mean=0,sd=0.0001)
        
        
        Nt<-Nt1-Nt1*rep(dV,each=patches)+rep(dV,each=patches)*disp_matrix%*%Nt1 #dispersal
        
        zt1<-(disp_matrix%*%(zt*Nt1*rep(dV,each=patches))+zt*Nt1*(1-rep(dV,each=patches)))/Nt #gene flow
        zt1[is.na(zt1)]<-zt[is.na(zt1)]
        zt<-zt1
        
        z_store[,,i]<-z
        
        Nt[Nt<10^-2]<-0
        N_store[,,i]<-N
        N<-Nt
        z<-zt
      }
      
      #matplot(t(N_store[40,,]), type='l', lty=1)
      #matplot(t(z_store[40,,]), type='l', lty=1)
      
      # Trait.df<-gather(data.frame(time=1:Tmax,t(apply(z_store,3,colMeans))),key = Species,value = Trait,X1:X20)
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
      # # ggsave(filename = paste("./figures/Trait distributions, V = ",V,", d = ",d,".pdf",sep=""),width = 8,height = 6)
      # # 
      # # 
      # ggplot(N.df,aes(x=Species,y=patch,size=N,fill=Env_trait))+
      #   geom_point(pch=21)+
      #   #scale_fill_gradientn(colors=rev(brewer.pal(9,name = "YlGnBu")))+
      #   scale_fill_viridis()+
      #   removeGrid()+
      #   theme_bw()+
      #   ylim(25,50)
      # # ggsave(filename = paste("./figures/Abundance and traits spatial, V = ",V,", d = ",d,".pdf",sep=""),width = 8,height = 6)
      # # 
      
      clim_ana<-c("all","no_analogue","analogue")
      for(ca in clim_ana){
        if(ca == "all"){
          patch_select<-1:patches
        }
        if(ca == "no_analogue"){
          patch_select<-Temp > max(Temp_initial)
        }
        if(ca == "analogue"){
          patch_select<-Temp < max(Temp_initial) & Temp > min(Temp_initial)
        }
        
        z_present<-z[patch_select,]
        z_present[N[patch_select,]==0]<-NA
        z_present[is.infinite(z_present)]<-NA
        
        z_prechange<-z_store[patch_select,,Burn_in]
        z_prechange[N_store[patch_select,,Burn_in]==0]<-NA
        z_prechange[is.infinite(z_prechange)]<-NA
        
        mean_z_change<-rep(NA,species)
        for(j in 1:species){
          if(sum(N[patch_select,j])>0){
            mean_z_change[j]<-abs(weighted.mean(z_present[,j],N[patch_select,j])-weighted.mean(z_prechange[,j],N_store[patch_select,j,Burn_in]))
          }
        }
        
        hold<-apply(z_present,2,sd,na.rm=TRUE)
        
        #network dissimilarity####
        Ints<-matrix(1,species,species)
        diag(Ints)<-0
        
        colnames(Ints)<-rownames(Ints)<-paste("s",1:species)
        
        cut_value<-0.5
        
        nets_pre<-apply(N_store[patch_select,,Burn_in],1,function(x){
          Int_strength<-abs(Ints*rep(x,each=species))
          Int_strength[x==0,]<-0
          Int_strength_cut<-quantile(Int_strength[Int_strength>0],cut_value)#mean(Int_strength[Int_strength>0])
          Int_strength[Int_strength<Int_strength_cut]<-0
          Ints2<-1*Int_strength>0
          hold.df<-t(data.frame(Ints2[x>0,x>0]))
          net1<-graph.adjacency(hold.df)
          return(net1) 
        })
        
        nets_post<-apply(N[patch_select,],1,function(x){
          Int_strength<-abs(Ints*rep(x,each=species))
          Int_strength[x==0,]<-0
          Int_strength_cut<-quantile(Int_strength[Int_strength>0],cut_value)#mean(Int_strength[Int_strength>0])
          Int_strength[Int_strength<Int_strength_cut]<-0
          Ints2<-1*Int_strength>0
          hold.df<-t(data.frame(Ints2[x>0,x>0]))
          net1<-graph.adjacency(hold.df)
          return(net1) 
        })
        
        regWeb_pre<-metaweb(nets_pre)
        regWeb_post<-metaweb(nets_post)
        
        #network_betaplot(regWeb_post,regWeb_pre)
        
        NetInds<-data.frame(GenInd2(get.adjacency(regWeb_post,sparse = F)))/data.frame(GenInd2(get.adjacency(regWeb_pre,sparse = F)))                             
        
        
        response.data1<-data.frame(Dispersal=d,
                                   Adapt_potential=V,
                                   Rep=r,
                                   Patches=ca,
                                   Value=c(mean(rowSums(N[patch_select,]>0))/mean(rowSums(N_store[patch_select,,Burn_in]>0)),
                                           sum(colSums(N[patch_select,])>0)/sum(colSums(N_store[patch_select,,Burn_in])>0),
                                           mean(rowSums(N[patch_select,]))/mean(rowSums(N_store[patch_select,,Burn_in])),
                                           mean(colSums(N[patch_select,]>0)[colSums(N[patch_select,]>0)>0])/mean(colSums(N_store[patch_select,,Burn_in]>0)[colSums(N_store[patch_select,,Burn_in]>0)>0]),
                                           mean(apply(z_present,2,sd,na.rm=TRUE),na.rm=TRUE)/mean(apply(z_prechange,2,sd,na.rm=TRUE),na.rm=TRUE),
                                           mean(mean_z_change, na.rm=TRUE),
                                           betalink2(regWeb_pre,regWeb_post,bf = B_jack_diss)$WN,
                                           NetInds$Ltot,
                                           NetInds$LD,
                                           NetInds$C,
                                           NetInds$Cbar
                                   ),
                                   Response=c("Local richness","Regional richness","Local biomass","Range size","Optima sd","Optima change","Network beta","Total links","Link density","Connectance","Compartmentalization"))
        
        if(d==dispV[1] & V==V_all[1] & r == 1 & ca=="all"){
          response.df<-response.data1
        } else {response.df<-rbind(response.df,response.data1)
        }
      }
    }
  }
}

response.df$Response<-factor(response.df$Response,levels=c("Local richness","Regional richness","Local biomass","Range size","Optima sd","Optima change","Network beta","Total links","Link density","Connectance","Compartmentalization"),ordered = TRUE)

response_means<-response.df %>% 
  group_by(Response,Dispersal,Adapt_potential,Patches) %>% 
  summarise(Mean=mean(Value, na.rm=T),Lower=quantile(Value,probs = 0.25,na.rm=T),Upper = quantile(Value,probs = 0.75,na.rm=T))

ggplot(filter(response_means,
              Response=="Local richness" |
                Response=="Regional richness" |
                Response=="Local biomass" |
                Response=="Range size")
       ,aes(x=Dispersal,y=Mean,group=Adapt_potential, color=as.character(Adapt_potential)))+
  #scale_color_viridis(trans="log",breaks=V_all)+
  geom_point()+
  geom_line()+
  facet_grid(Response~Patches,scales = "free_y")+
  scale_x_log10()+
  theme_bw()+
  removeGrid()
ggsave(filename = "./figures/Changing environment/Biodiversity function.pdf",width = 13,height =8 )

ggplot(filter(response_means,
              Response=="Optima sd" |
                Response=="Optima change")
       ,aes(x=Dispersal,y=Mean,group=Adapt_potential, color=as.character(Adapt_potential)))+
  #scale_color_viridis(trans="log",breaks=V_all)+
  geom_point()+
  geom_line()+
  facet_grid(Response~Patches,scales = "free_y")+
  scale_x_log10()+
  theme_bw()+
  removeGrid()
ggsave(filename = "./figures/Changing environment/Adaptation.pdf",width = 13,height =8 )

ggplot(filter(response_means,
              Response=="Network beta" |
                Response=="Total links" |
                Response=="Link density" |
                Response=="Connectance" |
                Response=="Compartmentalization")
       ,aes(x=Dispersal,y=Mean,group=Adapt_potential, color=Adapt_potential))+
  scale_color_viridis(trans="log",breaks=V_all)+
  geom_point()+
  geom_line()+
  facet_grid(Response~Patches,scales = "free_y")+
  scale_x_log10()+
  theme_bw()+
  removeGrid()
ggsave(filename = "./figures/Changing environment/Network change.pdf",width = 13,height =8 )


save(response.df,file = "./workspace/Evolving MC - change.RData")
