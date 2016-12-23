#packages####
library(dplyr)
library(ggplot2)
library(ggExtra)
library(igraph)
library(tidyr)
library(viridis)
library(betalink)
library(NetIndices)

source("./functions/EM env change functions.r")

#variables to contrast####
reps<-5
V_all<-c(0.001,0.01,0.1,1,10) #additive genetic variation in thermal optimum
dispV<-c(0.00001,0.0001,0.001,0.01,0.1,0.5)


#simulation parameters####
species<-40
patches<-50

nplants<-species*0.5
nherb<-species*0.3
npred<-species*0.2

#Environmental change####
burn_in<-5000
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


for(r in 1:reps){
  
  b11=-0.1
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
  
  for(disp in dispV){
    dV<-rnorm(species,mean = disp,sd=disp*0.25)
    for(V in V_all){
      print(paste("rep = ",r,", V = ",V,", d = ",disp,sep=""))
      vV<-c(rnorm(nplants,mean = V,sd=V*0.25),rnorm(nherb,mean = V*0.5,sd=V*0.5*0.25),rnorm(npred,mean = V*0.25,sd=V*0.25*0.25))
      
      N<-matrix(c(rep(5,nplants),rep(3,nherb),rep(1,npred)),species,patches)
      Zsave<-Nsave<-array(NA,dim=c(species,patches,Tmax))
      zmat<-matrix(z,species,patches)
      
      for(l in 1:Tmax){
        Temp<-Temp_I+ChangeV[l]
        A<-Env_perform2(env = Temp,z=zmat,sig_p = sig_p)
        
        g<-exp(C3+BB%*%N+A)
        
        Nt1<-g*N
        
        #change in trait z
        z_change<-z_up<-(exp(C3+Env_perform2(env = Temp,z = zmat+0.01,sig_p = sig_p)+BB%*%N)-g)
        z_down<-(exp(C3+Env_perform2(env = Temp,z = zmat-0.01,sig_p = sig_p)+BB%*%N)-g)
        
        z_up_down<-(z_up>0)*1
        z_up_down[z_down>0]<-((z_down>0)*-1)[z_down>0]
        
        z_change[z_down>z_up]<-z_down[z_down>z_up]
        z_change[z_change<0]<-0
        
        zt<-zmat+(Nt1*vV*z_change*z_up_down)
        
        Nt<-Nt1-Nt1*dV+dV*Nt1%*%disp_matrix
        
        zt1<-(dV*(zt*Nt1)%*%disp_matrix+(Nt1*zt*(1-dV)))/Nt
        zt1[is.na(zt1)]<-zt[is.na(zt1)]
        zt<-zt1
        
        Nt[Nt<10^-3]<-0
        N<-Nt
        Nsave[,,l]<-N
        zmat<-zt
        Zsave[,,l]<-zmat
      }
      
      finalCom.df<-data.frame(species=1:species,patch=rep(1:patches,each=species),N=c(N),z=c(zmat))
      
      ggplot(filter(finalCom.df,N>0),aes(x=species,y=patch,color=z,size=N))+
        geom_point()+
        scale_color_viridis()
      
      clim_ana<-c("all","no_analogue","analogue")
      for(ca in clim_ana){
        if(ca == "all"){
          patch_select<-1:patches
        }
        if(ca == "no_analogue"){
          patch_select<-Temp > max(Temp_I)
        }
        if(ca == "analogue"){
          patch_select<-Temp < max(Temp_I) & Temp > min(Temp_I)
        }
        
        z_present<-zmat[,patch_select]
        z_present[N[,patch_select]==0]<-NA
        z_present[is.infinite(z_present)]<-NA
        
        z_prechange<-Zsave[,patch_select,burn_in]
        z_prechange[Nsave[,patch_select,burn_in]==0]<-NA
        z_prechange[is.infinite(z_prechange)]<-NA
        
        mean_z_change<-rep(NA,species)
        for(j in 1:species){
          if(sum(N[j,patch_select])>0){
            mean_z_change[j]<-abs(weighted.mean(z_present[j,],N[j,patch_select])-weighted.mean(z_prechange[j,],Nsave[j,patch_select,burn_in]))
          }
        }
        
        hold<-apply(z_present,1,sd,na.rm=TRUE)
        
        #network dissimilarity####
        Ints<-matrix(1,species,species)
        diag(Ints)<-0
        
        colnames(Ints)<-rownames(Ints)<-paste("s",1:species)
        
        cut_value<-0.5
        
        nets_pre<-apply(Nsave[,patch_select,burn_in],2,function(x){
          Int_strength<-abs(Ints*rep(x,each=species))
          Int_strength[x==0,]<-0
          Int_strength_cut<-quantile(Int_strength[Int_strength>0],cut_value)#mean(Int_strength[Int_strength>0])
          Int_strength[Int_strength<Int_strength_cut]<-0
          Ints2<-1*Int_strength>0
          hold.df<-t(data.frame(Ints2[x>0,x>0]))
          net1<-graph.adjacency(hold.df)
          return(net1) 
        })
        
        nets_post<-apply(N[,patch_select],2,function(x){
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
        
        
        response.data1<-data.frame(Dispersal=disp,
                                   Adapt_potential=V,
                                   Rep=r,
                                   Patches=ca,
                                   Value=c(mean(colSums(N[,patch_select]>0))/mean(colSums(Nsave[,patch_select,burn_in]>0)),
                                           sum(rowSums(N[,patch_select])>0)/sum(rowSums(Nsave[,patch_select,burn_in])>0),
                                           mean(colSums(N[,patch_select]))/mean(colSums(Nsave[,patch_select,burn_in])),
                                           mean(rowSums(N[,patch_select]>0)[rowSums(N[,patch_select]>0)>0])/mean(rowSums(Nsave[,patch_select,burn_in]>0)[rowSums(Nsave[,patch_select,burn_in]>0)>0]),
                                           mean(apply(z_present,1,sd,na.rm=TRUE),na.rm=TRUE)/mean(apply(z_prechange,1,sd,na.rm=TRUE),na.rm=TRUE),
                                           mean(mean_z_change, na.rm=TRUE),
                                           betalink2(regWeb_pre,regWeb_post,bf = B_jack_diss)$WN,
                                           NetInds$Ltot,
                                           NetInds$LD,
                                           NetInds$C,
                                           NetInds$Cbar
                                   ),
                                   Response=c("Local richness","Regional richness","Local biomass","Range size","Optima sd","Optima change","Network beta","Total links","Link density","Connectance","Compartmentalization"))
        
        if(disp==dispV[1] & V==V_all[1] & r == 1 & ca=="all"){
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
       ,aes(x=Dispersal,y=Mean,group=Adapt_potential, color=as.character(Adapt_potential),fill=as.character(Adapt_potential)))+
  #scale_color_viridis(trans="log",breaks=V_all)+
  geom_ribbon(aes(ymin=Lower,ymax=Upper),alpha=0.2,color=NA)+
  geom_point()+
  geom_line()+
  facet_grid(Response~Patches,scales = "free_y")+
  scale_x_log10()+
  theme_bw()+
  removeGrid()
#ggsave(filename = "./figures/Changing environment/Biodiversity function.pdf",width = 13,height =8 )

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
#ggsave(filename = "./figures/Changing environment/Adaptation.pdf",width = 13,height =8 )

ggplot(filter(response_means,
              Response=="Network beta" |
                Response=="Total links" |
                Response=="Link density" |
                Response=="Connectance" |
                Response=="Compartmentalization")
       ,aes(x=Dispersal,y=Mean,group=Adapt_potential, color=as.character(Adapt_potential)))+
  #scale_color_viridis(trans="log",breaks=V_all)+
  geom_point()+
  geom_line()+
  facet_grid(Response~Patches,scales = "free_y")+
  scale_x_log10()+
  theme_bw()+
  removeGrid()
#ggsave(filename = "./figures/Changing environment/Network change.pdf",width = 13,height =8 )


#save(response.df,file = "./workspace/Evolving MC - change.RData")



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
