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


#functions####
B_jack_diss<-function(pm){
  with(pm, {
    (b+c)/(a+b+c)
  })}

B_jack_diss_gains<-function(pm){
  with(pm, {
    (c)/(a+b+c)
  })}

B_jack_diss_loss<-function(pm){
  with(pm, {
    (b)/(a+b+c)
  })}

betalink2<-function (n1, n2, bf = B01) 
{
  v1 <- igraph::V(n1)$name
  v2 <- igraph::V(n2)$name
  vs <- v1[v1 %in% v2]
  beta_S <- bf(betapart(v1, v2))
  e1 <- plyr::aaply(igraph::get.edgelist(n1), 1, function(x) stringr::str_c(x, 
                                                                            collapse = "--", paste = "_"))
  e2 <- plyr::aaply(igraph::get.edgelist(n2), 1, function(x) stringr::str_c(x, 
                                                                            collapse = "--", paste = "_"))
  beta_WN <- bf(betapart(e1, e2))
  if (length(vs) >= 2) {
    sn1 <- igraph::induced.subgraph(n1, which(igraph::V(n1)$name %in% 
                                                vs))
    sn2 <- igraph::induced.subgraph(n2, which(igraph::V(n2)$name %in% 
                                                vs))
    se1 <- plyr::aaply(igraph::get.edgelist(sn1), 1, function(x) stringr::str_c(x, 
                                                                                collapse = "--", paste = "_"))
    se2 <- plyr::aaply(igraph::get.edgelist(sn2), 1, function(x) stringr::str_c(x, 
                                                                                collapse = "--", paste = "_"))
    beta_OS <- bf(betapart(se1, se2))
    beta_ST <- beta_WN - beta_OS
  }
  else {
    beta_OS <- NaN
    beta_ST <- NaN
  }
  return(list(S = beta_S, OS = beta_OS, WN = beta_WN, ST = beta_ST))
}

GenInd2<-function (Flow = NULL, Tij = t(Flow), Import = NULL, Export = NULL, 
                   tol = 0) 
{
  if(length(Flow)==0){
    list(N = 0, T.. = 0, TST = 0, Lint = 0, 
         Ltot = 0, LD = 0, C = 0, Tijbar = 0, 
         TSTbar = 0, Cbar = 0)
  } else{
    N <- InternalNetwork(Tij, Import, Export)
    RateComp <- N$FlowToC - N$FlowFromC
    ncTij <- ncol(Tij)
    nrTij <- nrow(Tij)
    ncomp <- ncol(N$Tint)
    compNames <- rownames(N$Tint)
    intlinks <- length(which(N$Tint > tol))
    links <- length(which(Tij > tol))
    LD <- links/ncomp
    ExportSum <- sum(N$FlowTo[N$export])
    ImportSum <- sum(N$FlowFrom[N$import])
    Throughflow <- sum(N$Tint) + ImportSum - sum(RateComp[RateComp < 
                                                            0])
    Throughput <- sum(Tij)
    Avthrflow <- Throughflow/ncomp
    Connectance <- intlinks/ncomp/(ncomp - 1)
    Avlinkweight <- Throughput/links
    linkmat <- N$Tint
    linkmat[linkmat > 0] <- 1
    Cij <- matrix(nrow = ncomp, ncol = ncomp, 0)
    for (i in 1:ncomp) {
      int_i <- union(which(linkmat[i, ] > 0), which(linkmat[, 
                                                            i] > 0))
      for (j in 1:ncomp) {
        int_j <- union(which(linkmat[j, ] > 0), which(linkmat[, 
                                                              j] > 0))
        sect <- intersect(int_i, int_j)
        uni <- union(int_i, int_j)
        Cij[i, j] <- length(sect)/length(uni)
      }
    }
    Compart <- (sum(Cij) - ncomp)/ncomp/(ncomp - 1)
    list(N = ncomp, T.. = Throughput, TST = Throughflow, Lint = intlinks, 
         Ltot = links, LD = LD, C = Connectance, Tijbar = Avlinkweight, 
         TSTbar = Avthrflow, Cbar = Compart)
  }}
environment(GenInd2) <- environment(GenInd)

#simulation model####
reps<-3

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

V_all<-c(0.00001,0.001,0.1,1)#c(0.00001,0.0001,0.001,0.01,0.1,1,10) #additive genetic variation in thermal optimum
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
      
      matplot(t(N_store[40,,]), type='l', lty=1)
      matplot(t(z_store[40,,]), type='l', lty=1)
      
      
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
      
      ggplot(N.df,aes(x=Species,y=Env_trait,size=N,fill=Env_trait))+
        geom_point(pch=21)+
        #scale_fill_gradientn(colors=brewer.pal(9,name = "BuGn"))+
        scale_fill_viridis(option = "D")+
        geom_point(aes(x=Species,y=z_initial),col="red",shape=4,size=2)+
        removeGrid()+
        theme_bw()
      # ggsave(filename = paste("./figures/Trait distributions, V = ",V,", d = ",d,".pdf",sep=""),width = 8,height = 6)
      # 
      # 
      ggplot(N.df,aes(x=Species,y=patch,size=N,fill=Env_trait))+
        geom_point(pch=21)+
        #scale_fill_gradientn(colors=rev(brewer.pal(9,name = "YlGnBu")))+
        scale_fill_viridis()+
        removeGrid()+
        theme_bw()+
        ylim(25,50)
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
      
      #network dissimilarity####
      Ints<-matrix(1,species,species)
      diag(Ints)<-0
      
      colnames(Ints)<-rownames(Ints)<-paste("s",1:species)
      
      cut_value<-0
      
      nets_pre<-apply(N_store[,,Burn_in],1,function(x){
        Int_strength<-abs(Ints*rep(x,each=species))
        Int_strength[x==0,]<-0
        Int_strength_cut<-quantile(Int_strength[Int_strength>0],cut_value)#mean(Int_strength[Int_strength>0])
        Int_strength[Int_strength<Int_strength_cut]<-0
        Ints2<-1*Int_strength>0
        hold.df<-t(data.frame(Ints2[x>0,x>0]))
        net1<-graph.adjacency(hold.df)
        return(net1) 
      })
      
      nets_post<-apply(N,1,function(x){
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
                                 Value=c(mean(rowSums(N>0))/mean(rowSums(N_store[,,Burn_in]>0)),
                                         sum(colSums(N)>0)/sum(colSums(N_store[,,Burn_in])>0),
                                         mean(rowSums(N))/mean(rowSums(N_store[,,Burn_in])),
                                         mean(colSums(N>0)[colSums(N>0)>0])/mean(colSums(N_store[,,Burn_in]>0)[colSums(N_store[,,Burn_in]>0)>0]),
                                         mean(apply(z_present,2,sd,na.rm=TRUE),na.rm=TRUE)/mean(apply(z_prechange,2,sd,na.rm=TRUE),na.rm=TRUE),
                                         mean(mean_z_change, na.rm=TRUE),
                                         betalink2(regWeb_pre,regWeb_post,bf = B_jack_diss)$WN,
                                         NetInds$Ltot,
                                         NetInds$LD,
                                         NetInds$C,
                                         NetInds$Cbar
                                         ),
                                 Response=c("Local richness","Regional richness","Local biomass","Range size","Optima sd","Optima change","Network beta","Total links","Link density","Connectance","Compartmentalization"))
      if(d==dispV[1] & V==V_all[1] & r == 1){
        response.df<-response.data1
      } else {response.df<-rbind(response.df,response.data1)}
    }
  }
}

response.df$Response<-factor(response.df$Response,levels=c("Local richness","Regional richness","Local biomass","Range size","Optima sd","Optima change","Network beta","Total links","Link density","Connectance","Compartmentalization"),ordered = TRUE)

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
