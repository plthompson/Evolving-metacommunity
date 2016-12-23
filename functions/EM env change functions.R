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

T_perform<-function(Temp,z,zmax,sig_p){
  Tmat<-matrix(Temp,length(Temp),species) 
  wT<-exp(-((Tmat-z)/2*sig_p)^2)
  #wT2<-1-((Tmat-z)/(z-z+zmax))^2
  #wT[Tmat>=z]<-wT2[Tmat>=z]
  wT[wT<0]<-0
  return(wT)
}

Env_perform<-function(env,z,zmax=NA,sig_p){
  Tmat<-matrix(rep(env,each=species),species,length(env)) 
  wT<-exp(-((Tmat-z)/2*rep(sig_p,each=patches))^2)
  #wT2<-1-((Tmat-z)/(z-z+zmax))^2
  #wT[Tmat>=z]<-wT2[Tmat>=z]
  wT[wT<0]<-0
  wT<-wT-1
  return(wT)
}

Env_perform2<-function(env,z,zmax=NA,sig_p){
  Tmat<-matrix(rep(env,each=species),species,length(env)) 
  wT<-exp(-((Tmat-z)/2*rep(sig_p,each=patches))^2)
  #wT2<-1-((Tmat-z)/(z-z+zmax))^2
  #wT[Tmat>=z]<-wT2[Tmat>=z]
  wT[wT<0]<-0
  wT<-wT-1
  return(wT)
}

calc_net_change<-function(N,N1,zmat,z1){
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
    trophic_select<-c("all","plant","herbivore","predator")
    for(troph in trophic_select){
      if(troph == "all"){
        sp_select<-1:species
      } else {
        sp_select<-trophicV==troph
      }
      N_select<-N[sp_select,patch_select]
      N1_select<-N1[sp_select,patch_select]
      z_select<-zmat[sp_select,patch_select]
      z1_select<-z1[sp_select,patch_select]
      
      #remove z that is not present
      z_select[N_select==0]<-NA
      z_select[is.infinite(z_select)]<-NA
      z1_select[N1_select==0]<-NA
      z1_select[is.infinite(z_select)]<-NA
      
      #calculate mean change in z weighted by abundance
      mean_z_change<-mean(do.call(c,lapply(1:nrow(N_select),function(i) abs(weighted.mean(z_select[i,],N_select[i,])-weighted.mean(z1_select[i,],N1_select[i,])))),na.rm=T)
      
      #calculate change in z sd - for species that persist
      mean_z_sd_change<-mean(apply(z_select,1,sd,na.rm=T)[rowSums(N_select)>0],na.rm=T)/mean(apply(z1_select,1,sd,na.rm=T)[rowSums(N_select)>0],na.rm=T)
      
      #network dissimilarity####
      if(troph=="all"){
        Ints<-matrix(1,species,species)
        diag(Ints)<-0
        
        colnames(Ints)<-rownames(Ints)<-paste(trophicV,1:species)
        
        cut_value<-0.5
        
        nets_pre<-apply(N1_select,2,function(x){
          Int_strength<-abs(Ints*rep(x,each=species))
          Int_strength[x==0,]<-0
          Int_strength_cut<-quantile(Int_strength[Int_strength>0],cut_value)
          Int_strength[Int_strength<Int_strength_cut]<-0
          #Int_strength[colSums(Int_strength)==0,]<-0 # to remove the species that have no reciprocal effect
          Ints2<-1*Int_strength>0
          hold.df<-t(data.frame(Ints2[x>0,x>0]))
          net1<-graph.adjacency(hold.df)
          return(net1) 
        })
        
        nets_post<-apply(N[,patch_select],2,function(x){
          Int_strength<-abs(Ints*rep(x,each=species))
          Int_strength[x==0,]<-0
          Int_strength_cut<-quantile(Int_strength[Int_strength>0],cut_value)
          Int_strength[Int_strength<Int_strength_cut]<-0
          #Int_strength[colSums(Int_strength)==0,]<-0 # to remove the species that have no reciprocal effect
          Ints2<-1*Int_strength>0
          hold.df<-t(data.frame(Ints2[x>0,x>0]))
          net1<-graph.adjacency(hold.df)
          return(net1) 
        })
        
        regWeb_pre<-metaweb(nets_pre)
        regWeb_post<-metaweb(nets_post)
        
        #network_betaplot(regWeb_post,regWeb_pre)
        
        NetInds<-data.frame(GenInd2(get.adjacency(regWeb_post,sparse = F)))/data.frame(GenInd2(get.adjacency(regWeb_pre,sparse = F)))                             
        
        NetInds$NetDiss<-1-betalink2(regWeb_pre,regWeb_post,bf = B_jack_diss)$WN
        
      } else {
        NetInds<-data.frame(N=NA,T..=NA,TST=NA,Lint=NA,Ltot=NA,LD=NA,C=NA,Tijbar=NA,TSTbar=NA,Cbar=NA,NetDiss=NA)
      }
      
      #calculate change in richness, biomass, range size - range size only considers species that persist
      response.data1<-data.frame(Value=c(mean(colSums(N_select>0))/mean(colSums(N1_select>0)), #local richness
                                         mean(sum(rowSums(N_select)>0)/sum(rowSums(N1_select)>0)), #regional richness
                                         mean(colSums(N_select)/colSums(N1_select),na.rm=T), #local biomass
                                         mean((rowSums(N_select>0)/rowSums(N1_select>0))[rowSums(N_select)>0],na.rm=T),#range size
                                         mean_z_change,
                                         mean_z_sd_change,
                                         NetInds$NetDiss,
                                         NetInds$Ltot,
                                         NetInds$LD,
                                         NetInds$C,
                                         NetInds$Cbar),
                                 Response=c("Local S","Regional S","Local biomass","Range size","Optima change","Optima sd","Network simmilarity","Total links","Link density","Connectance","Compartmentalization"),
                                 Trophic=troph,
                                 Dispersal=disp,
                                 Genetic_variation=V,
                                 Rep=r,
                                 Patches=ca)
      if(ca=="all" & troph=="all"){
        response.df<-response.data1
      } else {response.df<-rbind(response.df,response.data1)
      }
    }
  }
  return(response.df)
}