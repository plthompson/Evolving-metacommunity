#individual based LV model
library(dplyr)
library(ggplot2)
library(tidyr)
library(igraph)
library(NetIndices)
library(vegan)


GenInd2<-function (Flow = NULL, Tij = t(Flow), Import = NULL, Export = NULL,  tol = 0) {
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

EM_IBM<-function(species = 80, patches = patches, mutation_r = 0.01, disp = 0.01, type = "co-exist", r = 0.5, changeTime = changeTime, intersp_mut_var = 0.25, intersp_disp_var = 0.25) {
  
  calc_net_change<-function(N,N1,zmat,z1,BB,trophicV){
    #network dissimilarity####
    Ints<-BB
    diag(Ints)<-0
    colnames(Ints)<-rownames(Ints)<-paste(trophicV,1:species)
    Ints[lower.tri(Ints)]<-0
    Ints[,trophicV=="plants"]<-0
    
    nets_pre<-apply(N1,1,function(x){
      Int_strength<-abs(Ints*rep(x,each=species))
      Int_strength[x==0,]<-0
      hold.df<-t(data.frame(Int_strength[x>0,x>0]))
      net1<-graph.adjacency(hold.df,weighted = T)
      return(net1) 
    })
    
    nets_post<-apply(N,1,function(x){
      Int_strength<-abs(Ints*rep(x,each=species))
      Int_strength[x==0,]<-0
      hold.df<-t(data.frame(Int_strength[x>0,x>0]))
      net1<-graph.adjacency(hold.df,weighted = T)
      return(net1) 
    })
    
    initialInts<-t(matrix(c(Ints)*rep(t(N1),each=species),ncol=dim(N1)[1]))
    initialInts_sub<-initialInts[,colSums(initialInts)!=0]
    
    Temp<-(Environment$environment+changet*changeTime)
    Temp_I<-(Environment$environment)
    
    BC_dist.df<-data.frame()
    for(i in 1:dim(N)[1]){
      comInts<-c(Ints*rep(N[i,],each=species))
      comInts<-comInts[colSums(initialInts)!=0]
      
      BC_dist<-as.matrix(vegdist(abs(rbind(comInts,initialInts_sub)),method="bray",binary = F))[1,-1]
      BC_dist.df<-rbind(BC_dist.df,data.frame(BC_dist=1-min(BC_dist),Patch=i,Closest_patch=which(BC_dist==min(BC_dist)),Analogue=Temp[i] <= max(Temp_I))[1,])
    }
    
    NetInds<-data.frame()
    for(i in BC_dist.df$Patch){
      if(length(E(nets_pre[[filter(BC_dist.df,Patch == i)$Closest_patch[1]]]))>0){
        if(length(E(nets_post[[i]]))>0){
          NetInds<-rbind(NetInds,data.frame(GenInd2(get.adjacency(nets_post[[i]],attr = "weight",sparse = F)))/data.frame(GenInd2(get.adjacency(nets_pre[[filter(BC_dist.df,Patch == i)$Closest_patch[1]]],attr = "weight",sparse = F))))
        } else {
          NetInds<-rbind(NetInds,0/data.frame(GenInd2(get.adjacency(nets_pre[[filter(BC_dist.df,Patch == i)$Closest_patch[1]]],attr = "weight",sparse = F))))
        }
      } else{
        NetInds<-rbind(NetInds,NA)
      }
    }
    
    BC_dist.df<-bind_cols(BC_dist.df,NetInds)
    
    Net_dis.df<-BC_dist.df %>% 
      group_by(Analogue) %>% 
      mutate(Temp_diff=Temp[Patch]-Temp_I[Closest_patch]) %>%
      dplyr::select(-Patch,-Closest_patch) %>% 
      summarise_all(funs(mean(.,na.rm=T))) %>% 
      mutate(Patches=c("no_analogue","analogue")) %>% 
      dplyr::select(-Analogue, -Lint)
    
    names(Net_dis.df)<-c("Network_similarity","Nodes","System_throughput", "System_throughflow","Links","Link_density","Connectance","Average_link_weight","Compartment_throughflow","Compartmentalization","Temp_diff","Patches")
    Net_dis.df<-Net_dis.df %>% 
      gather(key = Response,value = "Value",Network_similarity:Temp_diff) %>% 
      mutate(Dispersal=disp,Genetic_variation=mut,Rep=r, Trophic="all")
    
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
      if(type != "trophic"){
        trophic_select <-c("all")
      }
      for(troph in trophic_select){
        if(troph == "all"){
          sp_select<-1:species
        } else {
          sp_select<-trophicV==troph
        }
        N_select<-N[patch_select,sp_select]
        N1_select<-N1[patch_select,sp_select]
        z_select<-zmat[patch_select,sp_select]
        z1_select<-z1[patch_select,sp_select]
        
        #remove z that is not present
        z_select[N_select==0]<-NA
        z_select[is.infinite(z_select)]<-NA
        z1_select[N1_select==0]<-NA
        z1_select[is.infinite(z_select)]<-NA
        
        #calculate mean change in z weighted by abundance
        mean_z_change<-mean(do.call(c,lapply(1:ncol(N_select),function(i) abs(weighted.mean(z_select[,i],N_select[,i])-weighted.mean(z1_select[,i],N1_select[,i])))),na.rm=T)
        
        #calculate change in z sd - for species that persist
        mean_z_sd_change<-mean(apply(z_select,2,sd,na.rm=T)[colSums(N_select)>0],na.rm=T)/mean(apply(z1_select,2,sd,na.rm=T)[colSums(N_select)>0],na.rm=T)
        
        #calculate change in richness, biomass, range size - range size only considers species that persist
        response.data1<-data.frame(Value=c(mean(rowSums(N_select>0))/mean(rowSums(N1_select>0)), #local richness
                                           mean(sum(colSums(N_select)>0)/sum(colSums(N1_select)>0)), #regional richness
                                           mean(rowSums(N_select)/rowSums(N1_select),na.rm=T), #local biomass
                                           mean((colSums(N_select>0)/colSums(N1_select>0))[colSums(N_select)>0],na.rm=T),#range size
                                           mean_z_change,
                                           mean_z_sd_change),
                                   Response=c("Local S","Regional S","Local biomass","Range size","Optima change","Optima sd"),
                                   Trophic=troph,
                                   Dispersal=disp,
                                   Genetic_variation=mut,
                                   Rep=r,
                                   Patches=ca)
        if(ca=="all" & troph=="all"){
          response.df<-response.data1
        } else {response.df<-rbind(response.df,response.data1)
        }
      }
    }
    response.df<-bind_rows(response.df,Net_dis.df)
    return(response.df)
  }
  
  
  burnIn<-10000
  burnOut<-5000
  Tmax<-burnIn+changeTime+burnOut
  changeMag<-patches/4
  changet<-changeMag/(changeTime-1)
  
  Env_perform<-function(env,z,zmax=NA,sig_p){
    wT<-exp(-((env-z)/2*sig_p)^2)
    wT[wT<0]<-0
    wT<-wT-1
    return(wT)
  }
  
  Environment<-data.frame(patch = 1:patches, environment = c(1:(1+patches/2),(patches/2):2))
  if(mutation_r == "vary"){
    Species_traits<-data.frame(species = 1:species, z = seq(min(Environment$environment),max(Environment$environment),length = species), sig_p = 0.5,r = r, dispersal = rnorm(n = species,mean = disp,sd = disp*intersp_disp_var), mut_r = runif(n = species, min = 0, max = 0.01), offspring = 0, type = "plant")
    } else {
    Species_traits<-data.frame(species = 1:species, z = seq(min(Environment$environment),max(Environment$environment),length = species), sig_p = 0.5,r = r, dispersal = rnorm(n = species,mean = disp,sd = disp*intersp_disp_var), mut_r = rnorm(n = species,mean = as.numeric(mutation_r),sd = as.numeric(mutation_r)*intersp_mut_var), offspring = 0, type = "plant")
  }
  
  Species_traits$dispersal[Species_traits$dispersal>1]<-1
  Species_traits$dispersal[Species_traits$dispersal<0]<-0
  Species_traits$mut_r[Species_traits$mut_r<0]<-0
  
  int_scaler<-0.01
  
  B<-matrix(runif(n = species*species,-0.5,-0),nrow = species, ncol = species)*int_scaler*r #stable co-existence
  diag(B)<- -1*int_scaler*r
  
  if(type == "priority"){
    B<-matrix(runif(n = species*species,-0.5,-0),nrow = species, ncol = species)*int_scaler*r #stable co-existence
    high_comp<-sample(1:species,size = round(species*0.4),replace = FALSE)
    high_comp2<-sample(high_comp,size = length(high_comp),replace=FALSE)
    for(hc in 1:length(high_comp)){
      B[high_comp[hc],high_comp2[hc]]<- runif(min = -1.22,max = -1.18,n = 1) * int_scaler * r
      B[high_comp2[hc],high_comp[hc]]<- runif(min = -1.22,max = -1.18,n = 1) * int_scaler * r
    }
    diag(B)<- -1*int_scaler*r
  }
  
  if(type == "noInt"){
    B<-matrix(0,nrow=species,ncol=species)*int_scaler*r 
    diag(B)<- -1*int_scaler*r
  }
  
  if(type == "fw"){
    nplants<- species*0.5
    nherb<- species*0.3
    npred<- species*0.2
    
    
    plantV<-1:nplants
    herbV<-(nplants+1):(nplants+nherb)
    predV<-(species-npred+1):(species)
    trophicV<-factor(c(rep("plant",nplants),rep("herbivore",nherb),rep("predator",npred)),levels=c("plant","herbivore","predator"),ordered = T)
    
    Species_traits<-data.frame(species = 1:species, z = runif(n = species, min = min(Environment$environment), max = max(Environment$environment)), sig_p = 0.5,r = r, dispersal = rnorm(n = species,mean = disp,sd = disp*0.25), mut_r = rnorm(n = species,mean = mutation_r,sd = mutation_r*0.25), offspring = 0, type = trophicV)
    Species_traits$dispersal[Species_traits$dispersal>1]<-1
    Species_traits$dispersal[Species_traits$dispersal<0]<-0    
    Species_traits$r[Species_traits$type == "herbivore"] <- - 0.5 #need to give no intrinsic rate of growth for herbivores and predators
    Species_traits$r[Species_traits$type == "predator"] <- - 0.5 #need to give no intrinsic rate of growth for herbivores and predators
    
    b11 <--0.9
    b12 <- -1.5
    b21 <- 0.5
    b23 <- - 0.8
    b32 <- 0.5
    bdiag1 <- -1
    bdiag2 <- -0.75
    
    #tritrophic BB Matrix####
    B11 <- matrix(matrix(rnorm(nplants*nplants,mean = -0.7,sd=0.2),nrow=nplants,ncol=nplants)*int_scaler*r,nplants,nplants)
    B12 <- b12*matrix(runif(nplants*nherb),nplants,nherb)*int_scaler*r
    B13 <- matrix(0,nplants,npred)
    B21 <- b21*matrix(runif(nherb*nplants),nherb,nplants)*int_scaler*r
    B22 <- matrix(0,nherb,nherb)
    B23 <- b23*matrix(runif(nherb*npred),nherb,npred)*int_scaler*r
    B31 <- matrix(0,npred,nplants)
    B32 <- b32*matrix(runif(npred*nherb),npred,nherb)*int_scaler*r
    B33 <- matrix(0,npred,npred)
    B <- rbind(cbind(B11 ,B12, B13),cbind(B21,B22, B23),cbind(B31, B32, B33))
    diag(B) <- bdiag1*int_scaler*r
    diag(B[(nplants+1):species,(nplants+1):species])<-bdiag2*int_scaler*r
    #diag(B[(nplants+nherb+1):species,(nplants+nherb+1):species]) <- bdiag2*int_scaler*r
  }
  
  #N<-matrix(sample(size = species*patches,x = 5,replace = T), nrow=patches, ncol=species)
  
  N<-sapply(X = 1:species,FUN = function(X){
    hold<-rep(0,patches)
    hold[sample(1:patches,size = 5,replace=F)]<-5
    return(hold)})
  N<-matrix(c(N),nrow = patches, ncol = species)
  
  
  species_init<-c(sapply(1:species,FUN = function(x) {
    rep(x,colSums(N)[x])
  }))
  
  patch_init<-c()
  for(s in 1:species){
    for(p in 1:patches){
      patch_init<-c(patch_init,rep(p,N[p,s]))
    }
  }
  
  ind.df<-data.frame(individual = 1:sum(N), patch = patch_init, species = species_init)
  ind.df<-left_join(ind.df,Environment)
  ind.df<-left_join(ind.df,Species_traits)
  
  Nsave<-data.frame()
  hold.df<-data.frame(patch = 1:patches,species = rep(1:species,each = patches))
  sampleV<-seq(1100,Tmax, by=100)
  init_col<-seq(100,1000,by=100)
  
  pb <- txtProgressBar(min = 0, max = Tmax, style = 3)
  for(i in 1:Tmax){
    setTxtProgressBar(pb, i)
    if(i %in% sampleV){
      N.df<-ind.df %>% 
        group_by(patch, species, type) %>% 
        summarise(N = n(),upper_z = quantile(z,probs = 0.75),lower_z=quantile(z,probs = 0.25), z=mean(z), Environment = mean(environment), Offspring = sum(offspring),Mean_offspring = mean(offspring)) %>% 
        mutate(time = i)
      
      if(length(unique(N.df$patch))<patches){
        blank.df<-data.frame(patch = 1:patches, Environment = Environment$environment,species = 1, N= 0, upper_z =NA,lower_z = NA,Offspring = NA, Mean_offspring = NA,time = i)
        N.df<-bind_rows(N.df,blank.df %>% 
                          filter(!patch %in% unique(N.df$patch)))
      }
      
      Nsave<-bind_rows(Nsave,N.df)
    }
    
    N.mat<-left_join(hold.df,ind.df %>% 
                       group_by(patch, species) %>% 
                       summarise(N = n()), by = c("patch","species")) %>% 
      spread(key = species,value = N) %>% 
      data.matrix()
    N.mat<-N.mat[,-1]
    
    N.mat[is.na(N.mat)]<-0
    
    if(i == burnIn){
      N1<-N.mat
      
      z.mat<-left_join(hold.df,ind.df %>% 
                         group_by(patch, species) %>% 
                         summarise(z = mean(z)), by = c("patch","species")) %>% 
        spread(key = species,value = z) %>% 
        data.matrix()
      z.mat<-z.mat[,-1]
      
      z1<-z.mat
    }
    
    if(i>burnIn & i<burnIn+changeTime){
      Environment$environment<-Environment$environment+changet
      ind.df$environment<-ind.df$environment+changet
    }
    
    if(i %in% init_col) {
      col.df<-data.frame(individual = (max(ind.df$individual))+1:(species*patches), patch = 1:patches, species = rep(1:species, each = patches))
      col.df<-left_join(col.df,Environment, by = "patch")
      col.df<-left_join(col.df,Species_traits, by = "species")
      ind.df<-bind_rows(ind.df,col.df)
    }
    
    
    ind.df$env_effect<-Env_perform(env = ind.df$environment,z = ind.df$z,sig_p = ind.df$sig_p)*4*r#-abs(ind.df$environment-ind.df$z)*2
    
    ind.df <- left_join(ind.df,data.frame(patch = 1:patches, species = rep(1:species ,each = patches),ints = c(N.mat%*%t(B))), by = c("patch", "species"))
    
    ind.df$offspring<-rpois(n = nrow(ind.df),lambda = exp(ind.df$r+ind.df$ints+ind.df$env_effect))
    
    if(sum(ind.df$offspring)>0){
      reproducers<-ind.df[ind.df$offspring>0,]
      reproducers$parent<-reproducers$individual
      reproducers$individual<-1:nrow(reproducers)
      
      parents<-data.frame(individual = unlist(lapply(unique(reproducers$offspring),FUN = function(x){
        rep(reproducers$parent[reproducers$offspring==x], each = x)
      })))
      
      
      ind.df2<-left_join(parents, ind.df, by = "individual")
      names(ind.df2)[1]<-"parent"
      ind.df2$individual<-(max(ind.df$individual)+1):(max(ind.df$individual)+nrow(ind.df2))
      

      
      if(i>2000){
        ind.df2$z<-rnorm(n = nrow(ind.df2), mean = ind.df2$z, sd = ind.df2$mut_r)
        # mutation<-rnorm(n = nrow(ind.df2),mean = mutation_r, sd = mutation_r*intersp_mut_var)
        # mutation[mutation<0]<-0
        # ind.df2$z<-rnorm(n = nrow(ind.df2), mean = ind.df2$z, sd = mutation)
      } 
      
      #dispersal
      ind.df2$parent<-NULL
      #dispersers<-rbinom(n = nrow(ind.df2),size = 1,prob = ind.df2$dispersal)
      dispersal<-rnorm(n = nrow(ind.df2),mean = disp, sd = disp*intersp_disp_var)
      dispersal[dispersal<0]<-0 ; dispersal[dispersal>1]<-1
      dispersers<-rbinom(n = nrow(ind.df2),size = 1,prob = dispersal)
      
      ind.df2$patch[dispersers]<-ind.df2$patch[dispersers]+sample(c(-1,1),size =sum(dispersers),replace = TRUE)
      ind.df2$patch[ind.df2$patch<1]<-patches
      ind.df2$patch[ind.df2$patch>patches]<-1
      ind.df2$environment<-NULL
      ind.df2<-left_join(ind.df2,Environment, by = "patch")
      
    } else {
      ind.df2<-data.frame()
    }
    
    
    ind.df$survive <- rbinom(n = nrow(ind.df),size = 1,prob = 0.5)
    
    ind.df<-ind.df[ind.df$survive==1,]
    
    ind.df<-bind_rows(ind.df,ind.df2)
    ind.df$ints<-NULL
  }
  close(pb)
  
    N<-N.mat
    
    z.mat<-left_join(hold.df,ind.df %>% 
                       group_by(patch, species) %>% 
                       summarise(z = mean(z)), by = c("patch","species")) %>% 
      spread(key = species,value = z) %>% 
      data.matrix()
    z.mat<-z.mat[,-1]
    
    z<-z.mat
  
  net_change.df<-calc_net_change(N = N,N1 = N1,zmat = z,z1 = z1,BB = B, trophicV = rep("p", species))
  
  
  return(list(Nsave=Nsave,net_change.df=net_change.df))
}

patches<-30
changeTime<-15000

dispV<-c(0,0.00005,0.001,0.01)#c(0.00001,0.00005,0.0001,0.0005,0.001,0.01,0.1)
mutationV<-c(0,0.0025,0.005,0.0075,0.01,0.02, "vary")#c(0,0.001,0.005,0.01,0.02)
results.df<-data.frame()
for(rep in 1:5){
  for(disp in dispV){
    for(mut in mutationV){
      print(paste("Rep - ",rep,"; Disp - ", disp,"; mut - ", mut, sep= ""))
      
      fileConn<-file("./output_r1.txt")
      writeLines(paste("Rep - ",rep,"; Disp - ", disp,"; mut - ", mut, sep= ""), fileConn)
      close(fileConn)
      
      Nsave1<-EM_IBM(mutation_r = mut, disp = disp, patches = patches, changeTime = changeTime, type = "priority")
      
      Nsave<-Nsave1$Nsave
      
      ggplot(filter(Nsave),aes(x=patch,y=time,color=z))+
        geom_point(size = 0.5)+
        facet_wrap(~species)+
        scale_color_viridis()+
        geom_hline(yintercept = c(10000,25000))+
        removeGrid()
      ggsave(paste("disp = ",disp, " mut = ", mut,".png", sep = ""),width=11,height = 8)
      
      results.df<-bind_rows(results.df,Nsave1$net_change.df)
    }
  }
}

save(results.df, file = "./EM_results.RData")

