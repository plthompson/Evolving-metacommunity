#individual based LV model
library(dplyr)
library(ggplot2)
library(tidyr)
#library(viridis)
#library(vegan)

EM_IBM<-function(species = 80, patches = patches, mutation_r = 0.01, disp = 0.01, type = "co-exist", r = 0.5, changeTime = changeTime) {
  
  burnIn<-10000
  burnOut<-1000
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
  Species_traits<-data.frame(species = 1:species, z = seq(min(Environment$environment),max(Environment$environment),length = species), sig_p = 0.5,r = r, dispersal = rnorm(n = species,mean = disp,sd = disp*0.25), mut_r = rnorm(n = species,mean = mutation_r,sd = mutation_r*0.25), offspring = 0, type = "plant")
  Species_traits$dispersal[Species_traits$dispersal>1]<-1
  Species_traits$dispersal[Species_traits$dispersal<0]<-0
  
  int_scaler<-0.01
  
  B<-matrix(runif(n = species*species,-0.5,-0),nrow = species, ncol = species)*int_scaler*r #stable co-existence
  diag(B)<- -1*int_scaler*r
  
  if(type == "priority"){
    B<-matrix(rnorm(species*species,mean = -0.7,sd=0.2),nrow=species,ncol=species)*int_scaler*r #priority effects
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
      
      ind.df2$z<-rnorm(n = nrow(ind.df2), mean = ind.df2$z, sd = ind.df2$mut_r)
      
      #dispersal
      ind.df2$parent<-NULL
      dispersers<-rbinom(n = nrow(ind.df2),size = 1,prob = ind.df2$dispersal)
      
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
  
  
  return(Nsave)
}

patches<-30
changeTime<-3000

dispV<-c(0.0001,0.0005,0.001,0.01,0.1)
mutationV<-c(0,0.01,0.02,0.03,0.04,0.05,0.07)
results.df<-data.frame()
for(rep in 1:5){
  for(disp in dispV){
    for(mut in mutationV){
      
      Nsave<-EM_IBM(mutation_r = mut, disp = disp, patches = patches, changeTime = changeTime, type = "priority")
      
      analogue<-Nsave %>% 
        filter(time == max(Nsave$time)) %>%
        select(patch,Environment) %>% 
        mutate(analogue = Environment<=(1+patches/2)) %>%
        select(-Environment) %>% 
        group_by(patch) %>% 
        slice(1) %>% 
        ungroup()
      
      Nsave<-left_join(Nsave,analogue)
      
      Local<-Nsave %>% 
        filter(time==10000 | time == max(Nsave$time)) %>% 
        group_by(patch,time, analogue) %>% 
        summarise(S = sum(N>0), N = sum(N)) %>% 
        ungroup() %>% 
        group_by(time, analogue) %>% 
        summarise(S = mean(S), N = mean(N)) %>% 
        ungroup() %>% 
        group_by(analogue) %>% 
        summarise(S_ratio = last(S)/first(S), N_ratio = last(N)/first(N)) %>% 
        mutate(scale = "Local")
      
      Regional<-Nsave %>% 
        filter(time==10000 | time == max(Nsave$time)) %>% 
        group_by(time, analogue, species) %>% 
        summarise(N = sum(N)) %>% 
        ungroup() %>% 
        group_by(time, analogue) %>% 
        summarise(S = sum(N>0), N = sum(N)) %>% 
        ungroup() %>% 
        group_by(analogue) %>% 
        summarise(S_ratio = last(S)/first(S), N_ratio = last(N)/first(N)) %>% 
        mutate(scale = "Regional")
      
      z.results<-Nsave %>% 
        filter(time==10000 | time == max(Nsave$time)) %>% 
        complete(time,nesting(patch, species),fill = list(N = 0)) %>%
        select(-analogue) %>% 
        left_join(analogue) %>% 
        group_by(species, analogue) %>% 
        mutate(persist = last(N)>0) %>% 
        filter(persist == TRUE) %>% 
        ungroup() %>% 
        group_by(species) %>% 
        mutate(initial.z = weighted.mean(z, w = N)) %>% 
        group_by(species,time, analogue, initial.z) %>%
        summarise(z.sd = sd(z, na.rm = TRUE), z = weighted.mean(z,w = N)) %>% 
        filter(time == max(Nsave$time)) %>% 
        ungroup() %>% 
        group_by(species, analogue) %>% 
        summarise(z = mean(z-initial.z, na.rm = TRUE)) %>% 
        ungroup() %>% 
        group_by(analogue) %>% 
        summarise(z_change = mean(z, na.rm = TRUE)) %>% 
        mutate(scale = "Regional")
      
      results.hold<-bind_rows(Local,Regional)
      results.hold<- left_join(results.hold, z.results)
      results.hold$dispersal <- disp
      results.hold$rep <- rep
      results.hold$mutation_rate <- mut
      
      results.df<-bind_rows(results.df,results.hold)
    }
  }
}

save(results.df, file = "./EM_results_priority.RData")

