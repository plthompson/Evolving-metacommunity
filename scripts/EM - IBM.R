#individual based LV model
library(dplyr)
library(ggplot2)
library(tidyr)
#library(viridis)
#library(vegan)

EM_IBM<-function(species = 80, patches = 50, mutation_r = 0.1, disp = 0.01, type = "co-exist", r = 0.5, changeTime = 2000) {
  
  burnIn<-1000
  Tmax<-burnIn+changeTime+burnIn
  changeMag<-patches/4
  changet<-changeMag/(changeTime-1)
  
  Env_perform<-function(env,z,zmax=NA,sig_p){
    wT<-exp(-((env-z)/2*sig_p)^2)
    wT[wT<0]<-0
    wT<-wT-1
    return(wT)
  }
  
  Environment<-data.frame(patch = 1:patches, environment = c(1:(1+patches/2),(patches/2):2))
  Species_traits<-data.frame(species = 1:species, z = seq(min(Environment$environment),max(Environment$environment),length = species), sig_p = 0.5,r = r, dispersal = disp, offspring = 0)
  
  B<-matrix(runif(n = species*species,-0.5,-0),nrow = species, ncol = species)*0.02*r #stable co-existence
  
  if(type == "priority"){
    B<-matrix(rnorm(species*species,mean = -0.9,sd=0.2),nrow=species,ncol=species)*0.02*r #priority effects
  }
  
  if(type == "noInt"){
    B<-matrix(0,nrow=species,ncol=species)*0.02*r #priority effects
  }
  
  diag(B)<- -1*0.02
  
  N<-matrix(sample(size = species*patches,x = 5,replace = T), nrow=patches, ncol=species)
  
  species_init<-unlist(sapply(1:species,FUN = function(x) {
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
  
  pb <- txtProgressBar(min = 0, max = Tmax, style = 3)
  for(i in 1:Tmax){
    setTxtProgressBar(pb, i)
    N.df<-ind.df %>% 
      group_by(patch, species) %>% 
      summarise(N = n(),upper_z = quantile(z,probs = 0.75),lower_z=quantile(z,probs = 0.25), z=mean(z), Environment = mean(environment), Offspring = sum(offspring),Mean_offspring = mean(offspring)) %>% 
      mutate(time = i)
    
    if(length(unique(N.df$patch))<patches){
      blank.df<-data.frame(patch = 1:patches, Environment = Environment$environment,species = 1, N= 0, upper_z =NA,lower_z = NA,Offspring = NA, Mean_offspring = NA,time = i)
      N.df<-bind_rows(N.df,blank.df %>% 
                        filter(!patch %in% unique(N.df$patch)))
    }
    
    Nsave<-bind_rows(Nsave,N.df)
    
    N.mat<-left_join(hold.df,ind.df %>% 
                       group_by(patch, species) %>% 
                       summarise(N = n()), by = c("patch","species")) %>% 
      spread(key = species,value = N) %>% 
      select(-patch) %>% 
      data.matrix()
    N.mat[is.na(N.mat)]<-0
    
    if(i>burnIn & i<burnIn+changeTime){
      Environment$environment<-Environment$environment+changet
      ind.df$environment<-ind.df$environment+changet
    }
    
    ind.df$env_effect<-Env_perform(env = ind.df$environment,z = ind.df$z,sig_p = ind.df$sig_p)*4*r#-abs(ind.df$environment-ind.df$z)*2
    
    ind.df <- merge(ind.df,data.frame(patch = 1:patches, species = rep(1:species ,each = patches),ints = c(N.mat%*%B)))
    
    ind.df<-ind.df %>% 
      group_by(individual) %>% 
      mutate(lambda = exp(r+ints+env_effect),offspring = rpois(n = 1,lambda = lambda))
    
    if(sum(ind.df$offspring)>0){
      reproducers<-ind.df %>% 
        filter(offspring>0)
      reproducers$parent<-reproducers$individual
      reproducers$individual<-1:nrow(reproducers)
      parents<-data.frame(individual = unlist(lapply(reproducers$individual,FUN = function(x) {
        rep(reproducers$parent[x], reproducers$offspring[x])
      })))
      
      
      ind.df2<-merge(parents, ind.df)
      names(ind.df2)[1]<-"parent"
      ind.df2$individual<-(max(ind.df$individual)+1):(max(ind.df$individual)+nrow(ind.df2))
      
      ind.df2$z<-rnorm(n = nrow(ind.df2), mean = ind.df2$z, sd = mutation_r)
      
      #dispersal
      ind.df2$parent<-NULL
      dispersers<-rbinom(n = nrow(ind.df2),size = 1,prob = disp)
      
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
  
  sampleV<-seq(100,Tmax, by = 100)
  output<-filter(Nsave, time %in% sampleV)
  
  
  return(output)
}

dispV<-c(0.0001,0.001,0.01,0.1,1)
mutationV<-c(0,0.01,0.05,0.1,0.5)
results.df<-data.frame()
for(rep in 1:5){
  for(disp in dispV){
    for(mut in mutationV){
      
      Nsave<-EM_IBM(mutation_r = mut, disp = disp)
      
      N_pre<-Nsave %>% 
        filter(time == 1000)
      
      N_post<-Nsave %>% 
        filter(time == 4000)
      
      analogue<-Nsave %>% 
        filter(time == 4000) %>%
        select(patch,Environment) %>% 
        mutate(analogue = Environment<=25) %>%
        select(-Environment) %>% 
        group_by(patch) %>% 
        slice(1) %>% 
        ungroup()
      
      Nsave<-left_join(Nsave,analogue)
      
      Local<-Nsave %>% 
        filter(time==1000 | time == 4000) %>% 
        group_by(patch,time, analogue) %>% 
        summarise(S = sum(N>0), N = sum(N)) %>% 
        ungroup() %>% 
        group_by(time, analogue) %>% 
        summarise(S = mean(S), N = mean(N)) %>% 
        mutate(Scale = "Local")
      
      Regional<-Nsave %>% 
        filter(time==1000 | time == 4000) %>% 
        group_by(time, analogue, species) %>% 
        summarise(N = sum(N)) %>% 
        ungroup() %>% 
        group_by(time, analogue) %>% 
        summarise(S = sum(N>0), N = sum(N)) %>% 
        mutate(Scale = "Regional")
      
      results.hold<-bind_rows(Local,Regional)
      results.hold$dispersal <- disp
      results.hold$rep <- rep
      results.hold$mutation_rate <- mut
      
      results.df<-bind_rows(results.df,results.hold)
    }
  }
}

save(results.df, file = "./EM_results.RData")

