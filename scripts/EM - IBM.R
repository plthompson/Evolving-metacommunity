#individual based LV model
library(tidyverse)
library(viridis)
library(vegan)

EM_IBM<-function(species = 20, patches = 40, mutation_r = 0.1, disp = 0.01, type = "co-exist") {
  
  burnIn<-500
  changeTime<-500
  Tmax<-burnIn+changeTime
  changeMag<-10
  changet<-changeMag/(changeTime-1)
  
  Env_perform<-function(env,z,zmax=NA,sig_p){
    wT<-exp(-((env-z)/2*sig_p)^2)
    wT[wT<0]<-0
    wT<-wT-1
    return(wT)
  }
  
  Environment<-data.frame(patch = 1:patches, environment = c(1:(1+patches/2),(patches/2):2))
  Species_traits<-data.frame(species = 1:species, z = seq(min(Environment$environment),max(Environment$environment),length = species), sig_p = 0.5,r = 5, dispersal = disp)
  
  B<-matrix(runif(n = species*species,-0.75,-0.25),nrow = species, ncol = species)*0.1 #stable co-existence
  
  if(type == "priority"){
    B<-matrix(rnorm(species*species,mean = -0.9,sd=0.2),nrow=species,ncol=species)*0.1 #priority effects
  }
  
  diag(B)<- -1*0.1
  
  N<-matrix(sample(size = species*patches,x = 5:60,replace = T), nrow=patches, ncol=species)
  
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
      summarise(N = n(),upper_z = quantile(z,probs = 0.75),lower_z=quantile(z,probs = 0.25), z=mean(z)) %>% 
      mutate(time = i)
    
    Nsave<-bind_rows(Nsave,N.df)
    
    N.mat<-left_join(hold.df,ind.df %>% 
                       group_by(patch, species) %>% 
                       summarise(N = n()), by = c("patch","species")) %>% 
      spread(key = species,value = N) %>% 
      select(-patch) %>% 
      data.matrix()
    N.mat[is.na(N.mat)]<-0
    
    if(i>burnIn){
      Environment$environment<-Environment$environment+changet
      ind.df$environment<-ind.df$environment+changet
    }
    
    ind.df$env_effect<-Env_perform(env = ind.df$environment,z = ind.df$z,sig_p = ind.df$sig_p)*5#-abs(ind.df$environment-ind.df$z)*2
    
    ind.df <- merge(ind.df,data.frame(patch = 1:patches, species = rep(1:species ,each = patches),ints = c(N.mat%*%B)))
    
    ind.df<-ind.df %>% 
      group_by(individual) %>% 
      mutate(lambda = exp(r+ints+env_effect),offspring = rpois(n = 1,lambda = lambda), survive = rbinom(n = 1,size = 1,prob = lambda))
    
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
    
    
    ind.df$survive <- rbinom(n = nrow(ind.df),size = 1,prob = 0.95)
    
    ind.df<-ind.df[ind.df$survive==1,]
    
    ind.df<-bind_rows(ind.df,ind.df2)
    ind.df$ints<-NULL
  }
  close(pb)
  output<-filter(Nsave, time %in% round(seq(1,Tmax, length = 100)))
  return(output)
}

Nsave<-EM_IBM()

ggplot(Nsave,aes(x = time, y = z, color=as.factor(species), fill=as.factor(species), group = species))+
  geom_vline(xintercept = 500, lty = 2)+
  geom_ribbon(aes(ymin=lower_z,ymax = upper_z),alpha = 0.4,col=NA)+
  geom_line()+
  facet_wrap(~patch)+
  theme_bw()

ggplot(Nsave,aes(x = time, y = z, color=as.factor(patch), fill=as.factor(patch), group = patch))+
  geom_vline(xintercept = 500, lty = 2)+
  geom_ribbon(aes(ymin=lower_z,ymax = upper_z),alpha = 0.4,col=NA)+
  geom_line()+
  facet_wrap(~species)+
  theme_bw()

ggplot(Nsave,aes(x = time, y = N, color=as.factor(species), group = species))+
  geom_vline(xintercept = 500, lty = 2)+
  geom_line()+
  scale_color_viridis(discrete = TRUE)+
  facet_wrap(~patch)

Nsave %>%
  group_by(species,time) %>% 
  summarise(Regional_abundance = sum(N)) %>% 
  ungroup() %>% 
  group_by(time) %>% 
  summarise(Regional_S = n()) %>% 
  ggplot(aes(x=time, y = Regional_S))+
  geom_line()+
  theme_bw()


data.wide<-Nsave %>%
  filter(time %in% round(seq(100,Tmax, length = 50))) %>% 
  select(patch,species,time,N) %>% 
  spread(key = species,value = N,fill = 0) %>%
  select(-patch,-time)

data.wide<-data.wide[,colSums(data.wide)>0]
data.wide<-data.wide[rowSums(data.wide)>0,]

nmds<-data.wide %>%
  decostand(method = "hellinger") %>% 
  vegdist(method = "euclidean") %>% 
  metaMDS(trymax = 1)

data.scores <- as.data.frame(scores(nmds))

data.scores<-bind_cols(Nsave %>%
                         filter(time %in% round(seq(100,Tmax, length = 50))) %>% 
                         select(patch, species, time, N) %>% 
                         spread(key = species, value = N, fill = 0) %>% 
                         select(patch, time), data.scores)

ggplot(filter(data.scores, time<burnIn),aes(x=NMDS1,y = NMDS2, group=as.factor(patch), color = as.factor(patch)))+
  geom_point()+
  geom_path()+
  scale_color_viridis(discrete = TRUE)+
  theme_bw()

ggplot(filter(data.scores, time>burnIn),aes(x=NMDS1,y = NMDS2, group=as.factor(patch), color = as.factor(patch)))+
  geom_point()+
  geom_path()+
  scale_color_viridis(discrete = TRUE)+
  theme_bw()


ggplot(filter(data.scores, time>100),aes(x=NMDS1,y = NMDS2, group=as.factor(patch), color = time))+
  geom_point()+
  geom_path()+
  scale_color_viridis()+
  theme_bw()

