#individual based LV model
library(tidyverse)

species<-8
patches<-9

mutation_r<-0
disp<-0
Tmax<-200

Env_perform<-function(env,z,zmax=NA,sig_p){
  Tmat<-matrix(rep(env,each=species),species,length(env)) 
  wT<-exp(-((Tmat-z)/2*rep(sig_p,each=patches))^2)
  #wT2<-1-((Tmat-z)/(z-z+zmax))^2
  #wT[Tmat>=z]<-wT2[Tmat>=z]
  wT[wT<0]<-0
  wT<-wT-1
  return(wT)
}

B<-matrix(runif(n = species*species,-0.75,-0.25),nrow = species, ncol = species)*0.1 #stable co-existence
diag(B)<- -1*0.1

B<-matrix(rnorm(species*species,mean = -0.9,sd=0.2),nrow=species,ncol=species)*0.1 #priority effects
diag(B)<- -1*0.1

Environment<-data.frame(patch = 1:patches, environment = 20)

N<-matrix(sample(size = species*patches,x = 1:100,replace = T), nrow=patches, ncol=species)

species_init<-unlist(sapply(1:species,FUN = function(x) {
  rep(x,colSums(N)[x])
}))

patch_init<-c()
for(s in 1:species){
  for(p in 1:patches){
    patch_init<-c(patch_init,rep(p,N[p,s]))
  }
}

ind.df<-data.frame(individual = 1:sum(N), patch = patch_init, species = species_init, z = 20, disp_k = 0.001, r=5)
ind.df<-merge(ind.df,Environment)

Nsave<-data.frame()
hold.df<-data.frame(patch = 1:patches,species = rep(1:species,each = patches))

pb <- txtProgressBar(min = 0, max = Tmax, style = 3)
for(i in 1:Tmax){
  setTxtProgressBar(pb, i)
  N.df<-ind.df %>% 
    group_by(patch, species) %>% 
    summarise(N = n(), z=mean(z)) %>% 
    mutate(time = i)
  
  Nsave<-bind_rows(Nsave,N.df)
  
N.mat<-left_join(hold.df,ind.df %>% 
    group_by(patch, species) %>% 
    summarise(N = n()), by = c("patch","species")) %>% 
    spread(key = species,value = N) %>% 
  select(-patch) %>% 
  data.matrix()
N.mat[is.na(N.mat)]<-0

ind.df$env_effect<--abs(ind.df$environment-ind.df$z)*2

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
} else {
  ind.df2<-data.frame()
}

#dispersal
dispersers<-rbinom(n = nrow(ind.df2),size = 1,prob = disp)
ind.df2$patch[dispersers]<-sample(1:patches,size =sum(dispersers),replace = TRUE)
  

ind.df$survive <- rbinom(n = nrow(ind.df),size = 1,prob = 0.95)

ind.df<-ind.df[ind.df$survive==1,]

ind.df<-bind_rows(ind.df,ind.df2)
ind.df$ints<-NULL
}
close(pb)


ggplot(Nsave,aes(x = time, y = z, color=as.factor(species), group = species))+
  geom_line()+
  facet_wrap(~patch)

ggplot(Nsave,aes(x = time, y = N, color=as.factor(species), group = species))+
  geom_line()+
  facet_wrap(~patch)


