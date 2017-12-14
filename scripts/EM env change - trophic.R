#packages####
library(ggExtra)
library(igraph)
library(tidyverse)
library(viridis)
library(betalink)
library(NetIndices)
library(vegan)

source("./functions/EM env change functions.r")

#variables to contrast####
reps<-10
V_all<-c(0.001,0.01,0.1,1,10,100) #additive genetic variation in thermal optimum
dispV<-c(0.00001,0.0001,0.001,0.01,0.1,0.5)


#simulation parameters####
species<-40
patches<-50

nplants<-species*0.5
nherb<-species*0.3
npred<-species*0.2

#Environmental change####
eco_burn<-5000
evo_burn<-10000
burn_in<-eco_burn+evo_burn
change<-7000
change_mag<-7

ChangeV<-c(rep(0,burn_in),seq(0,change_mag-1,length=change))
Tmax<-length(ChangeV)

#Species interactions####
weight=1/80*3

plantV<-1:nplants
herbV<-(nplants+1):(nplants+nherb)
predV<-(species-npred+1):(species)
trophicV<-factor(c(rep("plant",nplants),rep("herbivore",nherb),rep("predator",npred)),levels=c("plant","herbivore","predator"),ordered = T)

results.df<-data.frame()
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
  
  C3<-c(rep(0.2,nplants),rep(0,species-nplants))
  
  #Environment
  Temp_I<-c(1:(1+patches/2),(patches/2):2)
  
  #Environmental optima
  z<-runif(n = species,min = min(Temp_I),max=max(Temp_I))
  z<-c(seq(min(Temp_I),max(Temp_I),length=nplants),seq(min(Temp_I),max(Temp_I),length=nherb),seq(min(Temp_I),max(Temp_I),length=npred))
  
  sig_p<-0.15
  
  #dispersal####
  # edges<-rep(1:patches,each=2)
  # edges<-c(edges[-1],edges[1])
  # graph<-make_graph(edges, directed = FALSE)
  # dist_mat<-distances(graph)
  # 
  # dd<-0.5#rnorm(n = species,mean=0.5,sd=0.1)
  # disp_matrix<-exp(-dd*dist_mat)
  # diag(disp_matrix)<-0
  # disp_matrix<-disp_matrix/rowSums(disp_matrix)
  # 
  # ddV<-rnorm(species,mean = dd,sd = 0.5*0.1)
  # 
  # disp_array<-array(NA,dim=c(patches,patches,species))
  # for(s in 1:species){
  #   disp_matrix<-exp(-ddV[s]*dist_mat)
  #   diag(disp_matrix)<-0
  #   disp_matrix<-disp_matrix/rowSums(disp_matrix)
  #   disp_array[,,s]<-disp_matrix
  # }
  # disp_array_initial<-disp_array
  # #matplot(disp_array[,25,],type='l')
  
  disp_matrix<-matrix(0,nrow = patches,ncol = patches)
  for(p in 1:(patches-1)){
    disp_matrix[p,p+1]<-0.5
    disp_matrix[p+1,p]<-0.5
  }
  disp_matrix[1,patches]<-0.5
  disp_matrix[patches,1]<-0.5
  
  for(disp in dispV){
    print(paste("Rep = ",r,", Disp = ",disp,sep=""))
    dV<-rnorm(species,mean = disp,sd=disp*0.25)
    output<-do.call(rbind,lapply(1:length(V_all),FUN = EMC_fun,no_troph_evo = TRUE))
    results.df<-rbind(results.df,output)
  }
}

results.df$Response<-factor(results.df$Response,levels=c("Local S","Regional S","Local biomass","Range size","Optima sd","Optima change","Temp_diff","Network_similarity","Nodes","Links","Link_density","Connectance","Average_link_weight","Compartmentalization","System_throughput","System_throughflow","Compartment_throughflow"),ordered = TRUE)

save(results.df,file = "./workspace/Evolving MC - change no TL evo.RData")

response_means<-results.df %>% 
  group_by(Response,Dispersal,Genetic_variation,Patches,Trophic) %>% 
  summarise(Mean=median(Value, na.rm=T),Lower=quantile(Value,probs = 0.25,na.rm=T),Upper = quantile(Value,probs = 0.75,na.rm=T))

#main figures####
ggplot(filter(response_means,
              Response=="Regional S",
              Trophic== "all",
              Patches!="all"),
       aes(x=Dispersal,y=Mean,group=Genetic_variation, color=as.character(Genetic_variation),fill=as.character(Genetic_variation)))+
  #scale_color_viridis(trans="log",breaks=V_all)+
  geom_ribbon(aes(ymin=Lower,ymax=Upper),alpha=0.2,color=NA)+
  geom_point()+
  scale_color_discrete(name="Genetic variation")+
  geom_line()+
  scale_fill_discrete(guide=FALSE)+
  facet_wrap(~Patches)+
  scale_x_log10(breaks=c(1e-5,1e-4,1e-3,1e-2,1e-1,1))+
  theme_bw()+
  ylab("Regional species richness")+
  removeGrid()+
  theme(strip.background = element_rect(colour="white", fill="white"))
ggsave("./figures/Regional diversity.png",width=7,height=3.5)

ggplot(filter(response_means,
              Response=="Network_similarity",
              Trophic== "all"),
       aes(x=Dispersal,y=Mean,group=Genetic_variation, color=as.character(Genetic_variation),fill=as.character(Genetic_variation)))+
  #scale_color_viridis(trans="log",breaks=V_all)+
  geom_ribbon(aes(ymin=Lower,ymax=Upper),alpha=0.2,color=NA)+
  geom_point()+
  scale_color_discrete(name="Genetic variation")+
  geom_line()+
  scale_fill_discrete(guide=FALSE)+
  facet_wrap(~Patches,scales = "free_y")+
  scale_x_log10(breaks=c(1e-5,1e-4,1e-3,1e-2,1e-1,1))+
  theme_bw()+
  ylab("Network similarity")+
  removeGrid()+
  theme(strip.background = element_rect(colour="white", fill="white"))
ggsave("./figures/Community reorganization.png",width=7,height=3.5)

ggplot(filter(response_means,
              Response=="Link_density" |
                Response=="Connectance" |
                Response=="Compartmentalization"|
                Response=="Average_link_weight",
              Trophic== "all"),
       aes(x=Dispersal,y=Mean,group=Genetic_variation, color=as.character(Genetic_variation),fill=as.character(Genetic_variation)))+
  scale_color_discrete(name="Genetic variation")+
  scale_fill_discrete(guide=FALSE)+
  geom_ribbon(aes(ymin=Lower,ymax=Upper),alpha=0.2,color=NA)+
  geom_point()+
  geom_line()+
  facet_grid(Patches~Response)+
  scale_x_log10(breaks=c(1e-5,1e-4,1e-3,1e-2,1e-1,1))+
  theme_bw()+
  removeGrid()+
  theme(strip.background = element_rect(colour="white", fill="white"))
ggsave(filename = "./figures/Network change.png",width = 12,height =6)

ggplot(filter(response_means,
              Response=="System_throughput" |
                Response=="System_throughflow" |
                Response=="Compartment_throughflow" ,
              Trophic== "all")
       ,aes(x=Dispersal,y=Mean,group=Genetic_variation, color=as.character(Genetic_variation),fill=as.character(Genetic_variation)))+
  #scale_color_viridis(trans="log",breaks=V_all)+
  geom_ribbon(aes(ymin=Lower,ymax=Upper),alpha=0.2,color=NA)+
  geom_point()+
  geom_line()+
  facet_grid(Response~Patches,scales = "free_y")+
  scale_x_log10()+
  theme_bw()+
  removeGrid()


ggplot(filter(response_means,
              Response=="Regional S",
              Trophic== "predator",
              Patches!="all"),
       aes(x=Dispersal,y=Mean,group=Genetic_variation, color=as.character(Genetic_variation),fill=as.character(Genetic_variation)))+
  #scale_color_viridis(trans="log",breaks=V_all)+
  geom_ribbon(aes(ymin=Lower,ymax=Upper),alpha=0.2,color=NA)+
  geom_point()+
  scale_color_discrete(name="Genetic variation")+
  geom_line()+
  scale_fill_discrete(guide=FALSE)+
  facet_wrap(~Patches)+
  scale_x_log10(breaks=c(1e-5,1e-4,1e-3,1e-2,1e-1,1))+
  scale_y_continuous(breaks=seq(0,1,by=0.25))+
  theme_bw()+
  ylab("Regional predator richness")+
  removeGrid()+
  theme(strip.background = element_rect(colour="white", fill="white"))
ggsave("./figures/Predator diversity.png",width=7,height=3.5)



#old figures####

ggplot(filter(response_means,
              Response=="Local S" |
                Response=="Regional S" |
                Response=="Local biomass" |
                Response=="Range size",
              Trophic=="all")
       ,aes(x=Dispersal,y=Mean,group=Genetic_variation, color=as.character(Genetic_variation),fill=as.character(Genetic_variation)))+
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
              Response=="Local S" |
                Response=="Regional S" |
                Response=="Local biomass" |
                Response=="Range size",
              Trophic=="plant")
       ,aes(x=Dispersal,y=Mean,group=Genetic_variation, color=as.character(Genetic_variation),fill=as.character(Genetic_variation)))+
  #scale_color_viridis(trans="log",breaks=V_all)+
  geom_ribbon(aes(ymin=Lower,ymax=Upper),alpha=0.2,color=NA)+
  geom_point()+
  geom_line()+
  facet_grid(Response~Patches,scales = "free_y")+
  scale_x_log10()+
  theme_bw()+
  removeGrid()

ggplot(filter(response_means,
              Response=="Local S" |
                Response=="Regional S" ,
              Patches !="all",
              Trophic=="predator")
       ,aes(x=Dispersal,y=Mean,group=Genetic_variation, color=as.character(Genetic_variation),fill=as.character(Genetic_variation)))+
  #scale_color_viridis(trans="log",breaks=V_all)+
  geom_ribbon(aes(ymin=Lower,ymax=Upper),alpha=0.2,color=NA)+
  geom_point()+
  geom_line()+
  facet_grid(Response~Patches,scales = "free_y")+
  scale_x_log10()+
  theme_bw()+
  removeGrid()

ggplot(filter(response_means,
              Response=="Local S" |
                Response=="Regional S" |
                Response=="Local biomass" |
                Response=="Range size",
              Trophic=="herbivore")
       ,aes(x=Dispersal,y=Mean,group=Genetic_variation, color=as.character(Genetic_variation),fill=as.character(Genetic_variation)))+
  #scale_color_viridis(trans="log",breaks=V_all)+
  geom_ribbon(aes(ymin=Lower,ymax=Upper),alpha=0.2,color=NA)+
  geom_point()+
  geom_line()+
  facet_grid(Response~Patches,scales = "free_y")+
  scale_x_log10()+
  theme_bw()+
  removeGrid()

ggplot(filter(response_means,
              Response=="Local S" |
                Response=="Regional S" |
                Response=="Local biomass" |
                Response=="Range size",
              Trophic=="predator")
       ,aes(x=Dispersal,y=Mean,group=Genetic_variation, color=as.character(Genetic_variation),fill=as.character(Genetic_variation)))+
  #scale_color_viridis(trans="log",breaks=V_all)+
  geom_ribbon(aes(ymin=Lower,ymax=Upper),alpha=0.2,color=NA)+
  geom_point()+
  geom_line()+
  facet_grid(Response~Patches,scales = "free_y")+
  scale_x_log10()+
  theme_bw()+
  removeGrid()

ggplot(filter(response_means,
              Response=="Optima sd" |
                Response=="Optima change" | 
              Response =="Temp_diff",
              Trophic=="all")
       ,aes(x=Dispersal,y=Mean,group=Genetic_variation, color=as.character(Genetic_variation)))+
  #scale_color_viridis(trans="log",breaks=V_all)+
  geom_point()+
  geom_line()+
  facet_grid(Response~Patches,scales = "free_y")+
  scale_x_log10()+
  theme_bw()+
  removeGrid()
#ggsave(filename = "./figures/Changing environment/Adaptation.pdf",width = 13,height =8 )

ggplot(filter(response_means,
              Response=="Optima sd" |
                Response=="Optima change",
              Trophic=="plant")
       ,aes(x=Dispersal,y=Mean,group=Genetic_variation, color=as.character(Genetic_variation)))+
  #scale_color_viridis(trans="log",breaks=V_all)+
  geom_point()+
  geom_line()+
  facet_grid(Response~Patches,scales = "free_y")+
  scale_x_log10()+
  theme_bw()+
  removeGrid()

ggplot(filter(response_means,
              Response=="Optima sd" |
                Response=="Optima change",
              Trophic=="herbivore")
       ,aes(x=Dispersal,y=Mean,group=Genetic_variation, color=as.character(Genetic_variation)))+
  #scale_color_viridis(trans="log",breaks=V_all)+
  geom_point()+
  geom_line()+
  facet_grid(Response~Patches,scales = "free_y")+
  scale_x_log10()+
  theme_bw()+
  removeGrid()

ggplot(filter(response_means,
              Response=="Optima sd" |
                Response=="Optima change",
              Trophic=="predator")
       ,aes(x=Dispersal,y=Mean,group=Genetic_variation, color=as.character(Genetic_variation)))+
  #scale_color_viridis(trans="log",breaks=V_all)+
  geom_point()+
  geom_line()+
  facet_grid(Response~Patches,scales = "free_y")+
  scale_x_log10()+
  theme_bw()+
  removeGrid()

ggplot(filter(response_means,
                Response=="Link_density" |
                Response=="Connectance" |
                Response=="Compartmentalization"|
                Response=="Average_link_weight",
              Trophic== "all")
       ,aes(x=Dispersal,y=Mean,group=Genetic_variation, color=as.character(Genetic_variation),fill=as.character(Genetic_variation)))+
  #scale_color_viridis(trans="log",breaks=V_all)+
  geom_ribbon(aes(ymin=Lower,ymax=Upper),alpha=0.2,color=NA)+
  geom_point()+
  geom_line()+
  facet_grid(Response~Patches,scales = "free_y")+
  scale_x_log10()+
  theme_bw()+
  removeGrid()
#ggsave(filename = "./figures/Changing environment/Network change.pdf",width = 13,height =8 )

ggplot(filter(response_means,
              Response=="System_throughput" |
                Response=="System_throughflow" |
                Response=="Compartment_throughflow" ,
              Trophic== "all")
       ,aes(x=Dispersal,y=Mean,group=Genetic_variation, color=as.character(Genetic_variation),fill=as.character(Genetic_variation)))+
  #scale_color_viridis(trans="log",breaks=V_all)+
  geom_ribbon(aes(ymin=Lower,ymax=Upper),alpha=0.2,color=NA)+
  geom_point()+
  geom_line()+
  facet_grid(Response~Patches,scales = "free_y")+
  scale_x_log10()+
  theme_bw()+
  removeGrid()



