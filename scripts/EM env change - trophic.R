Env_perform<-function(env,z,zmax,sig_p){
  Tmat<-matrix(env,length(env),species) 
  wT<-exp(-((Tmat-z)/2*rep(sig_p,each=patches))^2)
  #wT2<-1-((Tmat-z)/(z-z+zmax))^2
  #wT[Tmat>=z]<-wT2[Tmat>=z]
  wT[wT<0]<-0
  return(wT)
}

reps<-3

patches<-50

Burn_in<-20000
Change_time<-40000
Burn_out<-5000
Tmax<-Burn_in+Change_time+Burn_out

#specify the environment
env_min<-0
env_max_i<-0.8
env_max<-1.2
Temp_initial<-c(seq(env_min,env_max_i,length=(1+patches/2)),rev(seq(env_min,env_max_i,length=(1+patches/2))[-c(1,(1+patches/2))]))
Temp_changeV<-c(rep(0,Burn_in),seq(0,env_max-env_max_i,length=Change_time),rep(env_max-env_max_i,Burn_out))

#create food webs####
species<-40

nplants<-species*0.5
nherb<-species*0.3
npred<-species*0.2

trophV<-c(rep("p", nplants),rep('h',nherb),rep("c",npred))

b11=-0.1
b12=-0.3
b21=0.2
b23=-0.15
#b32=0.15
b32=0.12
bdiag1=-.008
bdiag2=0

weight<-(1/species)*3

B11 <- b11*matrix(runif(nplants*nplants),nplants,nplants)
B12 <- b12*matrix(runif(nplants*nherb),nplants,nherb)
B13 <- matrix(0,nplants,npred)
B21 <- b21*matrix(runif(nherb*nplants),nherb,nplants)
B22 <- matrix(0,nherb,nherb)
B23 <- b23*matrix(runif(nherb*npred),nherb,npred)
B31 <- matrix(0,npred,nplants)
B32 <- b32*matrix(runif(npred*nherb),npred,nherb)
B33 <- matrix(0,npred,npred)
BB  <- rbind(cbind(B11 ,B12, B13),cbind(B21,B22, B23),cbind(B31, B32, B33))
BB <- weight*BB
diag(BB) <- bdiag1
diag(BB[(nplants+nherb+1):species,(nplants+nherb+1):species]) <- bdiag2

colnames(BB)<-rownames(BB)<-c(paste("p",1:nplants),paste('h',1:nherb),paste("c",1:npred))

C<-c(rep(0.04,nplants),rep(0,species-nplants))


#dispersal####
edges<-rep(1:patches,each=2)
edges<-c(edges[-1],edges[1])
graph<-make_graph(edges, directed = FALSE)
dist_mat<-distances(graph)

dd<-0.5#rnorm(n = species,mean=0.5,sd=0.1)
disp_matrix<-exp(-dd*dist_mat)
diag(disp_matrix)<-0
disp_matrix<-disp_matrix/rowSums(disp_matrix)

ddV<-rnorm(species,mean = dd,sd = 0.1)

disp_array<-array(NA,dim=c(patches,patches,species))
for(s in 1:species){
  disp_matrix<-exp(-ddV[s]*dist_mat)
  diag(disp_matrix)<-0
  disp_matrix<-disp_matrix/rowSums(disp_matrix)
  disp_array[,,s]<-disp_matrix
}
disp_array_initial<-disp_array


tr_mc_func<-function(disp=0.001,troph_gen=1, trophic_env_specificity=F,trophic_dispersal_rate=F,vary_kernel=F, vary_kernel_byTrophic=F,active_disp=FALSE,selective_error=0,envChange=F,V=0.01){
  #variation in genetic variation
  V_species<-rnorm(n = species,mean = V,sd=V*0.1)
  
  #variation in dispersal
  dispV<-rnorm(n = species,mean = disp,sd = disp*0.1)
  if(trophic_dispersal_rate==T){
    dispV<-c(sort(dispV)[trophV=="p"][sample(1:nplants,size = nplants,replace=F)],
             sort(dispV)[trophV=="h"][sample(1:nherb,size = nherb,replace=F)],
             sort(dispV)[trophV=="c"][sample(1:npred,size = npred,replace=F)])
  }
  
  #environmental optima####
  z<-c(seq(from = env_min,to = env_max_i,length=nplants),
       seq(from = env_min,to = env_max_i,length=nherb),
       seq(from = env_min,to = env_max_i,length=npred)) #environmental optima
  
  sig_p<-rep(3,species) #rise in performance as enviromental match increases
  if(trophic_env_specificity==T){
    #sig_p<-c(rep(5,nplants),rep(3,nherb),rep(1,npred))
    sig_p<-c(rep(3,nplants),rep(2,nherb),rep(1,npred))
  }
  
  A<-Env_perform(env = Temp_initial,z,max(C),sig_p)*0.2-0.2
  
  #trophic feeding specificity
  BB1<-BB
  for(j in (nplants+1):species){
    int_select<-rbinom(species,size = 1,prob = troph_gen)
    BB1[,j]<-BB1[,j]*int_select
    BB1[j,]<-BB1[j,]*int_select
  }
  
  #initial conditions for metacommunity####
  N<-matrix(rep(c(rep(4,nplants),rep(0.1,nherb),rep(0.05,npred)),each=patches),patches,species)
  
  Nstore<-array(NA,dim=c(patches,species,Tmax))
  
  #simulation####
  if(active_disp==FALSE){
    for(i in 1:Tmax){
      Temp<-Temp_initial+Temp_changeV[i]
      A<-Env_perform(env = Temp,z = z,sig_p = sig_p)*0.2-0.2
      
      response_matrix<-exp(N%*%t(BB1)+A+rep(C,each=patches))
      Nt1<-N*response_matrix
      
      #change in trait z
      z_change<-z_up<-exp(N%*%t(BB1)+Env_perform(env = Temp,z+0.01,max(C),sig_p)*0.2-0.2+rep(C,each=patches))-
        response_matrix
      
      z_down<-exp(N%*%t(BB1)+Env_perform(env = Temp,z-0.01,max(C),sig_p)*0.2-0.2+rep(C,each=patches))-
        response_matrix
      
      z_up_down<-(z_up>0)*1
      z_up_down[z_down>0]<-((z_down>0)*-1)[z_down>0]
      
      z_change[z_down>z_up]<-z_down[z_down>z_up]
      z_change[z_change<0]<-0
      
      zt<-z+(rep(V_species,each=patches)*z_change*z_up_down)#+rnorm(species*patches,mean=0,sd=0.0001)
      
      
      #dispersal
      Nt<-Nt1-Nt1*rep(dispV,each=patches)+rep(dispV,each=patches)*disp_matrix%*%Nt1 #dispersal
      
      #gene flow
      geneflow_hold<-disp_matrix%*%(zt*Nt1*rep(dispV,each=patches))

      zt1<-(geneflow_hold+zt*Nt1*(1-rep(dispV,each=patches)))/Nt #gene flow
      zt1[is.na(zt1)]<-zt[is.na(zt1)]
      zt<-zt1
      
      Nt[Nt<10^-2]<-0
      Nstore[,,i]<-N
      N<-Nt
      z<-zt
    }
    matplot(t(Nstore[2,,]), type='l', col=as.numeric(as.factor(trophV)))
    } else{
    for(i in 1:Tmax){
      Temp<-Temp_initial+Temp_changeV[i]
      A<-Env_perform(env = Temp,z = z,sig_p = sig_p)*0.2-0.2
      
      response_matrix<-exp(N%*%t(BB1)+A+rep(C,each=patches))
      Nt1<-N*response_matrix
      
      #change in trait z
      z_change<-z_up<-exp(N%*%t(BB1)+Env_perform(env = Temp,z+0.01,max(C),sig_p)*0.2-0.2+rep(C,each=patches))-
        response_matrix
      
      z_down<-exp(N%*%t(BB1)+Env_perform(env = Temp,z-0.01,max(C),sig_p)*0.2-0.2+rep(C,each=patches))-
        response_matrix
      
      z_up_down<-(z_up>0)*1
      z_up_down[z_down>0]<-((z_down>0)*-1)[z_down>0]
      
      z_change[z_down>z_up]<-z_down[z_down>z_up]
      z_change[z_change<0]<-0
      
      zt<-z+(rep(V_species,each=patches)*z_change*z_up_down)#+rnorm(species*patches,mean=0,sd=0.0001)
      
      
      #calculate dispersal array
      disp_array1<-lapply((1:species)[trophV!="p"],function(j) {
        hold<-(response_matrix[,j]-min(response_matrix[,j]))/max(response_matrix[,j]-min(response_matrix[,j]))*disp_array_initial[,,j]
        hold1<-hold/colSums(hold)
      })
      disp_array[,,trophV!="p"]<-array(unlist(disp_array1),dim = c(patches,patches,species-nplants))
      
      #dispersal
      disp_hold<-do.call(cbind,lapply(1:species,function(j) disp_array[,,j]%*%Nt1[,j]))
      Nt<-Nt1-Nt1*rep(dispV,each=patches)+rep(dispV,each=patches)*disp_hold #dispersal
      
      #gene flow
      geneflow_hold<-do.call(cbind,lapply(1:species,function(j) disp_array[,,j]%*%(zt*Nt1*rep(dispV,each=patches))[,j]))
      
      zt1<-(geneflow_hold+zt*Nt1*(1-rep(dispV,each=patches)))/Nt #gene flow
      zt1[is.na(zt1)]<-zt[is.na(zt1)]
      zt<-zt1
      
      Nt[Nt<10^-2]<-0
      Nstore[,,i]<-N
      N<-Nt
      z<-zt
    }
  }
}

P_source<-sum(round((exp(rep(C,each=patches)+N%*%t(BB1)+A))[,trophV=="p"][N[,trophV=="p"]>0],digits = 2)>=1)/sum(N[,trophV=="p"]>0)
P_sink<-sum(round((exp(rep(C,each=patches)+N%*%t(BB1)+A))[,trophV=="p"][N[,trophV=="p"]>0],digits = 2)<1)/sum(N[,trophV=="p"]>0)

H_source<-sum(round((exp(rep(C,each=patches)+N%*%t(BB1)+A))[,trophV=="h"][N[,trophV=="h"]>0],digits = 2)>=1)/sum(N[,trophV=="h"]>0)
H_sink<-sum(round((exp(rep(C,each=patches)+N%*%t(BB1)+A))[,trophV=="h"][N[,trophV=="h"]>0],digits = 2)<1)/sum(N[,trophV=="h"]>0)

C_source<-sum(round((exp(rep(C,each=patches)+N%*%t(BB1)+A))[,trophV=="c"][N[,trophV=="c"]>0],digits = 2)>=1)/sum(N[,trophV=="c"]>0)
C_sink<-sum(round((exp(rep(C,each=patches)+N%*%t(BB1)+A))[,trophV=="c"][N[,trophV=="c"]>0],digits = 2)<1)/sum(N[,trophV=="c"]>0)

Source.df<-data.frame(Per_source=c(P_source,H_source,C_source,P_sink,H_sink,C_sink),T_level=c("Plants","Herbivores","Carnivores"),Source_sink=rep(c("Source","Sink"),each=3))
Source.df$T_level<-factor(Source.df$T_level,levels =c("Plants","Herbivores","Carnivores"),ordered=T)


#make food web networks####
Ints<-BB1
Ints[1:nplants,1:nplants]<-0
Ints[upper.tri(Ints)]<-0
diag(Ints)<-0
cut_value<-0
nets<-apply(N,1,function(x){
  Int_strength<-abs(Ints*rep(x,each=length(x)))
  Int_strength[x==0,]<-0
  Int_strength_cut<-quantile(Int_strength[Int_strength>0],cut_value)#mean(Int_strength[Int_strength>0])
  Int_strength[Int_strength<Int_strength_cut]<-0
  Ints2<-1*Int_strength>0
  hold.df<-t(data.frame(Ints2[x>0,x>0]))
  if(sum(x>0)<2){
    hold.df<-matrix(FALSE,1,1,dimnames = list(colnames(Int_strength)[x>0],colnames(Int_strength)[x>0]))
  }
  if(sum(x>0)!=0){
    net1<-graph.adjacency(hold.df)
    return(net1) 
  } 
})


return(list(Source.df,nets,N))
}