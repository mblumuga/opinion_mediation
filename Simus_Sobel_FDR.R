## Install REMP
#library(devtools)
#require(HIMA)
nsimu<-1
require(multilevel)
require(ggplot2)
nb_indiv<-200
mysd<-2
for (i in 1:nsimu)
{
  print(i)
  
  pmed<-.5
  nb_med<-5000
  EtoM<-500
  MtoS<-500
  env<-rnorm(nb_indiv)
  menv<-matrix(env,nrow=nb_indiv,ncol=nb_med,byrow=F)
  intercept<-rnorm(nb_med)
  mintercept<-matrix(intercept,nrow=nb_indiv,ncol=nb_med,byrow=T)
  effect<-rep(0,nb_med)
  ss<-sample(1:nb_med,size=EtoM,replace=F)
  effect[ss]<-rnorm(EtoM)
  meffect<-matrix(effect,nrow=nb_indiv,ncol=nb_med,byrow=T)
  residuals<-rnorm(nb_med*nb_indiv,sd=mysd)
  mres<-matrix(residuals,nrow=nb_indiv,ncol=nb_med,byrow=T)
  
  epiG<-mintercept+menv*meffect+mres
  
  effect2<-rnorm(MtoS)
  #ss2<-c(sample(ss,size=MtoS/2,replace=F),sample((1:nb_med)[-ss],size=MtoS/2,replace=F))
  ss2<-c(sample(ss,size=MtoS,replace=FALSE))
  
  truess<-intersect(ss,ss2)
  
  meffect2<-matrix(rnorm(MtoS),nrow=nb_indiv,ncol=MtoS,byrow=T)
  
  alpha0<-sum(effect2*effect[ss2])*(1-pmed)/pmed
  Y<-rnorm(1)+rowSums(meffect2*epiG[,ss2])+alpha0*env+rnorm(nb_indiv,sd=mysd)
  z_scores <- sapply(1:nb_med,FUN=function(x){sobel(env,epiG[,x],Y)$z.value})
  my_sobel <- 2*pnorm(abs(z_scores),lower.tail=F)
  
  ff<-fdrtool::fdrtool(z_scores)
  dat<-data.frame(pval=c(my_sobel,ff$pval),corr=c(rep("P-values",nb_med),rep("P-values after adjustment using fdrtools",nb_med)))
  # Change the legend position to "top" 
  # (possible values: "left","top", "right", "bottom")
  
  
  ggplot(dat, aes(x=pval, fill=corr)) +
    geom_histogram(position="identity",binwidth=0.02,alpha=0.5)+
    theme_grey(base_size = 15)+
    theme(legend.position="top")+
    theme(legend.title = element_blank())+
  scale_x_continuous(name="P-values") +
    scale_y_continuous(name="Counts")
  ggsave(filename="2hists.jpg")
  
  
  +#+
    #theme(legend.text=element_text(size=10))+
    #theme(axis.text=element_text(size=30)+
   

    #scale_fill_discrete(breaks=c("correction","no correction"),
            #              labels=c("Control", "Treatment 1", "Treatment 2"))
}
