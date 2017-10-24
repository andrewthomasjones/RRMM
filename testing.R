normal_densities<-function(fit, k, data){
  y<-array(NA, c(k,fit$g+1))
  xseq<-seq(min(data), max(data), length.out=k)
  y[,1]<-xseq
  
  for(i in 2:(fit$g+1)){
    y[,i]<-fit$pi[i-1]*dnorm(xseq, mean=fit$mu[i-1], sd=sqrt((fit$sigma[i-1])))
    
  }
  colnames(y)<-c("x", paste0('y', seq(1:fit$g)))
  return(y)
}

plot_normal_densities<-function(fit, data, k, bins,name){
  data1<-as.data.frame(normal_densities(fit,k, data))
  data2<-data.frame(data=data)
  p=ggplot(data=data1)+geom_histogram(data=data2, aes(y=..density..,x=data), fill='white', colour='grey', bins=bins)+geom_line(aes(x=x,  y=y1))+geom_line(aes(x=x,  y=y2))+geom_line(data=data1, aes(x=x,  y=y1+y2))
  p=p+scale_colour_discrete(guide=FALSE)+theme_bw()+xlab("z-scores")+ylab("Density")+ggtitle(name)
  print(p)
}

plot_normal_densities_3<-function(fit, data, k, bins,name){
  data1<-as.data.frame(normal_densities(fit,k, data))
  data2<-data.frame(data=data)
  p=ggplot(data=data1)+geom_histogram(data=data2, aes(y=..density..,x=data), fill='white', colour='grey', bins=bins)+geom_line(aes(x=x,  y=y1))+geom_line(aes(x=x,  y=y2))+geom_line(aes(x=x,  y=y3))+geom_line(data=data1, aes(x=x,  y=y1+y2+y3))
  p=p+scale_colour_discrete(guide=FALSE)+theme_bw()+xlab("z-scores")+ylab("Density")+ggtitle(name)
  print(p)
}


normal_densities_a<-function(fit, k){
  fit$G<-length(fit$lambda)
  y<-array(NA, c(k,fit$G+1))
  xseq<-seq(min(fit$x), max(fit$x), length.out=k)
  y[,1]<-xseq
  for(i in 2:(fit$G+1)){
    y[,i]<-fit$lambda[i-1]*dnorm(xseq, mean=fit$mu[i-1], sd=fit$sigma[i-1])
    
  }
  
  colnames(y)<-c("x", paste0('y', seq(1:fit$G)))
  return(y)
}

plot_normal_densities_a<-function(fit, data, k, bins,name){
  data1<-as.data.frame(normal_densities_a(fit,k))
  data2<-data.frame(data=data)
  p=ggplot(data=data1)+geom_histogram(data=data2, aes(y=..density..,x=data), fill='white', colour='grey', bins=bins)+geom_line(aes(x=x,  y=y1))+geom_line(aes(x=x,  y=y2))+geom_line(data=data1, aes(x=x,  y=y1+y2))
  p=p+scale_colour_discrete(guide=FALSE)+theme_bw()+xlab("z-scores")+ylab("Density")+ggtitle(name)
  print(p)
}


hiv1<-MMEst(hivdata, 3, 0.00001, 1000)
plot_normal_densities_3(hiv1, hivdata, 1000, 50, "Robust Norm HIV")

hed1<-MMEst(hedenfalk$z, 2, 0.00001, 2000)
plot_normal_densities(hed1, hedenfalk$z, 1000, 50, "Robust Norm Hedenfalk")

col1<-MMEst(Colon$z, 2, 0.00001, 2000)
plot_normal_densities(col1, Colon$z, 1000, 50, "Robust Norm Colon")


hiv_fit_norm <- normalmixEM(x=hivdata, lambda =pi_0(hivdata,3), k=2, maxrestarts=100)
hed_fit_norm <- normalmixEM(x=hedenfalk$z,k=2, lambda = pi_0(hedenfalk$z,0.2),  maxrestarts=100) 
col_fit_norm <- normalmixEM(x=Colon$z, lambda =pi_0(Colon$z,0.2), k=2, maxrestarts=100)
plot_normal_densities_a(hiv_fit_norm,hivdata, 1000, 50, "Hiv Emp. Normal")
plot_normal_densities_a(col_fit_norm,Colon$z, 1000, 50, "Colon Normal")
plot_normal_densities_a(hed_fit_norm,hedenfalk$z, 1000, 50, "Hedenfalk Normal")

