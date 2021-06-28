##################################################
###Function: Initialization Parameters       #####
##################################################
# Inputs:
# u: response matrix, a N by J numerical matrix
# domain: K, numerical value of number of latent dimension
# indic: factor loading indicator matrix, a J by K numerical matrix, 
#       reflecting each item loads on which domain
#################################################
##### Outputs a list of initialized parameters                                                             
##### a0: item discrimination parameters, a J by K matrix                                                           
##### b0: item difficulty parameters, vector of length J
##### eta0: variational parameters(\etla(\xi) in the paper), a N by J matrix  
##### eps0: variational parameters(\xi in the paper), a N by J matrix                          
##### Sigma: covariance matrix for domains, a K by K matrix
##################################################
library(Rcpp)
#some rcpp functions can handle missing data
src<-'#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector ide(NumericMatrix u, NumericVector scale, NumericVector p, NumericVector q, NumericVector y,
double s){
int item=u.ncol();
NumericVector r(item);
for (int i=0; i<item; ++i){
NumericVector u1=na_omit(u(_,i));
NumericVector scale1=scale[!is_na(u(_,i))];
NumericVector x1= scale1[u1==1];
NumericVector x2= scale1[u1==0];
double md=(mean(x1)-mean(x2));
if(Rcpp::traits::is_nan<REALSXP>(md)){
double md=0;
//}
//if(Rf_isNull(x2)){
//Rcout << "x2 is NULL." << std::endl;
//x2.fill(0);
}
r[i]=md/s*p[i]*q[i]/y[i];
}
return r;
}
'
sourceCpp(code=src)
#some functions
#identification: u=response
identify<-function(u){
  scale=rowSums(u,na.rm = T) #num of item responsed correctly per examinee
  p=apply(u,2,mean,na.rm=T) #the frequency per item
  p=replace(p,p==0,0.001)
  p=replace(p,p==1,0.999)
  q=1-p
  y=qnorm(p,0,1) #inverse of the standard normal CDF
  y=dnorm(y,0,1) #the density function
  s=sd(scale)
  r=NULL
  #r<-ide(u,scale,p,q,y,s)
  for (i in 1:dim(u)[2]) {
    u1=u[!is.na(u[,i]),i]
    scale1=scale[!is.na(u[,i])]
    x1=scale1[u1==1]
    x2=scale1[u1==0]
    if(identical(x1,numeric(0))){
      x1=0
    }
    if(identical(x2,numeric(0))){
      x2=0
    }
    r[i]=(mean(x1)-mean(x2))/s*p[i]*q[i]/y[i]
  }
  return(r)
}
#initialization
init<-function(u,domain,indic){
  r=identify(u)
  person=dim(u)[1]
  r[r>0.9]=0.9
  r[r<0]=abs(r[r<0][1])
  r[r==0]=0.0001
  a0=t(rep(1,domain)%o%(r/sqrt(1-r^2)))*indic
  a0=replace(a0,a0>4,4)
  b0=-qnorm(colSums(u,na.rm=T)/person,0,1)/r
  b0[b0>4]=4
  b0[b0<(-4)]=-4
  Sigma = diag(domain)
  theta=matrix(rnorm(person*domain,0,1),nrow=person)
  #person*item
  xi=array(1,person)%*%t(b0)-theta%*%t(a0)
  eta0=(exp(xi)/(1+exp(xi))-0.5)/(2*xi)
  eta0[is.na(u)]=NA
  eps0=xi
  return(list(a0,b0,eta0,eps0,Sigma))
}