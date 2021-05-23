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
src<-"#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector ide(NumericMatrix u, NumericVector scale, NumericVector p, NumericVector q, NumericVector y,
double s){
int item=u.ncol();
NumericVector r(item);
for (int i=0; i<item; ++i){
NumericVector x1= scale[u(_,i)==1];
NumericVector x2= scale[u(_,i)==0];
r[i]=(mean(x1)-mean(x2))/s*p[i]*q[i]/y[i];
}
return r;
}
"
sourceCpp(code=src)

#identification: make sure the initialized values are identified
identify<-function(u){
  scale=rowSums(u) #num of item responsed correctly per examinee
  p=apply(u,2,mean) #the frequency per item
  q=1-p
  y=qnorm(p,0,1) #inverse of the standard normal CDF
  y=dnorm(y,0,1) #the density function
  s=sd(scale)
  r<-ide(u,scale,p,q,y,s)
  return(r)
}
#initialization
init<-function(u,domain,indic){
  r=identify(u)
  person=dim(u)[1]
  r[r>0.9]=0.9
  r[r<0]=abs(r[r<0][1])
  a0=t(rep(1,domain)%o%(r/sqrt(1-r^2)))*indic
  b0=-qnorm(colSums(u)/person,0,1)/r
  Sigma = diag(domain)
  theta=matrix(rnorm(person*domain,0,1),nrow=person)
  #person*item
  xi=array(1,person)%*%t(b0)-theta%*%t(a0)
  eta0=(exp(xi)/(1+exp(xi))-0.5)/(2*xi)    
  eps0=xi
  return(list(a0,b0,eta0,eps0,Sigma))
}