##################################################
#####     Function: GVEM_M2PL_Rotation       #####
##################################################
# Inputs:
# u: response matrix, a N by J numerical matrix
# domain: K, numerical value of number of latent dimension
#################################################
##### Outputs a list of initialized parameters                                                             
##### ra: item discrimination parameters, a J by K matrix                                                           
##### rb: item difficulty parameters, vector of length J
##### reta: variational parameters(\eta(\xi) in the paper), a N by J matrix  
##### reps: variational parameters(\xi in the paper), a N by J matrix                          
##### rsigma: covariance matrix for domains, a K by K matrix
##### mu_i: mean parameter for each person, a K by N matrix
##### sig_i: covariance matrix for each person, a K by K by N array
##### n: the number of iterations
##### rk: factor loadings after the transformation, a J by K matrix
##### Q_mat: Q-matrix to indicate the loading structure,  a J by K matrix
##### GIC: numerical value of GIC
##################################################

library(Rcpp)
library(Matrix)
library(psych)
library(gtools)
#e step
e_efa2pl<-'
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
List eefa2(const arma::mat&u, const int& domain, const int& person, 
const int& item, const arma::mat& eta, const arma::mat& new_a,const arma::vec& new_b){
arma::mat MU=mat(domain,person,fill::zeros);
arma::cube SIGMA=zeros(domain,domain,person);
arma::mat Spart=mat(domain,domain,fill::zeros);
arma::mat xi=mat(person,item,fill::zeros);
for(int i=0; i<person; ++i){
arma::mat sigma_part=arma::mat(domain,domain,arma::fill::zeros);
//mu_part.zeros();
arma::vec mu_part= arma::zeros(domain);
for(int j=0; j<item; ++j){
sigma_part=sigma_part+eta(i,j)*trans(new_a.row(j))*new_a.row(j);
mu_part=mu_part+trans((2*eta(i,j)*new_b(j)+u(i,j)-0.5)*new_a.row(j));
}
arma::mat Sigma= eye(domain,domain);
arma::mat sigmahat=solve((Sigma+2*sigma_part),eye(domain,domain));
arma::vec muhat=sigmahat*mu_part;
SIGMA.slice(i)=sigmahat;
MU.col(i)=muhat;
arma::mat apro=new_a*(sigmahat+muhat*trans(muhat))*trans(new_a);
xi.row(i)=trans(sqrt(square(new_b)-2*new_b%(new_a*muhat)+apro.diag()));
//Spart=Spart+sigmahat+muhat*trans(muhat);
}
return List::create(Named("SIGMA") = SIGMA,
Named("MU") = MU,Named("xi") = xi);
}
'
sourceCpp(code=e_efa2pl)

#update a
a_efa2pl<-'
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
arma::mat aefa2(const arma::mat&u,const int& domain,const int& person, const int& item, const arma::mat& eta,
const arma::vec& new_b,const arma::cube& SIGMA, const arma::mat& MU){
arma::mat new_a=mat(item,domain,fill::zeros);
for(int j=0; j<item; ++j){
arma::vec a_nu= zeros(domain);
arma::mat a_de= zeros(domain,domain);
for(int i=0; i<person; ++i){
arma::mat sigma=SIGMA.slice(i);
arma::vec mu=MU.col(i);
a_de=a_de+eta(i,j)*sigma+eta(i,j)*(mu*trans(mu));
a_nu=a_nu+(u(i,j)-0.5+2*new_b(j)*eta(i,j))*mu;
}
new_a.row(j)=trans(solve(a_de,eye(domain,domain))*a_nu/2);
}
return new_a;
}
'
sourceCpp(code=a_efa2pl)

#change the sign of a and sigma
rt<-function(A,Sig){
  domain=dim(A)[2]
  #change the sign
  sign_reversed = (-1)^(colSums(A)<0)
  A=A%*%diag(sign_reversed)
  sign_mat=sign_reversed%*%t(sign_reversed)
  Sig=sign_mat*Sig
  return(list(ra=A,rsigma=Sig))
}

#lower-bound
lb_2pl<-'
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
double lb2pl(const arma::mat&u,const arma::mat& xi, 
const arma::mat& sig,const arma::mat& new_a,
const arma::vec& new_b,const arma::cube& SIGMA, const arma::mat& MU){
int person=u.n_rows;
int item=u.n_cols;
int domain=sig.n_cols;
arma::mat eta =(exp(xi)/(1+exp(xi))-0.5)/(2*xi); 
double sum = 0; 
double sum2 = 0;
for(int i=0; i<person; ++i){
for(int j=0; j<item; ++j){
sum = sum + log(exp(xi(i,j))/(1+exp(xi(i,j))))+(0.5-u(i,j))*new_b(j)+(u(i,j)-0.5)*as_scalar(new_a.row(j)*MU.col(i));
sum = sum - 0.5*xi(i,j);
double apro=as_scalar(new_a.row(j)*(SIGMA.slice(i)+MU.col(i)*trans(MU.col(i)))*trans(new_a.row(j)));
sum = sum - eta(i,j)*(new_b(j)*new_b(j)-2*new_b(j)*as_scalar(new_a.row(j)*MU.col(i))+apro - xi(i,j)*xi(i,j));
}
sum2=sum2+trace(solve(sig,eye(domain,domain))*(SIGMA.slice(i)+MU.col(i)*trans(MU.col(i))));
}
double lb=sum - 0.5*sum2 + person/2*log(det(solve(sig,eye(domain,domain))));
return lb;
}
'
sourceCpp(code=lb_2pl)

#the main function
vem_2PLEFA <- function(u, domain,updateProgress=NULL) {
  person=dim(u)[1]
  item=dim(u)[2]
  converge = 1
  
  #initialization
  r=identify(u)
  person=dim(u)[1]
  r[r>0.9]=0.9
  r[r<0]=abs(r[r<0][1])
  new_a=t(rep(1,domain)%o%(r/sqrt(1-r^2)))
  new_b=-qnorm(colSums(u)/person,0,1)/r
  Sigma = diag(domain)
  theta=matrix(rnorm(person*domain,0,1),nrow=person)
  #person*item
  xi=array(1,person)%*%t(new_b)-theta%*%t(new_a)
  eta=(exp(xi)/(1+exp(xi))-0.5)/(2*xi)
  
  
  # initial=init(u,domain,indic)
  MU=matrix(0,nrow=domain,ncol=person)
  SIGMA=array(0,dim=c(domain,domain,person))
  n = 0
  # new_a = initial[[1]] * indic
  # new_b = initial[[2]]
  # eta= initial[[3]]
  # xi=initial[[4]]
  # Sigma=initial[[5]]
  # par_Sigma = Sigma
  while(converge==1){
    if (is.function(updateProgress)) {
      text <- paste0("n=:", n)
      updateProgress(detail = text)
    }
    #update mu and sigma for each person
    rs1<-eefa2(u,domain, person, item, eta, new_a, new_b)
    xi=rs1$xi
    SIGMA=rs1$SIGMA
    MU=rs1$MU
    eta=(exp(xi)/(1+exp(xi))-0.5)/(2*xi)
    eta=replace(eta,abs(xi)< 0.01,0.125)
    #Sigma=Spart/person
    #d_temp=sqrt(diag(diag(Sigma)))
    #Sigma=solve(d_temp)%*%Sigma%*%solve(d_temp)
    
    #update b
    par_b=new_b
    b_part = t(new_a%*%MU)
    new_b=colSums(0.5-u+2*eta*b_part)/colSums(2*eta)
    par_a=new_a
    #update a
    new_a=aefa2(u,domain,person,item,eta,new_b,SIGMA,MU)
    if (norm(as.vector(new_a)-as.vector(par_a),type="2")+norm(new_b-par_b,type="2") <0.0001){
      converge=0
    }
    n=n+1
  }
  #rotate factors: promax
  rotated=Promax(new_a)
  new_a=rotated$loadings[1:item,]
  #cut-off points:0.2
  #new_a=replace(new_a,abs(new_a)< 0.2,0)
  Sigma=rotated$Phi
  #change the sign
  rt1<-rt(new_a,Sigma)
  new_a=rt1$ra
  Sigma=rt1$rsigma
  #calculate the factor loadings
  rk<-matrix(NA,item,domain)
  for (j in 1:item) {
    ra1=new_a[j,]
    B=diag(domain)+ra1%*%t(ra1)*Sigma
    rk[j,]=t(solve(chol(B))%*%ra1)
  }
  Q_mat<-(abs(rk)>0.3)*1
  #gic
  lbound=lb2pl(u,xi,Sigma,new_a,new_b,SIGMA,MU)
  gic=log(log(person))*log(person)*sum(Q_mat+item) - 2*lbound
  return(list(ra=new_a,rb=new_b,reta = eta,reps=xi,rsigma = Sigma,mu_i = MU,sig_i = SIGMA,n=n,rk=rk,Q_mat=Q_mat,GIC=gic))
}