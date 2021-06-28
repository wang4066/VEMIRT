##################################################
#####     Function: GVEM_M2PL_CFA            #####
##################################################
# Inputs:
# u: response matrix, a N by J numerical matrix
# domain: K, numerical value of number of latent dimension
# indic: factor loading indicator matrix, a J by K numerical matrix, 
#       reflecting each item loads on which domain
#################################################
##### Outputs a list of updated parameters                                                             
##### ra: item discrimination parameters, a J by K matrix                                                           
##### rb: item difficulty parameters, vector of length J
##### reta: variational parameters(\eta(\xi) in the paper), a N by J matrix  
##### reps: variational parameters(\xi in the paper), a N by J matrix                          
##### rsigma: covariance matrix for domains, a K by K matrix
##### mu_i: mean parameter for each person, a K by N matrix
##### sig_i: covariance matrix for each person, a K by K by N array
##### n: the number of iterations
##### Q_mat: Q-matrix to indicate the loading structure,  a J by K matrix
##### GIC,AIC,BIC : model fit index
##################################################
library(Rcpp)
library(Matrix)
library(psych)
library(gtools)
#e step
e_cfa2pl<-'
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
List ecfa2(const arma::mat&u,const arma::mat& Sigma, const int& domain, const int& person, 
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
if(NumericVector::is_na(eta(i,j))){
continue;}
sigma_part=sigma_part+eta(i,j)*trans(new_a.row(j))*new_a.row(j);
mu_part=mu_part+trans((2*eta(i,j)*new_b(j)+u(i,j)-0.5)*new_a.row(j));
}
arma::mat sigmahat=solve((solve(Sigma,eye(domain,domain))+2*sigma_part),eye(domain,domain));
arma::vec muhat=sigmahat*mu_part;
SIGMA.slice(i)=sigmahat;
MU.col(i)=muhat;
arma::mat apro=new_a*(sigmahat+muhat*trans(muhat))*trans(new_a);
xi.row(i)=trans(sqrt(square(new_b)-2*new_b%(new_a*muhat)+apro.diag()));
Spart=Spart+sigmahat+muhat*trans(muhat);
}
return List::create(Named("Spart") = Spart,Named("SIGMA") = SIGMA,
Named("MU") = MU,Named("xi") = xi);
}
'
sourceCpp(code=e_cfa2pl)

#update a
a_cfa2pl<-'
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
arma::mat acfa2(const arma::mat&u,const arma::mat&indic,const int& person, 
const int& item, const int& domain, const arma::mat& eta,
const arma::vec& new_b,const arma::cube& SIGMA, const arma::mat& MU,const arma::mat& new_a1){
arma::mat new_a=new_a1;
for(int j=0; j<item; ++j){
arma::rowvec Ind=indic.row(j);
arma::uvec iind=find(Ind==1);
int di=std::count(Ind.begin(),Ind.end(), 1);
arma::vec a_nu= zeros(di);
arma::mat a_de= zeros(di,di);
for(int i=0; i<person; ++i){
arma::mat sigma=SIGMA.slice(i);
sigma=sigma.submat(iind,iind);
arma::vec mu=MU.col(i);
mu=mu.elem(iind);
if(NumericVector::is_na(eta(i,j))){
continue;}
a_de=a_de+eta(i,j)*sigma+eta(i,j)*(mu*trans(mu));
a_nu=a_nu+(u(i,j)-0.5+2*new_b(j)*eta(i,j))*mu;
}
arma::uvec id(1);
id.at(0)=j;
//new_a.submat(id,iind)=trans(inv(a_de)*a_nu/2);
new_a.submat(id,iind)=trans(solve(a_de,eye(di,di))*a_nu/2);
}
return new_a;
}
'
sourceCpp(code=a_cfa2pl)

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
if(NumericVector::is_na(eta(i,j))){
continue;}
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
vem_2PLCFA <- function(u,domain, indic,updateProgress=NULL) {
  person=dim(u)[1]
  item=dim(u)[2]
  #initialization
  initial=init(u,domain,indic)
  new_a = initial[[1]] * indic
  new_b = initial[[2]]
  eta= initial[[3]]
  xi=initial[[4]]
  Sigma=initial[[5]]
  converge = 0
  MU=matrix(0,nrow=domain,ncol=person)
  SIGMA=array(0,dim=c(domain,domain,person))
  n = 0
  par_Sigma = Sigma
  while(converge==0 && rankMatrix(Sigma) == domain && n < 5500){
    if (is.function(updateProgress)) {
      text <- paste0("n=:", n)
      updateProgress(detail = text)
    }
    par_Sigma = Sigma
    #update MU, SIGMA, sigma, eta
    rs1<-ecfa2(u,Sigma,domain, person, item, eta, new_a, new_b)
    xi=rs1$xi
    xi[is.na(u)]=NA
    Spart=rs1$Spart
    SIGMA=rs1$SIGMA
    MU=rs1$MU
    eta=(exp(xi)/(1+exp(xi))-0.5)/(2*xi)
    eta=replace(eta,abs(eta)< 0.01,0.125)
    Sigma=Spart/person
    d_temp=sqrt(diag(diag(Sigma)))
    Sigma=solve(d_temp)%*%Sigma%*%solve(d_temp)
    
    #update b
    par_b=new_b
    b_part = t(new_a%*%MU)
    new_b=colSums(0.5-u+2*eta*b_part,na.rm = T)/colSums(2*eta,na.rm = T)
    par_a=new_a
    #update a
    new_a=acfa2(u, indic, person, item, domain, eta, new_b, SIGMA, MU,new_a)
    if (norm(as.vector(new_a)-as.vector(par_a),type="2")+norm(new_b-par_b,type="2")+
        norm(as.vector(Sigma)-as.vector(par_Sigma),type="2") <0.0001){
      converge=1
    }
    n=n+1
  }
  if (rankMatrix(Sigma) < domain){
    rsigma = par_Sigma
  }else{
    rsigma = Sigma
  }
  #new_a=replace(new_a,new_a< 0,0)
  new_a=replace(new_a,abs(new_a)< 0.001,0)
  Q_mat = (new_a != 0)*1
  #gic
  lbound=lb2pl(u,xi,Sigma,new_a,new_b,SIGMA,MU)
  gic=log(log(person))*log(person)*sum(Q_mat+item) - 2*lbound
  #AIC, BIC
  bic = log(person)*sum(Q_mat+item) - 2*lbound
  aic = 2*sum(Q_mat+item) -2*lbound
  return(list(ra=new_a,rb=new_b,reta = eta,reps=xi,rsigma = Sigma,
              mu_i = MU,sig_i = SIGMA,n=n,Q_mat=Q_mat,GIC=gic,AIC=aic,
              BIC=bic))
}