##################################################
#####     Function: GVEM_M3PL_Rotation       #####
##################################################
# Inputs:
# u: response matrix, a N by J numerical matrix
# domain: K, numerical value of number of latent dimension
# samp: numerical value of subsample for each iteration
# forgetrate: forget rate for the stochastic algorithm
# mu_b, sigma2_b: the mean and variance parameters, prior distribution of b parameters
# Alpha, Beta: the alpha and beta parameters, prior distribution of g parameters
#################################################
##### Outputs a list of initialized parameters                                                             
##### ra: item discrimination parameters, a J by K matrix                                                           
##### rb: item difficulty parameters, vector of length J
##### rc: item guessing parameters, vector of length J
##### rs: variational parameters(s in the paper), a N by J matrix  
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
library(testit)
#e step
e_efa3pl<-'
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
List eefa3(const arma::mat&u,const int& domain, const arma::vec& id_1,
const int& item, const arma::mat& eta, const arma::mat& new_s,const arma::mat& new_a,const arma::vec& new_b,
const arma::cube& SIGMA1,const arma::mat& MU1,const arma::mat& xi1){
arma::mat MU=MU1;
arma::cube SIGMA=SIGMA1;
arma::mat xi=xi1;
for(int i=0; i<id_1.n_elem; ++i){
int k=id_1(i)-1;
arma::mat sigma_part=arma::mat(domain,domain,arma::fill::zeros);
//mu_part.zeros();
arma::vec mu_part= arma::zeros(domain);
for(int j=0; j<item; ++j){
sigma_part=sigma_part+eta(k,j)*(1-u(k,j)+new_s(k,j)*u(k,j))*trans(new_a.row(j))*new_a.row(j);
mu_part=mu_part+trans((2*eta(k,j)*new_b(j)+u(k,j)-0.5)*(1-u(k,j)+new_s(k,j)*u(k,j))*new_a.row(j));
}
arma::mat Sigma= eye(domain,domain);
arma::mat sigmahat=solve((Sigma+2*sigma_part),eye(domain,domain));
arma::vec muhat=sigmahat*mu_part;
SIGMA.slice(k)=sigmahat;
MU.col(k)=muhat;
arma::mat mukk=sigmahat+muhat*trans(muhat);
arma::mat muk=new_a*mukk*trans(new_a);
xi.row(k)=trans(sqrt(square(new_b)-2*new_b%(new_a*muhat)+muk.diag()));
}
return List::create(Named("SIGMA") = SIGMA,
Named("MU") = MU,Named("xi") = xi);
}
'
sourceCpp(code=e_efa3pl)

#update s
s_3pl<-'
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
arma::mat s3pl(const int& item, const arma::vec& id_1,const arma::mat& s,const arma::mat& eta,
const arma::mat& new_a,const arma::vec& s_num,const arma::vec& new_b,const arma::mat&u,const arma::cube& SIGMA,
const arma::mat& MU,const arma::mat& xi){
arma::mat new_s=s;
for(int i=0; i<id_1.n_elem; ++i){
int k=id_1(i)-1;
for(int j=0; j<item; ++j){
double s_denom = 0;
if (u(k,j) == 0){
new_s(k,j) = 1;
}else{
s_denom = s_denom + log(exp(xi(k,j))/(1+exp(xi(k,j))))+(0.5-u(k,j))*new_b(j)+(u(k,j)-0.5)*as_scalar(new_a.row(j)*MU.col(k))-0.5*xi(k,j);
arma::mat mid=SIGMA.slice(k)+MU.col(k)*trans(MU.col(k));
double aa=as_scalar(new_a.row(j)*mid*trans(new_a.row(j)));
s_denom = s_denom - eta(k,j)*as_scalar(new_b(j)*new_b(j)-2*new_b(j)*new_a.row(j)*MU.col(k)+aa-xi(k,j)*xi(k,j));
s_denom = exp(-s_denom);
new_s(k,j) = s_num(j)/(s_num(j)+s_denom); 
}
}
}
return new_s;
}
'
sourceCpp(code=s_3pl)

#update eta
eta_3pl<-'
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
arma::mat eta3pl(const int& item, const arma::vec& id_1,const arma::mat& ea,const arma::mat& xi){
arma::mat eta=ea;
for(int i=0; i<id_1.n_elem; ++i){
int k=id_1(i)-1;
for(int j=0; j<item; ++j){
if(xi(k,j)<0.01){
eta(k,j)=0.125;
}
}
}
return eta;
}
'
sourceCpp(code=eta_3pl)

#update b
b_3pl<-'
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
List b3pl(const int& item, const arma::vec& id_1,const arma::mat& new_s,const arma::mat& eta,
const arma::mat& u, const arma::mat& new_a,const double& dec_st,const arma::mat& MU,
const arma::vec& prev_b_num,const arma::vec& prev_b_denom,const double& mu_b,const double& sigma2_b){
arma::vec b_denom= arma::zeros(item);
arma::vec b_num= arma::zeros(item);
for(int j=0; j<item; ++j){
for(int i=0; i<id_1.n_elem; ++i){
int k=id_1(i)-1;
b_denom(j) = b_denom(j) + 2*(1-u(k,j)+new_s(k,j)*u(k,j))*eta(k,j);
b_num(j) = b_num(j) + (1-u(k,j)+new_s(k,j)*u(k,j))*(0.5-u(k,j)+2*eta(k,j)*as_scalar(new_a.row(j)*MU.col(k)));
}
}
b_denom = b_denom + 1/sigma2_b;
b_num = b_num + mu_b/sigma2_b;
arma::vec new_b=(dec_st*b_num+(1-dec_st)*prev_b_num)/(dec_st*b_denom+(1-dec_st)*prev_b_denom);
return List::create(Named("new_b") = new_b,Named("b_denom") = b_denom,Named("b_num") = b_num);
}'
sourceCpp(code=b_3pl)

#update a
a_efa3pl<-'
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
List aefa3(const arma::mat&u,const arma::vec& id_1,const int& item,const int& domain, const arma::mat& eta, 
const arma::mat& a,const arma::vec& new_b,const arma::mat& new_s,const arma::cube& SIGMA, const arma::mat& MU,const arma::mat& prev_a_num,
const arma::cube& prev_a_denom, const double& dec_st){
arma::mat new_a=a;
arma::mat a_num=zeros(domain,item); 
arma::cube a_denom=zeros(domain,domain,item);
for(int j=0; j<item; ++j){
arma::mat a_denom_sub=a_denom.slice(j);
for(int k=0; k<id_1.n_elem; ++k){
int i=id_1(k)-1;
arma::mat sigma=SIGMA.slice(i);
arma::vec mu=MU.col(i);
double a1=1-u(i,j)+new_s(i,j)*u(i,j);
a_denom_sub=a_denom_sub+a1*eta(i,j)*(sigma+mu*trans(mu));
a_num.col(j)=a_num.col(j)+a1*(u(i,j)-0.5+2*new_b(j)*eta(i,j))*mu; 
}
if(rank(a_denom_sub) < domain){
break;
}
a_denom.slice(j)=a_denom_sub;
arma::mat prev_a_denom_sub=prev_a_denom.slice(j);
arma::mat a2=solve((dec_st*a_denom_sub+(1-dec_st)*prev_a_denom_sub),eye(domain,domain));
new_a.row(j)=trans(a2*(dec_st*a_num.col(j)+(1-dec_st)*prev_a_num.col(j))/2);
}
return List::create(Named("new_a") = new_a,Named("a_denom") = a_denom,Named("a_num") = a_num);
}
'
sourceCpp(code=a_efa3pl)

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
lb_3pl<-'
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
double lb3pl(const arma::mat&u,const arma::mat& xi, const arma::mat& s, const arma::vec& id_1,
const arma::mat& new_a,const arma::vec& new_c,const arma::mat& sig,
const arma::vec& new_b,const arma::cube& SIGMA, const arma::mat& MU,const double& Alpha,
const double& Beta,const double& mu_b,const double& sigma2_b){
int item=u.n_cols;
int domain=new_a.n_cols;
arma::mat eta =(exp(xi)/(1+exp(xi))-0.5)/(2*xi); 
double sum1 = 0;
double sum2 = 0;
double sum3 = 0;
for(int k=0; k<id_1.n_elem; ++k){
int i=id_1(k)-1;
for(int j=0; j<item; ++j){
double a1=1-u(i,j)+s(i,j)*u(i,j);
double a2=log(exp(xi(i,j))/(1+exp(xi(i,j))))+(0.5-u(i,j))*new_b(j);
double a3=(u(i,j)-0.5)*as_scalar(new_a.row(j)*MU.col(i));
double apro=as_scalar(new_a.row(j)*(SIGMA.slice(i)+MU.col(i)*trans(MU.col(i)))*trans(new_a.row(j)));
double a4=eta(i,j)*(new_b(j)*new_b(j)-2*new_b(j)*as_scalar(new_a.row(j)*MU.col(i))+apro - xi(i,j)*xi(i,j));
sum1 = sum1 + a1*(a2+a3-0.5*xi(i,j)-a4);
sum3=sum3 + a1*log(1-new_c(j)) + u(i,j)*(1-s(i,j))*log(new_c(j));
}
sum2=sum2+trace(solve(sig,eye(domain,domain))*(SIGMA.slice(i)+MU.col(i)*trans(MU.col(i))));
}
double sum_b = as_scalar(-1/(2*sigma2_b)*trans(new_b-mu_b)*(new_b-mu_b));
double sum_c = sum((Alpha-1)*log(new_c)+(Beta-1)*log(1-new_c));
double lb=sum1 - 0.5*sum2 + id_1.n_elem/2*log(det(eye(domain,domain)))+sum3+ sum_b + sum_c;
return lb;
}
'
sourceCpp(code=lb_3pl)

#the main function
sgvem_3PLEFA <- function(u,domain, samp,forgetrate,mu_b,sigma2_b,Alpha,Beta,updateProgress=NULL) {
  person=dim(u)[1]
  item=dim(u)[2]
  converge = 1
  aveg_p=colMeans(u==1)
  inite=person
  #initialization
  r=identify(u)
  person=dim(u)[1]
  r[r>0.9]=0.9
  r[r<0]=abs(r[r<0][1])
  new_a=t(rep(1,domain)%o%(r/sqrt(1-r^2)))
  new_b=-qnorm(colSums(u)/person,0,1)/r
  new_c=rep(0.1,item)
  theta=matrix(rnorm(person*domain,0,1),nrow=person)
  #person*item
  xi=array(1,person)%*%t(new_b)-theta%*%t(new_a)
  eta=(exp(xi)/(1+exp(xi))-0.5)/(2*xi)
  new_s=matrix(NA,nrow=person,ncol=item)
  for (i in 1:person) {
    for (j in 1:item) {
      new_s[i,j]<-ifelse(u[i,j]==0,1,runif(1,aveg_p[j],1))
    }
  }
  
  SIGMA=array(0,dim=c(domain,domain,person))
  for (i in 1:person) {
    SIGMA[,,i]=diag(domain)
  }
  MU=matrix(0,nrow=domain,ncol=person)
  prev_a_num = matrix(0,nrow=domain,ncol=item)
  prev_a_denom = array(0,dim=c(domain,domain,item))
  prev_b_num = matrix(0,nrow=item,ncol=1)
  prev_b_denom = matrix(0,nrow=item,ncol=1)
  prev_c_num = matrix(0,nrow=item,ncol=1)
  prev_c_denom = 0
  
  n = 1
  #lb=NULL
  prev_lb = 0
  nlb<-n20<-NULL
  while(converge==1 &&  n<6000){
    if (is.function(updateProgress)) {
      text <- paste0("n=:", n)
      updateProgress(detail = text)
    }
    if(n==1){
      dec_st=1
      id_1 = sample(person,inite)
    }else{
      #dec_st = (n+1)^(-forgetrate)*0.25
      dec_st = (n+1)^(-forgetrate)
      id_1 = sample(person,samp)
    }
    s_num = (1-new_c)/new_c
    #par_Sigma=Sigma
    #update MU and SIGMA,xi,Spart
    rs1<-eefa3(u, domain, id_1, item, eta, new_s, new_a, new_b, SIGMA, MU, xi)
    
    SIGMA=rs1$SIGMA
    MU=rs1$MU
    xi=rs1$xi
    
    #update s
    new_s=s3pl(item, id_1, new_s, eta, new_a, s_num, new_b, u, SIGMA, MU, xi)
    
    #update eta
    eta=(exp(xi)/(1+exp(xi))-0.5)/(2*xi)
    eta=eta3pl(item,id_1,eta,abs(xi))
    
    
    #update b
    par_b=new_b
    
    rs2=b3pl(item,id_1,new_s,eta,u,new_a,dec_st,MU,prev_b_num,prev_b_denom,mu_b,sigma2_b)
    new_b=rs2$new_b
    b_denom= rs2$b_denom
    b_num = rs2$b_num
    #update a
    par_a=new_a
    rs3=aefa3(u, id_1, item, domain, eta, new_a, new_b, new_s, SIGMA, MU, prev_a_num, prev_a_denom, dec_st)
    new_a=rs3$new_a
    a_denom= rs3$a_denom
    a_num = rs3$a_num
    #update c
    par_c = new_c;
    c_num = colSums(u[id_1,]*(1-new_s[id_1,]))+Alpha - 1
    c_denom = length(id_1) + Alpha + Beta - 2
    new_c = (dec_st*c_num+(1-dec_st)*prev_c_num)/(dec_st*c_denom+(1-dec_st)*prev_c_denom)
    
    #save old sums for SAEM updates
    prev_a_num = (1-dec_st)*prev_a_num + dec_st*a_num
    prev_a_denom = (1-dec_st)*prev_a_denom + dec_st*a_denom
    
    prev_b_num = (1-dec_st)*prev_b_num + dec_st*b_num
    prev_b_denom = (1-dec_st)*prev_b_denom + dec_st*b_denom 
    
    prev_c_num = (1-dec_st)*prev_c_num + dec_st*c_num
    prev_c_denom = (1-dec_st)*prev_c_denom + dec_st*c_denom
    
    prev_lb=(1-dec_st)*prev_lb+ dec_st*lb3pl(u, xi, new_s, id_1,new_a, new_c,diag(domain), new_b, SIGMA, MU,
                                          Alpha, Beta, mu_b, sigma2_b)
    nlb<-append(nlb,prev_lb)
    
    if(n%%20==0 && n>20 && has_error(Promax(new_a),silent=T)==0){
      n20=append(n20,abs(mean(nlb[(n-39):(n-20)])-mean(nlb[(n-19):n])))
      if(abs(mean(nlb[(n-39):(n-20)])-mean(nlb[(n-19):n]))<0.006 && sum(is.na(Promax(new_a)$loadings[1:item,]))==0){
        converge=0
      }
    }
    n=n+1
  }
  #rotate factors: promax
  rotated=Promax(new_a)
  new_a=rotated$loadings[1:item,]
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
  lbound=lb3pl(u,xi,new_s,person,new_a, new_c, Sigma, new_b, SIGMA, MU,Alpha, Beta, mu_b, sigma2_b)
  gic=log(log(person))*log(person)*sum(Q_mat+item) - 2*lbound
  return(list(ra=new_a,rb=new_b,rc=new_c,rs = new_s,
              reta = eta,reps=xi,rsigma = Sigma,mu_i = MU,
              sig_i = SIGMA,n=n,Q_mat=Q_mat,GIC=gic,rk=rk))
}