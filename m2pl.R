library(shiny)
library(psych)
library(Rcpp)
library(testit)
library(Matrix)
library(gtools)
####initialization fuctions#####
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

#####cfa_2p functions###
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

#####efa_2pl functions######
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
arma::mat Sigma= eye(domain,domain);
arma::mat sigmahat=solve((solve(Sigma,eye(domain,domain))+2*sigma_part),eye(domain,domain));
arma::vec muhat=sigmahat*mu_part;
SIGMA.slice(i)=sigmahat;
MU.col(i)=muhat;
arma::mat apro=new_a*(sigmahat+muhat*trans(muhat))*trans(new_a);
xi.row(i)=trans(sqrt(square(new_b)-2*new_b%(new_a*muhat)+apro.diag()));
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
if(NumericVector::is_na(eta(i,j))){
continue;}
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
  new_a=replace(new_a,new_a>4,4)
  new_b=-qnorm(colSums(u,na.rm=T)/person,0,1)/r
  new_b[new_b>4]=4
  new_b[new_b<(-4)]=-4
  Sigma = diag(domain)
  theta=matrix(rnorm(person*domain,0,1),nrow=person)
  #person*item
  xi=array(1,person)%*%t(new_b)-theta%*%t(new_a)
  eta=(exp(xi)/(1+exp(xi))-0.5)/(2*xi)
  eta[is.na(u)]=NA
  
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
    xi[is.na(u)]=NA
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
    new_b=colSums(0.5-u+2*eta*b_part,na.rm = T)/colSums(2*eta,na.rm = T)
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
  #AIC, BIC
  bic = log(person)*sum(Q_mat+item) - 2*lbound
  aic = 2*sum(Q_mat+item) -2*lbound
  return(list(ra=new_a,rb=new_b,reta = eta,reps=xi,rsigma = Sigma,mu_i = MU,
              sig_i = SIGMA,n=n,rk=rk,Q_mat=Q_mat,GIC=gic,AIC=aic,
              BIC=bic))
}

#####lasso_c1_2pl functions######
#update a without penalty
na_lc12pl<-'
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
arma::mat nalc12pl(const arma::mat&u,const arma::mat&indic,const arma::vec& nopenalty_col,
const int& person, const arma::mat& eta,const arma::vec& new_b,const arma::cube& SIGMA, const arma::mat& MU){
int item=indic.n_rows;
int domain=indic.n_cols;
arma::mat new_a=zeros(item,domain);
for(int j=0; j<nopenalty_col.n_elem; ++j){
int k=nopenalty_col(j)-1;
arma::rowvec Ind=indic.row(k);
arma::uvec iind=find(Ind==1);
arma::vec a_nu= zeros(iind.n_elem);
arma::mat a_de= zeros(iind.n_elem,iind.n_elem);
for(int i=0; i<person; ++i){
arma::mat sigma=SIGMA.slice(i);
sigma=sigma.submat(iind,iind);
arma::vec mu=MU.col(i);
mu=mu.elem(iind);
a_de=a_de+eta(i,k)*sigma+eta(i,k)*(mu*trans(mu));
a_nu=a_nu+(u(i,k)-0.5+2*new_b(k)*eta(i,k))*mu;
}
arma::uvec id(1);
id.at(0)=k;
new_a.submat(id,iind)=trans(solve(a_de,eye(iind.n_elem,iind.n_elem))*a_nu/2);
}
return new_a;
}
'
sourceCpp(code=na_lc12pl)

#update a with penalty
pa_lc12pl<-'
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
arma::mat palc12pl(const arma::mat&u,const int& domain,const int& person,const double& lbd, const arma::vec& sdf, const arma::mat& eta,
const arma::mat& a,
const arma::vec& new_b,const arma::cube& SIGMA, const arma::mat& MU){
arma::mat new_a=a;
for(int j=0; j<sdf.n_elem; ++j){
int m=sdf(j)-1;
for(int k=0; k<domain; ++k){
double delta = 0; 
double deriv2 = 0;
for(int i=0; i<person; ++i){
arma::mat sigma=SIGMA.slice(i);
arma::vec mu=MU.col(i);
arma::mat sigmumu=sigma+mu*trans(mu);
arma::uvec A=regspace<uvec>(0, 1,k-1);
arma::uvec B=regspace<uvec>(k+1, 1, domain-1);
arma::uvec ex_k=join_vert(A,B);
arma::uvec id(1);
id.at(0)=m;
arma::uvec id2(1);
id2.at(0)=k;
delta = delta + (u(i,m)-0.5)*mu(k)+2*new_b(m)*eta(i,m)*mu(k)-2*eta(i,m)*dot(trans(new_a.submat(id,ex_k)),sigmumu.submat(ex_k,id2));
deriv2 = deriv2 + 2*eta(i,m)*sigmumu(k,k);
}
if(abs(delta)>lbd){
double S=sign(delta)*(abs(delta) - lbd);
new_a(m,k) = S/deriv2;
}else{
new_a(m,k) = 0;
}
}
}
return new_a;
}
'
sourceCpp(code=pa_lc12pl)


#lasso with constraint 1 function
vem_2PLEFA_L1_const1 <- function(u,new_a,new_b,eta,xi,Sigma, domain,lbd,indic,nopenalty_col,updateProgress=NULL) {
  person=dim(u)[1]
  item=dim(u)[2]
  converge = 0
  is_singular = 0
  
  MU=matrix(0,nrow=domain,ncol=person)
  SIGMA=array(0,dim=c(domain,domain,person))
  n = 0
  par_Sigma = Sigma
  #nv<-NULL
  while(converge==0 && rankMatrix(Sigma) == domain && n < 1000){
    if (is.function(updateProgress)) {
      text <- paste0("n=:", n)
      updateProgress(detail = text)
    }
    par_Sigma = Sigma
    #update MU, SIGMA, Sigma
    rs1<-ecfa2(u,Sigma,domain, person, item, eta, new_a, new_b)
    xi=rs1$xi
    SIGMA=rs1$SIGMA
    MU=rs1$MU
    Spart=rs1$Spart
    Sigma=Spart/person
    d_temp=sqrt(diag(diag(Sigma)))
    Sigma=solve(d_temp)%*%Sigma%*%solve(d_temp)
    eta=(exp(xi)/(1+exp(xi))-0.5)/(2*xi)
    eta=replace(eta,abs(xi)< 0.01,0.125)
    
    
    #update b
    par_b=new_b
    b_part = t(new_a%*%MU)
    new_b=colSums(0.5-u+2*eta*b_part)/colSums(2*eta)
    
    par_a=new_a
    #update a
    new_a1=nalc12pl(u, indic, nopenalty_col, person, eta, new_b, SIGMA, MU)
    new_a1[-nopenalty_col,]=new_a[-nopenalty_col,]
    #L1-penalty
    sdf=setdiff(1:item,nopenalty_col)
    new_a=palc12pl(u, domain, person, lbd, sdf, eta, new_a1, new_b, SIGMA, MU)
    #par_a=new_a2
    #nv<-append(nv,norm(as.vector(new_a)-as.vector(par_a),type="2")+norm(new_b-par_b,type="2")+
    #             norm(as.vector(Sigma)-as.vector(par_Sigma),type="2"))
    if (norm(as.vector(new_a)-as.vector(par_a),type="2")+norm(new_b-par_b,type="2")+
        norm(as.vector(Sigma)-as.vector(par_Sigma),type="2")<0.001){
      converge=1
    }
    n=n+1
  }
  if (rankMatrix(Sigma) < domain){
    rsigma = par_Sigma
    is_singular = 1
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
  return(list(ra=new_a,rb=new_b,reta = eta,reps=xi,rsigma = Sigma,mu_i = MU,sig_i = SIGMA,n=n,
              Q_mat=Q_mat,GIC=gic,AIC=aic,
              BIC=bic))
}

#main function: choose optimal lambda
vem_2PLEFA_L1_const1_all<-function(u,domain,indic,updateProgress=NULL){
  lbd=seq(2,20,2)
  person=dim(u)[1]
  item=dim(u)[2]
  nopenalty_col=which(rowSums(indic)==1)
  #initialization
  initial=init(u,domain,indic)
  new_a = initial[[1]]
  new_b = initial[[2]]
  eta= initial[[3]]
  xi=initial[[4]]
  Sigma=initial[[5]]
  rl<-vector("list",length(lbd))
  gic<-NULL
  for(j in 1:length(lbd)){
    if (is.function(updateProgress)) {
      text <- paste0("j=:", j)
      updateProgress(detail = text)
    }
    r0=vem_2PLEFA_L1_const1(u,new_a,new_b,eta,xi,Sigma,domain,lbd[j],indic,nopenalty_col,updateProgress=NULL) 
    rl [[j]]=vem_2PLCFA(u,domain, r0$Q_mat,updateProgress=NULL)
    lbound=lb2pl(u,rl[[j]]$reps,rl[[j]]$rsigma,rl[[j]]$ra,
                 rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i)
    gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat+item) - 2*lbound
    new_a = r0$ra
    new_b = r0$rb
    eta= r0$reta
    xi=r0$reps
    Sigma=r0$rsigma
  }
  id=which.min(gic)
  #adjust the range of lambda
  if(id==1){
    lbd=seq(0.1,2,length.out=6)
    rl<-vector("list",length(lbd))
    gic<-NULL
    for(j in 1:length(lbd)){
      r0=vem_2PLEFA_L1_const1(u,new_a,new_b,eta,xi,Sigma,domain,lbd[j],indic,nopenalty_col,updateProgress=NULL) 
      rl [[j]]=vem_2PLCFA(u,domain, r0$Q_mat,updateProgress=NULL)
      lbound=lb2pl(u,rl[[j]]$reps,rl[[j]]$rsigma,rl[[j]]$ra,
                   rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i)
      gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat+item) - 2*lbound
      new_a = r0$ra
      new_b = r0$rb
      eta= r0$reta
      xi=r0$reps
      Sigma=r0$rsigma
    }}
  else if(id==10){
    lbd=seq(20,40,length.out=6)
    rl<-vector("list",length(lbd))
    gic<-NULL
    for(j in 1:length(lbd)){
      r0=vem_2PLEFA_L1_const1(u,new_a,new_b,eta,xi,Sigma,domain,lbd[j],indic,nopenalty_col,updateProgress=NULL) 
      rl [[j]]=vem_2PLCFA(u,domain, r0$Q_mat,updateProgress=NULL)
      lbound=lb2pl(u,rl[[j]]$reps,rl[[j]]$rsigma,rl[[j]]$ra,
                   rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i)
      gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat+item) - 2*lbound
      new_a = r0$ra
      new_b = r0$rb
      eta= r0$reta
      xi=r0$reps
      Sigma=r0$rsigma
    }
  }
  id=which.min(gic)
  rs=rl[[id]]
  rs$lbd=lbd[id]
  rs$id=id
  return(rs)
}

#####lasso_c2_2pl functions######
#update a without penalty
na_lc22pl<-'
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
arma::mat nalc22pl(const arma::mat&u,const int& domain,
const arma::mat& a,const arma::vec& nopenalty_col,
const arma::vec& lastone,
const int& person, const arma::mat& eta,const arma::vec& new_b,const arma::cube& SIGMA, const arma::mat& MU){
arma::mat new_a=a;
for(int j=0; j<nopenalty_col.n_elem; ++j){
int k=nopenalty_col(j)-1;
int l=lastone(j)-1;
double a_nu= 0;
double a_de= 0;
for(int i=0; i<person; ++i){
double sigma=SIGMA(l,l,i);
double mu=MU(l,i);
a_de=a_de+eta(i,k)*sigma+eta(i,k)*(mu*mu);
a_nu=a_nu+(u(i,k)-0.5+2*new_b(k)*eta(i,k))*mu;
}
new_a(k,l)=1/a_de*a_nu/2;
}
return new_a;
}
'
sourceCpp(code=na_lc22pl)

#update a with penalty: 1:domain, off-diagonal
pa_lc22pl<-'
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
arma::mat palc22pl(const arma::mat&u,const arma::mat& a,
const arma::vec& nopenalty_col,
const arma::vec& lastone,const double& lbd,
const int& person, const arma::mat& eta,const arma::vec& new_b,const arma::cube& SIGMA, const arma::mat& MU){
arma::mat new_a1=a;
int domain=new_a1.n_cols;
for(int j=0; j<nopenalty_col.n_elem; ++j){
int n=nopenalty_col(j)-1;
int l=lastone(j)-1;
arma::uvec A=regspace<uvec>(0, 1,l-1);
arma::uvec B=regspace<uvec>(l+1, 1, domain-1);
arma::uvec sdf=join_vert(A,B);
for(int k=0; k<sdf.n_elem; ++k){
double delta = 0; 
double deriv2 = 0;
int m=sdf(k);
for(int i=0; i<person; ++i){
arma::mat sigma=SIGMA.slice(i);
arma::vec mu=MU.col(i);
arma::mat sigmumu=sigma+mu*trans(mu);
arma::uvec A1=regspace<uvec>(0, 1,m-1);
arma::uvec B1=regspace<uvec>(m+1, 1, domain-1);
arma::uvec ex_k=join_vert(A1,B1);
arma::uvec id(1);
id.at(0)=n;
arma::uvec id2(1);
id2.at(0)=m;
delta = delta + (u(i,n)-0.5)*mu(m)+2*new_b(n)*eta(i,n)*mu(m)-2*eta(i,n)*dot(trans(new_a1.submat(id,ex_k)),sigmumu.submat(ex_k,id2));
deriv2 = deriv2 + 2*eta(i,n)*sigmumu(m,m);
}
if(abs(delta)>lbd){
double S=sign(delta)*(abs(delta) - lbd);
new_a1(n,m) = S/deriv2;
}else{
new_a1(n,m) = 0;
}
}
}
return new_a1;
}
'
sourceCpp(code=pa_lc22pl)

#update a with penalty: domain+1:item
pa_lc22pl1<-'
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
arma::mat palc22pl1(const arma::mat&u,const int& domain,const int& item,const int& person,const double& lbd, const arma::mat& eta, arma::mat& new_a,
const arma::vec& new_b,const arma::cube& SIGMA, const arma::mat& MU,const arma::vec& pc){
for(int n=0; n<pc.n_elem; ++n){
int j=pc(n)-1;
for(int k=0; k<domain; ++k){
double delta = 0; 
double deriv2 = 0;
for(int i=0; i<person; ++i){
arma::mat sigma=SIGMA.slice(i);
arma::vec mu=MU.col(i);
arma::mat sigmumu=sigma+mu*trans(mu);
arma::uvec A=regspace<uvec>(0, 1,k-1);
arma::uvec B=regspace<uvec>(k+1, 1, domain-1);
arma::uvec ex_k=join_vert(A,B);
arma::uvec id(1);
id.at(0)=j;
arma::uvec id2(1);
id2.at(0)=k;
delta = delta + (u(i,j)-0.5)*mu(k)+2*new_b(j)*eta(i,j)*mu(k)-2*eta(i,j)*dot(trans(new_a.submat(id,ex_k)),sigmumu.submat(ex_k,id2));
deriv2 = deriv2 + 2*eta(i,j)*sigmumu(k,k);
}
if(abs(delta)>lbd){
double S=sign(delta)*(abs(delta) - lbd);
new_a(j,k) = S/deriv2;
}else{
new_a(j,k) = 0;
}
}
}
return new_a;
}
'
sourceCpp(code=pa_lc22pl1)

#lasso with constraint 2 function
vem_2PLEFA_L1_const2 <- function(u,new_a,new_b,eta,xi,Sigma,  domain,lbd,
                                 indic,nopenalty_col,updateProgress=NULL) {
  person=dim(u)[1]
  item=dim(u)[2]
  converge = 0
  is_singular = 0
  
  MU=matrix(0,nrow=domain,ncol=person)
  SIGMA=array(0,dim=c(domain,domain,person))
  n = 0
  
  par_Sigma = Sigma
  while(converge==0 && rankMatrix(Sigma) == domain && n < 1000){
    if (is.function(updateProgress)) {
      text <- paste0("n=:", n)
      updateProgress(detail = text)
    }
    par_Sigma = Sigma
    #update Sigma, MU, SIGMA
    rs1<-ecfa2(u,Sigma,domain, person, item, eta, new_a, new_b)
    xi=rs1$xi
    SIGMA=rs1$SIGMA
    MU=rs1$MU
    Spart=rs1$Spart
    Sigma=Spart/person
    d_temp=sqrt(diag(diag(Sigma)))
    Sigma=solve(d_temp)%*%Sigma%*%solve(d_temp)
    eta=(exp(xi)/(1+exp(xi))-0.5)/(2*xi)
    eta=replace(eta,abs(xi)< 0.01,0.125)
    
    
    #update b
    par_b=new_b
    b_part = t(new_a%*%MU)
    new_b=colSums(0.5-u+2*eta*b_part)/colSums(2*eta)
    
    par_a=new_a
    #update a
    #find the last one for each item by using indicator matrix
    lastone=apply(indic[nopenalty_col,], 1, function(x) tail(which(x!=0),1))
    new_a=nalc22pl(u, domain,new_a,nopenalty_col,lastone, person, eta, new_b, SIGMA, MU)
    #L1-penalty: off-diagnoal
    new_a=palc22pl(u, new_a,nopenalty_col,lastone, lbd, person, eta, new_b, SIGMA, MU)
    #upper-tiangular should be zero
    new_a=replace(new_a,indic==0,0)
    #domain+1:item
    #find penaly columns
    pc=setdiff(1:item,nopenalty_col)
    new_a=palc22pl1(u, domain, item, person, lbd, eta, new_a, new_b, SIGMA, MU,pc)
    #new_a=replace(new_a,new_a< 0,0)
    #par_a=new_a2
    if (norm(as.vector(new_a)-as.vector(par_a),type="2")+norm(new_b-par_b,type="2")+
        norm(as.vector(Sigma)-as.vector(par_Sigma),type="2")<0.001){
      converge=1
    }
    n=n+1
  }
  if (rankMatrix(Sigma) < domain){
    rsigma = par_Sigma
    is_singular = 1
  }else{
    rsigma = Sigma
  }
  new_a=replace(new_a,abs(new_a)< 0.001,0)
  Q_mat = (new_a != 0)*1
  #gic
  lbound=lb2pl(u,xi,Sigma,new_a,new_b,SIGMA,MU)
  gic=log(log(person))*log(person)*sum(Q_mat+item) - 2*lbound
  #AIC, BIC
  bic = log(person)*sum(Q_mat+item) - 2*lbound
  aic = 2*sum(Q_mat+item) -2*lbound
  return(list(ra=new_a,rb=new_b,reta = eta,reps=xi,rsigma = Sigma,mu_i = MU,sig_i = SIGMA,n=n,
              Q_mat=Q_mat,GIC=gic,AIC=aic,
              BIC=bic))
}

#main function: choose optimal lambda
vem_2PLEFA_L1_const2_all<-function(u,domain,indic,non_pen,updateProgress=NULL){
  lbd=seq(2,20,2)
  person=dim(u)[1]
  item=dim(u)[2]
  nopenalty_col=c(which(rowSums(indic)<domain),non_pen)
  #initialization
  initial=init(u,domain,indic)
  new_a = initial[[1]]
  new_b = initial[[2]]
  eta= initial[[3]]
  xi=initial[[4]]
  Sigma=initial[[5]]
  rl<-vector("list",length(lbd))
  gic<-NULL
  for(j in 1:length(lbd)){
    if (is.function(updateProgress)) {
      text <- paste0("j=:", j)
      updateProgress(detail = text)
    }
    r0=vem_2PLEFA_L1_const2(u,new_a,new_b,eta,xi,Sigma,  domain,lbd[j],indic,nopenalty_col,updateProgress=NULL) 
    rl [[j]]=vem_2PLCFA(u,domain, r0$Q_mat,updateProgress=NULL)
    lbound=lb2pl(u,rl[[j]]$reps,rl[[j]]$rsigma,rl[[j]]$ra,
                 rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i)
    gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat+item) - 2*lbound
    new_a = r0$ra
    new_b = r0$rb
    eta= r0$reta
    xi=r0$reps
    Sigma=r0$rsigma
  }
  id=which.min(gic)
  #adjust the range of lambda
  if(id==1){
    lbd=seq(0.1,2,length.out=6)
    rl<-vector("list",length(lbd))
    gic<-NULL
    for(j in 1:length(lbd)){
      r0=vem_2PLEFA_L1_const2(u,new_a,new_b,eta,xi,Sigma,  domain,lbd[j],indic,nopenalty_col,updateProgress=NULL) 
      rl [[j]]=vem_2PLCFA(u,domain, r0$Q_mat,updateProgress=NULL)
      lbound=lb2pl(u,rl[[j]]$reps,rl[[j]]$rsigma,rl[[j]]$ra,
                   rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i)
      gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat+item) - 2*lbound
      new_a = r0$ra
      new_b = r0$rb
      eta= r0$reta
      xi=r0$reps
      Sigma=r0$rsigma
    }}
  else if(id==10){
    lbd=seq(20,40,length.out=6)
    rl<-vector("list",length(lbd))
    gic<-NULL
    for(j in 1:length(lbd)){
      r0=vem_2PLEFA_L1_const2(u,new_a,new_b,eta,xi,Sigma,  domain,lbd[j],indic,nopenalty_col,updateProgress=NULL) 
      rl [[j]]=vem_2PLCFA(u,domain, r0$Q_mat,updateProgress=NULL)
      lbound=lb2pl(u,rl[[j]]$reps,rl[[j]]$rsigma,rl[[j]]$ra,
                   rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i)
      gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat+item) - 2*lbound
      new_a = r0$ra
      new_b = r0$rb
      eta= r0$reta
      xi=r0$reps
      Sigma=r0$rsigma
    }
  }
  id=which.min(gic)
  rs=rl[[id]]
  rs$lbd=lbd[id]
  rs$id=id
  return(rs)
}

#####adaptive_lasso_c1_2pl functions######
#update penalty a
pa_al2pl<-'
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
arma::mat paal2pl(const arma::mat&u,const int& domain,const int& person,const double& lbd, const arma::vec& sdf, const arma::mat& eta,
const arma::mat& a,
const arma::vec& new_b,const arma::cube& SIGMA, const arma::mat& MU,const arma::mat& weights){
arma::mat new_a=a;
for(int j=0; j<sdf.n_elem; ++j){
int m=sdf(j)-1;
for(int k=0; k<domain; ++k){
double delta = 0; 
double deriv2 = 0;
for(int i=0; i<person; ++i){
arma::mat sigma=SIGMA.slice(i);
arma::vec mu=MU.col(i);
arma::mat sigmumu=sigma+mu*trans(mu);
arma::uvec A=regspace<uvec>(0, 1,k-1);
arma::uvec B=regspace<uvec>(k+1, 1, domain-1);
arma::uvec ex_k=join_vert(A,B);
arma::uvec id(1);
id.at(0)=m;
arma::uvec id2(1);
id2.at(0)=k;
delta = delta + (u(i,m)-0.5)*mu(k)+2*new_b(m)*eta(i,m)*mu(k)-2*eta(i,m)*dot(trans(new_a.submat(id,ex_k)),sigmumu.submat(ex_k,id2));
deriv2 = deriv2 + 2*eta(i,m)*sigmumu(k,k);
}
if(abs(delta)>lbd*(1/weights(m,k))){
double S=sign(delta)*(abs(delta) - lbd*(1/weights(m,k)));
new_a(m,k) = S/deriv2;
}else{
new_a(m,k) = 0;
}
}
}
return new_a;
}
'
sourceCpp(code=pa_al2pl)


#adaptive lasso with constraint 1 function
vem_2PLEFA_adaptive_const1 <- function(u,new_a,new_b,eta,xi,Sigma, domain,lbd,indic,nopenalty_col,weights,updateProgress=NULL) {
  person=dim(u)[1]
  item=dim(u)[2]
  converge = 1
  #is_singular = 0
  MU=matrix(0,nrow=domain,ncol=person)
  SIGMA=array(0,dim=c(domain,domain,person))
  n = 0
  par_Sigma = Sigma
  while(converge==1 && rankMatrix(Sigma) == domain && n<3000){
    if (is.function(updateProgress)) {
      text <- paste0("n=:", n)
      updateProgress(detail = text)
    }
    #update Sigma, MU, SIGMA
    par_Sigma = Sigma
    rs1<-ecfa2(u,Sigma,domain, person, item, eta, new_a, new_b)
    xi=rs1$xi
    SIGMA=rs1$SIGMA
    MU=rs1$MU
    Spart=rs1$Spart
    Sigma=Spart/person
    d_temp=sqrt(diag(diag(Sigma)))
    Sigma=solve(d_temp)%*%Sigma%*%solve(d_temp)
    eta=(exp(xi)/(1+exp(xi))-0.5)/(2*xi)
    eta=replace(eta,abs(xi)< 0.01,0.125)
    
    
    #update b
    par_b=new_b
    b_part = t(new_a%*%MU)
    new_b=colSums(0.5-u+2*eta*b_part)/colSums(2*eta)
    
    par_a=new_a
    #update a
    new_a1=nalc12pl(u, indic, nopenalty_col, person, eta, new_b, SIGMA, MU)
    new_a1[-nopenalty_col,]=new_a[-nopenalty_col,]
    #adaptive lasso penalty
    sdf=setdiff(1:item,nopenalty_col)
    new_a=paal2pl(u, domain, person, lbd, sdf, eta, new_a1, new_b, SIGMA, MU,weights)
    #par_a=new_a2
    if (norm(as.vector(new_a)-as.vector(par_a),type="2")+norm(new_b-par_b,type="2")<0.0001){
      converge=0
    }
    n=n+1
  }
  if (rankMatrix(Sigma) < domain){
    rsigma = par_Sigma
    is_singular = 1
  }else{
    rsigma = Sigma
  }
  #new_a=replace(new_a,new_a< 0,0)
  new_a=replace(new_a,abs(new_a)< 0.001,0)
  Q_mat = (new_a != 0)*1
  return(list(ra=new_a,rb=new_b,reta = eta,reps=xi,rsigma = Sigma,mu_i = MU,sig_i = SIGMA,n=n,
              Q_mat=Q_mat))
}

#main function: choose optimal lambda
vem_2PLEFA_adaptive_const1_all<-function(u,domain,indic,gamma,updateProgress=NULL){
  lbd=seq(2,20,2)
  nopenalty_col=which(rowSums(indic)==1)
  person=dim(u)[1]
  item=dim(u)[2]
  #initialization
  initial=init(u,domain,indic)
  new_a = initial[[1]]
  new_b = initial[[2]]
  eta= initial[[3]]
  xi=initial[[4]]
  Sigma=initial[[5]]
  #create weights: use CFA
  wa=vem_2PLCFA(u,domain,indic,updateProgress=NULL)$ra
  weights = abs(wa)^gamma+1e-05
  rl<-vector("list",length(lbd))
  gic<-NULL
  for(j in 1:length(lbd)){
    if (is.function(updateProgress)) {
      text <- paste0("j=:", j)
      updateProgress(detail = text)
    }
    rl [[j]]=vem_2PLEFA_adaptive_const1(u,new_a,new_b,eta,xi,Sigma, domain,lbd[j],indic,nopenalty_col,weights,updateProgress=NULL)
    lbound=lb2pl(u,rl[[j]]$reps,rl[[j]]$rsigma,rl[[j]]$ra,
                 rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i)
    gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat+item) - 2*lbound
    new_a = rl[[j]]$ra
    new_b = rl[[j]]$rb
    eta= rl[[j]]$reta
    xi=rl[[j]]$reps
    Sigma=rl[[j]]$rsigma
  }
  id=which.min(gic)
  #adjust the range of lambda
  if(id==1){
    lbd=seq(0.1,2,length.out=6)
    rl<-vector("list",length(lbd))
    gic<-NULL
    for(j in 1:length(lbd)){
      rl [[j]]=vem_2PLEFA_adaptive_const1(u,new_a,new_b,eta,xi,Sigma, domain,lbd[j],indic,nopenalty_col,weights,updateProgress=NULL)
      lbound=lb2pl(u,rl[[j]]$reps,rl[[j]]$rsigma,rl[[j]]$ra,
                   rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i)
      gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat+item) - 2*lbound
      new_a = rl[[j]]$ra
      new_b = rl[[j]]$rb
      eta= rl[[j]]$reta
      xi=rl[[j]]$reps
      Sigma=rl[[j]]$rsigma
    }}
  else if(id==10){
    lbd=seq(20,40,length.out=6)
    rl<-vector("list",length(lbd))
    gic<-NULL
    for(j in 1:length(lbd)){
      rl [[j]]=vem_2PLEFA_adaptive_const1(u,new_a,new_b,eta,xi,Sigma, domain,lbd[j],indic,nopenalty_col,weights,updateProgress=NULL)
      lbound=lb2pl(u,rl[[j]]$reps,rl[[j]]$rsigma,rl[[j]]$ra,
                   rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i)
      gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat+item) - 2*lbound
      new_a = rl[[j]]$ra
      new_b = rl[[j]]$rb
      eta= rl[[j]]$reta
      xi=rl[[j]]$reps
      Sigma=rl[[j]]$rsigma
    }
  }
  id=which.min(gic)
  rs=vem_2PLCFA(u,domain, rl[[id]]$Q_mat,updateProgress=NULL)
  rs$lbd=lbd[id]
  rs$id=id
  return(rs)
}

#####adaptive_lasso_c2_2pl functions######
#update a with penalty: 1:domain, off-diagonal
pa_alc22pl<-'
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
arma::mat paalc22pl(const arma::mat&u,const arma::mat& a,
const arma::mat&indic,const arma::vec& nopenalty_col,
const arma::vec& lastone,const double& lbd,
const int& person, const arma::mat& eta,const arma::vec& new_b,const arma::cube& SIGMA, const arma::mat& MU,
const arma::mat&weights){
arma::mat new_a1=a;
int domain=new_a1.n_cols;
for(int j=0; j<nopenalty_col.n_elem; ++j){
int n=nopenalty_col(j)-1;
int l=lastone(j)-1;
arma::uvec A=regspace<uvec>(0, 1,l-1);
arma::uvec B=regspace<uvec>(l+1, 1, domain-1);
arma::uvec sdf=join_vert(A,B);
for(int k=0; k<sdf.n_elem; ++k){
double delta = 0; 
double deriv2 = 0;
int m=sdf(k);
for(int i=0; i<person; ++i){
arma::mat sigma=SIGMA.slice(i);
arma::vec mu=MU.col(i);
arma::mat sigmumu=sigma+mu*trans(mu);
arma::uvec A1=regspace<uvec>(0, 1,m-1);
arma::uvec B1=regspace<uvec>(m+1, 1, domain-1);
arma::uvec ex_k=join_vert(A1,B1);
arma::uvec id(1);
id.at(0)=n;
arma::uvec id2(1);
id2.at(0)=m;
delta = delta + (u(i,n)-0.5)*mu(m)+2*new_b(n)*eta(i,n)*mu(m)-2*eta(i,n)*dot(trans(new_a1.submat(id,ex_k)),sigmumu.submat(ex_k,id2));
deriv2 = deriv2 + 2*eta(i,n)*sigmumu(m,m);
}
if(abs(delta)>lbd*(1/weights(n,m))){
double S=sign(delta)*(abs(delta) - lbd*(1/weights(n,m)));
new_a1(n,m) = S/deriv2;
}else{
new_a1(n,m) = 0;
}
}
}
return new_a1;
}
'
sourceCpp(code=pa_alc22pl)

#update a with penalty: domain+1:item
pa_alc22pl1<-'
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
arma::mat paalc22pl1(const arma::mat&u,const int& domain,const int& item,const int& person,const double& lbd, const arma::mat& eta, arma::mat& new_a,
const arma::vec& new_b,const arma::cube& SIGMA, const arma::mat& MU,const arma::vec& pc,const arma::mat&weights){
for(int n=0; n<pc.n_elem; ++n){
int j=pc(n)-1;
for(int k=0; k<domain; ++k){
double delta = 0; 
double deriv2 = 0;
for(int i=0; i<person; ++i){
arma::mat sigma=SIGMA.slice(i);
arma::vec mu=MU.col(i);
arma::mat sigmumu=sigma+mu*trans(mu);
arma::uvec A=regspace<uvec>(0, 1,k-1);
arma::uvec B=regspace<uvec>(k+1, 1, domain-1);
arma::uvec ex_k=join_vert(A,B);
arma::uvec id(1);
id.at(0)=j;
arma::uvec id2(1);
id2.at(0)=k;
delta = delta + (u(i,j)-0.5)*mu(k)+2*new_b(j)*eta(i,j)*mu(k)-2*eta(i,j)*dot(trans(new_a.submat(id,ex_k)),sigmumu.submat(ex_k,id2));
deriv2 = deriv2 + 2*eta(i,j)*sigmumu(k,k);
}
if(abs(delta)>lbd*(1/weights(j,k))){
double S=sign(delta)*(abs(delta) - lbd*(1/weights(j,k)));
new_a(j,k) = S/deriv2;
}else{
new_a(j,k) = 0;
}
}
}
return new_a;
}
'
sourceCpp(code=pa_alc22pl1)



#adaptive lasso with constraint 2 function
vem_2PLEFA_adaptive_const2 <- function(u, new_a,new_b,eta,xi,Sigma, domain,lbd,
                                       indic,nopenalty_col,weights,updateProgress=NULL) {
  person=dim(u)[1]
  item=dim(u)[2]
  converge = 1
  
  
  MU=matrix(0,nrow=domain,ncol=person)
  SIGMA=array(0,dim=c(domain,domain,person))
  n = 0
  
  par_Sigma = Sigma
  while(converge==1 && rankMatrix(Sigma) == domain && n<5000){
    if (is.function(updateProgress)) {
      text <- paste0("n=:", n)
      updateProgress(detail = text)
    }
    par_Sigma = Sigma
    #update MU, SIGMA, Sigma
    rs1<-ecfa2(u,Sigma,domain, person, item, eta, new_a, new_b)
    xi=rs1$xi
    SIGMA=rs1$SIGMA
    MU=rs1$MU
    Spart=rs1$Spart
    Sigma=Spart/person
    d_temp=sqrt(diag(diag(Sigma)))
    Sigma=solve(d_temp)%*%Sigma%*%solve(d_temp)
    eta=(exp(xi)/(1+exp(xi))-0.5)/(2*xi)
    eta=replace(eta,abs(xi)< 0.01,0.125)
    
    
    #update b
    par_b=new_b
    b_part = t(new_a%*%MU)
    new_b=colSums(0.5-u+2*eta*b_part)/colSums(2*eta)
    
    par_a=new_a
    #update a
    #find the last one for each item by using indicator matrix
    lastone=apply(indic[nopenalty_col,], 1, function(x) tail(which(x!=0),1))
    new_a=nalc22pl(u, domain,new_a,nopenalty_col,lastone, person, eta, new_b, SIGMA, MU)
    #L1-penalty: off-diagnoal
    new_a=paalc22pl(u, new_a,indic,nopenalty_col,lastone, lbd, person, eta, new_b, SIGMA, MU,weights)
    #upper-tiangular should be zero
    new_a=replace(new_a,indic==0,0)
    #domain+1:item
    #find penaly columns
    pc=setdiff(1:item,nopenalty_col)
    new_a=paalc22pl1(u, domain, item, person, lbd, eta, new_a, new_b, SIGMA, MU,pc,weights)
    
    #par_a=new_a2
    if (norm(as.vector(new_a)-as.vector(par_a),type="2")+norm(new_b-par_b,type="2")+
        norm(as.vector(Sigma)-as.vector(par_Sigma),type="2")<0.001){
      converge=0
    }
    n=n+1
  }
  if (rankMatrix(Sigma) < domain){
    rsigma = par_Sigma
    #is_singular = 1
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
  return(list(ra=new_a,rb=new_b,reta = eta,reps=xi,rsigma = Sigma,mu_i = MU,sig_i = SIGMA,n=n,
              Q_mat=Q_mat,GIC=gic,AIC=aic,
              BIC=bic))
}

#main function: choose optimal lambda
vem_2PLEFA_adaptive_const2_all<-function(u,domain,gamma,indic,non_pen,updateProgress=NULL){
  lbd=seq(2,20,2)
  person=dim(u)[1]
  item=dim(u)[2]
  nopenalty_col=c(which(rowSums(indic)<domain),non_pen)
  #initialization
  initial=init(u,domain,indic)
  new_a = initial[[1]]
  new_b = initial[[2]]
  eta= initial[[3]]
  xi=initial[[4]]
  Sigma=initial[[5]]
  #create weights: use CFA
  wa=vem_2PLCFA(u,domain,indic,updateProgress=NULL)$ra
  weights = abs(wa)^gamma+1e-05
  rl<-vector("list",length(lbd))
  gic<-NULL
  for(j in 1:length(lbd)){
    if (is.function(updateProgress)) {
      text <- paste0("j=:", j)
      updateProgress(detail = text)
    }
    rl [[j]]=vem_2PLEFA_adaptive_const2(u,new_a,new_b,eta,xi,Sigma, domain,lbd[j],
                                        indic,nopenalty_col,weights,updateProgress=NULL)
    lbound=lb2pl(u,rl[[j]]$reps,rl[[j]]$rsigma,rl[[j]]$ra,
                 rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i)
    gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat+item) - 2*lbound
    new_a = rl[[j]]$ra
    new_b = rl[[j]]$rb
    eta= rl[[j]]$reta
    xi=rl[[j]]$reps
    Sigma=rl[[j]]$rsigma
  }
  id=which.min(gic)
  #adjust the range of lambda
  if(id==1){
    lbd=seq(0.1,2,length.out=6)
    rl<-vector("list",length(lbd))
    gic<-NULL
    for(j in 1:length(lbd)){
      rl [[j]]=vem_2PLEFA_adaptive_const2(u,new_a,new_b,eta,xi,Sigma, domain,lbd[j],indic,nopenalty_col,weights,updateProgress=NULL)
      lbound=lb2pl(u,rl[[j]]$reps,rl[[j]]$rsigma,rl[[j]]$ra,
                   rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i)
      gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat+item) - 2*lbound
      new_a = rl[[j]]$ra
      new_b = rl[[j]]$rb
      eta= rl[[j]]$reta
      xi=rl[[j]]$reps
      Sigma=rl[[j]]$rsigma
    }}
  else if(id==10){
    lbd=seq(20,40,length.out=6)
    rl<-vector("list",length(lbd))
    gic<-NULL
    for(j in 1:length(lbd)){
      rl [[j]]=vem_2PLEFA_adaptive_const2(u,new_a,new_b,eta,xi,Sigma, domain,lbd[j],indic,nopenalty_col,weights,updateProgress=NULL)
      lbound=lb2pl(u,rl[[j]]$reps,rl[[j]]$rsigma,rl[[j]]$ra,
                   rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i)
      gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat+item) - 2*lbound
      new_a = rl[[j]]$ra
      new_b = rl[[j]]$rb
      eta= rl[[j]]$reta
      xi=rl[[j]]$reps
      Sigma=rl[[j]]$rsigma
    }
  }
  id=which.min(gic)
  rs=vem_2PLCFA(u,domain, rl[[id]]$Q_mat,updateProgress=NULL)
  rs$lbd=lbd[id]
  rs$id=id
  return(rs)
}

#####shiny app######
options(shiny.maxRequestSize = 30*1024^2)
# Define UI for data upload app ----
ui <- navbarPage("VEMIRT-2PL Shiny App",
                 
                 # App title ----
                 tabPanel("EFA",
                          
                          # Sidebar layout with input and output definitions ----
                          sidebarLayout(
                            
                            # Sidebar panel for inputs ----
                            sidebarPanel(
                              
                              # Input: Select a file ----
                              fileInput("file1", "Choose Data CSV File",
                                        multiple = TRUE,
                                        accept = c("text/csv",
                                                   "text/comma-separated-values,text/plain",
                                                   ".csv")),
                              
                              # Horizontal line ----
                              tags$hr(),
                              
                              # Input: Checkbox if file has header ----
                              checkboxInput("header1", "Header", TRUE),
                              
                              # Input: Select separator ----
                              radioButtons("sep1", "Separator",
                                           choices = c(Comma = ",",
                                                       Semicolon = ";",
                                                       Tab = "\t"),
                                           selected = ","),
                              
                              
                              # Horizontal line ----
                              tags$hr(),
                              
                              # Input: Select number of rows to display ----
                              radioButtons("disp1", "Display",
                                           choices = c(Head = "head",
                                                       All = "all"),
                                           selected = "head"),
                              
                              # Horizontal line ----
                              tags$hr(),
                              #Input: Select EFA methods ----
                              selectInput("method", 
                                          label = "Choose an EFA method",
                                          choices = list("Rotation"='Rot', 
                                                         "Adaptive Lasso"='AL',
                                                         "Lasso"="La"),
                                          selected = "Rot"),
                              
                              
                              #Choose Constraint 1
                              conditionalPanel(
                                condition = "input.method == 'AL' || input.method == 'La'",
                                radioButtons("cons", "Constraint",
                                             choices = c("Constraint 1" = "c1",
                                                         "Constraint 2" = "c2"),
                                             selected = "c1"),
                                helpText("Note: For constraint 1, you need to set a K by K submatrix of the indicator 
                                         matrix to be an identity matrix. For constraint 2, you need to set a K by K submatrix of the indicator 
                                         matrix to be a triangular matrix. ") 
                              ),
                              
                              
                              conditionalPanel(
                                condition = "input.method != 'Rot'",
                                # Input: Select the indicator dile ----
                                fileInput("file5", "Choose Loading Indicator CSV File",
                                          multiple = TRUE,
                                          accept = c("text/csv",
                                                     "text/comma-separated-values,text/plain",
                                                     ".csv")),
                                
                                helpText("Note:  The loading indicator file should be a binary Item by Factor matirx with headers. 
                                         0 refers to the item is irrelevant with this factor, 1 otherwise.")
                              ),
                              
                              #input the last number
                              conditionalPanel(
                                condition = "(input.method == 'AL' || input.method == 'La') && input.cons == 'c2' ",
                                numericInput("np", "Item Index", min=1, max=1000, value=5),
                                
                                helpText("Note:  Input the index of an item which load on all domains and satisfies with constraint 2")
                              ),
                              #input the domain
                              tags$hr(),
                              numericInput("dom", "domain", min=1, max=10, value=5),
                              
                              #input the gamma
                              conditionalPanel(
                                condition = "input.method == 'AL' ",
                                numericInput("gm", "Gamma", min=0.1, max=5, value=2),
                                
                                helpText("Note:  Previous literature recommends 0.5, 1, 2 as the gamma value")
                              ),
                              
                              # Horizontal line ----
                              tags$hr(),
                              actionButton("go1", "Run"),
                              
                              
                              # Horizontal line ----
                              tags$hr(),
                              radioButtons("checkGroup1", "Download Results",
                                           choices = list("All Results" = "all",
                                                          "Item Parameters" = "item",
                                                          "Covariance Matrix" = "cov"),
                                           selected = "all"),
                              # Button
                              downloadButton("downloadData", "Download Results")
                              
                            ),
                            
                            # Main panel for displaying outputs ----
                            mainPanel(
                              
                              uiOutput("tab"),
                              uiOutput("tab2"),
                              
                              h2("Data"),
                              # Output: Data file ----
                              tableOutput("contents1"),
                              
                              #warning
                              h1(textOutput("warn")),
                              
                              h2("Item Parameter Results"),
                              tableOutput("par1"),
                              
                              h2("Covariance Matrix"),
                              tableOutput("cov1"),
                              
                              h2("Q Matrix"),
                              tableOutput("fl1"),
                              
                              h2("Model Fit"),
                              textOutput("gic1"),
                              textOutput("aic1"),
                              textOutput("bic1"),
                              
                            )
                            
                          )
                 ),
                 tabPanel("CFA",
                          # Sidebar layout with input and output definitions ----
                          sidebarLayout(
                            
                            # Sidebar panel for inputs ----
                            sidebarPanel(
                              
                              # Input: Select a file ----
                              fileInput("file3", "Choose Data CSV File",
                                        multiple = TRUE,
                                        accept = c("text/csv",
                                                   "text/comma-separated-values,text/plain",
                                                   ".csv")),
                              
                              # Horizontal line ----
                              tags$hr(),
                              
                              # Input: Checkbox if file has header ----
                              checkboxInput("header3", "Header", TRUE),
                              
                              # Input: Select separator ----
                              radioButtons("sep3", "Separator",
                                           choices = c(Comma = ",",
                                                       Semicolon = ";",
                                                       Tab = "\t"),
                                           selected = ","),
                              
                              
                              # Horizontal line ----
                              tags$hr(),
                              
                              # # Input: Select number of rows to display ----
                              radioButtons("disp3", "Display",
                                           choices = c(Head = "head",
                                                       All = "all"),
                                           selected = "head"),
                              
                              #Horizontal line ----
                              tags$hr(),
                              
                              # Input: Select the indicator dile ----
                              fileInput("file4", "Choose Loading Indicator CSV File",
                                        multiple = TRUE,
                                        accept = c("text/csv",
                                                   "text/comma-separated-values,text/plain",
                                                   ".csv")),
                              
                              helpText("Note: The loading indicator file should be a binary Item by Factor matirx. 
               0 refers to the item is irrelevant with this factor, 1 otherwise."),
                              # Horizontal line ----
                              tags$hr(),
                              
                              # Input: Checkbox if file has header ----
                              checkboxInput("header4", "Header", TRUE),
                              
                              # Input: Select separator ----
                              radioButtons("sep4", "Separator",
                                           choices = c(Comma = ",",
                                                       Semicolon = ";",
                                                       Tab = "\t"),
                                           selected = ","),
                              
                              
                              # Horizontal line ----
                              tags$hr(),
                              
                              #Input: Select number of rows to display ----
                              radioButtons("disp4", "Display",
                                           choices = c(Head = "head",
                                                       All = "all"),
                                           selected = "head"),
                            
                              
                              
                              # Horizontal line ----
                              tags$hr(),
                              actionButton("go3", "Run"),
                              
                              tags$hr(),
                              radioButtons("checkGroup2", "Download Results",
                                           choices = list("All Results" = "all",
                                                          "Item Parameters" = "item",
                                                          "Covariance Matrix" = "cov"),
                                           selected = "all"),
                              #Button
                              downloadButton("downloadDataall2", "Download Results")
                              #downloadButton("downloadDataall2", "Download All Results")
                              # downloadButton("downloadDataitem2", "Download Item Parameter Results"),
                              # downloadButton("downloadDatasigma2", "Download Covariance Results")
                              
                              
                              
                            ),
                            
                            # Main panel for displaying outputs ----
                            mainPanel(
                              
                              h2("Data"),
                              # Output: Data file ----
                              tableOutput("contents3"),
                              
                              
                              h2("Loading Indicator Matrix"),
                              #Output: indicator matrix,
                              tableOutput("contents4"),
                              
                              
                              h2("Item Parameter Results"),
                              tableOutput("par2"),
                              
                              h2("Covariance Matrix"),
                              tableOutput("cov2"),
                              
                              h2("Model Fit"),
                              textOutput("gic2"),
                              textOutput("aic2"),
                              textOutput("bic2"),
                              
                              
                            )
                            
                          )
                 ))

# Define server logic to read selected file ----
server <- function(input, output,session) {
  url <- a("Parallel Analysis", href="https://www.google.com/")
  output$tab <- renderUI({
    tagList("Other Functions:", url)
  })
  url2 <- a("M3PL estimation", href="https://www.google.com/",)
  output$tab2 <- renderUI({
    tagList(url2)
  })
  #EFA
  u<-reactive({req(input$file1)
    df <- data.matrix(read.csv(input$file1$datapath,
                               header = input$header1,
                               sep = input$sep1))
    return(df)
  })
  output$contents1 <- renderTable({
    if(input$disp1 == "head") {
      return(u()[1:6,])
    }
    else {
      return(u())
    }
    
  })
  domain<-reactive({req(input$dom)})
  
  result0<-reactive({
    # Create a Progress object
    progress <- shiny::Progress$new()
    progress$set(message = "Iteration times", value = 0)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())
    updateProgress <- function(value = NULL, detail = NULL) {
      if (is.null(value)) {
        value <- progress$getValue()
        value <- value + (progress$getMax() - value) / 10
      }
      progress$set(value = value, detail = detail)
    }
    if (input$method == "Rot") {
      result<-vem_2PLEFA(u(),domain(),updateProgress)}
    else if(input$method == "AL" && input$cons=="c1" ){
      indic<-reactive({req(input$file5)
        df2 <- data.matrix(read.csv(input$file5$datapath))
        return(df2)
      })
      gamma=reactive({req(input$gm)})
      result<-vem_2PLEFA_adaptive_const1_all(u(),domain(),indic(),gamma(),updateProgress)}
    else if(input$method == "AL" && input$cons=="c2"){
      indic<-reactive({req(input$file5)
        df2 <- data.matrix(read.csv(input$file5$datapath))
        return(df2)
      })
      gamma=reactive({req(input$gm)})
      non_pen=reactive({req(input$np)})
      result<-vem_2PLEFA_adaptive_const2_all(u(),domain(),gamma(),indic(),non_pen(),updateProgress)}
    else if(input$method == "La" && input$cons=="c1"){
      indic<-reactive({req(input$file5)
        df2 <- data.matrix(read.csv(input$file5$datapath))
        return(df2)
      })
      result<-vem_2PLEFA_L1_const1_all(u(),domain(),indic(),updateProgress)}
    else if(input$method == "La" && input$cons=="c2"){
      indic<-reactive({req(input$file5)
        df2 <- data.matrix(read.csv(input$file5$datapath))
        return(df2)
      })
      non_pen=reactive({req(input$np)})
      result<-vem_2PLEFA_L1_const2_all(u(),domain(),indic(),non_pen(),updateProgress)}
    return(result)
  }
  )
  
  output$warn<-renderText({
    input$go1
    isolate({
      if(input$method != "Rot" && (result0 ()$lbd == 0.1 || result0 ()$lbd == 40)){
        return("Warning: The optimal penalty parameter may be out of range, a different gamma value is suggested")
      }else{
        return(NULL)
      }})
  })
  
  output$par1<-renderTable({
    input$go1
    isolate({
        m<-cbind(result0()$ra,result0()$rb)
        colnames(m)<-c(paste("a",1:domain(),sep=""),"b")
      return(m)
    })
    
  })
  
  
  output$cov1<-renderTable({
    input$go1
    isolate({
      m<-result0()$rsigma
      rownames(m)<-colnames(m)<-paste("dimension",1:domain(),sep="")
      return(m)
    })
    
  })
  output$fl1<-renderTable({
    input$go1
    isolate({
      m<-result0()$Q_mat
      colnames(m)<-paste("dimension",1:input$dom,sep="")
      #rownames(m)<-paste("Item",1:dim(m)[1],sep="")
      return(m)
    })
    
  })
  
  output$gic1<-renderText({
    input$go1
    isolate({
      paste("GIC =", round(result0()$GIC,2))
    })
  })
  
  output$aic1<-renderText({
    input$go1
    isolate({
      paste("AIC =", round(result0()$AIC,2))
    })
  })
  
  output$bic1<-renderText({
    input$go1
    isolate({
      paste("BIC =", round(result0()$BIC,2))
    })
  })
  #Downloadable csv of selected dataset ----
  
  output$downloadData <- downloadHandler(
    filename = function() {
      if(input$checkGroup1=="all"){
        paste("GVEMEFA",input$checkGroup1 ,"results.rds",sep="")}else{
          paste("GVEMEFA",input$checkGroup1 ,"results.csv",sep="")}
      
    },
    content = function(file) {
      if(input$checkGroup1=="all"){
        saveRDS(result0(), file)}else if(input$checkGroup1=="item"){
            m<-cbind(result0()$ra,result0()$rb)
            colnames(m)<-c(paste("a",1:domain(),sep=""),"b")
          write.csv(m,file,row.names = F)
        }else if(input$checkGroup1=="cov"){
          m<-result0()$rsigma
          rownames(m)<-colnames(m)<-paste("dimension",1:domain(),sep="")
          write.csv(m,file,row.names = F)
        }
    }
  )
  
  #CFA output
  u1<-reactive({req(input$file3)
    df <- data.matrix(read.csv(input$file3$datapath,
                               header = input$header3,
                               sep = input$sep3))
    return(df)
  })
  output$contents3 <- renderTable({
    if(input$disp3 == "head") {
      return(u1()[1:6,])
    }
    else {
      return(u1())
    }
    
  })
  
  indic1<-reactive({req(input$file4)
    df2 <- data.matrix(read.csv(input$file4$datapath,
                                header = input$header4,
                                sep = input$sep4))
    return(df2)
  })
  
  output$contents4 <- renderTable({
    
    if(input$disp4 == "head") {
      return(indic1()[1:6,])
    }
    else {
      return(indic1())
    }
    
  })
  domain1<-reactive(dim(indic1())[2])
  
  result<-reactive({
    # Create a Progress object
    progress <- shiny::Progress$new()
    progress$set(message = "Iteration times", value = 0)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())
    updateProgress <- function(value = NULL, detail = NULL) {
      if (is.null(value)) {
        value <- progress$getValue()
        value <- value + (progress$getMax() - value) / 10
      }
      progress$set(value = value, detail = detail)
    }
      result<-vem_2PLCFA(u1(),domain1(),indic1(),updateProgress) 
    return(result)
  })
  
  output$par2<-renderTable({
    input$go3
    isolate({
        m<-cbind(result()$ra,result()$rb)
        colnames(m)<-c(paste("a",1:domain1(),sep=""),"b")
      return(m)
    })
    
  })
  
  output$cov2<-renderTable({
    input$go3
    isolate({
      m<-result()$rsigma
      rownames(m)<-colnames(m)<-paste("dimension",1:domain1(),sep="")
      return(m)
    })
    
  })
  
  output$gic2<-renderText({
    input$go3
    isolate({
      paste("GIC =", round(result()$GIC,2))
    })
  })
  
  output$aic2<-renderText({
    input$go3
    isolate({
      paste("AIC =", round(result()$AIC,2))
    })
  })
  
  output$bic2<-renderText({
    input$go3
    isolate({
      paste("BIC =", round(result()$BIC,2))
    })
  })
  #Downloadable csv of selected dataset ----
  output$downloadDataall2 <- downloadHandler(
    filename = function() {
      if(input$checkGroup2=="all"){
        paste("GVEMCFA",input$checkGroup2 ,"results.rds",sep="")}else{
          paste("GVEMCFA",input$checkGroup2 ,"results.csv",sep="")}
      
    },
    content = function(file) {
      if(input$checkGroup2=="all"){
        saveRDS(result(), file)}else if(input$checkGroup2=="item"){
            m<-cbind(result()$ra,result()$rb)
            colnames(m)<-c(paste("a",1:domain1(),sep=""),"b")
          write.csv(m,file,row.names = F)
        }else if(input$checkGroup2=="cov"){
          m<-result()$rsigma
          rownames(m)<-colnames(m)<-paste("dimension",1:domain1(),sep="")
          write.csv(m,file,row.names = F)
        }
    }
  )
}
# Run the app ----
shinyApp(ui, server)
