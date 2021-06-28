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

#####cfa_3pl functions######
#e step
e_cfa3pl<-'
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
List ecfa3pl(const arma::mat&u,const arma::mat& Sigma, const int& domain, const arma::vec& id_1,
const int& item, const arma::mat& eta, const arma::mat& new_s,const arma::mat& new_a,const arma::vec& new_b,
const arma::cube& SIGMA1,const arma::mat& MU1,const arma::mat& xi1){
arma::mat MU=MU1;
arma::cube SIGMA=SIGMA1;
arma::mat Spart=mat(domain,domain,fill::zeros);
arma::mat xi=xi1;
for(int i=0; i<id_1.n_elem; ++i){
int k=id_1(i)-1;
arma::mat sigma_part=arma::mat(domain,domain,arma::fill::zeros);
arma::vec mu_part= arma::zeros(domain);
for(int j=0; j<item; ++j){
sigma_part=sigma_part+eta(k,j)*(1-u(k,j)+new_s(k,j)*u(k,j))*trans(new_a.row(j))*new_a.row(j);
mu_part=mu_part+trans((2*eta(k,j)*new_b(j)+u(k,j)-0.5)*(1-u(k,j)+new_s(k,j)*u(k,j))*new_a.row(j));
}
arma::mat sigmahat=solve((solve(Sigma,eye(domain,domain))+2*sigma_part),eye(domain,domain));
arma::vec muhat=sigmahat*mu_part;
SIGMA.slice(k)=sigmahat;
MU.col(k)=muhat;
arma::mat mukk=sigmahat+muhat*trans(muhat);
arma::mat muk=new_a*mukk*trans(new_a);
xi.row(k)=trans(sqrt(square(new_b)-2*new_b%(new_a*muhat)+muk.diag()));
Spart=Spart+sigmahat+muhat*trans(muhat);
}
return List::create(Named("SIGMA") = SIGMA,Named("MU") = MU,Named("xi") = xi,Named("Spart") = Spart);
}
'
sourceCpp(code=e_cfa3pl)

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
a_cfa3pl<-'
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
List acfa3(const arma::mat&u,const arma::mat&indic,const arma::vec& id_1,const int& item,const int& domain, const arma::mat& eta, 
         const arma::mat& a,const arma::vec& new_b,const arma::mat& new_s,const arma::cube& SIGMA, const arma::mat& MU,const arma::mat& prev_a_num,
         const arma::cube& prev_a_denom, const double& dec_st){
  arma::mat new_a=a;
  arma::mat a_num=zeros(domain,item); 
  arma::cube a_denom=zeros(domain,domain,item);
  for(int j=0; j<item; ++j){
    arma::rowvec Ind=indic.row(j);
    arma::uvec iind=find(Ind==1);
    arma::mat a_denom_sub=a_denom.slice(j);
    arma::vec a_num_sub=a_num.col(j);
    arma::vec a_num_sub1=a_num_sub.elem(iind);
    for(int i=0; i<id_1.n_elem; ++i){
      int k=id_1(i)-1;
      arma::mat sigma=SIGMA.slice(k);
      sigma=sigma.submat(iind,iind);
      arma::vec mu=MU.col(k);
      mu=mu.elem(iind);
      double a1=1-u(k,j)+new_s(k,j)*u(k,j);
      a_denom_sub.submat(iind,iind)=a_denom_sub.submat(iind,iind)+a1*eta(k,j)*(sigma+mu*trans(mu));
      a_num_sub1=a_num_sub1+a1*(u(k,j)-0.5+2*new_b(j)*eta(k,j))*mu; 
    }
    a_denom.slice(j)=a_denom_sub;
    a_num_sub.elem(iind)=a_num_sub1;
    a_num.col(j)=a_num_sub;
    arma::mat prev_a_denom_sub=prev_a_denom.slice(j);
    prev_a_denom_sub=prev_a_denom_sub.submat(iind,iind);
    arma::mat a2=solve((dec_st*a_denom_sub.submat(iind,iind)+(1-dec_st)*prev_a_denom_sub),eye(iind.n_elem,iind.n_elem));
    arma::vec prev_a_num_sub=prev_a_num.col(j);
    prev_a_num_sub=prev_a_num_sub.elem(iind);
    arma::uvec id(1);
    id.at(0)=j;
    new_a.submat(id,iind)=trans(a2*(dec_st*a_num_sub1+(1-dec_st)*prev_a_num_sub)/2);
  }
  return List::create(Named("new_a") = new_a,Named("a_denom") = a_denom,Named("a_num") = a_num);
}
'
sourceCpp(code=a_cfa3pl)


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
saem_3PLCFA <- function(u, domain, indic,samp,forgetrate,mu_b,sigma2_b,Alpha,
                        Beta,updateProgress=NULL) {
  person=dim(u)[1]
  item=dim(u)[2]
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
  prev_lb = 0
  converge = 1
  n = 1
  #initialization
  new_c=rep(0.15,item)
  new_s=matrix(0,nrow=person,ncol=item)
  aveg_p=colMeans(u==1)
  for (i in 1:person) {
    for (j in 1:item) {
      if(u[i,j]==0){
        new_s[i,j]=1;
      }else{
        new_s[i,j]=runif(1,aveg_p[j],1)
      }
    }
  }
  initial=init(u,domain,indic)
  new_a = initial[[1]]
  new_b = initial[[2]]
  eta= initial[[3]]
  Sigma=initial[[5]]
  xi=eta
  inite=person
  prev_lb = 0
  nlb<-n20<-NULL
  while(converge==1 && rankMatrix(Sigma) == domain && n<6000){
    if (is.function(updateProgress)) {
      text <- paste0("n=:", n)
      updateProgress(detail = text)
    }
    if(n==1){
      #if(n<21){
      dec_st=1
      id_1 = sample(person,inite)
    }else{
      dec_st = (n+1)^(-forgetrate)*0.25
      id_1 = sample(person,samp)
    }
    s_num = (1-new_c)/new_c
    par_Sigma=Sigma
    #update MU and SIGMA,xi,Spart
    rs1<-ecfa3pl(u, Sigma, domain, id_1, item, eta, new_s, new_a, new_b, SIGMA, MU, xi)
    SIGMA=rs1$SIGMA
    MU=rs1$MU
    xi=rs1$xi
    Spart=rs1$Spart
    #update s
    new_s=s3pl(item,id_1,new_s,eta,new_a,s_num,new_b,u,SIGMA,MU,xi)
    #update eta
    eta=(exp(xi)/(1+exp(xi))-0.5)/(2*xi)
    eta=eta3pl(item,id_1,eta,abs(xi))
    #update Sigma
    Sigma=Spart/length(id_1)
    d_temp=sqrt(diag(diag(Sigma)))
    Sigma=solve(d_temp)%*%Sigma%*%solve(d_temp)
    
    #update b
    par_b=new_b
    rs2=b3pl(item,id_1,new_s,eta,u,new_a,dec_st,MU,prev_b_num,prev_b_denom,mu_b,sigma2_b)
    new_b=rs2$new_b
    b_denom= rs2$b_denom
    b_num = rs2$b_num
    
    #update a
    par_a=new_a
    rs3=acfa3(u, indic, id_1,item, domain, eta, new_a, 
              new_b, new_s, SIGMA, MU, prev_a_num, prev_a_denom, dec_st)
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
    
    #check the lower bound
    old=prev_lb
    prev_lb=(1-dec_st)*prev_lb+ dec_st*lb3pl(u, xi, new_s, id_1, 
                                             new_a, new_c, Sigma, new_b, SIGMA, MU, 
                                             Alpha, Beta, mu_b, sigma2_b)
    nlb<-append(nlb,prev_lb)
    
    #calculate the mean lower bound difference for every 20 iterations
    if(n%%20==0 && n>20){
      n20=append(n20,abs(mean(nlb[(n-39):(n-20)])-mean(nlb[(n-19):n])))
      if(abs(mean(nlb[(n-39):(n-20)])-mean(nlb[(n-19):n]))<0.006){
        #if (nrm1<0.32){  
        converge=0
      }
    }
    n=n+1
  }
  if (rankMatrix(Sigma) < domain){
    rsigma = par_Sigma
  }else{
    rsigma = Sigma
  }
  new_a=replace(new_a,abs(new_a)< 0.001,0)
  Q_mat = (new_a != 0)*1
  #gic
  lbound=lb3pl(u,xi,new_s,person,new_a, new_c, Sigma, new_b, SIGMA, MU,Alpha, Beta, mu_b, sigma2_b)
  gic=log(log(person))*log(person)*sum(Q_mat+item) - 2*lbound
  #AIC, BIC
  bic = log(person)*sum(Q_mat+item) - 2*lbound
  aic = 2*sum(Q_mat+item) -2*lbound
  return(list(ra=new_a,rb=new_b,rc=new_c,rs = new_s,
              reta = eta,reps=xi,rsigma = rsigma,mu_i = MU,sig_i = SIGMA,n=n,
              Q_mat=Q_mat,GIC=gic,AIC=aic,
              BIC=bic))
}


#####efa_3pl functions######
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
  new_a=replace(new_a,new_a>4,4)
  new_b=-qnorm(colSums(u,na.rm=T)/person,0,1)/r
  new_b[new_b>4]=4
  new_b[new_b<(-4)]=-4
  new_c=rep(0.1,item)
  theta=matrix(rnorm(person*domain,0,1),nrow=person)
  #person*item
  xi=array(1,person)%*%t(new_b)-theta%*%t(new_a)
  eta=(exp(xi)/(1+exp(xi))-0.5)/(2*xi)
  eta[is.na(u)]=NA
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
  #AIC, BIC
  bic = log(person)*sum(Q_mat+item) - 2*lbound
  aic = 2*sum(Q_mat+item) -2*lbound
  return(list(ra=new_a,rb=new_b,rc=new_c,rs = new_s,
              reta = eta,reps=xi,rsigma = Sigma,mu_i = MU,
              sig_i = SIGMA,n=n,Q_mat=Q_mat,GIC=gic,rk=rk,AIC=aic,
              BIC=bic))
}

#####lasso_c1_3pl functions######
#update a without penalty
na_lc13pl<-'
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
List nalc13pl(const arma::mat&u,const arma::mat&indic,const arma::vec& nopenalty_col,const arma::mat& eta,const arma::mat& s,
const arma::vec& new_b,const arma::cube& SIGMA, const arma::mat& MU,
const arma::mat& a,const arma::vec& id_1,const arma::mat& prev_a_num,const arma::cube& prev_a_denom, 
const double& dec_st){
int domain=indic.n_cols;
int item=u.n_cols;
arma::mat new_a=a;
arma::mat a_num=mat(domain,item,fill::zeros);
arma::cube a_denom=zeros(domain,domain,item);
for(int j=0; j<nopenalty_col.n_elem; ++j){
int k=nopenalty_col(j)-1;
arma::rowvec Ind=indic.row(k);
arma::uvec iind=find(Ind==1);
arma::mat a_denom_sub=a_denom.slice(j);
arma::uvec id2(1);
id2(0)=k;
for(int i=0; i<id_1.n_elem; ++i){
int l=id_1(i)-1;
double a1=1-u(l,k)+s(l,j)*u(l,k);
arma::mat sigma=SIGMA.slice(l);
arma::uvec id(1);
id(0)=l;
arma::mat a2=sigma.submat(iind,iind)+MU.submat(iind,id)*trans(MU.submat(iind,id));
a_denom_sub.submat(iind,iind)=a_denom_sub.submat(iind,iind)+a1*eta(l,k)*a2;
a_num.submat(iind,id2)=a_num.submat(iind,id2)+a1*(u(l,k)-0.5+2*new_b(k)*eta(l,k))*MU(iind,id);  
}
arma::mat prev_a_denom_sub=prev_a_denom.slice(k);
arma::mat a3=dec_st*a_denom_sub.submat(iind,iind) + (1-dec_st)*prev_a_denom_sub.submat(iind,iind);
arma::mat a4=(dec_st*a_num.submat(iind,id2)+(1-dec_st)*prev_a_num.submat(iind,id2))/2;
new_a.submat(id2,iind)=trans(solve(a3,eye(iind.n_elem,iind.n_elem))*a4);
a_denom.slice(k)=a_denom_sub;
}
return List::create(Named("new_a") = new_a,Named("a_denom") = a_denom,Named("a_num") = a_num);
}
'
sourceCpp(code=na_lc13pl)

#update a with penalty
pa_lc13pl<-'
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
List palc13pl(const arma::mat&u,const int& domain,const int& item,const arma::vec& id_1,const double& lbd, 
const arma::mat& eta,const arma::mat& s,const arma::mat& a,const arma::vec& new_b,const arma::cube& SIGMA, 
const arma::mat& MU,const double& dec_st,const arma::mat& prev_delta,const arma::mat& prev_deriv2,const arma::vec& sdf){
arma::mat new_a=a;
arma::mat delta=mat(item,domain,fill::zeros);
arma::mat deriv2=mat(item,domain,fill::zeros);
for(int n=0; n<sdf.n_elem; ++n){
int j=sdf(n)-1;
for(int k=0; k<domain; ++k){
for(int i=0; i<id_1.n_elem; ++i){
int l=id_1(i)-1;
arma::mat sigma=SIGMA.slice(l);
arma::vec mu=MU.col(l);
arma::mat sigmumu=sigma+mu*trans(mu);
arma::uvec A=regspace<uvec>(0, 1,k-1);
arma::uvec B=regspace<uvec>(k+1, 1, domain-1);
arma::uvec ex_k=join_vert(A,B);
double a1=1-u(l,j)+s(l,j)*u(l,j);
arma::uvec id(1);
id.at(0)=j;
arma::uvec id2(1);
id2.at(0)=k;
delta(j,k) = delta(j,k) + a1*((u(l,j)-0.5)*mu(k)+2*new_b(j)*eta(l,j)*mu(k)-2*eta(l,j)*dot(trans(new_a.submat(id,ex_k)),sigmumu.submat(ex_k,id2)));
deriv2(j,k) = deriv2(j,k) + a1*2*eta(l,j)*sigmumu(k,k);
}
double avg_delta = dec_st*delta(j,k) + (1-dec_st)*prev_delta(j,k);
if(abs(avg_delta)>lbd){
double S=sign(avg_delta)*(abs(avg_delta) - lbd);
new_a(j,k) = S/(dec_st*deriv2(j,k) + (1-dec_st)*prev_deriv2(j,k));
}else{
new_a(j,k) = 0;
}
}
}
return List::create(Named("new_a") = new_a,Named("delta") = delta,Named("deriv2") = deriv2);
}
'
sourceCpp(code=pa_lc13pl)

#lasso with constraint 1 function
sgvem_3PLEFA_L1_const1 <- function(u,new_a,new_b,new_c,new_s,eta,xi,Sigma, domain,inite,samp,forgetrate,lbd,indic,
                                   mu_b,sigma2_b,Alpha,Beta,nopenalty_col,updateProgress=NULL) {
  person=dim(u)[1]
  item=dim(u)[2]
  converge = 1
  MU=matrix(0,nrow=domain,ncol=person)
  SIGMA=array(diag(domain),dim=c(domain,domain,person))
  xi=eta
  prev_a_num = matrix(0,nrow=domain,ncol=item)  
  prev_a_denom = array(0,dim=c(domain,domain,item)) 
  prev_delta = matrix(0,nrow=item,ncol=domain)
  prev_deriv2 = matrix(0,nrow=item,ncol=domain);
  prev_b_num = matrix(0,nrow=item,ncol=1)
  prev_b_denom = matrix(0,nrow=item,ncol=1)
  prev_c_num = matrix(0,nrow=item,ncol=1) 
  prev_c_denom = 0
  n = 1
  par_Sigma = Sigma
  nlb<-n20<-NULL
  prev_lb = 0
  while(converge==1 && rankMatrix(Sigma) == domain && n < 8000){
    if (is.function(updateProgress)) {
      text <- paste0("n=:", n)
      updateProgress(detail = text)
    }
    if(n==1){
      dec_st=1
      id_1 = sample(person,inite)
    }else{
      dec_st = (n+1)^(-forgetrate)*0.25
      #dec_st = (n+1)^(-forgetrate)
      id_1 = sample(person,samp)
    }
    s_num = (1-new_c)/new_c
    par_Sigma = Sigma
    #update MU and SIGMA,xi,Spart
    rs1<-ecfa3pl(u, Sigma, domain, id_1, item, eta, new_s, new_a, new_b, SIGMA, MU, xi)
    SIGMA=rs1$SIGMA
    MU=rs1$MU
    xi=rs1$xi
    Spart=rs1$Spart
    
    #update s
    new_s=s3pl(item,id_1,new_s,eta,new_a,s_num,new_b,u,SIGMA,MU,xi)
    #update eta
    eta=(exp(xi)/(1+exp(xi))-0.5)/(2*xi)
    eta=eta3pl(item,id_1,eta,abs(xi))
    #update Sigma
    Sigma=Spart/length(id_1)
    d_temp=sqrt(diag(diag(Sigma)))
    Sigma=solve(d_temp)%*%Sigma%*%solve(d_temp)
    
    #update b
    par_b=new_b
    rs2=b3pl(item,id_1,new_s,eta,u,new_a,dec_st,MU,prev_b_num,prev_b_denom,mu_b,sigma2_b)
    new_b=rs2$new_b
    b_denom= rs2$b_denom
    b_num = rs2$b_num
    
    
    par_a=new_a
    
    #update a
    rs3=nalc13pl(u, indic, nopenalty_col, eta, new_s, new_b, SIGMA, MU, new_a, id_1, prev_a_num, prev_a_denom, dec_st)
    new_a1=rs3$new_a
    new_a1[-nopenalty_col,]=new_a[-nopenalty_col,]
    a_denom=rs3$a_denom
    a_num=rs3$a_num
    #L1-penalty
    sdf=setdiff(1:item,nopenalty_col)
    rs4=palc13pl(u, domain, item, id_1, lbd, eta, new_s, new_a1, new_b, SIGMA, MU, dec_st, prev_delta, prev_deriv2,sdf)
    new_a=rs4$new_a
    delta=rs4$delta
    deriv2=rs4$deriv2
    
    #if all zeros in some rows, recovery to the previous a
    id_allzero=which(rowSums(new_a)==0)
    if(sum(rowSums(new_a)==0)!=0){
      new_a[id_allzero,] = par_a[id_allzero,] 
      break
    }
    
    
    #update c
    # par_c = new_c;
    # c_num = u[id_1,]*(1-new_s[id_1,])
    # new_c = (dec_st*colSums(c_num)+(1-dec_st)*prev_c_num)/(dec_st*length(id_1)+(1-dec_st)*prev_c_denom)  
    par_c = new_c;
    c_num = colSums(u[id_1,]*(1-new_s[id_1,]))+Alpha - 1
    c_denom = length(id_1) + Alpha + Beta - 2
    new_c = (dec_st*c_num+(1-dec_st)*prev_c_num)/(dec_st*c_denom+(1-dec_st)*prev_c_denom)
    
    #save old sums for SAEM updates
    prev_a_num = (1-dec_st)*prev_a_num + dec_st*a_num
    prev_a_denom = (1-dec_st)*prev_a_denom + dec_st*a_denom
    
    prev_b_num = (1-dec_st)*prev_b_num + dec_st*b_num
    prev_b_denom = (1-dec_st)*prev_b_denom + dec_st*b_denom 
    
    # prev_c_num = (1-dec_st)*prev_c_num + dec_st*colSums(c_num)
    # prev_c_denom = (1-dec_st)*prev_c_denom + dec_st*length(id_1)
    prev_c_num = (1-dec_st)*prev_c_num + dec_st*c_num
    prev_c_denom = (1-dec_st)*prev_c_denom + dec_st*c_denom
    
    prev_delta = (1-dec_st)*prev_delta + dec_st*delta
    prev_deriv2 = (1-dec_st)*prev_deriv2 + dec_st*deriv2
    
    #par_a=new_a2
    #nv<-append(nv,norm(as.vector(new_a)-as.vector(par_a),type="2")+norm(new_b-par_b,type="2")+
    #             norm(as.vector(Sigma)-as.vector(par_Sigma),type="2"))
    prev_lb=(1-dec_st)*prev_lb+ dec_st*lb3pl(u, xi, new_s, id_1, 
                                             new_a, new_c, Sigma, new_b, SIGMA, MU, 
                                             Alpha, Beta, mu_b, sigma2_b)
    nlb<-append(nlb,prev_lb)
    
    if(n%%20==0 && n>20){
      n20=append(n20,abs(mean(nlb[(n-39):(n-20)])-mean(nlb[(n-19):n])))
      if(abs(mean(nlb[(n-39):(n-20)])-mean(nlb[(n-19):n]))<0.006){
        converge=0
      }
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
  lbound=lb3pl(u,xi,new_s,c(1:person),new_a,new_c,Sigma, new_b, SIGMA, MU, 
               Alpha, Beta, mu_b, sigma2_b)
  gic=log(log(person))*log(person)*sum(Q_mat+item) - 2*lbound
  #AIC, BIC
  bic = log(person)*sum(Q_mat+item) - 2*lbound
  aic = 2*sum(Q_mat+item) -2*lbound
  return(list(ra=new_a,rb=new_b,rc=new_c,rs = new_s,reta = eta,reps=xi,rsigma = rsigma,rs = new_s,
              mu_i = MU,sig_i = SIGMA,n=n,Q_mat=Q_mat,n20=n20,GIC=gic,AIC=aic,BIC=bic))
}

#main function: choose optimal lambda
sgvem_3PLEFA_L1_const1_all<-function(u,domain,indic,samp,forgetrate,
                                     mu_b,sigma2_b,Alpha,Beta,updateProgress=NULL){
  lbd=seq(2,20,2)
  inite=dim(u)[1]
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
  new_c=rep(0.15,item)
  new_s=matrix(0,nrow=person,ncol=item)
  aveg_p=colMeans(u==1)
  for (i in 1:person) {
    for (j in 1:item) {
      if(u[i,j]==0){
        new_s[i,j]=1;
      }else{
        new_s[i,j]=runif(1,aveg_p[j],1)
      }
    }
  }
  rl<-vector("list",length(lbd))
  gic<-NULL
  for(j in 1:length(lbd)){
    if (is.function(updateProgress)) {
      text <- paste0("j=:", j)
      updateProgress(detail = text)
    }
    r0=sgvem_3PLEFA_L1_const1(u,new_a,new_b,new_c,new_s,eta,xi,Sigma, domain,
                              inite,samp,forgetrate,lbd[j],indic,
                              mu_b,sigma2_b,Alpha,Beta,nopenalty_col,updateProgress=NULL)
    rl [[j]]=saem_3PLCFA(u,domain, r0$Q_mat,samp,forgetrate,mu_b,sigma2_b,Alpha,
                         Beta,updateProgress=NULL)
    lbound=lb3pl(u,rl[[j]]$reps,rl[[j]]$rs,c(1:person),rl[[j]]$ra,rl[[j]]$rc,rl[[j]]$rsigma,
                 rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i,Alpha, Beta, mu_b, sigma2_b)
    gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat+item) - 2*lbound
    new_a = r0$ra
    new_b = r0$rb
    new_c = r0$rc
    new_s = r0$rs
    eta= r0$reta
    xi=r0$reps
  }
  id=which.min(gic)
  #adjust the range of lambda
  if(id==1){
    lbd=seq(0.1,2,length.out=6)
    rl<-vector("list",length(lbd))
    gic<-NULL
    for(j in 1:length(lbd)){
      r0=sgvem_3PLEFA_L1_const1(u,new_a,new_b,new_c,new_s,eta,xi,Sigma, domain,
                                inite,samp,forgetrate,lbd[j],indic,
                                mu_b,sigma2_b,Alpha,Beta,nopenalty_col,updateProgress=NULL)
      rl [[j]]=saem_3PLCFA(u,domain, r0$Q_mat,samp,forgetrate,mu_b,sigma2_b,Alpha,
                           Beta,updateProgress=NULL)
      lbound=lb3pl(u,rl[[j]]$reps,rl[[j]]$rs,c(1:person),rl[[j]]$ra,rl[[j]]$rc,rl[[j]]$rsigma,
                   rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i,Alpha, Beta, mu_b, sigma2_b)
      gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat+item) - 2*lbound
      new_a = r0$ra
      new_b = r0$rb
      new_c = r0$rc
      new_s = r0$rs
      eta= r0$reta
      xi=r0$reps
    }}
  else if(id==10){
    lbd=seq(20,40,length.out=6)
    rl<-vector("list",length(lbd))
    gic<-NULL
    for(j in 1:length(lbd)){
      r0=sgvem_3PLEFA_L1_const1(u,new_a,new_b,new_c,new_s,eta,xi,Sigma, domain,
                                inite,samp,forgetrate,lbd[j],indic,
                                mu_b,sigma2_b,Alpha,Beta,nopenalty_col,updateProgress=NULL)
      rl [[j]]=saem_3PLCFA(u,domain, r0$Q_mat,samp,forgetrate,mu_b,sigma2_b,Alpha,
                           Beta,updateProgress=NULL)
      lbound=lb3pl(u,rl[[j]]$reps,rl[[j]]$rs,c(1:person),rl[[j]]$ra,rl[[j]]$rc,rl[[j]]$rsigma,
                   rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i,Alpha, Beta, mu_b, sigma2_b)
      gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat+item) - 2*lbound
      new_a = r0$ra
      new_b = r0$rb
      new_c = r0$rc
      new_s = r0$rs
      eta= r0$reta
      xi=r0$reps
    }
  }
  id=which.min(gic)
  rs=rl[[id]]
  rs$lbd=lbd[id]
  rs$id=id
  return(rs)
}


#####lasso_c2_3pl functions######
#non_penalty a
na_lc23pl<-'
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
List nalc23pl(const arma::mat&u,const int& domain,const int& item,
const arma::mat& a,const arma::vec& nopenalty_col,
const arma::vec& lastone,const arma::mat& s, const arma::mat& eta,const arma::vec& new_b,
const arma::cube& SIGMA, const arma::mat& MU,const arma::vec& id_1,const arma::mat& prev_a_num,const arma::cube& prev_a_denom, 
const double& dec_st){
arma::mat new_a=a;
arma::mat a_num=mat(domain,item,fill::zeros);
arma::cube a_denom=zeros(domain,domain,item);
for(int j=0; j<nopenalty_col.n_elem; ++j){
int k=nopenalty_col(j)-1;
int l=lastone(j)-1;
for(int i=0; i<id_1.n_elem; ++i){
int n=id_1(i)-1;
double sigma=SIGMA(l,l,n);
double mu=MU(l,n);
double a1=1-u(n,k)+s(n,k)*u(n,k);
a_denom(l,l,k)=a_denom(l,l,k)+a1*eta(n,k)*sigma+eta(n,k)*(mu*mu);
a_num(l,k)=a_num(l,k)+a1*(u(n,k)-0.5+2*new_b(k)*eta(n,k))*mu;
}
new_a(k,l)=1/(dec_st*a_denom(l,l,k)+(1-dec_st)*prev_a_denom(l,l,k))*(dec_st*a_num(l,k)+(1-dec_st)*prev_a_num(l,k))/2;

}
return List::create(Named("new_a") = new_a,Named("a_denom") = a_denom,Named("a_num") = a_num);
}
'
sourceCpp(code=na_lc23pl)

#update a with penalty: 1:domain, off-diagonal
pa_lc23pl<-'
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
List palc23pl(const arma::mat&u,const arma::mat& a,
const arma::vec& nopenalty_col,const arma::mat& s,
const arma::vec& lastone,const double& lbd,
const arma::vec& id_1, const arma::mat& eta,const arma::vec& new_b,const arma::cube& SIGMA, const arma::mat& MU,
const double& dec_st,const arma::mat& prev_delta,const arma::mat& prev_deriv2){
arma::mat new_a1=a;
int domain=new_a1.n_cols;
int item=u.n_cols;
arma::mat delta=mat(item,domain,fill::zeros);
arma::mat deriv2=mat(item,domain,fill::zeros);
for(int j=0; j<nopenalty_col.n_elem; ++j){
int n=nopenalty_col(j)-1;
int l=lastone(j)-1;
arma::uvec A=regspace<uvec>(0, 1,l-1);
arma::uvec B=regspace<uvec>(l+1, 1, domain-1);
arma::uvec sdf=join_vert(A,B);
for(int k=0; k<sdf.n_elem; ++k){
int m=sdf(k);
for(int p=0; p<id_1.n_elem; ++p){
int i=id_1(p)-1;
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
double a2=1-u(i,n)+s(i,n)*u(i,n);
delta(n,m) = delta(n,m) + a2*(u(i,n)-0.5)*mu(m)+2*new_b(n)*eta(i,n)*mu(m)-2*eta(i,n)*dot(trans(new_a1.submat(id,ex_k)),sigmumu.submat(ex_k,id2));
deriv2(n,m) = deriv2(n,m) + 2*a2*eta(i,n)*sigmumu(m,m);
}
double avg_delta = dec_st*delta(n,m) + (1-dec_st)*prev_delta(n,m);
if(abs(avg_delta)>lbd){
double S=sign(avg_delta)*(abs(avg_delta) - lbd);
new_a1(n,m) = S/(dec_st*deriv2(n,m) + (1-dec_st)*prev_deriv2(n,m));
}else{
new_a1(n,m) = 0;
}
}
}
return List::create(Named("new_a") = new_a1,Named("delta") = delta,Named("deriv2") = deriv2);
}
'
sourceCpp(code=pa_lc23pl)

#update a with penalty: domain+1:item
pa_lc23pl1<-'
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
List palc23pl1(const arma::mat&u,const arma::mat& a,const arma::mat& s,const double& lbd,
const arma::vec& id_1, const arma::mat& eta,const arma::vec& new_b,const arma::cube& SIGMA, const arma::mat& MU,
const double& dec_st,const arma::mat& prev_delta,const arma::mat& prev_deriv2,const arma::vec& pc, 
const arma::mat& delta1, const arma::mat& deriv22){
arma::mat new_a1=a;
int domain=new_a1.n_cols;
arma::mat delta=delta1;
arma::mat deriv2=deriv22;
for(int j=0; j<pc.n_elem; ++j){
int n=pc(j)-1;
for(int m=0; m<domain; ++m){
for(int p=0; p<id_1.n_elem; ++p){
      int i=id_1(p)-1;
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
      double a2=1-u(i,n)+s(i,n)*u(i,n);
      delta(n,m) = delta(n,m) + a2*(u(i,n)-0.5)*mu(m)+2*new_b(n)*eta(i,n)*mu(m)-2*eta(i,n)*dot(trans(new_a1.submat(id,ex_k)),sigmumu.submat(ex_k,id2));
      deriv2(n,m) = deriv2(n,m) + 2*a2*eta(i,n)*sigmumu(m,m);
    }
    double avg_delta = dec_st*delta(n,m) + (1-dec_st)*prev_delta(n,m);
    if(abs(avg_delta)>lbd){
      double S=sign(avg_delta)*(abs(avg_delta) - lbd);
      new_a1(n,m) = S/(dec_st*deriv2(n,m) + (1-dec_st)*prev_deriv2(n,m));
    }else{
      new_a1(n,m) = 0;
    }
}
}
return List::create(Named("new_a") = new_a1,Named("delta") = delta,Named("deriv2") = deriv2);
}
'
sourceCpp(code=pa_lc23pl1)


#lasso with constraint 2 function
sgvem_3PLEFA_lasso_const2 <- function(u,new_a,new_b,new_c,new_s,eta,xi,Sigma,inite,samp,forgetrate,domain,lbd,
                                      mu_b,sigma2_b,Alpha,Beta,indic,nopenalty_col,updateProgress=NULL) {
  person=dim(u)[1]
  item=dim(u)[2]
  converge = 1
  MU=matrix(0,nrow=domain,ncol=person)
  SIGMA=array(diag(domain),dim=c(domain,domain,person))
  #xi=eta
  prev_a_num = matrix(0,nrow=domain,ncol=item)  
  prev_a_denom = array(0,dim=c(domain,domain,item)) 
  prev_delta = matrix(0,nrow=item,ncol=domain)
  prev_deriv2 = matrix(0,nrow=item,ncol=domain);
  prev_b_num = matrix(0,nrow=item,ncol=1)
  prev_b_denom = matrix(0,nrow=item,ncol=1)
  prev_c_num = matrix(0,nrow=item,ncol=1) 
  prev_c_denom = 0
  n = 1
  par_Sigma = Sigma
  #lb=NULL
  #nv<-NULL
  nlb<-n20<-NULL
  prev_lb = 0
  while(converge==1 && rankMatrix(Sigma) == domain && n < 8000){
    if (is.function(updateProgress)) {
      text <- paste0("n=:", n)
      updateProgress(detail = text)
    }
    if(n==1){
      dec_st=1
      id_1 = sample(person,inite)
    }else{
      dec_st = (n+1)^(-forgetrate)*0.25
      #dec_st = (n+1)^(-forgetrate)
      id_1 = sample(person,samp)
    }
    s_num = (1-new_c)/new_c
    par_Sigma = Sigma
    #update MU and SIGMA,xi,Spart
    rs1<-ecfa3pl(u, Sigma, domain, id_1, item, eta, new_s, new_a, new_b, SIGMA, MU, xi)
    SIGMA=rs1$SIGMA
    MU=rs1$MU
    xi=rs1$xi
    Spart=rs1$Spart
    
    #update s
    new_s=s3pl(item,id_1,new_s,eta,new_a,s_num,new_b,u,SIGMA,MU,xi)
    #update eta
    eta=(exp(xi)/(1+exp(xi))-0.5)/(2*xi)
    eta=eta3pl(item,id_1,eta,abs(xi))
    
    Sigma=Spart/length(id_1)
    d_temp=sqrt(diag(diag(Sigma)))
    Sigma=solve(d_temp)%*%Sigma%*%solve(d_temp)
    
    #update b
    par_b=new_b
    rs2=b3pl(item,id_1,new_s,eta,u,new_a,dec_st,MU,prev_b_num,prev_b_denom,mu_b,sigma2_b)
    new_b=rs2$new_b
    b_denom= rs2$b_denom
    b_num = rs2$b_num
    
    
    par_a=new_a
    
    #update a
    #find the last one for each item by using indicator matrix
    lastone=apply(indic[nopenalty_col,], 1, function(x) tail(which(x!=0),1))
    rs3=nalc23pl(u, domain,item,new_a,nopenalty_col,lastone,new_s,eta, 
                 new_b, SIGMA, MU, id_1, prev_a_num, prev_a_denom, dec_st)
    new_a=rs3$new_a
    a_denom=rs3$a_denom
    a_num=rs3$a_num
    #L1-penalty: off-diagnoal
    rs4=palc23pl(u, new_a, nopenalty_col, new_s, lastone, lbd, id_1, eta, new_b, 
                 SIGMA, MU, dec_st, prev_delta, prev_deriv2)
    new_a=rs4$new_a
    delta=rs4$delta
    deriv2=rs4$deriv2
    #upper-tiangular should be zero
    new_a=replace(new_a,indic==0,0)
    #domain+1:item
    #find penaly columns
    pc=setdiff(1:item,nopenalty_col)
    rs5=palc23pl1(u, new_a, new_s, lbd, id_1, eta, new_b, SIGMA, MU, dec_st, 
                  prev_delta, prev_deriv2, pc, delta, deriv2)
    new_a=rs5$new_a
    delta=rs5$delta
    deriv2=rs5$deriv2
    
    #new_a=replace(new_a,new_a< -1,0)
    
    # id_allzero=which(rowSums(new_a)==0)
    # new_a[id_allzero,] = par_a[id_allzero,]
    id_allzero=which(rowSums(new_a)==0)
    if(sum(rowSums(new_a)==0)!=0){
      new_a[id_allzero,] = par_a[id_allzero,]
      break
    }
    
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
    
    prev_delta = (1-dec_st)*prev_delta + dec_st*delta
    prev_deriv2 = (1-dec_st)*prev_deriv2 + dec_st*deriv2
    
    
    prev_lb=(1-dec_st)*prev_lb+ dec_st*lb3pl(u, xi, new_s, id_1, 
                                             new_a, new_c, Sigma, new_b, SIGMA, MU, 
                                             Alpha, Beta, mu_b, sigma2_b)
    nlb<-append(nlb,prev_lb)
    
    if(n%%20==0 && n>20){
      n20=append(n20,abs(mean(nlb[(n-39):(n-20)])-mean(nlb[(n-19):n])))
      if(abs(mean(nlb[(n-39):(n-20)])-mean(nlb[(n-19):n]))<0.006){
        converge=0
      }
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
  #new_a=replace(new_a,abs(new_a)< 0.001,0)
  Q_mat = (new_a != 0)*1
  lbound=lb3pl(u,xi,new_s,c(1:person),new_a,new_c,Sigma, new_b, SIGMA, MU, 
               Alpha, Beta, mu_b, sigma2_b)
  gic=log(log(person))*log(person)*sum(Q_mat+item) - 2*lbound
  #AIC, BIC
  bic = log(person)*sum(Q_mat+item) - 2*lbound
  aic = 2*sum(Q_mat+item) -2*lbound
  return(list(ra=new_a,rb=new_b,rc=new_c,rs = new_s,reta = eta,reps=xi,rsigma = rsigma,rs = new_s,
              mu_i = MU,sig_i = SIGMA,n=n,Q_mat=Q_mat,n20=n20,GIC=gic,AIC=aic,
              BIC=bic))
}

#main function: choose optimal lambda
sgvem_3PLEFA_L1_const2_all<-function(u,domain,samp,forgetrate,
                                     mu_b,sigma2_b,Alpha,Beta,indic,non_pen,updateProgress=NULL){
  lbd=seq(2,20,2)
  inite=dim(u)[1]
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
  new_c=rep(0.15,item)
  new_s=matrix(0,nrow=person,ncol=item)
  aveg_p=colMeans(u==1)
  for (i in 1:person) {
    for (j in 1:item) {
      if(u[i,j]==0){
        new_s[i,j]=1;
      }else{
        new_s[i,j]=runif(1,aveg_p[j],1)
      }
    }
  }
  rl<-vector("list",length(lbd))
  gic<-NULL
  for(j in 1:length(lbd)){
    if (is.function(updateProgress)) {
      text <- paste0("j=:", j)
      updateProgress(detail = text)
    }
    r0=sgvem_3PLEFA_lasso_const2(u,new_a,new_b,new_c,new_s,eta,xi,Sigma,inite,samp,
                                 forgetrate,domain,lbd[j],mu_b,sigma2_b,Alpha,Beta,indic,nopenalty_col,updateProgress=NULL) 
    rl [[j]]=saem_3PLCFA(u,domain, r0$Q_mat,samp,forgetrate,mu_b,sigma2_b,Alpha,
                         Beta,updateProgress=NULL)
    lbound=lb3pl(u,rl[[j]]$reps,rl[[j]]$rs,c(1:person),rl[[j]]$ra,rl[[j]]$rc,rl[[j]]$rsigma,
                 rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i,Alpha, Beta, mu_b, sigma2_b)
    gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat+item) - 2*lbound
    new_a = r0$ra
    new_b = r0$rb
    new_c = r0$rc
    new_s = r0$rs
    eta= r0$reta
    xi=r0$reps
  }
  id=which.min(gic)
  #adjust the range of lambda
  if(id==1){
    lbd=seq(0.1,2,length.out=6)
    rl<-vector("list",length(lbd))
    gic<-NULL
    for(j in 1:length(lbd)){
      r0=sgvem_3PLEFA_lasso_const2(u,new_a,new_b,new_c,new_s,eta,xi,Sigma,inite,samp,
                                   forgetrate,domain,lbd[j],mu_b,sigma2_b,Alpha,Beta,indic,nopenalty_col,updateProgress=NULL) 
      rl [[j]]=saem_3PLCFA(u,domain, r0$Q_mat,samp,forgetrate,mu_b,sigma2_b,Alpha,
                           Beta,updateProgress=NULL)
      lbound=lb3pl(u,rl[[j]]$reps,rl[[j]]$rs,c(1:person),rl[[j]]$ra,rl[[j]]$rc,rl[[j]]$rsigma,
                   rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i,Alpha, Beta, mu_b, sigma2_b)
      gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat+item) - 2*lbound
      new_a = r0$ra
      new_b = r0$rb
      new_c = r0$rc
      new_s = r0$rs
      eta= r0$reta
      xi=r0$reps
    }}
  else if(id==10){
    lbd=seq(20,40,length.out=6)
    rl<-vector("list",length(lbd))
    gic<-NULL
    for(j in 1:length(lbd)){
      r0=sgvem_3PLEFA_lasso_const2(u,new_a,new_b,new_c,new_s,eta,xi,Sigma,inite,samp,
                                   forgetrate,domain,lbd[j],mu_b,sigma2_b,Alpha,Beta,
                                   indic,nopenalty_col,updateProgress=NULL) 
      rl [[j]]=saem_3PLCFA(u,domain, r0$Q_mat,samp,forgetrate,mu_b,sigma2_b,Alpha,
                           Beta,updateProgress=NULL)
      lbound=lb3pl(u,rl[[j]]$reps,rl[[j]]$rs,c(1:person),rl[[j]]$ra,rl[[j]]$rc,rl[[j]]$rsigma,
                   rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i,Alpha, Beta, mu_b, sigma2_b)
      gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat+item) - 2*lbound
      new_a = r0$ra
      new_b = r0$rb
      new_c = r0$rc
      new_s = r0$rs
      eta= r0$reta
      xi=r0$reps
    }
  }
  id=which.min(gic)
  rs=rl[[id]]
  rs$lbd=lbd[id]
  rs$id=id
  return(rs)
}


#####adaptive_lasso_c1_3pl functions######
#update penalty a
pa_al3pl<-'
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
List paal3pl(const arma::mat&u,const int& domain,const int& item,const arma::vec& id_1,const double& lbd, 
const arma::mat& eta,const arma::mat& s,const arma::mat& a,const arma::vec& new_b,const arma::cube& SIGMA, 
const arma::mat& MU,const double& dec_st,const arma::mat& prev_delta,const arma::mat& prev_deriv2,
const arma::mat& weights,const arma::vec& sdf){
arma::mat new_a=a;
arma::mat delta=mat(item,domain,fill::zeros);
arma::mat deriv2=mat(item,domain,fill::zeros);
for(int n=0; n<sdf.n_elem; ++n){
int j=sdf(n)-1;
for(int k=0; k<domain; ++k){
for(int i=0; i<id_1.n_elem; ++i){
int l=id_1(i)-1;
arma::mat sigma=SIGMA.slice(l);
arma::vec mu=MU.col(l);
arma::mat sigmumu=sigma+mu*trans(mu);
arma::uvec A=regspace<uvec>(0, 1,k-1);
arma::uvec B=regspace<uvec>(k+1, 1, domain-1);
arma::uvec ex_k=join_vert(A,B);
double a1=1-u(l,j)+s(l,j)*u(l,j);
arma::uvec id(1);
id.at(0)=j;
arma::uvec id2(1);
id2.at(0)=k;
delta(j,k) = delta(j,k) + a1*((u(l,j)-0.5)*mu(k)+2*new_b(j)*eta(l,j)*mu(k)-2*eta(l,j)*dot(trans(new_a.submat(id,ex_k)),sigmumu.submat(ex_k,id2)));
deriv2(j,k) = deriv2(j,k) + a1*2*eta(l,j)*sigmumu(k,k);
}
double avg_delta = dec_st*delta(j,k) + (1-dec_st)*prev_delta(j,k);
if(abs(avg_delta)>lbd*(1/weights(j,k))){
double S=sign(avg_delta)*(abs(avg_delta) - lbd*(1/weights(j,k)));
new_a(j,k) = S/(dec_st*deriv2(j,k) + (1-dec_st)*prev_deriv2(j,k));
}else{
new_a(j,k) = 0;
}
}
}
return List::create(Named("new_a") = new_a,Named("delta") = delta,Named("deriv2") = deriv2);
}
'
sourceCpp(code=pa_al3pl)


#adaptive lasso with constraint 1 function
sgvem_3PLEFA_adapt_const1 <- function(u,new_a,new_b,new_c,new_s,eta,xi,Sigma,inite,samp,forgetrate,domain,lbd,indic,
                                      mu_b,sigma2_b,Alpha,Beta,weights,nopenalty_col,updateProgress=NULL) {
  person=dim(u)[1]
  item=dim(u)[2]
  converge = 1
  # new_s=matrix(0,nrow=person,ncol=item)
  # aveg_p=colMeans(u==1)
  # for (i in 1:person) {
  #   for (j in 1:item) {
  #     if(u[i,j]==0){
  #       new_s[i,j]=1;
  #     }else{
  #       new_s[i,j]=runif(1,aveg_p[j],1)
  #     }
  #   }
  # }
  MU=matrix(0,nrow=domain,ncol=person)
  SIGMA=array(diag(domain),dim=c(domain,domain,person))
  #xi=eta
  prev_a_num = matrix(0,nrow=domain,ncol=item)  
  prev_a_denom = array(0,dim=c(domain,domain,item)) 
  prev_delta = matrix(0,nrow=item,ncol=domain)
  prev_deriv2 = matrix(0,nrow=item,ncol=domain);
  prev_b_num = matrix(0,nrow=item,ncol=1)
  prev_b_denom = matrix(0,nrow=item,ncol=1)
  prev_c_num = matrix(0,nrow=item,ncol=1) 
  prev_c_denom = 0
  n = 1
  par_Sigma = Sigma
  #lb=NULL
  #nv<-NULL
  nlb<-n20<-NULL
  prev_lb = 0
  while(converge==1 && rankMatrix(Sigma) == domain && n < 8000){
    if (is.function(updateProgress)) {
      text <- paste0("n=:", n)
      updateProgress(detail = text)
    }
    if(n==1){
      dec_st=1
      id_1 = sample(person,inite)
    }else{
      dec_st = (n+1)^(-forgetrate)*0.25
      #dec_st = (n+1)^(-forgetrate)
      id_1 = sample(person,samp)
    }
    s_num = (1-new_c)/new_c
    par_Sigma = Sigma
    #update MU and SIGMA,xi,Spart
    rs1<-ecfa3pl(u, Sigma, domain, id_1, item, eta, new_s, new_a, new_b, SIGMA, MU, xi)
    SIGMA=rs1$SIGMA
    MU=rs1$MU
    xi=rs1$xi
    Spart=rs1$Spart
    
    #update s
    new_s=s3pl(item,id_1,new_s,eta,new_a,s_num,new_b,u,SIGMA,MU,xi)
    #update eta
    eta=(exp(xi)/(1+exp(xi))-0.5)/(2*xi)
    eta=eta3pl(item,id_1,eta,abs(xi))
    
    Sigma=Spart/length(id_1)
    d_temp=sqrt(diag(diag(Sigma)))
    Sigma=solve(d_temp)%*%Sigma%*%solve(d_temp)
    
    #update b
    par_b=new_b
    rs2=b3pl(item,id_1,new_s,eta,u,new_a,dec_st,MU,prev_b_num,prev_b_denom,mu_b,sigma2_b)
    new_b=rs2$new_b
    b_denom= rs2$b_denom
    b_num = rs2$b_num
    
    
    par_a=new_a
    
    #update a
    rs3=nalc13pl(u, indic, nopenalty_col, eta, new_s, new_b, SIGMA, MU, new_a, id_1, prev_a_num, prev_a_denom, dec_st)
    new_a1=rs3$new_a
    new_a1[-nopenalty_col,]=new_a[-nopenalty_col,]
    a_denom=rs3$a_denom
    a_num=rs3$a_num
    #L1-penalty
    sdf=setdiff(1:item,nopenalty_col)
    rs4=paal3pl(u,domain, item, id_1, lbd, eta, new_s, new_a1, new_b, SIGMA, MU, dec_st, 
                prev_delta, prev_deriv2,weights,sdf)
    new_a=rs4$new_a
    delta=rs4$delta
    deriv2=rs4$deriv2
    
    #new_a=replace(new_a,new_a< -1,0)
    
    id_allzero=which(rowSums(new_a)==0)
    if(sum(rowSums(new_a)==0)!=0){
      new_a[id_allzero,] = par_a[id_allzero,]
      break
    }
    #new_a[id_allzero,] = par_a[id_allzero,]
    
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
    
    prev_delta = (1-dec_st)*prev_delta + dec_st*delta
    prev_deriv2 = (1-dec_st)*prev_deriv2 + dec_st*deriv2
    
    
    prev_lb=(1-dec_st)*prev_lb+ dec_st*lb3pl(u, xi, new_s, id_1, 
                                             new_a, new_c, Sigma, new_b, SIGMA, MU, 
                                             Alpha, Beta, mu_b, sigma2_b)
    nlb<-append(nlb,prev_lb)
    
    if(n%%20==0 && n>20){
      n20=append(n20,abs(mean(nlb[(n-39):(n-20)])-mean(nlb[(n-19):n])))
      if(abs(mean(nlb[(n-39):(n-20)])-mean(nlb[(n-19):n]))<0.006){
        converge=0
      }
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
  lbound=lb3pl(u,xi,new_s,c(1:person),new_a,new_c,Sigma, new_b, SIGMA, MU, 
               Alpha, Beta, mu_b, sigma2_b)
  gic=log(log(person))*log(person)*sum(Q_mat+item) - 2*lbound
  #AIC, BIC
  bic = log(person)*sum(Q_mat+item) - 2*lbound
  aic = 2*sum(Q_mat+item) -2*lbound
  return(list(ra=new_a,rb=new_b,rc=new_c,rs = new_s,reta = eta,reps=xi,rsigma = Sigma,rs = new_s,
              mu_i = MU,sig_i = SIGMA,n=n,Q_mat=Q_mat,n20=n20,GIC=gic,AIC=aic,BIC=bic))
}


#main function: choose optimal lambda
sgvem_3PLEFA_adaptive_const1_all<-function(u,domain,indic,samp,forgetrate,mu_b,sigma2_b,Alpha,Beta,gamma,updateProgress=NULL){
  inite=dim(u)[1]
  person=dim(u)[1]
  item=dim(u)[2]
  lbd=seq(2,20,2)
  nopenalty_col=which(rowSums(indic)==1)
  #initialization
  initial=init(u,domain,indic)
  new_a = initial[[1]]
  new_b = initial[[2]]
  eta= initial[[3]]
  xi=initial[[4]]
  Sigma=initial[[5]]
  new_c=rep(0.15,item)
  new_s=matrix(0,nrow=person,ncol=item)
  aveg_p=colMeans(u==1)
  for (i in 1:person) {
    for (j in 1:item) {
      if(u[i,j]==0){
        new_s[i,j]=1;
      }else{
        new_s[i,j]=runif(1,aveg_p[j],1)
      }
    }
  }
  #create weights: use CFA
  wa=saem_3PLCFA(u,domain,indic,samp,forgetrate,
                 mu_b,sigma2_b,Alpha,Beta,updateProgress=NULL)$ra
  weights = abs(wa)^gamma+1e-05
  rl<-vector("list",length(lbd))
  gic<-NULL
  for(j in 1:length(lbd)){
    if (is.function(updateProgress)) {
      text <- paste0("j=:", j)
      updateProgress(detail = text)
    }
    rl [[j]]=sgvem_3PLEFA_adapt_const1(u,new_a,new_b,new_c,new_s,eta,xi,Sigma,
                                       inite,samp,forgetrate,domain,lbd[j],indic,
                                       mu_b,sigma2_b,Alpha,Beta,weights,nopenalty_col,updateProgress=NULL)
    lbound=lb3pl(u,rl[[j]]$reps,rl[[j]]$rs,c(1:person),rl[[j]]$ra,rl[[j]]$rc,rl[[j]]$rsigma,
                 rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i,Alpha, Beta, mu_b, sigma2_b)
    gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat+item) - 2*lbound
    new_a = rl[[j]]$ra
    new_b = rl[[j]]$rb
    new_c = rl[[j]]$rc
    new_s = rl[[j]]$rs
    eta= rl[[j]]$reta
    xi=rl[[j]]$reps
  }
  id=which.min(gic)
  #adjust the range of lambda
  if(id==1){
    lbd=seq(0.1,2,length.out=6)
    rl<-vector("list",length(lbd))
    gic<-NULL
    for(j in 1:length(lbd)){
      rl [[j]]=sgvem_3PLEFA_adapt_const1(u,new_a,new_b,new_c,new_s,eta,xi,Sigma,
                                         inite,samp,forgetrate,domain,lbd[j],indic,
                                         mu_b,sigma2_b,Alpha,Beta,weights,nopenalty_col,updateProgress=NULL)
      lbound=lb3pl(u,rl[[j]]$reps,rl[[j]]$rs,c(1:person),rl[[j]]$ra,rl[[j]]$rc,rl[[j]]$rsigma,
                   rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i,Alpha, Beta, mu_b, sigma2_b)
      gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat+item) - 2*lbound
      new_a = rl[[j]]$ra
      new_b = rl[[j]]$rb
      new_c = rl[[j]]$rc
      new_s = rl[[j]]$rs
      eta= rl[[j]]$reta
      xi=rl[[j]]$reps
    }}
  else if(id==10){
    lbd=seq(20,40,length.out=6)
    rl<-vector("list",length(lbd))
    gic<-NULL
    for(j in 1:length(lbd)){
      rl [[j]]=sgvem_3PLEFA_adapt_const1(u,new_a,new_b,new_c,new_s,eta,xi,Sigma,
                                         inite,samp,forgetrate,domain,lbd[j],indic,
                                         mu_b,sigma2_b,Alpha,Beta,weights,nopenalty_col,updateProgress=NULL)
      lbound=lb3pl(u,rl[[j]]$reps,rl[[j]]$rs,c(1:person),rl[[j]]$ra,rl[[j]]$rc,rl[[j]]$rsigma,
                   rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i,Alpha, Beta, mu_b, sigma2_b)
      gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat+item) - 2*lbound
      new_a = rl[[j]]$ra
      new_b = rl[[j]]$rb
      new_c = rl[[j]]$rc
      new_s = rl[[j]]$rs
      eta= rl[[j]]$reta
      xi=rl[[j]]$reps
    }
  }
  id=which.min(gic)
  rs=saem_3PLCFA(u,domain, rl[[id]]$Q_mat,samp,forgetrate,mu_b,sigma2_b,Alpha,
                 Beta,updateProgress=NULL)
  rs$lbd=lbd[id]
  rs$id=id
  return(rs)
}


#####adaptive_lasso_c2_3pl functions######
#update a with penalty: 1:domain, off-diagonal
pa_alc23pl<-'
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
List paalc23pl(const arma::mat&u,const arma::mat& a,
const arma::vec& nopenalty_col,const arma::mat& s,
const arma::vec& lastone,const double& lbd,
const arma::vec& id_1, const arma::mat& eta,const arma::vec& new_b,const arma::cube& SIGMA, const arma::mat& MU,
const double& dec_st,const arma::mat& prev_delta,const arma::mat& prev_deriv2,const arma::mat&weights){
arma::mat new_a1=a;
int domain=new_a1.n_cols;
int item=u.n_cols;
arma::mat delta=mat(item,domain,fill::zeros);
arma::mat deriv2=mat(item,domain,fill::zeros);
for(int j=0; j<nopenalty_col.n_elem; ++j){
int n=nopenalty_col(j)-1;
int l=lastone(j)-1;
arma::uvec A=regspace<uvec>(0, 1,l-1);
arma::uvec B=regspace<uvec>(l+1, 1, domain-1);
arma::uvec sdf=join_vert(A,B);
for(int k=0; k<sdf.n_elem; ++k){
int m=sdf(k);
for(int p=0; p<id_1.n_elem; ++p){
int i=id_1(p)-1;
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
double a2=1-u(i,n)+s(i,n)*u(i,n);
delta(n,m) = delta(n,m) + a2*(u(i,n)-0.5)*mu(m)+2*new_b(n)*eta(i,n)*mu(m)-2*eta(i,n)*dot(trans(new_a1.submat(id,ex_k)),sigmumu.submat(ex_k,id2));
deriv2(n,m) = deriv2(n,m) + 2*a2*eta(i,n)*sigmumu(m,m);
}
double avg_delta = dec_st*delta(n,m) + (1-dec_st)*prev_delta(n,m);
if(avg_delta>lbd*(1/weights(n,m))){
double S=sign(avg_delta)*(abs(avg_delta) - lbd*(1/weights(n,m)));
new_a1(n,m) = S/(dec_st*deriv2(n,m) + (1-dec_st)*prev_deriv2(n,m));
}else{
new_a1(n,m) = 0;
}
}
}
return List::create(Named("new_a") = new_a1,Named("delta") = delta,Named("deriv2") = deriv2);
}
'
sourceCpp(code=pa_alc23pl)

#update a with penalty: domain+1:item
pa_alc23pl1<-'
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
List paalc23pl1(const arma::mat&u,const arma::mat& a,const arma::mat& s,const double& lbd,
const arma::vec& id_1, const arma::mat& eta,const arma::vec& new_b,const arma::cube& SIGMA, const arma::mat& MU,
const double& dec_st,const arma::mat& prev_delta,const arma::mat& prev_deriv2,const arma::vec& pc, 
const arma::mat& delta1, const arma::mat& deriv22,const arma::mat&weights){
arma::mat new_a1=a;
int domain=new_a1.n_cols;
arma::mat delta=delta1;
arma::mat deriv2=deriv22;
for(int j=0; j<pc.n_elem; ++j){
int n=pc(j)-1;
for(int m=0; m<domain; ++m){
for(int p=0; p<id_1.n_elem; ++p){
      int i=id_1(p)-1;
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
      double a2=1-u(i,n)+s(i,n)*u(i,n);
      delta(n,m) = delta(n,m) + a2*(u(i,n)-0.5)*mu(m)+2*new_b(n)*eta(i,n)*mu(m)-2*eta(i,n)*dot(trans(new_a1.submat(id,ex_k)),sigmumu.submat(ex_k,id2));
      deriv2(n,m) = deriv2(n,m) + 2*a2*eta(i,n)*sigmumu(m,m);
    }
    double avg_delta = dec_st*delta(n,m) + (1-dec_st)*prev_delta(n,m);
    if(avg_delta>lbd*(1/weights(n,m))){
    double S=sign(avg_delta)*(abs(avg_delta) - lbd*(1/weights(n,m)));
      new_a1(n,m) = S/(dec_st*deriv2(n,m) + (1-dec_st)*prev_deriv2(n,m));
    }else{
      new_a1(n,m) = 0;
    }
}
}
return List::create(Named("new_a") = new_a1,Named("delta") = delta,Named("deriv2") = deriv2);
}
'
sourceCpp(code=pa_alc23pl1)



#adaptive lasso with constraint 2 function
sgvem_3PLEFA_adapt_const2 <- function(u,new_a,new_b,new_c,new_s,eta,xi,Sigma,inite,samp,forgetrate,domain,lbd,
                                      mu_b,sigma2_b,Alpha,Beta,weights,indic,nopenalty_col,updateProgress=NULL) {
  person=dim(u)[1]
  item=dim(u)[2]
  converge = 1
  MU=matrix(0,nrow=domain,ncol=person)
  SIGMA=array(diag(domain),dim=c(domain,domain,person))
  #xi=eta
  prev_a_num = matrix(0,nrow=domain,ncol=item)  
  prev_a_denom = array(0,dim=c(domain,domain,item)) 
  prev_delta = matrix(0,nrow=item,ncol=domain)
  prev_deriv2 = matrix(0,nrow=item,ncol=domain);
  prev_b_num = matrix(0,nrow=item,ncol=1)
  prev_b_denom = matrix(0,nrow=item,ncol=1)
  prev_c_num = matrix(0,nrow=item,ncol=1) 
  prev_c_denom = 0
  n = 1
  par_Sigma = Sigma
  #lb=NULL
  #nv<-NULL
  nlb<-n20<-NULL
  prev_lb = 0
  while(converge==1 && rankMatrix(Sigma) == domain && n < 8000){
    if (is.function(updateProgress)) {
      text <- paste0("n=:", n)
      updateProgress(detail = text)
    }
    if(n==1){
      dec_st=1
      id_1 = sample(person,inite)
    }else{
      dec_st = (n+1)^(-forgetrate)*0.25
      #dec_st = (n+1)^(-forgetrate)
      id_1 = sample(person,samp)
    }
    s_num = (1-new_c)/new_c
    par_Sigma = Sigma
    #update MU and SIGMA,xi,Spart
    rs1<-ecfa3pl(u, Sigma, domain, id_1, item, eta, new_s, new_a, new_b, SIGMA, MU, xi)
    SIGMA=rs1$SIGMA
    MU=rs1$MU
    xi=rs1$xi
    Spart=rs1$Spart
    
    #update s
    new_s=s3pl(item,id_1,new_s,eta,new_a,s_num,new_b,u,SIGMA,MU,xi)
    #update eta
    eta=(exp(xi)/(1+exp(xi))-0.5)/(2*xi)
    eta=eta3pl(item,id_1,eta,abs(xi))
    
    Sigma=Spart/length(id_1)
    d_temp=sqrt(diag(diag(Sigma)))
    Sigma=solve(d_temp)%*%Sigma%*%solve(d_temp)
    
    #update b
    par_b=new_b
    rs2=b3pl(item,id_1,new_s,eta,u,new_a,dec_st,MU,prev_b_num,prev_b_denom,mu_b,sigma2_b)
    new_b=rs2$new_b
    b_denom= rs2$b_denom
    b_num = rs2$b_num
    
    
    par_a=new_a
    
    #update a
    #find the last one for each item by using indicator matrix
    lastone=apply(indic[nopenalty_col,], 1, function(x) tail(which(x!=0),1))
    rs3=nalc23pl(u, domain,item,new_a,nopenalty_col,lastone,new_s,eta, 
                 new_b, SIGMA, MU, id_1, prev_a_num, prev_a_denom, dec_st)
    new_a=rs3$new_a
    a_denom=rs3$a_denom
    a_num=rs3$a_num
    #L1-penalty: off-diagnoal
    rs4=paalc23pl(u, new_a, nopenalty_col, new_s, lastone, lbd, id_1, eta, new_b, 
                  SIGMA, MU, dec_st, prev_delta, prev_deriv2,weights)
    new_a=rs4$new_a
    delta=rs4$delta
    deriv2=rs4$deriv2
    
    #upper-tiangular should be zero
    new_a=replace(new_a,indic==0,0)
    #domain+1:item
    #find penaly columns
    pc=setdiff(1:item,nopenalty_col)
    rs5=paalc23pl1(u, new_a, new_s, lbd, id_1, eta, new_b, SIGMA, MU, dec_st, 
                   prev_delta, prev_deriv2, pc, delta, deriv2,weights)
    new_a=rs5$new_a
    delta=rs5$delta
    deriv2=rs5$deriv2
    
    #new_a=replace(new_a,new_a< -1,0)
    
    id_allzero=which(rowSums(new_a)==0)
    if(sum(rowSums(new_a)==0)!=0){
      new_a[id_allzero,] = par_a[id_allzero,] 
      break
    }
    
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
    
    prev_delta = (1-dec_st)*prev_delta + dec_st*delta
    prev_deriv2 = (1-dec_st)*prev_deriv2 + dec_st*deriv2
    
    
    prev_lb=(1-dec_st)*prev_lb+ dec_st*lb3pl(u, xi, new_s, id_1, 
                                             new_a, new_c, Sigma, new_b, SIGMA, MU, 
                                             Alpha, Beta, mu_b, sigma2_b)
    nlb<-append(nlb,prev_lb)
    
    if(n%%20==0 && n>20){
      n20=append(n20,abs(mean(nlb[(n-39):(n-20)])-mean(nlb[(n-19):n])))
      if(abs(mean(nlb[(n-39):(n-20)])-mean(nlb[(n-19):n]))<0.006){
        converge=0
      }
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
  lbound=lb3pl(u,xi,new_s,c(1:person),new_a,new_c,Sigma, new_b, SIGMA, MU, 
               Alpha, Beta, mu_b, sigma2_b)
  gic=log(log(person))*log(person)*sum(Q_mat+item) - 2*lbound
  #AIC, BIC
  bic = log(person)*sum(Q_mat+item) - 2*lbound
  aic = 2*sum(Q_mat+item) -2*lbound
  return(list(ra=new_a,rb=new_b,rc=new_c,rs = new_s,reta = eta,reps=xi,rsigma = Sigma,rs = new_s,
              mu_i = MU,sig_i = SIGMA,n=n,Q_mat=Q_mat,n20=n20,GIC=gic,AIC=aic,
              BIC=bic))
}


#main function: choose optimal lambda
sgvem_3PLEFA_adaptive_const2_all<-function(u,domain,samp,forgetrate,
                                           mu_b,sigma2_b,Alpha,Beta,indic,non_pen,gamma,updateProgress=NULL){
  lbd=seq(2,20,2)
  inite=dim(u)[1]
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
  new_c=rep(0.15,item)
  new_s=matrix(0,nrow=person,ncol=item)
  aveg_p=colMeans(u==1)
  for (i in 1:person) {
    for (j in 1:item) {
      if(u[i,j]==0){
        new_s[i,j]=1;
      }else{
        new_s[i,j]=runif(1,aveg_p[j],1)
      }
    }
  }
  #create weights: use CFA
  wa=saem_3PLCFA(u,domain,indic,samp,forgetrate,
                 mu_b,sigma2_b,Alpha,Beta,updateProgress=NULL)$ra
  weights = abs(wa)^gamma+1e-05
  rl<-vector("list",length(lbd))
  gic<-NULL
  for(j in 1:length(lbd)){
    if (is.function(updateProgress)) {
      text <- paste0("j=:", j)
      updateProgress(detail = text)
    }
    rl [[j]]=sgvem_3PLEFA_adapt_const2(u,new_a,new_b,new_c,new_s,eta,xi,Sigma,inite,samp,
                                       forgetrate,domain,lbd[j],mu_b,sigma2_b,Alpha,Beta,
                                       weights,indic,nopenalty_col,updateProgress=NULL) 
    lbound=lb3pl(u,rl[[j]]$reps,rl[[j]]$rs,c(1:person),rl[[j]]$ra,rl[[j]]$rc,rl[[j]]$rsigma,
                 rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i,Alpha, Beta, mu_b, sigma2_b)
    gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat+item) - 2*lbound
    new_a = rl[[j]]$ra
    new_b = rl[[j]]$rb
    new_c = rl[[j]]$rc
    new_s = rl[[j]]$rs
    eta= rl[[j]]$reta
    xi=rl[[j]]$reps
  }
  id=which.min(gic)
  #adjust the range of lambda
  if(id==1){
    lbd=seq(0.1,2,length.out=6)
    rl<-vector("list",length(lbd))
    gic<-NULL
    for(j in 1:length(lbd)){
      rl [[j]]=sgvem_3PLEFA_adapt_const2(u,new_a,new_b,new_c,new_s,eta,xi,Sigma,inite,samp,
                                         forgetrate,domain,lbd[j],mu_b,sigma2_b,Alpha,Beta,
                                         weights,indic,nopenalty_col,updateProgress=NULL) 
      lbound=lb3pl(u,rl[[j]]$reps,rl[[j]]$rs,c(1:person),rl[[j]]$ra,rl[[j]]$rc,rl[[j]]$rsigma,
                   rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i,Alpha, Beta, mu_b, sigma2_b)
      gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat+item) - 2*lbound
      new_a = rl[[j]]$ra
      new_b = rl[[j]]$rb
      new_c = rl[[j]]$rc
      new_s = rl[[j]]$rs
      eta= rl[[j]]$reta
      xi=rl[[j]]$reps
    }}
  else if(id==10){
    lbd=seq(20,40,length.out=6)
    rl<-vector("list",length(lbd))
    gic<-NULL
    for(j in 1:length(lbd)){
      rl [[j]]=sgvem_3PLEFA_adapt_const2(u,new_a,new_b,new_c,new_s,eta,xi,Sigma,inite,samp,
                                         forgetrate,domain,lbd[j],mu_b,sigma2_b,Alpha,Beta,
                                         weights,indic,nopenalty_col,updateProgress=NULL) 
      lbound=lb3pl(u,rl[[j]]$reps,rl[[j]]$rs,c(1:person),rl[[j]]$ra,rl[[j]]$rc,rl[[j]]$rsigma,
                   rl[[j]]$rb,rl[[j]]$sig_i,rl[[j]]$mu_i,Alpha, Beta, mu_b, sigma2_b)
      gic[j]=log(log(person))*log(person)*sum(rl[[j]]$Q_mat+item) - 2*lbound
      new_a = rl[[j]]$ra
      new_b = rl[[j]]$rb
      new_c = rl[[j]]$rc
      new_s = rl[[j]]$rs
      eta= rl[[j]]$reta
      xi=rl[[j]]$reps
    }
  }
  id=which.min(gic)
  rs=saem_3PLCFA(u,domain, rl[[id]]$Q_mat,samp,forgetrate,mu_b,sigma2_b,Alpha,
                 Beta,updateProgress=NULL)
  rs$lbd=lbd[id]
  rs$id=id
  return(rs)
}



#####shiny app####
options(shiny.maxRequestSize = 30*1024^2)
# Define UI for data upload app ----
ui <- navbarPage("VEMIRT-3PL Shiny App",
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
                              
                              # Horizontal line ----
                              tags$hr(),
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
                              
                              tags$hr(),
                              sliderInput("fr1", "Step size", min = 0.5, max = 1.0, value = 0.51),
                              numericInput("sample1", "Subsample", min=1, max=1000, value=10),
                              numericInput("bmu2", "Prior Distribution for b (Mean)", min=-3, max=3, value=0),
                              numericInput("bsig2", "Prior Distribution for b (Variance)", min=0.1, max=20, value=1),
                              numericInput("alpha2", "Prior Distribution for g (Alpha)", min=0, max=50, value=2),
                              numericInput("beta2", "Prior Distribution for g (Beta)", min=0, max=50, value=5),
                              
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
                 #####CFA
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
                              
                              
                              tags$hr(),
                              sliderInput("fr2", "Step size", min = 0.5, max = 1.0, value = 0.6),
                              
                              numericInput("sample2", "Subsample", min=1, max=1000, value=10),
                              numericInput("bmu1", "Prior Distribution for b (Mean)", min=-3, max=3, value=0),
                              numericInput("bsig1", "Prior Distribution for b (Variance)", min=0.1, max=20, value=1),
                              numericInput("alpha1", "Prior Distribution for g (Alpha)", min=0, max=50, value=2),
                              numericInput("beta1", "Prior Distribution for g (Beta)", min=0, max=50, value=5),
                              
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
  url2 <- a("M2PL estimation", href="https://www.google.com/",)
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
      result<-sgvem_3PLEFA(u(),domain(),input$sample1,input$fr1,input$bmu2,
                           input$bsig2,input$alpha2,
                           input$beta2,updateProgress)}
    else if(input$method == "AL" && input$cons=="c1"){
      indic<-reactive({req(input$file5)
        df2 <- data.matrix(read.csv(input$file5$datapath))
        return(df2)
      })
      gamma=reactive({req(input$gm)})
      
      result<-sgvem_3PLEFA_adaptive_const1_all(u(),domain(),indic(),
                                               input$sample1,input$fr1,input$bmu2,
                                               input$bsig2,input$alpha2,
                                               input$beta2,gamma(),updateProgress)
    }
    else if(input$method == "AL" && input$cons=="c2"){
      indic<-reactive({req(input$file5)
        df2 <- data.matrix(read.csv(input$file5$datapath))
        return(df2)
      })
      gamma=reactive({req(input$gm)})
      non_pen=reactive({req(input$np)})
      result<-sgvem_3PLEFA_adaptive_const2_all(u(),domain(),input$sample1,input$fr1,input$bmu2,
                                               input$bsig2,input$alpha2,
                                               input$beta2,indic(),non_pen(),
                                               gamma(),updateProgress)}
    else if(input$method == "La" && input$cons=="c1"){
      indic<-reactive({req(input$file5)
        df2 <- data.matrix(read.csv(input$file5$datapath))
        return(df2)
      })
      result<-sgvem_3PLEFA_L1_const1_all(u(),domain(),indic(),input$sample1,input$fr1,input$bmu2,
                                         input$bsig2,input$alpha2,
                                         input$beta2,updateProgress)}
    else if(input$method == "La" && input$cons=="c2"){
      indic<-reactive({req(input$file5)
        df2 <- data.matrix(read.csv(input$file5$datapath))
        return(df2)
      })
      non_pen=reactive({req(input$np)})
      result<-sgvem_3PLEFA_L1_const2_all(u(),domain(),input$sample1,input$fr1,input$bmu2,
                                         input$bsig2,input$alpha2,
                                         input$beta2,indic(),non_pen(),updateProgress)
    }
    return(result)
  }
  )
  
  output$warn<-renderText({
    input$go1
    isolate({
      if(input$method != "Rot" && (result0 ()$lbd == 0.1 || result0 ()$lbd == 40)){
        return("Warning: The optimal penalty parameter may be out of range")
      }else{
        return(NULL)
      }})
  })
  
  output$par1<-renderTable({
    input$go1
    isolate({
      m<-cbind(result0()$ra,result0()$rb,result0()$rc)
      colnames(m)<-c(paste("a",1:domain(),sep=""),"b","g")
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
  # output$downloadData <- downloadHandler(
  #   # filename = function() {
  #   #   paste(input$mod, "GVEMresults.rds", sep = "")
  #   # },
  #   filename = "GVEMEFAresults.rds",
  #   content = function(file) {
  #     saveRDS(result0(), file)
  #   }
  # )
  
  output$downloadData <- downloadHandler(
    filename = function() {
      if(input$checkGroup1=="all"){
        paste("GVEMEFA",input$checkGroup1 ,"results.rds",sep="")}else{
          paste("GVEMEFA",input$checkGroup1 ,"results.csv",sep="")}
      
    },
    content = function(file) {
      if(input$checkGroup1=="all"){
        saveRDS(result0(), file)}else if(input$checkGroup1=="item"){
          m<-cbind(result0()$ra,result0()$rb,result0()$rc)
          colnames(m)<-c(paste("a",1:domain(),sep=""),"b","g") 
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
    result<-saem_3PLCFA(u1(),domain1(),indic1(),input$sample2,input$fr2,input$bmu1,
                        input$bsig1,input$alpha1,
                        input$beta1,updateProgress)
    return(result)
  })
  
  output$par2<-renderTable({
    input$go3
    isolate({
      m<-cbind(result()$ra,result()$rb,result()$rc)
      colnames(m)<-c(paste("a",1:domain1(),sep=""),"b","g") 
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