##################################################
# Function: GVEM_M2PL_Adapt Lasso + Constraint 1##
##################################################
# Inputs:
# u: response matrix, a N by J numerical matrix
# domain: K, numerical value of number of latent dimension
# indic: factor loading indicator matrix, a J by K numerical matrix, 
#       reflecting each item loads on which domain
# gamma: numerical value of adaptive lasso parameter
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
##### lbd: numerical value of penalty parameter lambda
##### id: numerical value of the position of lambda  
##################################################
library(Rcpp)
library(Matrix)
library(psych)
library(gtools)
source("cfa_2pl.r")

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