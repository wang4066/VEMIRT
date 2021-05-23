##################################################
# Function: GVEM_M3PL_Adapt Lasso + Constraint 2##
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
##### Q_mat: Q-matrix to indicate the loading structure,  a J by K matrix
##### GIC: numerical value of GIC
##### lbd: numerical value of penalty parameter lambda
##### id: numerical value of the position of lambda
##################################################

library(Rcpp)
library(Matrix)
library(psych)
library(gtools)
source("cfa_3pl.r")
#update a
al_c23pl<-'
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
List alc23pl(const arma::mat&u,const arma::mat& a,const double& lbd,const arma::mat&weights,
const arma::mat& eta,const arma::vec& new_b,const arma::cube& SIGMA, const arma::mat& s,
const arma::mat& MU,const arma::vec& id_1,arma::mat& prev_a_num,const arma::cube& prev_a_denom, 
const double& dec_st,const arma::mat& prev_delta,const arma::mat& prev_deriv2){
int domain=a.n_cols;
int item=a.n_rows;
arma::mat new_a=a;
arma::mat a_num=mat(domain,item,fill::zeros);
arma::cube a_denom=zeros(domain,domain,item);
arma::mat delta=mat(item,domain,fill::zeros);
arma::mat deriv2=mat(item,domain,fill::zeros);
for(int j=0; j<item; ++j){
for(int k=0; k<domain; ++k){
if(k==j){
for(int i=0; i<id_1.n_elem; ++i){
int l=id_1(i)-1;
double sigma=SIGMA(k,k,l);
double mu=MU(k,l);
double a1=1-u(l,j)+s(l,j)*u(l,j);
a_denom(k,k,j)=a_denom(k,k,j)+a1*eta(l,j)*(sigma+(mu*mu));
a_num(k,j)=a_num(k,j)+a1*(u(l,j)-0.5+2*new_b(j)*eta(l,j))*mu;
}
new_a(j,k)=1/(dec_st*a_denom(k,k,j)+(1-dec_st)*prev_a_denom(k,k,j))*(dec_st*a_num(k,j)+(1-dec_st)*prev_a_num(k,j))/2;
}else if(k>j){
new_a(j,k)=0;
}else{
for(int i=0; i<id_1.n_elem; ++i){
int l=id_1(i)-1;
arma::mat sigma=SIGMA.slice(l);
arma::vec mu=MU.col(l);
arma::mat sigmumu=sigma+mu*trans(mu);
arma::uvec A=regspace<uvec>(0, 1,k-1);
arma::uvec B=regspace<uvec>(k+1, 1, domain-1);
arma::uvec ex_k=join_vert(A,B);
arma::uvec id(1);
id.at(0)=j;
arma::uvec id2(1);
id2.at(0)=k;
double a2=1-u(l,j)+s(l,j)*u(l,j);
delta(j,k) = delta(j,k) + a2*((u(l,j)-0.5)*mu(k)+2*new_b(j)*eta(l,j)*mu(k)-2*eta(l,j)*dot(trans(new_a.submat(id,ex_k)),sigmumu.submat(ex_k,id2)));
deriv2(j,k) = deriv2(j,k) + a2*2*eta(l,j)*sigmumu(k,k);
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
}
return List::create(Named("new_a") = new_a,Named("a_denom") = a_denom,Named("a_num") = a_num,
Named("delta") = delta,Named("deriv2") = deriv2);
}
'
sourceCpp(code=al_c23pl)



#adaptive lasso with constraint 2 function
sgvem_3PLEFA_adapt_const2 <- function(u,new_a,new_b,new_c,new_s,eta,xi,Sigma,inite,samp,forgetrate,domain,lbd,
                                      mu_b,sigma2_b,Alpha,Beta,weights,updateProgress=NULL) {
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
    rs3=alc23pl(u, new_a, lbd, weights, eta, new_b, SIGMA, new_s, MU, id_1, 
             prev_a_num, prev_a_denom, dec_st, prev_delta, prev_deriv2)
    new_a=rs3$new_a
    a_denom=rs3$a_denom
    a_num=rs3$a_num
    delta=rs3$delta
    deriv2=rs3$deriv2
    
    
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
  return(list(ra=new_a,rb=new_b,rc=new_c,rs = new_s,reta = eta,reps=xi,rsigma = Sigma,rs = new_s,
              mu_i = MU,sig_i = SIGMA,n=n,Q_mat=Q_mat,n20=n20))
}


#main function: choose optimal lambda
sgvem_3PLEFA_adaptive_const2_all<-function(u,domain,samp,forgetrate,
                                           mu_b,sigma2_b,Alpha,Beta,updateProgress=NULL){
  gamma=2
  lbd=seq(2,20,2)
  inite=dim(u)[1]
  person=dim(u)[1]
  item=dim(u)[2]
  #create a indicator matrix satisfying the constraint 2
  a1=diag(domain)
  a1[lower.tri(a1,diag=FALSE)]=1
  indic=rbind(a1,matrix(1,nrow=(dim(u)[2]-domain),ncol=domain))
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
                                       weights,updateProgress=NULL) 
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
                                         weights,updateProgress=NULL) 
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
                                         weights,updateProgress=NULL) 
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
