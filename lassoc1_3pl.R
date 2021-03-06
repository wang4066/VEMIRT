##################################################
#   Function: GVEM_M3PL_Lasso + Constraint 1    ##
##################################################
# Inputs:
# u: response matrix, a N by J numerical matrix
# domain: K, numerical value of number of latent dimension
# indic: factor loading indicator matrix, a J by K numerical matrix, 
#       reflecting each item loads on which domain
# samp: numerical value of subsample for each iteration
# forgetrate: forget rate for the stochastic algorithm
# mu_b, sigma2_b: the mean and variance parameters, prior distribution of b parameters
# Alpha, Beta: the alpha and beta parameters, prior distribution of g parameters
#################################################
##### Outputs a list of updated parameters                                                             
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
##### GIC,AIC,BIC : model fit index
##### lbd: numerical value of penalty parameter lambda
##### id: numerical value of the position of lambda
##################################################
library(Rcpp)
library(Matrix)
library(psych)
library(gtools)
source("cfa_3pl.r")

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
