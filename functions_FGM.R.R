#####################
##### Libraries #####
#####################
library(MASS)
library(invgamma)
library(Directional)
library(movMF)
library(Rfast)
library(truncnorm)

#####################
##### Functions #####
#####################
# Function to generate samples from Dirichlet distribution
rDR=function(no,prob)
{
  generated_DR=sapply(1:no,function(w){comp=sapply(prob,function(u){return(rgamma(1,shape=u,rate=1))});
  comp=comp/sum(comp); return(comp)})
  return(generated_DR)
}

# Function to calculate norm of a vector
norm_vec=function(a)
{
  return(sqrt(sum(a^2)))
}

# Function to compute density of vMF distribution 
dvmf=function(x,mu,tau) # x and mu are vector of same dimension, and tau is a scalar
{
  D=length(mu)
  mBf=besselI(x=tau,nu=((D/2)-1), expon.scaled = FALSE)
  Cp_tau=mBf*((tau/(2*pi))^((D/2)-1))/(2*pi)
  dvMF=Cp_tau*(exp(tau*(t(mu)%*%x)))
  return(dvMF)
}

# Newton correction to MLE of tau
newton_corr=function(k_hat,D,R_bar)
{
  Apk=besselI(k_hat,nu=((D/2)-1),expon.scaled=FALSE)
  k1_hat=k_hat-((Apk-R_bar)/(1-Apk^2-((D-1)/k_hat)*Apk))
  return(k1_hat)
}


ln_vmf_c=function(tau,d) 
{
  res=(d/2-1)*log(tau)-log(besselI(tau,nu=(d/2-1),expon.scaled = TRUE)) # exp(-\tau)*I_{(d/2-1)}(\tau)
  return(res)
}

# This function calculates FG-kernel density given the parameters 
spca_kern=function(x,cntr,rd,sigma.sq,mu,tau)
{ 
  d=length(x); z=(x-cntr)
  dinom=norm_vec(tau*mu+(rd*z/sigma.sq))
  ln_h_x_theta=ln_vmf_c(tau,d)-d*log(2*pi*sigma.sq)/2-ln_vmf_c(dinom,d)-((norm_vec(z))^2+rd^2-2*sigma.sq*dinom)/(2*sigma.sq)-tau
  return(exp(ln_h_x_theta))
}

# Main function to calculate FG mixture results 
FG_mixture=function(x,K,M,sigma_hat,Iter,tau_def,alpha,mm)
{
  #@param: mm: Additional parameters indicating the number of extra labels assigned (m in Neal Alg 8)
  #@param: x:data; K:#spheres; M:#vMF karnels; Iter:#iterations; tau_def:initial value of tau;
  
  # Set parameters
  sig_frac=0 #@param: sig_frac: a proportion of iterations below which we do not update sigma;  
  samp_IS=15   #@param: samp_IS: minimum # samples required in a class to update tau in independent sampler;
  Iter2=2000 #@param: Iter2:#iterations for MH update of tau; sigma_hat:the maximum permissible value of sigma;
  sigma_1=NA #@param: sigma_1:detault value of sigma before sig_frac updates;
  
  if(missing(Iter)){Iter=1000;}
  if(missing(M)){M=5;}
  if(missing(K)){K=floor(nrow(x)/M);}
  if(missing(alpha)){alpha=1}
  if(missing(samp_IS)){samp_IS=15}
  if(missing(sig_frac)){sig_frac=0}  
  if(missing(Iter2)){Iter2=2000;}
  if(missing(mm)){mm=10;}
  if(missing(tau_def)){tau_def=1;}
  n=nrow(x); D=ncol(x);
  x=as.matrix(x)
  # estimating sigma  
  if(is.na(sigma_1))
  {
    library(FNN)
    KNN=get.knn(x,k=2)
    MEAN=apply((KNN$nn.dist^2),1,mean)
    sigma_1=2*mean((KNN$nn.dist[,1])^2)
  }
  if(missing(sigma_hat))
  {
    library(FNN)
    KNN=get.knn(x,k=2)
    MEAN=apply((KNN$nn.dist^2),1,mean)
    sigma_hat=10*mean(MEAN)
  }
  
  #####################
  ##### Initialize ####
  #####################
  # 1. sigma
  sig1=1; sig2=0.01
  sigma.sq=rinvgamma(1,shape=sig1,rate=sig2)
  sigma.sq=sigma_1  
  # 2. Weights 
  if(missing(K)){K=30}
  if(missing(M)){M=5}
  w_vec=rep(1/K,K) # hyperpriors on lambda
  Lmbd=rDR(1,w_vec)  # Lambda: weights of the von-Mises & hyperpriors on pi
  Pi=pi_vec=matrix(data=NA,ncol=M,nrow=K)
  for(i in 1:K)
  {
    pi_vec[i,]=rep(1/M,M)
    Pi[i,]=rDR(1,pi_vec[i,]) # Pi : weights of the circles
  }
  # 3. Centers and redii for the spheres (k-circles)
  cntr=matrix(data=NA,nrow=D,ncol=K) # ith column gives center  of ith circle
  rd=NULL
  x.quant=quantile(x[,1],c(1:K)/K) 
  quant=sapply(1:K,function(k){return(which((x[,1]>(x.quant[k]-.15))&(x[,1]<(x.quant[k]+.15)))[1])}) 
  quant[missing(quant)]=sample(1:nrow(x),length(which(missing(quant)==TRUE)),replace=FALSE )
  mu_cntr=x[quant,]
  sgm_cntr=rep(1,K)
  a_rd=rep(1,K)
  b_rd=rep(5,K)
  for(i in 1:K)
  {
    cntr[,i]=mvrnorm(1,mu_cntr[i,],sgm_cntr[i]*diag(D))
    rd[i]=rtruncnorm(1,a=0,b=Inf,mean=a_rd[i],sd=b_rd[i])
  }
  rand=sample(1:K,K,replace=FALSE)
  cntr=cntr[,rand]
  # 4. VonMises parameters
  m_mu=array(data=(-1/sqrt(D)),c(D,K,M)) # hyperparameters for the centers of vMF
  k_mu=0.1
  mu=array(data=NA,c(D,K,M))   # centers of the different vMF distribution
  tau=matrix(data=0.1,nrow=K,ncol=M) # radii of the different vMF distribution (fixed at 0.5)
  for(i in 1:K)
  {
    for(j in 1:M)
    {
      mu[,i,j]=rvmf(1,m_mu[,i,j],k_mu*tau[i,j])
    }
  }
  if(missing(tau_def)){tau_def=10}
  tau=matrix(data=tau_def,nrow=K,ncol=M)
  
  sigma_c=1; sigma_r=5; mu_r=1 
  k_hat=1;  # k_hat is 'b' in the paper
  ############################################
  ############ Latent variables ##############
  ############################################
  Y=matrix(data=NA,nrow=nrow(x),ncol=ncol(x)) # Y=vvT(y-c)/r; and x_i=y_i+e_i
  W=matrix(data=NA,nrow=nrow(x),ncol=K) # Indicator variable of Sphere levels
  Z=matrix(data=NA,nrow=nrow(x),ncol=M) # Indicator variable of vMF levels
  
  #$$$$$$$$$$$$$$$$$$$$$$$$$$ add values of W, Z, Y $$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  Y=rvmf(n,mu[,1,1],k=0.001)
  sigma.trace=NULL
  weights=matrix(data=NA,nrow=n,ncol=K)
  for(i in 1:n)
  {
    weights=sapply(1:K,function(u){return(Lmbd[u]*dmvnorm(x[i,],rd[u]*Y[i,]+cntr[,u],sigma=sigma.sq*diag(D)))})
    if(sum(weights)==0)
    {
      dst=apply(cntr,2,function(u){return(norm_vec(x[i,]-u))})
      W[i,]=0; W[i,which.min(dst)]=1;
    }
    else
    {
      if(sum(weights)==0)
      {
        W[i,]=rep(0,ncol(W))
        W[i,which.max(weights)]=1
      }
      else
      {
        weights=weights/sum(weights)
        W[i,]=rmultinom(n=1,size=1,prob = weights)  
      }
    }
  }
  drop=which(colSums(W)==0)
  if(length(drop)>0)
  {
    W=W[,-drop]; mu=mu[,-drop,]; tau=tau[-drop,]; cntr=cntr[,-drop]; rd=rd[-drop]; Pi=Pi[-drop,]; 
  }
  K=ncol(W)
  wgt_new=colSums(W)
  w_vec_new=wgt_new+1/K   
  Lmbd=rDR(1,w_vec_new)

  if(missing(mm)){mm=floor((n-K*M-1)/M)}
  ################################################################################
  for(i in 1:n)
  {
    l_use=which(W[i,]==1) 
    x_pr_sp=Y[i,]
    weights=sapply(1:M,function(m)
    {return(Pi[l_use,m]*dvmf(x_pr_sp,mu=mu[,l_use,m],tau=tau[l_use,m]))})
    if(sum(weights)==0)
    {
      Z[i,]=rep(0,ncol(Z))
      Z[i,which.max(weights)]=1
    }
    else
    {
      weights=weights/sum(weights)
      Z[i,]=rmultinom(n=1,size=1,prob = weights)
    }
  }
  for(l in 1:K)
  {
    idx=which(W[,l]==1)
    l_idx=length(idx)
    if(l_idx>0)
    {
      if(l_idx==1) {wgt_new=Z[idx,]}
      else
      {wgt_new=colSums(Z[idx,])}
      pi_vec_new=wgt_new+1/M    #pi_vec_new=wgt_new+pi_vec[l,]  
      Pi[l,]=rDR(1,pi_vec_new) # Pi : weights of the circles
    }
  } 
  
  ##############
  IDX=matrix(data=NA,nrow=n,ncol=2)
  IDX[,1]=sapply(1:n,function(u){which(W[u,]==1)})
  IDX[,2]=sapply(1:n,function(u){which(Z[u,]==1)})
  
  for(i in 1:n)
  {
    idx1_y=IDX[i,1]; idx2_y=IDX[i,2];
    c_i=cntr[,idx1_y]; r_i=rd[idx1_y];
    tau_i=tau[idx1_y,idx2_y]; mu_i=mu[,idx1_y,idx2_y];
    comp=(r_i*(x[i,]-c_i)/sigma.sq)+tau_i*mu_i
    tau_y=norm_vec(comp)
    mu_y=comp/tau_y
    Y[i,]=rvmf(1,mu_y,tau_y)
  }
  
  tau_pr_vec=seq(.1,20,0.25)  
  likelihood_chain=NULL
  ######################################
  
  ###########################################################################
  ############### List of extra paramters for NEal Algorithm 8 ##############
  ###########################################################################
  N_COMP=n
  CNTR_extra=matrix(data=(rnorm((N_COMP*D),mean=0,sd=sigma_c)),nrow=D)   
  RD_extra=rtruncnorm(N_COMP,a=0,b=Inf,mean=mu_r,sd=sigma_r)
  MU_extra=array(data=NA,c(D,N_COMP,M))   # centers of the different vMF distribution
  TAU_extra=matrix(data=tau_def,nrow=N_COMP,ncol=M)
  for(i in 1:N_COMP)
  {
    for(j in 1:M)
    {
      MU_extra[,i,j]=rvmf(1,m_mu[,1,j],k_mu*TAU_extra[i,j]) # mu is uniformly spread, but tau will be strict
    }
  }
  PI_extra=matrix(data=NA,N_COMP,M)
  for(i in 1:N_COMP)
  {
    PI_extra[i,]=as.vector(rDR(1,rep((1/M),M)) )
  }
  ######################
  ##### Iterations #####
  ######################
  if(missing(Iter)){Iter=1500} # no. of iterations
  for(iter in 1:Iter)
  {
    if((iter%%50)==0) 
    {   print(iter); print(date())   }
    
    ###################################################################################
    #print("update Lambda")
    for(i in 1:n)
    {
      idx_0=which(W[i,]==1)
      if(length(which(W[,idx_0]==1))==1)
      {
        k_neg=(ncol(W)-1);
        new_labs=(k_neg+2):(k_neg+mm); n_comp=length(new_labs)
        cntr_extra=matrix(data=NA,nrow=D,ncol=(n_comp+1)); rd_extra=NULL;
        mu_extra=array(data=NA,dim=c(D,(n_comp+1),M)); tau_extra=matrix(data=NA,nrow=(n_comp+1),ncol=M)
        Pi_extra=matrix(data=NA,nrow=(n_comp+1),ncol=M)
        cntr_extra[,1]=cntr[,idx_0]; rd_extra[1]=rd[idx_0]; 
        mu_extra[,1,]=mu[,idx_0,]; tau_extra[1,]=tau[idx_0,]; Pi_extra[1,]=Pi[idx_0,]; 
        rand1=sample(1:N_COMP,n_comp,replace=FALSE); cntr_extra[,2:(n_comp+1)]=CNTR_extra[,rand1] 
        rand1=sample(1:N_COMP,n_comp,replace=FALSE); rd_extra[2:(n_comp+1)]=RD_extra[rand1]
        rand1=sample(1:N_COMP,n_comp,replace=FALSE); mu_extra[,2:(n_comp+1),]=MU_extra[,rand1,]
        rand1=sample(1:N_COMP,n_comp,replace=FALSE); tau_extra[2:(n_comp+1),]=TAU_extra[rand1,]
        rand1=sample(1:N_COMP,n_comp,replace=FALSE); Pi_extra[2:(n_comp+1),]=PI_extra[rand1,]
        n_comp=n_comp+1
        
        W=as.matrix(W[,-idx_0]);
        K=ncol(W); 
        cntr=cntr[,-idx_0]; rd=rd[-idx_0]; mu.temp=mu[,-idx_0,]; tau=tau[-idx_0,]; Pi=Pi[-idx_0,]
        cntr=as.matrix(cntr);
        if(K==1)
        {
          tau=t(as.matrix(tau));
          Pi=t(as.matrix(Pi)); 
          mu.temp=array(data=NA,dim=c(D,1,M))
          mu.temp[,1,]=mu[,-idx_0,]
        }
        mu=mu.temp
        idx_1=NULL
        idx_1[1:K]=sapply(1:K,function(l){sum(W[,l])*sum(sapply(1:M,function(m){return(Pi[l,m]*spca_kern(x[i,],cntr[,l],rd[l],sigma.sq,mu[,l,m],tau[l,m]))}))})
        idx_1[(K+1):(K+n_comp)]=(alpha/mm)*sapply(1:n_comp,function(l)
          {sum(sapply(1:M,function(m){return(Pi_extra[l,m]*spca_kern(x[i,],cntr_extra[,l],rd_extra[l],sigma.sq,mu_extra[,l,m],tau_extra[l,m]))}))})
        if(sum(idx_1)==0)
        {
          weights=rep(0,length(idx_1))
          weights[which.max(idx_1)]=1
        }
        else
        {
          idx_1=idx_1/sum(idx_1)
          weights=rmultinom(n=1,size=1,prob=idx_1)
        }
        idx_2=which(weights==1)
        if(idx_2>K)
        { 
          temp=rep(0,n); temp[i]=1;
          W=cbind(W,temp)
          colnames(W)=NULL
          cntr=cbind(cntr,cntr_extra[,(idx_2-K)])
          rd=c(rd,rd_extra[(idx_2-K)])
          tau=rbind(tau,tau_extra[(idx_2-K),])
          mu.temp=array(data=NA,c(D,(K+1),M))
          mu.temp[,-(K+1),]=mu
          mu.temp[,(K+1),]=mu_extra[,(idx_2-K),]
          mu=mu.temp
          pi_vec_new=rep(1/M,M)
          Pi=rbind(Pi,t(rDR(1,pi_vec_new)))
          K=K+1;
        }
        else
        {
          W[i,idx_2]=1
        }
      }
      else
      {
        k_neg=ncol(W)
        new_labs=(k_neg+1):(k_neg+mm); n_comp=length(new_labs)
        rand1=sample(1:N_COMP,n_comp,replace=FALSE); cntr_extra=CNTR_extra[,rand1] 
        rand1=sample(1:N_COMP,n_comp,replace=FALSE); rd_extra=RD_extra[rand1]
        rand1=sample(1:N_COMP,n_comp,replace=FALSE); mu_extra=MU_extra[,rand1,]
        rand1=sample(1:N_COMP,n_comp,replace=FALSE); tau_extra=TAU_extra[rand1,]
        rand1=sample(1:N_COMP,n_comp,replace=FALSE); Pi_extra=PI_extra[rand1,]
        
        idx_1=NULL
        idx_1[1:K]=sapply(1:K,function(l){
        if(l==idx_0)
        {return((sum(W[,l])-1)*sum(sapply(1:M,function(m){return(Pi[l,m]*spca_kern(x[i,],cntr[,l],rd[l],sigma.sq,mu[,l,m],tau[l,m]))})))}
        else
        {return(sum(W[,l])*sum(sapply(1:M,function(m){return(Pi[l,m]*spca_kern(x[i,],cntr[,l],rd[l],sigma.sq,mu[,l,m],tau[l,m]))})))}
        })  
        if(sum(idx_1)==0)
        {
          weights=rep(0,length(idx_1))
          weights[which.max(idx_1)]=1
        }
        else
        {
          idx_1=idx_1/sum(idx_1)
          weights=rmultinom(n=1,size=1,prob =idx_1) 
        }
        idx_2=which(weights==1)
        if(idx_2>K)
        {
          W[i,idx_0]=0
          temp=rep(0,n); temp[i]=1;
          W=cbind(W,temp)
          colnames(W)=NULL;
          cntr=cbind(cntr,cntr_extra[,(idx_2-K)])
          rd=c(rd,rd_extra[(idx_2-K)])
          tau=rbind(tau,tau_extra[(idx_2-K),])
          mu.temp=array(data=NA,c(D,(K+1),M))
          mu.temp[,-(K+1),]=mu
          mu.temp[,(K+1),]=mu_extra[,(idx_2-K),]
          mu=mu.temp
          pi_vec_new=rep(1/M,M)
          Pi=rbind(Pi,t(rDR(1,pi_vec_new)))
          K=K+1;
        }        
        else
        {
          W[i,]=rep(0,ncol(W))
          W[i,idx_2]=1
        }
      }
    }
    wgt_new=colSums(W)
    Lmbd=rDR(1,wgt_new)

    
    ################################################################################
#    print("update Pi")
    for(i in 1:n)
    {
      l_use=which(W[i,]==1) ; 
       weights=sapply(1:M,function(m)
       {return(Pi[l_use,m]*spca_kern(x[i,],cntr[,l_use],rd[l_use],sigma.sq,mu[,l_use,m],tau[l_use,m])  )})
       if(sum(weights)==0)
       {
         weights=rep(0,length(weights))
         weights[which.max(weights)]=1
       }
       else
       {
         weights=weights/sum(weights)    
       }
     Z[i,]=rmultinom(n=1,size=1,prob = weights)
    }
    for(l in 1:K)
    {
      idx=which(W[,l]==1)
      l_idx=length(idx)
      if(l_idx>0)
      {
        if(l_idx==1) {wgt_new=Z[idx,]}
        else
        {wgt_new=colSums(Z[idx,])}
        pi_vec_new=wgt_new+pi_vec[l,]  # This is the updated lambda
        Pi[l,]=rDR(1,pi_vec_new) # Pi : weights of the circles
      }
    } 
    
    ###############################################################################
    # Update Y #
    IDX=matrix(data=NA,nrow=n,ncol=2)
    IDX[,1]=sapply(1:n,function(u){which(W[u,]==1)})
    IDX[,2]=sapply(1:n,function(u){which(Z[u,]==1)})
    
    for(i in 1:n)
    {
      idx1_y=IDX[i,1]; idx2_y=IDX[i,2];
      c_i=cntr[,idx1_y]; r_i=rd[idx1_y];
      tau_i=tau[idx1_y,idx2_y]; mu_i=mu[,idx1_y,idx2_y];
      comp=(r_i*(x[i,]-c_i)/sigma.sq)+tau_i*mu_i
      tau_y=norm_vec(comp)
      mu_y=comp/tau_y
      Y[i,]=rvmf(1,mu_y,tau_y)
    }
    
    ###############################################################################
    # update sigma
    if(iter>(sig_frac*Iter)) #500; do not update sigma at all
    {
      sig.shape=sig1+n*D/2;
      rate.add=0
      for(l in 1:K)
      {
        idx=which(IDX[,1]==l)
        l_idx=length(idx)
        if(l_idx>0)
        {
          rate.add=rate.add+sum(sapply(idx,function(u){
            return(t(x[u,]-rd[l]*Y[u,]-cntr[,l])%*%(x[u,]-rd[l]*Y[u,]-cntr[,l]))}))
        }
      }
      sig.rate=sig2+rate.add/2
      sigma.sq=rinvgamma(1,shape=sig.shape,rate=sig.rate)
      sigma.trace[iter]=sigma.sq  
      if(sigma.sq>sigma_hat){sigma.sq=sigma_hat} 
    }  
    
    ###############################################################################
    # update c
    #print("update C, r, mu")
    for(l in 1:K)
    {
      idx=which(IDX[,1]==l)
      l_idx=length(idx)
      cntr.var=((l_idx/sigma.sq)+(1/sigma_c))^(-1)
      if(l_idx==0) {cntr.mean=rep(0,D)}
      if(l_idx==1) {cntr.mean=(x[idx,]-rd[l]*Y[idx,])/sigma.sq}
      if(l_idx>1)
      {
        cntr.mean=colSums(x[idx,]-rd[l]*Y[idx,])/sigma.sq  
      }
      cntr[,l]=mvrnorm(1,(cntr.var*cntr.mean),cntr.var*diag(D))
      
      # update r
      if(l_idx>0)
      { #rd[l]=mean(sapply(idx,function(u){return(norm_vec(x[u,]-cntr[,l]))})) }
        rd_var=(((1/sigma_r)+(l_idx/sigma.sq)))^(-1)
        if(l_idx==0){rd_mean=mu_r/sigma_r}
        if(l_idx>0){
          rd_mean=mu_r/sigma_r + 
            sum(sapply(idx,function(ws){return(t(Y[ws,])%*%(x[ws,]-cntr[,l]))}))/sigma.sq }
        rd[l]=rtruncnorm(1,a=0,b=Inf,mean=(rd_var*rd_mean),sd=sqrt(rd_var))
      }
      
      # update mu
      if(l_idx>0)
      {
        for(m in 1:M)
        {
          #print("update MU")
          idx2=idx[which(IDX[idx,2]==m)]
          l_idx2=length(idx2)
          if(l_idx2>0)
          {
          # Updates on location parameter
            if(l_idx2==1)
            {
              comp=Y[idx2,]
            }
            else
            {
              comp=colSums(Y[idx2,])
            }
            comp1=(k_mu*m_mu[,1,m]+comp) #k_mu is the parameter b in the text
            location=tau[l,m]*comp1
            scale=norm_vec(location)
            mu[,l,m]=rvmf(1,location/scale,scale)
            
            # print("update Tau")
            # Estimating tau
            if(iter>2)
            {
              if(l_idx2<samp_IS) 
              {
                tau[l,m]=tau[l,m]	
              }
              else
              {
                R_bar=norm_vec(colSums(Y[idx2,]))/l_idx2
                k_hat=R_bar*(D-R_bar^2)/(1-R_bar^2)
                if(k_hat<200){
                  for(iter2 in 1:10)
                  {
                    k1_hat=newton_corr(k_hat,D,R_bar)
                    if(k1_hat>.2){k_hat=k1_hat}
                  }
                }
                if(k_hat>200){ k_hat=200 }
                a=1
                mu_0=rep((1/sqrt(D)),D)
                comp0=mu[,l,m]%*%(k_mu*mu_0+comp)
                count=0; 
                z=rgamma(1,shape=2,rate=2/k_hat) #abs(comp0)
                for(iter2 in 1:Iter2)
                {
                  z_new=rgamma(1,shape=2,rate=(2/k_hat))
                  log_fznew_fz=((l_idx2+a)*(((D/2-1))*(log(z_new)-log(z))+log(besselI(x=z,nu=((D/2)-1),expon.scaled = FALSE)) 
                                            - log(besselI(x=z_new,nu=((D/2)-1),expon.scaled = FALSE))   )+(z_new-z)*comp0)
                  log_qz_qznew=log(dgamma(z,shape=2,rate=(2/k_hat)))-log(dgamma(z_new,shape=2,rate=(2/k_hat)))
                  
                  accept_prob=min(exp(log_fznew_fz+log_qz_qznew),1)
                  if(rbinom(1,size=1,prob=accept_prob)==1)
                  {
                    z=z_new; count=count+1
                  }
                }
                tau[l,m]=z
              }
            }
          }
        } 
      }
    }
    tau[which(tau>100,arr.ind=TRUE)]=100
    
    ## Likelihood chain
    if(iter>0)
    {
      fx0=NULL
      Lmbd=sapply(1:K,function(u){sum(W[,u])})
      Lmbd=Lmbd/sum(Lmbd)
      for(i in 1:nrow(x))
      {
        x0=x[i,]
        fvmf=NULL
        for(l in 1:K)
        {
          fvmf[l]=sum(sapply(1:M,function(m){return(Pi[l,m]*spca_kern(x0,cntr[,l],rd[l],sigma.sq,mu[,l,m],tau[l,m]))}))
        }
        fx0[i]=t(fvmf)%*%Lmbd
      }
      likelihood_chain[iter]=sum(log(fx0))
    }
  }
  y=matrix(data=NA,nrow=nrow(Y),ncol=ncol(Y))
  for(i in 1:n)
  {
    idx=which(W[i,]==1)
    y[i,]=Y[i,]*rd[idx]+cntr[,idx]
  }
  par(mfrow=c(1,2))
  par(mar=c(4,4,4,4))
  plot(x,lwd=3,main="Data points",cex=.3)
  plot(y,cex=.3,lwd=3,main="Predicted values")
  Lmbd=sapply(1:K,function(u){sum(W[,u])})
  Lmbd=Lmbd/sum(Lmbd)
  Results=list(Lmbd,Pi,cntr,rd,mu,tau,sigma.trace,likelihood_chain)
  names(Results)=c("Sp_wgts","Kern_wgts","Centers","Radius","mu","tau","trace_sigma","likelihood_chain")
  return(Results)
}

# Given the FG-mixture results and a matrix of test samples, this function provides density of the test sample under FG-mixture approach
get_density_FGM=function(x1_te,Res1)
{
  x1_te=as.matrix(x1_te)
  cntr=Res1$Centers; rd=Res1$Radius; mu=Res1$mu; tau=Res1$tau; Lmbd=Res1$Sp_wgts; Pi=Res1$Kern_wgts;
  fx0=NULL ; K=nrow(tau); M=ncol(tau); D=ncol(x1_te); ss=Res1$trace_sigma[length(Res1$trace_sigma)];
  #  if(ss<0.01){ss=0.01} ##################### Ask about this step ##########
  for(i in 1:nrow(x1_te))
  {
    x0=x1_te[i,]
    fvmf=NULL
    for(l in 1:K)
    {
      fvmf[l]=sum(sapply(1:M,function(m){return(Pi[l,m]*spca_kern(x0,cntr[,l],rd[l],ss,mu[,l,m],tau[l,m]))}))
    }
    fx0[i]=t(fvmf)%*%Lmbd
  }
  return(fx0)
}  

# Given the FG-mixture results and a number N, this function generates N samples from estimated density, and also provides scatter plots with and withour error of the same
plot_density_FGM=function(N,Res1)
{
  cntr=Res1$Centers; rd=Res1$Radius; mu=Res1$mu; tau=Res1$tau; Lmbd=Res1$Sp_wgts; Pi=Res1$Kern_wgts; 
  K=nrow(tau); M=ncol(tau);  D=nrow(cntr); sigma.sq=Res1$trace_sigma[length(Res1$trace_sigma)]
  z=matrix(data=NA,nrow=N,ncol=D)
  sp_label=sapply(1:N,function(u){v=rmultinom(n=1,size=1,prob = Lmbd); return(which.max(v))})
  count=0
  for(ii in 1:K)
  {
    idx1=which(sp_label==ii); 
    if(length(idx1)>0)	
    {
      kr_label=sapply(1:length(idx1),function(u){v=rmultinom(n=1,size=1,prob=Pi[ii,]); return(which.max(v))})
      for(jj in 1:M)
      {
        idx2=idx1[which(kr_label==jj)]; 
        if(length(idx2)>0)
        {
          yy=rvmf(length(idx2),mu[,ii,jj],tau[ii,jj])
          if(length(idx2)==1)
          {
            z[idx2,]=cntr[,ii]+rd[ii]*yy
          }
          else
          {
            z[idx2,]=t(apply(yy,1,function(u){return(cntr[,ii]+rd[ii]*u)})) 
          }
        }
      }
    }
  }
  error=mvrnorm(N,rep(0,D),sigma.sq*diag(D))
  w=z+error
  if(D==2)
  {
    par(mfrow=c(1,2))
    plot(z[,1],z[,2],main="Without error",lwd=3,cex=.3)
    plot(w[,1],w[,2],main="Predictive",lwd=3,cex=.3)
  }
  if(D==3)
  {
    library(plot3D)
    scatter3D(z[,1],z[,2],z[,3],main="Without error",lwd=3,cex=.3)
    scatter3D(w[,1],w[,2],w[,3],main="Predictive",lwd=3,cex=.3)
  }
  return(list(z,w))
}
