###############################################################################################################
################################################# README ######################################################
###############################################################################################################

The R-file FG_functions.R contains all the functions required to estimate density using FG-mixture approach.
The four main functions are:

(1) spca_kern: This function calculates FG-kernel density at a given point. Inputs to this function
are: 
(i) x: point at which the density is to be estimated;
(ii) cntr (center);
(iii) rd (radius);
(iv) sigma.sq (sigma^2); 
(v) mu; 
(vi)tau.
The output is the FG density at point x given the parameters.  


(2) FG_mixture : This is the main function which estimates density given a data matrix. Inputs are 
the following parameters:
(i) x: data in nXd format
(ii) K: The initial number of spheres (default=#data-points/M); 
(iii) M: number of kernels within each sphere (default=5); 
(iv) Iter: number of iterations (default=1000); 
(v) sigma_hat: an upper-bound for sigma^2 (default sigma_hat is estimated from the data); 
(vi)tau_def : A default value for vMF precision parameter tau (default=1); 
(vii) alpha: The precision parameter of DP (default=1); 
(viii) mm: Number of empty spheres (see Neal Algorithm 8, default=10).  

The output of this function contains of list of variable as follows:
(a) Sp_wgts: A vector containing the sphere weights;
(b) Kern_wgts: A matrix containing kernel weights;
(c) Centers: A matrix containing centers of the spheres;
(d) Radius: A vector containing radii of the spheres;
(e) mu: A three-dimensional array containing centers of the vMF kernels;
(f) tau: A matrix containing precision parameters of the vMF kernels;
(g) trace_sigma: The chain of sigma^2 values;
(h) likelihood_chain: The likelihood chain. 


(3) get_density_FGM:  This function provides estimated density of a set of test points using FG-mixture
approach. The inputs are:
(a) x_te: A matrix of test data-points in nXd format;
(b) Res1: A list obtained as the output of the function FG_mixture.

The output is the vector of estimated densities using FG-kernel approach.

(4) plot_density_FGM: This function provides sample points from predictive density using FG-mixture
approach. The inputs are:
(a) N: number of sample points to be generated;
(b) Res1: A list obtained as the output of the function FG_mixture.

The output is a matrix of simulated sample points from the predictive density in NXd format. The scatter 
plot of the generated sample points with and without adding the error component (estimated sigma^2) is also shown.


##################################################################################################################

We demonstrate the use of the functions using the following example:  
## Generating data:
library(ider)
x_train=gendata(DataName = "SShape", n = 750, noise =0.01, seed=1)
x_test=gendata(DataName = "SShape", n = 300, noise =0.01, seed=100)
library(plot3D)
scatter3D(x_train[,1],x_train[,2],x_train[,3])

## Density estimation using FG-mixture
source('functions_FGM.R.R')
FGM=FG_mixture(x_train,M=3,Iter=1000,alpha=1,mm=10) # It takes time
## Checking the results:
FGM$FGM$Sp_wgts
FGM$Kern_wgts
plot(FGM$trace_sigma)
plot(FGM$likelihood_chain)
FGM$Centers
FGM$Radius


## Scatter plot of sample points from predictive density 
predicted_values=plot_density_FGM(N=nrow(x_train),FGM)
pred_without_error=predicted_values[[1]]
pred_with_error=predicted_values[[2]]
# compare
scatter3D(x_train[,1],x_train[,2],x_train[,3],col="black",
	cex.symbols = .1,lwd=3,theta=0,phi =0,xlim=c(-1.1,1.1),ylim=c(-1.1,1.1),zlim=c(-1.1,1.1))
scatter3D(pred_with_error[,1],pred_with_error[,2],pred_with_error[,3],col="black",
	cex.symbols = .1,lwd=3,theta=0,phi =0,xlim=c(-1.1,1.1),ylim=c(-1.1,1.1),zlim=c(-1.1,1.1))


## Get estimated density of test points:
FG_density=get_density_FGM(x_test,FGM)

##################################################################################################################