\name{JMdesign-package}
\alias{JMdesign-package}
\alias{JMdesign}
\docType{package}
\title{
Joint Modeling of Longitudinal and Survival Data - Power Calculation
}
\description{
R function for power determination in joint modeling of 
longitudinal and survival data with k-th order trajectories and unknown 
variance-covariance matrix Sigma_theta.
}
\details{
\tabular{ll}{
Package: \tab JMdesign\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2014-10-21\cr
License: \tab GPL-2 \cr
}
The package contains the R-function \code{powerLongSurv} to perform power 
calculations for joint modeling of longitudinal and survival data when trajectories are 
of k-th order and the variance-covariance matrix Sigma_theta is unknown. 
}
\author{
Emil A. Cornea, Liddy M. Chen, Bahjat F. Qaqish, Haitao Chu, and Joseph G. Ibrahim

Maintainer: Shannon T. Holloway <sthollow@ncsu.edu>
}
\references{
L. M. Chen, J. G. Ibrahim, and H. Chu. Sample size and power determination
in joint modeling of longitudinal and survival data.
Statist. Med. 2011, 30 2295-2309
}

\keyword{ package }
\seealso{
   \code{\link{powerLongSurv}}, \code{\link{powerLongSurv-class}}, \code{\link{show-methods}}
}
\examples{
## Example 1.
## **********
## Input elements of Sigma_theta in forumula 4.6;
SigmaTheta <- matrix(c(1.2,0.0,0.0,0.0,0.7,0.0,0.0,0.0,0.8),nrow=3,ncol=3)

N        <-  200; # Total sample size;
nevents  <-  140; # Number of events;
tmedian  <-  0.7; # median survival;
meantf   <-  1.4; # mean follow-up time;
beta     <-  0.2; # Effect of the trajectory;
alpha    <-  0.05;# Type-I Error (2-sided);
sigmae_2 <-  0.09; # measurement error;

## schedule of measurement;
t <- c(0.4, 0.8, 1.2, 1.6, 2) ; # maximum 2 year follow-up;

## Input estimated proportion subjects with 2,3,4,5,6 measurements;
## This is \xi in formula 4.6;
## The data is obtained from the simulated data for the calculation in table 2;
p <- c(0.3, 0.4, 0.15, 0.1, 0.05);

## Input the order of trajectories
ordtraj <- 1 ## linear trajectories

## Call function
## Linear Trajectories
pLSl <- powerLongSurv(N, nevents, tmedian, meantf, p, t, SigmaTheta,
                      sigmae_2, ordtraj, beta, alpha=0.05)
pLSl
show(pLSl)
unclass(pLSl)

## Constant Trajectories
powerLongSurv(N, nevents, tmedian, meantf, p, t, SigmaTheta, sigmae_2,
              ordtraj=0, beta, alpha=0.05)
	      
## Quadratic Trajectories
powerLongSurv(N, nevents, tmedian, meantf, p, t, SigmaTheta,sigmae_2,
              ordtraj=2, beta, alpha=0.05)

## ***************************************************************************

## Example 2.
## **********
## Input elements of Sigma_theta in forumula 4.6;
SigmaTheta <- matrix(c(1.2,0.0,0.0,0.0,0.7,0.0,0.0,0.0,0.8),nrow=3,ncol=3)

N        <-  200; # Total sample size;
nevents  <-  140; # Number of events;
tmedian  <-  0.7; # median survival;
meantf   <-  1.4; # mean follow-up time;
beta     <-  0.2; # Effect of the trajectory;
alpha    <-  0.05;# Type-I Error (2-sided);
sigmae_2 <-  0.09; # measurement error;

## schedule of measurement;
t <- c(0.4, 0.8, 1.2, 1.6);

## Input estimated proportion subjects with 2,3,4,5,6 measurements;
## This is \xi in formula 4.6;
## The data is obtained from the simulated data for the calculation in table 2;
p  <-  c(0.3, 0.4, 0.2, 0.1);

## Input the order of trajectories
ordtraj <- 2 ## quadratic trajectories

## Call function
## Quadratic Trajectories
pLSq <- powerLongSurv(N, nevents, tmedian, meantf, p, t, SigmaTheta,
                      sigmae_2, ordtraj, beta, alpha=0.05)
pLSq
show(pLSq)
unclass(pLSq)

## Constant Trajectories
powerLongSurv(N, nevents, tmedian, meantf, p, t, SigmaTheta, sigmae_2,
              ordtraj=0, beta, alpha=0.05)

## Linear Trajectories
powerLongSurv(N, nevents, tmedian, meantf, p, t, SigmaTheta, sigmae_2,
              ordtraj=1, beta, alpha=0.05)
}
