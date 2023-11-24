#' Power calculation in joint modeling of longitudinal and 
#'  survival data - k-th Order Trajectories and Unknown Sigma
#'
#' Compute the power in joint modeling of longitudinal and survival data 
#'   when the variance-covariance matrix Sigma_Theta is unknown and the 
#'   trajectories are order k.
#'  
#' The function computes power for a one-sided test, either
#'  \deqn{ H_0: \beta = 0 \quad \mbox{ and } \quad H_{1A}: \beta > 0}{% 
#'    H_0: \beta = 0  and  H_1A: \beta > 0}
#'  or 
#'  \deqn{ H_0: \beta = 0 \quad \mbox{ and } \quad H_{1B}: \beta < 0}{% 
#'    H_0: \beta = 0  and   H_1B: \beta < 0}
#'  with Type I error \eqn{\alpha}. The choice of the alternative is 
#'  determined by the sign of \eqn{\beta}. Negative values for \eqn{\beta} indicate that 
#'  the alternative hypothesis is \eqn{H_{1B}}, while \eqn{\beta \ge 0} indicates that it is 
#'  \eqn{H_{1A}}.
#'  
#'  It creates a \code{powerLongSurv} object.
#'
#' The function \code{powerLongSurv} is used to calculate the power in joint 
#' modeling of longitudinal and survival data.
#' 
#' @param N numeric specifying the total sample size; minimum 20.
#' @param nevents numeric specifying the number of events; at least 20 and at 
#'   most N.
#' @param tmedian numeric specifying the median survival time; positive
#' @param meantf numeric specifying the mean follow-up time; positive and no 
#'   greater than max(t). 
#' @param p numeric vector of estimated subject proportions with 2,3,... 
#'   measurements, respectively, zero proportions allowed.
#' @param t numeric vector of measurement times, distinct positive components; 
#'     same length as \code{p}.
#' @param SigmaTheta numeric matrix specifying the covariance matrix Sigma_Theta
#' @param sigmae_2 numeric specifying the measurement error; positive.
#' @param ordtraj integer specifying the order of trajectories, must be less 
#'   the order of Sigma_Theta
#' @param beta numeric specifying the effect of the trajectory; default value 0.
#' @param alpha numeric, strictly between 0.0 and 1.0, specifying the Type-I 
#'   Error (2-sided), default value 0.05. 
#' @param tol numeric, For floating point objects x and y, if |x-y| <= tol, x==y. 
#'   Passed to R function all.equal.
#'
#' @returns An object of S4 class \code{powerLongSurv},
#'  which has the following 12 components
#'    \item{title }{character string}
#'    \item{subtitle }{character string}
#'    \item{t }{numeric vector}
#'    \item{p }{numeric vector}
#'    \item{N }{integer}
#'    \item{nevents }{integer}
#'    \item{censr }{numeric}
#'    \item{tmedian }{numeric}
#'    \item{meantf }{numeric}
#'    \item{SigmaTheta }{numeric matrix}
#'    \item{ordtraj }{integer}
#'    \item{BSigma }{numeric matrix}
#'    \item{beta }{numeric}
#'    \item{alpha }{numeric}
#'    \item{power }{numeric}
#'  
#' @references
#'  L. M. Chen, J. G. Ibrahim, and H. Chu. Sample size and power determination
#'  in joint modeling of longitudinal and survival data.
#'  Statist. Med. 2011, 30 2295-2309
#'  
#' @author Emil A. Cornea, Liddy M. Chen, Bahjat F. Qaqish, Haitao Chu, and 
#'   Joseph G. Ibrahim
#'   
#' @seealso \code{\link{powerLongSurv-class}}, \code{\link{show-methods}}
#' 
#' @examples
#'  ## Example 1.
#'  ## **********
#'  ## Input elements of Sigma_theta in forumula 4.6;
#'  SigmaTheta <- matrix(c(1.2,0.0,0.0,0.0,0.7,0.0,0.0,0.0,0.8),nrow=3,ncol=3)
#'  
#'  N        <-  200; # Total sample size;
#'  nevents  <-  140; # Number of events;
#'  tmedian  <-  0.7; # median survival;
#'  meantf   <-  1.4; # mean follow-up time;
#'  beta     <-  0.2; # Effect of the trajectory;
#'  alpha    <-  0.05;# Type-I Error (2-sided);
#'  sigmae_2 <-  0.09; # measurement error;
#'  
#'  ## schedule of measurement;
#'  t <- c(0.4, 0.8, 1.2, 1.6, 2) ; # maximum 2 year follow-up;
#'  
#'  ## Input estimated proportion subjects with 2,3,4,5,6 measurements;
#'  ## This is \xi in formula 4.6;
#'  ## The data is obtained from the simulated data for the calculation in table 2;
#'  p <- c(0.3, 0.4, 0.15, 0.1, 0.05);
#'  
#'  ## Input the order of trajectories
#'  ordtraj <- 1 ## linear trajectories
#'  
#'  ## Call function
#'  ## Linear Trajectories
#'  pLSl <- powerLongSurv(N, nevents, tmedian, meantf, p, t, SigmaTheta,
#'                        sigmae_2, ordtraj, beta, alpha=0.05)
#'  pLSl
#'  show(pLSl)
#'  unclass(pLSl)
#'  
#'  ## Constant Trajectories
#'  powerLongSurv(N, nevents, tmedian, meantf, p, t, SigmaTheta, sigmae_2,
#'                ordtraj=0, beta, alpha=0.05)
#'  
#'  ## Quadratic Trajectories
#'  powerLongSurv(N, nevents, tmedian, meantf, p, t, SigmaTheta, sigmae_2,
#'                ordtraj=2, beta, alpha=0.05)
#'  
#'  ## ***************************************************************************
#'  
#'  ## Example 2.
#'  ## **********
#'  ## Input elements of Sigma_theta in forumula 4.6;
#'  SigmaTheta <- matrix(c(1.2,0.0,0.0,0.0,0.7,0.0,0.0,0.0,0.8),nrow=3,ncol=3)
#'  
#'  N        <-  200; # Total sample size;
#'  nevents  <-  140; # Number of events;
#'  tmedian  <-  0.7; # median survival;
#'  meantf   <-  1.4; # mean follow-up time;
#'  beta     <-  0.2; # Effect of the trajectory;
#'  alpha    <-  0.05;# Type-I Error (2-sided);
#'  sigmae_2 <-  0.09; # measurement error;
#'  
#'  ## schedule of measurement;
#'  t <- c(0.4, 0.8, 1.2, 1.6);
#'  
#'  ## Input estimated proportion subjects with 2,3,4,5,6 measurements;
#'  ## This is \xi in formula 4.6;
#'  ## The data is obtained from the simulated data for the calculation in table 2;
#'  p <- c(0.3, 0.4, 0.2, 0.1);
#'  
#'  ## Input the order of trajectories
#'  ordtraj <- 2 ## quadratic trajectories
#'  
#'  ## Call function
#'  ## Quadratic Trajectories
#'  pLSq <- powerLongSurv(N,nevents,tmedian,meantf,p,t,SigmaTheta,sigmae_2,ordtraj,beta, alpha = 0.05)
#'  pLSq
#'  show(pLSq)
#'  unclass(pLSq)
#'  
#'  ## Constant Trajectories
#'  powerLongSurv(N, nevents, tmedian, meantf, p, t, SigmaTheta, sigmae_2,
#'                ordtraj=0, beta, alpha=0.05)
#'  
#'  ## Linear Trajectories
#'  powerLongSurv(N, nevents, tmedian, meantf, p, t, SigmaTheta, sigmae_2,
#'                ordtraj=1, beta, alpha=0.05)
#'
#' @importFrom stats pnorm qnorm
#' @export
powerLongSurv <- function(N,
                          nevents,
                          tmedian,
                          meantf,
                          p,
                          t,
                          SigmaTheta,
                          sigmae_2,
                          ordtraj,
                          beta = 0, 
                          alpha = 0.05,
                          tol = 1.5e-8){
 

  errorTesting(N=N,
               nevents=nevents,
               tmedian=tmedian,
               meantf=meantf,
               p=p,
               t=t,
               SigmaTheta=SigmaTheta,
               sigmae_2=sigmae_2,
               ordtraj=ordtraj,
               beta=beta, 
               alpha=alpha,
               tol=tol)

 
  # Auxiliary computing functions.
  # *****************************

  ##Calculate integral(x^i eta exp(-eta x) dx) from 0 to a (0<=a<Infinity), 
  ##i=0,1,...,k.
  ##Returns a (k+1) vector.
  integrals <- function(k,a,eta){
                 if( k >= 1 ) {
                   return((1.0 - cumsum(c(1,cumprod((eta*a)/1:k)))*exp(-eta*a))*
                           c(1,cumprod(1:k/eta)) )
                 } else if( k == 0 ) {
                   return( c(1.0 - exp(-eta*a)) )
                 }
               }

  ## Given a matrix A, calculate the sums on diagonals parallel to 
  ## secondary diagonal. 
  diag2Sums <- function(A){
                 A <- as.matrix(A)
                 m <- dim(A)
                 dsums <- numeric(m[1]+m[2]-1)
                 for(j in 1:m[2]) { 
                   dsums[j:(m[1]+j-1)] <- dsums[j:(m[1]+j-1)] + A[,j]
                 }
                 return(dsums)
               }

  f <- function(BSigma, Rn, sigmae_2){
         n <- nrow(Rn)
         tRn <- t(Rn)
         vn <- diag(n)*sigmae_2 + Rn %*% BSigma %*% tRn
         wn <- solve(vn)
         S <- BSigma %*% tRn %*% wn %*% Rn %*% BSigma
         return(S)
       }


  # Main Calculations 
  # *****************
  title <- "Joint Modeling of Longitudinal and Survival Data"
  tnames <- c("Constant", "Linear", "Quadratic", "Cubic", 
              paste(ordtraj, "-th Order", sep="") )
  subtitle <- paste("Power Calculation for Unknown Sigma - ",
                   tnames[min(ordtraj+1,5)], " Trajectory", sep="")

  o <- order(t)
  t <- t[o]
  p <- p[o]
  tt <- c(0,t)

  k <- ordtraj
  BSigma <- SigmaTheta[1:(k+1),1:(k+1)]
    
  # The weighted average of Sigma_theta_hat
  n <- length(tt)
  R <- matrix(1, nrow=n, ncol=1)
  if( k >= 1 ){
    R0 <- matrix(tt, nrow=n, ncol=k, byrow=FALSE)
    R <- cbind(R, R0^rep(1:k,each=n))
  }

  X <- matrix(0, nrow=k+1, ncol=k+1)
  for(i in 2:n){
    Ri <- as.matrix(R[1:i,])
    X <- X + p[i-1] * f(BSigma, Ri, sigmae_2)
  }

  eta <- log(2.0)/tmedian
  censr <- 100.0 * (1.0-nevents/N)
  et <- integrals(k=2*k, a=meantf, eta=eta)
  et[1] <- 1.0
  et[-1] <- (N/nevents)*et[-1]
  zalpha <- stats::qnorm(p={1 - alpha})
  if(beta < 0.0) zalpha <- -zalpha
  mus <- sqrt(diag2Sums(A=X) %*% et) * beta
  prepower <- sqrt(nevents)*mus - qnorm(p={1-alpha/2})
  power <- stats::pnorm(q=prepower)
  if(beta < 0.0)  power <- 1.0 - power

  # output
  # ******
  pls <- new("powerLongSurv",
             title=title,
             subtitle=subtitle,
             t=t,
             p=p,
             N=as.integer(N),
             nevents=as.integer(nevents),
             censr=censr,
             tmedian=tmedian, 
             meantf=meantf,
             SigmaTheta=as.matrix(SigmaTheta), 
             ordtraj=as.integer(ordtraj),
             BSigma=as.matrix(BSigma),
             beta=beta,
             alpha=alpha,
             power=as.numeric(power))

  return(pls)

} ## End of "powerLongSurv()"

