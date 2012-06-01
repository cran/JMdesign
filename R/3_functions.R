
powerLongSurv<-function(N,nevents,tmedian,meantf,p,t,SigmaTheta,sigmae_2,ordtraj,beta=0, alpha = 0.05)
{
 
    # Error checking for the parameter entries
    # ****************************************

    error_lst<-list()
    error_lst[[1]]<-"'N' - total sample size, minimum is 20."
    error_lst[[2]]<-"'nevents' - number of events must be at least 20 and at most N."
    error_lst[[3]]<-"'tmedian' - median survival time must be positive."
    error_lst[[4]]<-"'meantf' - mean follow-up time must be positive and no greater than max(t)."
    error_lst[[5]]<-"'p' - subject proportions per measurements must be entered as a numeric vector."
    error_lst[[6]]<-"Subject proportions (in vector 'p') must be between 0.0 and 1.0; zero proportion allowed."
    error_lst[[7]]<-"Subject proportions (in vector 'p') must sum up to 1."
    error_lst[[8]]<-"'t' - measurement times must be entered as a numeric vector."
    error_lst[[9]]<-"Measurement times (in vector 't') must be distinct."
    error_lst[[10]]<-"length(t)>2, - at least three measurements."
    error_lst[[11]]<-"All measurement times (in vector 't') must be positive."
    error_lst[[12]]<-"Lengths of 'p' and 't' must agree."
    error_lst[[13]]<-"'SigmaTheta' must be a square matrix."
    error_lst[[14]]<-"'SigmaTheta' must be a symmetric matrix."
    error_lst[[15]]<-"'SigmaTheta' must be a positive-definite matrix."
    error_lst[[16]]<-"Measurement error, 'sigmae_2', must be positive."
    error_lst[[17]]<-"Thajectory order,'ordtraj', must be an integer positive (>=0) and at least one less than the order of the matrix 'SigmaTheta'."
    error_lst[[18]]<-"The significance level, 'alpha', must be strictly between 0.0 and 1.0." 
 
    # Error checking function
    error.powerLongSurv<-function(N,nevents,tmedian,meantf,p,t,SigmaTheta,sigmae_2,ordtraj,beta, alpha){

        if(N<20)  return(1)
        if(nevents<20 || nevents>N) return(2)
        if(tmedian<=0) return(3)
        if(meantf<=0 || meantf>max(t)) return(4)
        if(!is.numeric(p)) return(5)
        if(sum(p*(1-p)<=0)>0) return(6)
        if(sum(p)!=1) return(7)
        if(!is.numeric(t)) return(8)
        t<-as.vector(t)
        if(length(t)>length(unique(t))) return(9)
        if(length(t)<=2) return(10)
        o=order(t);
        if(t[o][1]<=0) return(11)
        if(length(p) != length(t)) return(12)
        m=dim(SigmaTheta)
        if(m[1] != m[2]) return(13)
        if(sum(SigmaTheta != t(SigmaTheta))>0) return(14)
        eig=eigen(SigmaTheta, symmetric=TRUE, only.values = TRUE)
        if(sum(eig$values<=0)>0) return(15)
        if(sigmae_2<=0) return(16)
        if(ordtraj>=length(eig$values)) return(17)
        if(alpha<=0 || alpha>=1) return(18)
        return(0)
    } # End of error.powerLongSurv

    rc<-error.powerLongSurv(N,nevents,tmedian,meantf,p,t,SigmaTheta,sigmae_2,
              ordtraj,beta, alpha)

    if(rc>=1)
        stop(error_lst[[rc]])

 
    # Auxiliary computing functions.
    # *****************************

    integrals <- function(k,a,eta){
        ##Calculate integral(x^i eta exp(-eta x) dx) from 0 to a (0<=a<Infinity), i=0,1,...,k.
        ##Returns a (k+1) vector.
        if(k>=1) integs=(1-cumsum(c(1,cumprod((eta*a)/1:k)))*exp(-eta*a))*c(1,cumprod(1:k/eta))
        if(k==0) integs=c(1-exp(-eta*a))
        return(integs)
    }

    diag2Sums <- function(A){
        ## Given a matrix A, calculate the sums on diagonals parallel to 
        ## secondary diagonal. 
        A=as.matrix(A)
        m=dim(A)
        dsums=0*1:(m[1]+m[2]-1)
        for(j in 1:m[2]) 
            dsums[j:(m[1]+j-1)]=dsums[j:(m[1]+j-1)]+A[,j]
        return(dsums)
    }

    f<-function(BSigma, Rn, sigmae_2){
        n<-nrow(Rn);
        vn<-diag(n)*sigmae_2 + Rn %*% BSigma %*% t(Rn);
        wn<-solve(vn);
        S <-BSigma %*% t(Rn) %*% wn %*% Rn %*% BSigma;
        return(S)
    }


    # Main Calculations 
    # *****************

    title="Joint Modeling of Longitudinal and Survival Data"
    if(ordtraj==0) subtitle="Power Calculation for Unknown Sigma - Constant Trajectory"
    if(ordtraj==1) subtitle="Power Calculation for Unknown Sigma - Linear Trajectory"
    if(ordtraj==2) subtitle="Power Calculation for Unknown Sigma - Quadratic Trajectory"
    if(ordtraj==3) subtitle="Power Calculation for Unknown Sigma - Cubic Trajectory"
    if(ordtraj >= 4) subtitle=paste("Power Calculation for Unknown Sigma - ", ordtraj, "-th Order Trajectory")

    o=order(t); t=t[o]; p=p[o]; # sort 't' internally.
    pp=c(0,p); tt=c(0,t);

    k=ordtraj
    BSigma=SigmaTheta[1:(k+1),1:(k+1)]
    
    # The weighted average of Sigma_theta_hat;
    n=length(tt)
    R=matrix(1,nrow=n,ncol=1)
    if(k>=1){
        R0=matrix(tt,nrow=n,ncol=k,byrow=FALSE)
        R=cbind(R,R0^rep(1:k,each=n))
    }
    X=matrix(0,nrow=k+1,ncol=k+1);
    for(i in 2:length(tt)){
        Ri=as.matrix(R[1:i,]);
        X=X+ pp[i] * f(BSigma, Ri, sigmae_2);
    }

    eta<-log(2)/tmedian; # 1 / E[T];
    censr<-100 * (1-nevents/N);
    et=integrals(2*k,meantf,eta)
    et[1]=1
    et[-1]=(N/nevents)*et[-1]
    zalpha <- qnorm(1 - alpha);
    if(beta < 0) zalpha <- -zalpha
    mus<-sqrt(diag2Sums(X)%*%et) * beta;
    prepower<-sqrt(nevents)*mus - qnorm(1-alpha/2);
    power<-pnorm(prepower);
    if(beta < 0)  power <- 1-power;

    # output
    # ******
    pls=new("powerLongSurv",title=title,subtitle=subtitle,t=t,p=p,N=as.integer(N),
              nevents=as.integer(nevents),censr=censr,tmedian=tmedian, meantf=meantf,
              SigmaTheta=as.matrix(SigmaTheta), ordtraj=as.integer(ordtraj),BSigma=as.matrix(BSigma),
              beta=beta,alpha=alpha,power=as.numeric(power));
    return(pls)

} ## End of "powerLongSurv()"

