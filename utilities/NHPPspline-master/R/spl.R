#' The spline-based arrival rate function
#'
#' This function returns the value of the spline at time t.
#' @param t A scalar or vector of times at which the function will be evaluated
#' @return The value of the spline function at t is returned
#' @examples
#' ## spline-based intensity defined on [0,10]
#' d<-3      #spline function degree
#' knots<-seq(-3,13,by=1)
#' n<-length(knots)-(d+1)
#' c<-runif(n,1,5)
#' t<-1.5
#' spl(t,knots,c)  # the spline-based intensity evaluated at time 1.5
#' t<-seq(0,10,by=0.01)
#' plot(t,spl(t,knots,c),type="l") # a plot of the spline-based intensity over the interval [0,10]

#' @export

## A spline-based method for modelling and generating a nonhomogeneous Poisson process

### Three functions
# 1. To fit a spline-based intensity function from NHPP arrival time observations
# 2. To generate arrivals from the resulting spline-based intensity function
# 3. The spline-based function

## 3.
# This function returns the value of the spline at time t.

## Inputs
# t - a scalar or vector of times at which the function will be evaluated

## Outputs
# The value of the spline function at t is returned

spl<-function(t,knots,c){
  return(f_b(t,knots,c,d=3))
}


##DESCRIPTION

#Given a vector x calculates f(x)
#For each x it finds the knot interval it lies in and uses this (knot_mu,knot_{mu+1})

f<-function(x,knots,c,d)
{
  kn<-length(knots)
  n<-length(c)
  if( x>=max(knots) | x<min(knots)){
    print("x outside of knot range")
    return(0)
  }
  else{
    mu_in<-max(which(knots<=x))   # the index of t_mu
    if(mu_in-d<=0){
      R<-rep(0,mu_in)   # sets R at the required length
      for(p in 1:mu_in){
        R[p]<-B_i_3(p,x,knots)
      }
      O<-c[1:mu_in]
      f<-R%*%O
      return(f)
    }else if(mu_in>n){
      R<-rep(0, (kn-mu_in))
      for(p in (mu_in-d):(mu_in-d-1+(kn-mu_in))){
        R[(d+1)-(mu_in-p)]<-B_i_3(p,x,knots)
      }
      O<-c[ (mu_in-d ): (mu_in-d-1+(kn-mu_in))]
      f<-R%*%O
      return(f)
    }
    else{
      O<-c[(mu_in-d):mu_in]
      for(k in d:1){
        R<-matrix(rep(0,k*(k+1)),nrow=k)
        for(i in 1:dim(R)[1])
        {
          R[i,i]<-(knots[mu_in+i]-x)/(knots[mu_in+i]-knots[mu_in+i-k])
          R[i,i+1]<- (x - knots[mu_in+i-k])/(knots[mu_in+i]-knots[mu_in+i-k])
        }
        O<-R%*%O
      }
      f<-O
      return(f)
    }
  }
}
f_b<-Vectorize(f,"x")
