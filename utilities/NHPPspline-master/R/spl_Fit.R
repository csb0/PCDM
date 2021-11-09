#' Fitting a spline function from cardinal B-splines
#'
#' This function fits a spline function from cardinal B-splines and returns the optimal spline coefficients and the knot points it was built upon
#' @param  Arrivals A vector of arrival times. If multiple days the days are appended and not sorted.
#' @param m A scalar number of days of arrival data
#' @param  n_arrs A vector containint the number of arrivals on each observed day
#' @param  Tstart The start of the period of observation
#' @param  Tend The end of the period of observation
#' @param  kn The number of knots on which to build the spline function
#' @param  cyclic A TRUE/FALSE value stating whether the spline function should be constrained to be cyclic or not
#' @param  c_init (optional) If provided this is a scalar that tell us initial constant value from which to start the search for the optimal spline coefficients.
#' @param  plot (optional) If provided this is says whether to output a pplot or not. If you do not want a plot do not use this parameter, if you do make pplot=0
#' @return Opt_Spl_Coeffs - The optimised spline coefficients
#' @return Knots - The knots on which the spline function was built
#' @examples
#' ## Fitting the spline-based intensity given arrivals from a cyclic sinusoidal
#' data(Arrs)
#' data(n_Arrs)
#' m<-length(n_Arrs)
#' Tstart<-0
#' Tend<-24
#' kn<-50
#' cyclic=TRUE
#' spl_Fit(as.numeric(Arrs),m,as.numeric(n_Arrs),Tstart,Tend,kn,cyclic,pplot=0)

#' @export

## A spline-based method for modelling and generating a nonhomogeneous Poisson process

### Three functions
# 1. To fit a spline-based intensity function from NHPP arrival time observations
# 2. To generate arrivals from the resulting spline-based intensity function
# 3. The spline-based function

## 1.
# Fits a spline function from cardinal B-splines
# Returns the optimal spline coefficients and the knot points it was built upon

## Prerequisits
# This function requires the packages psych, pracma, parallel, Matrix and optiSolve

## Inputs
# Arrivals - a vector of arrival times. If multiple days the days are appended and not sorted.
# m - a scalar number of days of arrival data
# n_arrs - a vector containint the number of arrivals on each observed day
# Tstart - The start of the period of observation
# Tend - The end of the period of observation
# kn - The number of knots on which to build the spline function
# cyclic - a TRUE/FALSE value stating whether the spline function should be constrained to be cyclic or not
## Optional parameters
# c_init - if provided this is a scalar that tell us initial constant value from which to start the search for the optimal spline coefficients.
# plot - if provided this is says whether to output a pplot or not. If you do not want a plot do not use this parameter, if you do make pplot=0

## Outputs
# Opt_Spl_Coeffs - the optimised spline coefficients
# Knots - the knots on which the spline function was built

spl_Fit <- function(Arrivals,m,n_arrs,Tstart,Tend,kn,cyclic,c_init=NULL,pplot=NULL)
{
    # setting up the knot sequence
    min_boundary<-Tstart
    max_boundary<-Tend
    inner_knots<-seq(min_boundary,max_boundary,by=((max_boundary-min_boundary)/(kn-7)))  #knots between 0 and Tend
    knots<-c(-inner_knots[4],-inner_knots[3],-inner_knots[2],inner_knots,Tend+(Tend-inner_knots[length(inner_knots)-1]),Tend+(Tend-inner_knots[length(inner_knots)-2]),Tend+(Tend-inner_knots[length(inner_knots)-3]))

    d<- 3            # the degree of the Spline function
    n<- kn-(d+1)     # the number of B-splines

    ##TIME SAVERS##
    # vector containting integral of all n B-splines on [Tstart,Tend]
    int_B<-rep(0,n)
    for(i in 1:n){int_B[i]<-integrate(B_spline_val_b,lower=0,upper=Tend,knots=knots,d=d,p=i,rel.tol=.Machine$double.eps^0.1)$value}
    ##TIME SAVERS##

    len <- length(Arrivals)
    arrs_unordered <- Arrivals
    all_arrs<-sort(arrs_unordered)
    m_i <- n_arrs

    Omeg<-Omega(d,n,knots,2,Tend)    # the integral of the squared second derivatives matrix

    ## TIME SAVERS ##
    # the columns of this matrix are the knot values and we give the value of the B_spline over all arrivals at
    B_spl<-B_spline_val_c(all_arrs,knots,d,seq(1,n,by=1))
    B_uo<-B_spline_val_c(arrs_unordered,knots,d,seq(1,n,by=1))
    day_index<-rep(1:m,m_i)
    mat<- cbind(B_uo,arrs_unordered,day_index)
    ## ##

    ##Trust region algorithm
    Delta_max<-5
    Delta<-1
    eta<-0.2   #accept step
    epsilon<-0.05
    max_it<-10000

    # warm start
    lam<-0
    if(is.null(c_init)==FALSE){  c_init<-rep(c_init,n)
    } else{c_init<-rep(mean(n_arrs)/Tend,n)}

    c_warm<-opt_c(c_init,Delta,Delta_max,eta,epsilon,max_it,lam,all_arrs,m,Omeg,knots,B_spl,int_B,d,Tend,cyclic)

    ## check we're moving towards the minimum RIC
    lam<-5.12
    c_eta<-opt_c(c_warm,Delta,Delta_max,eta,epsilon,max_it,lam,all_arrs,m,Omeg,knots,B_spl,int_B,d,Tend,cyclic)
    RIC_big<-RIC(lam,c_eta,knots,all_arrs,m,Omeg,Tend,m_i,arrs_unordered,B_spl,B_uo,int_B)
    c_eta<-opt_c(c_warm,Delta,Delta_max,eta,epsilon,max_it,lam/2,all_arrs,m,Omeg,knots,B_spl,int_B,d,Tend,cyclic)
    RIC_small<-RIC(lam/2,c_eta,knots,all_arrs,m,Omeg,Tend,m_i,arrs_unordered,B_spl,B_uo,int_B)

    while (RIC_big < RIC_small )
    {
      RIC_small<-RIC_big
      lam<-lam*2
      c_eta<-opt_c(c_warm,Delta,Delta_max,eta,epsilon,max_it,lam,all_arrs,m,Omeg,knots,B_spl,int_B,d,Tend,cyclic)
      RIC_big<-RIC(lam,c_eta,knots,all_arrs,m,Omeg,Tend,m_i,arrs_unordered,B_spl,B_uo,int_B)
    }

    ## now we are moving towards the minimum find the narrowed interval O in which penalty that minimises the RIC lies
    while (RIC_small < RIC_big)
    {
      RIC_big<-RIC_small
      lam<-lam/2
      c_eta<-opt_c(c_warm,Delta,Delta_max,eta,epsilon,max_it,lam/2,all_arrs,m,Omeg,knots,B_spl,int_B,d,Tend,cyclic)
      RIC_small<-RIC(lam/2,c_eta,knots,all_arrs,m,Omeg,Tend,m_i,arrs_unordered,B_spl,B_uo,int_B)
    }
    lam_low=lam     # lower bound of O
    lam_high=2*lam  # upper bound of O

    ## more intensive search of the interval O
    opt<-optimise(f=opt_pen, interval=c(lam_low,lam_high), c_init=c_warm,knots=knots,arrs_unordered=arrs_unordered,m=m,m_i=m_i,Omeg=Omeg,Tend=Tend,Delta=Delta,Delta_max=Delta_max,eta=eta,epsilon=epsilon,max_it=max_it,B=B_spl,B_uo=B_uo,int_B=int_B ,lower=lam_low,upper=lam_high,d=d,cyclic=cyclic,tol=1e4)
    lam_opt<-opt$minimum   # optimal penalty
    c_op<-opt_c(c_warm,Delta,Delta_max,eta,epsilon,max_it,lam_opt,all_arrs,m,Omeg,knots,B_spl,int_B,d,Tend,cyclic)    #optimal spline coefficients for the optimal penalty

    if(is.null(pplot)==FALSE){     ## optional plot of the spline function
      q<-seq(knots[d+1],knots[n+1],by=0.1)
      plot(q,f_b(q,knots,c_op,d),type="l",col="blue",xlim=c(0,Tend),ylim=c(0,15),lwd=1,ylab="Rate Function",xlab="Time, t")#,xlim=c(0,Tend))
    }

    OO<-list(c_op,knots) #as.num(t(output))
    names(OO) <-c("Opt_Spl_Coeffs", "Knots")
    return(OO)
}



## DESCRIPTION
# calculates the penalty on the NHPP log-likelihood

# the two derivatives multiplied together at time x
omega<-function(x,i,j,knots,d,r){
  out<-c()
  a<-rth_deriv_B_mu(x,knots,d,r,i)
  b<-rth_deriv_B_mu(x,knots,d,r,j)
  if( a!=0 & b!=0){
    out<-a*b
  }
  else{ out<-0 }
  return(out)
}
omega<-Vectorize(omega,"x")

# the integral of omega
Omega<- function(d,n,knots,r,Tend)
{
  d<-3
  Omg<-matrix(rep(0,n*n),nrow=n)
  for(i in 1:n){
    for(j in i:min((i+d),n)){
      ith<-i:(i+4)
      jth<-j:(j+4)
      intersct<-intersect(ith,jth)
      min<-knots[max(min(intersct),d+1)]
      max<-knots[min(max(intersct),n+1)]
      int<-integrate(omega,lower=min,upper=max,i,j,knots,d,r)$value
      if(abs(int)<0.000001){int<-0}
      Omg[i,j]<- int
      Omg[j,i]<- int
    }
  }
  return(Omg)
}

##DESCRIPTION

#for x in interval [t_mu,t_mu+1) this gives the non-zero B_splines that are >0 in that interval
### the max(knots) will not return a value due to the open condition on the end of the interval


non_zero_B_spline<-function(x,knots,d)
{
  kn<-length(knots)
  n<-kn-(d+1)
  if( x >= max(knots) || x < min(knots)){
    print("x outside of knot range")
    return(0)
  }
  else
  {
    mu_in<-max(which(knots<=x))      #gives the index of the knot interval x falls in
    if(mu_in-d<=0){
      R<-rep(0,mu_in)   # this sets R up to be the required length
      for(p in 1:mu_in){
        R[p]<-B_i_3(p,x,knots)
      }
      return(R)
    }else if(mu_in>n){
      R<-rep(0, (kn-mu_in))
      for(p in (mu_in-d):(mu_in-d-1+(kn-mu_in))){
        R[(d+1)-(mu_in-p)]<-B_i_3(p,x,knots)
      }
      return(R)
    }
    else{
      for(k in d:1){
        R<-matrix(rep(0,k*(k+1)),nrow=k)
        for(l in 1:k)
        {
          R[l,l]<-(knots[mu_in+l]-x)/(knots[mu_in+l]-knots[mu_in+l-k])
          R[l,l+1]<- (x - knots[mu_in+l-k])/(knots[mu_in+l]-knots[mu_in+l-k])
        }
        if(k < d)
        {
          R_last<-R%*%R_last
        }
        else{R_last<-R}
      }
      return(R_last)
    }
  }
}


B_spline_val<-function(x,knots,d,p)   #p is the knot we want to know about
{
  mu_in<-max(which(knots<=x))     #gives the index of the knot interval x falls in
  kn<-length(knots)
  n<-length(knots) - (d+1)
  if(p>mu_in | p<(mu_in-d)){
    return(0)
  }
  else{
    if(mu_in-d<=0){
      R<-rep(0,mu_in)   # this sets R up to be the required length
      for(q in 1:mu_in){
        R[q]<-B_i_3(q,x,knots)
      }
      return(R[p])
    }else if(mu_in > n){
      R<-rep(0,(kn-mu_in))
      for(q in (mu_in-d):(mu_in-d-1+(kn-mu_in))){
        R[(d+1)-(mu_in-q)]<-B_i_3(q,x,knots)
      }
      return(R[(d+1-(mu_in-p))])
    }
    else{
      for(k in d:1){
        R<-matrix(rep(0,k*(k+1)),nrow=k)
        for(l in 1:k)
        {
          R[l,l]<-(knots[mu_in+l]-x)/(knots[mu_in+l]-knots[mu_in+l-k])
          R[l,l+1]<- (x - knots[mu_in+l-k])/(knots[mu_in+l]-knots[mu_in+l-k])
        }
        if(k < d)
        {
          R_last<-R%*%R_last
        }
        else{R_last<-R}
      }
      q<-(d+1) - (mu_in-p)
      return(R_last[q])
    }
  }
}

B_spline_val_b <- Vectorize(B_spline_val, "x")
B_spline_val_c <- Vectorize(B_spline_val_b, "p")

## DESCRIPTION
# The four functions below are the recurrence formulae
# required for the construction of a cubic B-spline
# once the knot sequence of the spline function is
# fixed these functions construct the B-spline basis
# functions

#cubic
B_i_3<-function(i,x,knots){
  d<-3
  out<- (x-knots[i])/(knots[i+d] - knots[i])*B_i_2(i,x,knots) + (knots[i+d+1] - x)/(knots[i+d+1] -knots[i+1])*B_i_2(i+1,x,knots)
  return(out)
}
B_i_3_b<-Vectorize(B_i_3,"x") # an array x can now be used as the input and a vector of cubic B-spline outputs returned
B_i_3_c<-Vectorize(B_i_3,"i") # an array i can now be used as the input and a vector of cubic B-spline outputs returned

#quadratic
B_i_2<-function(i,x,knots){
  B_2<-c( B_i_1(i,x,knots)*((x-knots[i])/(knots[i+2]-knots[i])) + B_i_1(i+1,x,knots)*((knots[i+3]-x)/(knots[i+3]-knots[i+1])))#, B_1[2]*((x-knots[i+1])/(knots[i+3]-knots[i+1]))  +   B_1[3]*((knots[i+4]-x)/(knots[i+4]-knots[i+2]))      )
  return(B_2)
}
B_i_2_b<-Vectorize(B_i_2,"x")
B_i_2_c<-Vectorize(B_i_2,"i")

#linear
B_i_1<-function(i,x,knots){
  B_1<-c(B_i_0(i,x,knots)*((x-knots[i])/(knots[i+1]-knots[i])) + B_i_0(i+1,x,knots)*((knots[i+2]-x)/(knots[i+2]-knots[i+1]))) # , B_0[2]*((x-knots[i+1])/(knots[i+2]-knots[i+1])) + B_0[3]*((knots[i+3]-x)/(knots[i+3]-knots[i+2])), B_0[3]*((x-knots[i+2])/(knots[i+3]-knots[i+2])) + B_0[4]*((knots[i+4]-x)/(knots[i+4]-knots[i+3])))
  return(B_1)
}
B_i_1_b<-Vectorize(B_i_1,"x")
B_i_1_c<-Vectorize(B_i_1,"i")

#constant
B_i_0<-function(i,x,knots){
  if(x >= knots[(i)] && x < knots[i+1]){
    return(1)
  }else{
    return(0)
  }
}
B_i_0_b<-Vectorize(B_i_0,"x")
B_i_0_c<-Vectorize(B_i_0,"i")


### DESCRIPTION

# calculation of the derviatives of the spline function and B-spline functions

#calculates the 1st derivative of f(x) given x in (knot_mu,knot_{mu+1})
D_B_i_1<-function(i,x,knots){
  d<-3
  mu_in<-max(which(knots<=x))
  if ( i < (mu_in-d) | i > mu_in){
    return(0)
  }else{
    out<-d*( B_i_2(i,x,knots)/(knots[i+d]-knots[i]) - B_i_2(i+1,x,knots)/(knots[i+d+1]-knots[i+1]))
    return(out)
  }
}
D_B_i_1_b<-Vectorize(D_B_i_1,"x")

#calculates the 2nd derivative of f(x) given x in (knot_mu,knot_{mu+1})
D_B_i_2<-function(i,x,knots){
  d<-2
  mu_in<-max(which(knots<=x))
  if ( i < (mu_in-d-1) | i > mu_in){
    return(0)
  }else{
    out<-d*( (d-1)*(B_i_1(i,x,knots)/(knots[i+d-1]-knots[i]) - B_i_1(i+1,x,knots)/(knots[i+d]-knots[i+1]))  -  (d-1)*(B_i_1(i+1,x,knots)/(knots[i+d]-knots[i+1]) - B_i_1(i+2,x,knots)/(knots[i+d+1]-knots[i+2]))   )
    return(out)
  }
}
D_B_i_2_b<-Vectorize(D_B_i_2,"x")


rth_deriv<-function(x,knots,c,d,r)
{
  kn<-length(knots)
  n<- kn - (d+1)
  if( x>max(knots) | x<min(knots)){
    print("x outside of knot range")
    return(0)
  }
  else
  {
    mu_in<-max(which(knots<=x))
    if(mu_in-d<=0){
      if(r==2){
        R<-rep(0,mu_in)
        for(q in 1:mu_in){
          R[q]<-D_B_i_2(q,x,knots)
        }
        R<-R%*%c[1:mu_in]
        return(R)
      }
      else{
        R<-rep(0,mu_in)
        for(q in 1:mu_in){
          R[q]<-D_B_i_1(q,x,knots)
        }
        R<-R%*%c[1:mu_in]
        return(R)
      }
    }else if(mu_in > n){
      if(r==2){
        R<-rep(0,(kn-mu_in))
        for(p in (mu_in-d):(mu_in-d-1+(kn-mu_in))){
          R[(d+1)-(mu_in-p)]<-D_B_i_2(p,x,knots)
        }
        R<-R%*%c[(mu_in-d ): (mu_in-d-1+(kn-mu_in))]
        return(R)
      }
      if(r==1){
        R<-rep(0,(kn-mu_in))
        for(p in (mu_in-d):(mu_in-d-1+(kn-mu_in))){
          R[(d+1)-(mu_in-p)]<-D_B_i_1(q,x,knots)
        }
        R<-R%*%c[(mu_in-d ): (mu_in-d-1+(kn-mu_in))]
        return(R)
      }
    }
    else{
      O<-c[(mu_in-d):mu_in]
      R<-R_k(x,1,mu_in,knots)
      if((d-r)>=2){
        for(i in 2:(d-r)){
          R<-R%*%R_k(x,i,mu_in,knots)
        }
      }
      for(j in (d-r+1):d){
        R<-R%*%d_dx_R_k(j,mu_in,knots)
      }
      R<-(factorial(d)/factorial(d-r))*R%*%O
    }
    return(R)
  }
}
rth_deriv_b<-Vectorize(rth_deriv,"x")


rth_deriv_sq<-function(x,knots,c,d,r)
{
  kn<-length(knots)
  n<- kn - (d+1)
  if( x>max(knots) | x<min(knots)){
    print("x outside of knot range")
    return(0)
  }
  else
  {
    mu_in<-max(which(knots<=x))
    if(mu_in-d<=0){
      if(r==2){
        R<-rep(0,mu_in)
        for(q in 1:mu_in){
          R[q]<-D_B_i_2(q,x,knots)
        }
        R<-R%*%c[1:mu_in]
        return(R^2)
      }
      else{
        R<-rep(0,mu_in)
        for(q in 1:mu_in){
          R[q]<-D_B_i_1(q,x,knots)
        }
        R<-R%*%c[1:mu_in]
        return(R^2)
      }
    }else if(mu_in > n){
      if(r==2){
        R<-rep(0,(kn-mu_in))
        for(p in (mu_in-d):(mu_in-d-1+(kn-mu_in))){
          R[(d+1)-(mu_in-p)]<-D_B_i_2(p,x,knots)
        }
        R<-R%*%c[(mu_in-d ): (mu_in-d-1+(kn-mu_in))]
        return(R^2)
      }
      if(r==1){
        R<-rep(0,(kn-mu_in))
        for(p in (mu_in-d):(mu_in-d-1+(kn-mu_in))){
          R[(d+1)-(mu_in-p)]<-D_B_i_1(q,x,knots)
        }
        R<-R%*%c[(mu_in-d ): (mu_in-d-1+(kn-mu_in))]
        return(R^2)
      }
    }
    else{
      O<-c[(mu_in-d):mu_in]
      R<-R_k(x,1,mu_in,knots)
      if((d-r)>=2){
        for(i in 2:(d-r)){
          R<-R%*%R_k(x,i,mu_in,knots)
        }
      }
      for(j in (d-r+1):d){
        R<-R%*%d_dx_R_k(j,mu_in,knots)
      }
      R<-(factorial(d)/factorial(d-r))*R%*%O
    }
    return(R^2)
  }
}
rth_deriv_sq_b<- Vectorize(rth_deriv_sq, "x")


#this function gives the rth derivative of the jth B_spline at x
rth_deriv_B_mu<-function(x,knots,d,r,j)
{
  kn<-length(knots)
  n<- kn - (d+1)
  mu_in<-max(which(knots<=x))  #where is j relative to mu_in
  if(j<(mu_in-d) | j>min(mu_in,n) )
  {
    return(0)
  }
  else{
    if(r==1){
      return(D_B_i_1(j,x,knots))
    }
    else{
      return(D_B_i_2(j,x,knots))
    }
  }
}
rth_deriv_B_mu_b<-Vectorize(rth_deriv_B_mu,"x")
rth_deriv_B_mu_c<-Vectorize(rth_deriv_B_mu,"j")

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



# A second-order Taylor series approximation of the likelihood function
# the model of the likelihood - where H should be PD
m_p<-function(c,knots,all_arrs,lam,m,Omeg,Tend,g,H,p){
  m_p<- l_p(c,knots,all_arrs,lam,m,Omeg,Tend) + t(g)%*%p + 0.5*t(p)%*%H%*%p     #the negative likelihood  + the gradient + the Hessian
  return(m_p)
}

## DESCRIPTION
# calculates the RIC for a fixed penalty parameter lam

opt_pen<-function(lam,c_init,knots,arrs_unordered,m,m_i,Omeg,Tend,Delta,Delta_max,eta,epsilon,max_it,B,B_uo,int_B,d,cyclic){
  all_arrs<-sort(arrs_unordered)
  c<-opt_c(c_init,Delta,Delta_max,eta,epsilon,max_it,lam,all_arrs,m,Omeg,knots,B,int_B,d,Tend,cyclic)
  ric<-RIC(lam,c,knots,all_arrs,m,Omeg,Tend,m_i,arrs_unordered,B,B_uo,int_B)
  return(ric)
}


## DESCRIPTION
# calculates the RIC score given a penalty lam and the optimal spline coefficients

RIC<-function(lam,c,knots,all_arrs,m,Omeg,Tend,m_i,arrs_unordered,B,B_uo,int_B){
  n<-length(c)
  I<- I_calc(m_i,arrs_unordered,c,knots,lam,m,Omeg,B_uo,int_B)
  J<- -(1/m)*Hess(c,knots,all_arrs,lam,m,Omeg,B)
  return(-2*l_p(c,knots,all_arrs,lam,m,Omeg,Tend) + 2*tr(I%*%solve(J)))
}


## DESCRIPTION
# trust region optimisation (TRO) algorithm
# n - the dimension of the optimisation problem

# solves the TRO subproblem
# aiming to minimise -m(p)

Opt_alt<-function(B,g,n,Delta,c,knots,d,Tend,cyclic,all_arrs,lam,m,Omeg){
  yip<-0.000000001

  obj<- quadfun(-0.5*B,-g,0)#,-l_p(c,knots,all_arrs,lam,m,Omeg,Tend))
  qcon<-quadcon(spMatrix(n,n,seq(1,n,by=1),seq(1,n,by=1),rep(1,n)),a=rep(0,n),d=0,dir="<=",Delta,use=TRUE)

  if(cyclic == TRUE){
    A<-matrix(rep(0,length(c)*3),nrow=3)
    for(i in 1:length(c)){
      A[1,i]<-B_spline_val(0,knots,d,i)-B_spline_val(Tend-yip,knots,d,i)
      A[2,i]<-rth_deriv_B_mu(x=0,knots,d,1,i)-rth_deriv_B_mu(x=(Tend-yip),knots,d,1,i)
      A[3,i]<-rth_deriv_B_mu(x=0,knots,d,2,i)-rth_deriv_B_mu(x=(Tend-yip),knots,d,2,i)
    }
    rownames(A)<-c("1","2","3")
    dir=c("==","==","==")
    val = c(0,0,0)
    lcon<-lincon( A=A , dir = dir,val=val)
  }else{
    A<-matrix(rep(0,length(c)),nrow=1,byrow=T)
    rownames(A)<-c("1")
    lcon<-lincon(A,dir=rep("==",nrow(A)),val=0)
  }
  lbc <- lbcon(val = -c + rep(0.0001,length(g)))

  op<-cop(obj, max=FALSE, lb=lbc , lc=lcon, qc=qcon)

  X<-c
  ### here is the problem
  names(X)<-paste(1:length(c))
  result <- solvecop(op,solver="alabama",X=X,quiet=TRUE)
  return (result$x )
}

## DESCRIPTION

# Finds the optimal values of the spline coefficients
# Uses trust region optimisation to fiund optimum

opt_c<-function(c_init,Delta,Delta_max,eta,epsilon,max_it,lam,all_arrs,m,Omeg,knots,B,int_B,d,Tend,cyclic){
  c_k<-c_init
  for(k in 1:max_it)
  {
    H<- Hess(c_k,knots,all_arrs,lam,m,Omeg,B)
    g<- Grad(c_k,knots,all_arrs,lam,m,Omeg,B,int_B)

    p_k<- Opt_alt(H,g,length(c_k),Delta,c_k,knots,d,Tend,cyclic,all_arrs,lam,m,Omeg)    #takes in the Hessian and gradient of the likelihood
    l<- l_p(c_k,knots,all_arrs,lam,m,Omeg,Tend)
    rho_k<- ( -l + l_p(c_k+p_k,knots,all_arrs,lam,m,Omeg,Tend)) / (-l + m_p(c_k,knots,all_arrs,lam,m,Omeg,Tend,g,H,p_k))
    if(rho_k < 0.25)
    {
      Delta<-0.25*Delta
    }
    if(rho_k > 0.75 & t(p_k)%*%p_k >= Delta^2-0.001){
      Delta<- min(2*Delta,Delta_max)
    }
    else{
      Delta<-Delta
    }
    if(rho_k > eta){
      c_k <- c_k + p_k
    }
    if(t(p_k)%*%p_k< epsilon^2)
    {
      break
    }
  }
  c_opt<-c_k
  return(c_opt)
}

## This function takes a NHPP described by a spline and gives back it's gradient at a point

#c - is the vector of spline coefficients at this step
#all_arrs - the observed arrivals over
#m - days
#Omeg - describes the curvature of the B-spline
#lam - the penalty on the curvature

#the  gradient of the likelihood
Grad<-function(c,knots,all_arrs,lam,m,Omeg,B,int_B){
  int<-c()
  d<-3
  n<-length(c)
  f_lin<-f_b(all_arrs,knots,c,d)
  S_sum<-apply(B/f_lin,2,sum)
  Grad_p <- S_sum - m*int_B - 0.5*lam*t(c)%*%Omeg
  return(as.vector(Grad_p))
}

#### The Likelihood function

# the penalised NHPP likelihood that we will try to maximise
l_p<-function(c,knots,all_arrs,lam,m,Omeg,Tend){

  l <- sum(log(f_b(all_arrs,knots,c,3))) - m*integral(f_b,xmin=0,xmax=Tend,knots=knots,c=c,d=3) - 0.5*lam*t(c)%*%Omeg%*%c
  return(l)
}


l_inf<-function(c,knots,all_arrs,lam,m,Omeg,Tend,B,int_B){
  a<-which(c<0)
  print(a)
  if(length(a)>0){c[a]=0}
  value <- sum(log(f_b(all_arrs,knots,c,3))) - m*integral(f_b,xmin=0,xmax=Tend,knots=knots,c=c,d=3) - 0.5*lam*t(c)%*%Omeg%*%c
  gradient<-Grad(c,knots,all_arrs,lam,m,Omeg,B,int_B)
  hessian<-Hess(c,knots,all_arrs,lam,m,Omeg,B)
  print(value)
  print(c)
  l<-list(value,gradient,hessian)
  names(l)<-c("value","gradient","hessian")
  return(l)
}


I_calc<-function(m_i,arrs_unordered,c,knots,lam,m,Omeg,B_uo,int_B){
  n<-length(c)
  sum<-matrix(rep(0,n*n),nrow=n)
  cs<-c(0,cumsum(m_i))
  for(i in 1:m){
    arrs<-arrs_unordered[(1+cs[i]):cs[(i+1)]]
    G<-Grad(c,knots,arrs,lam,1,Omeg/m,B_uo[(1+cs[i]):cs[(i+1)],],int_B)     # finds the gradient of the arrivals from a single day
    sum<-sum + G%*%t(G)
  }
  return((1/m)*sum)
}


##DESCRIPTION

#calculates the Hessian of a NHPP likelihood

Hess<-function(c,knots,all_arrs,lam,m,Omeg,B){
  d<-3
  Hess_p<- -0.5*lam*Omeg
  n<-length(c)
  f_sq<-f_b(all_arrs,knots,c,d)^2
  for(i in 1:n){
    for(k in i:min((i+d),n)){
      B_ik<-0
      for(p in 1:length(all_arrs)){
        B_ik <- B_ik + (B[p,i]*B[p,k])/(f_sq[p])
      }
      Hess_p[i,k]<- Hess_p[i,k] - B_ik
      Hess_p[k,i]<- Hess_p[i,k]
    }
  }
  return(Hess_p)
}




