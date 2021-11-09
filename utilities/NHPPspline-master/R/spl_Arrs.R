#' Arrival time generation for a cubic spline function composed of cardinal B-splines
#'
#' This thinning-based algorithm returns arrivals from the cubic spline function arrival rate
#' @param knots The knots upon which the B-spline is built - output 2 of the spl_Fit function (spl_Fit$Knots)
#' @param c The spline-coefficients of the spline function - output 1 of the spl_Fit function (spl_Fit$Opt_Spl_Coeffs)
#' @param m The number of days of arrivals required
#' @return Arrivals - An unsorted vector containing the arrivals from day 1 to m one day after another
#' @examples
#' ## spline-based intensity defined on [0,10]
#' d<-3        #spline function degree
#' knots<-seq(-3,13,by=1)
#' n<-length(knots)-(d+1)
#' c<-runif(n,1,5)
#' m<-5
#' ## plot the spline-based function and the arrivals
#' x<-seq(0,10-0.01,by=0.01)
#' plot(x,spl(x,knots,c),lwd=0.5,ylim=c(0,5),type="l")
#' Arrs<-spl_Arrs(knots,c,m)
#' points(sort(Arrs),rep(0,length(Arrs)))

#' @export

## A spline-based method for modelling and generating a nonhomogeneous Poisson process

### Three functions
# 1. To fit a spline-based intensity function from NHPP arrival time observations
# 2. To generate arrivals from the resulting spline-based intensity function
# 3. The spline-based function


## 2.
# Arrival time generation for a cubic spline function composed of cardinal B-splines
# This is a thinning algorithm with a piecewise-linear majorising function

## Inputs
# knots - the knots upon which the B-spline is built - output 2 of the spl_Fit function (spl_Fit$Knots)
# c - the spline-coefficients of the spline function - output 1 of the spl_Fit function (spl_Fit$Opt_Spl_Coeffs)
# m - the number of days of arrivals required

## Outputs
# Arrivals - an unsorted vector containing the arrivals from day 1 to m one day after another

spl_Arrs<-function(knots,c,m){
  d<-3  # cubic
  n<-length(knots)-4
  Tend<-knots[n+1]
  min<-knots[d+1]
  max<-knots[length(knots)-d]

  Arrs<-c()
  eff<-c()     ###efficiency vector
  k_diff<-knots[2]-knots[1]   # how much we translate the arrivals by at each stage
  n<-length(c)
  maxB<-B_spline_val(knots[3],knots,3,1)
  for(i in 1:m){
    Arrs_m<-c()
    for(j in 1:n){
      ### we need to thin each B-spline in turn
      maj_arrs<-PL_maj_inversion(k_diff,maxB,knots[1:5],c[j]) ### arrivals from the B-spline majorising function
      if(length(maj_arrs)>0){
        p_thin<- rep(1,length(maj_arrs)) - B_i_3_b(1,maj_arrs,knots[1:5])*c[j] / PL_maj_b(k_diff,maxB,knots[1:5],c[j],maj_arrs)
        ind<- rbinom(length(maj_arrs), 1, p_thin)
        acc<-which(ind==0)
        Acc<-maj_arrs[acc]
        eff<-c(eff,length(Acc)/length(maj_arrs))
        Acc<-Acc + (j-1)*k_diff   #translation in x axis
        Acc<-Acc[Acc>0 & Acc<=Tend]
        if(length(Acc>0)){Arrs_m<-c(Arrs_m,Acc)
        }
      }
      Arrs_m<-sort(Arrs_m)
    }
    Arrs<-c(Arrs,Arrs_m)
  }
  A<-which(Arrs<=max & Arrs>=min) ### only returns arrivals in interval of interest
  return(Arrs[A])
}

##DESCRIPTION

# generating arrivals from a piecewise linear function
# to be used in the thinning algorithm for generating arrivals from the spline-based intensity

PL_maj<-function(k_diff,maxB,locknot_1,c_j,t){
  #since the B-splines are cubic we have 5 knots in locknot_1
  m<-maxB/k_diff
  c1<- -m*knots[1]
  c2<- m*knots[5]
  if( t <= knots[2] ){
    out<- (m*t + c1)*c_j
  }else if(t>=knots[4]){
    out<- (-m*t + c2)*c_j
  }else{
    out<-maxB*c_j      #flat region in middle
  }
  return(out)
}
PL_maj_b<-Vectorize(PL_maj,"t")

## inversion to generate a single arrival from the cubic B-spline given the last arrival at time t_curr
PL_maj_single_arr<-function(k_diff,maxB,locknot_1,coeff,t_curr){
  lam<-c(0,coeff*maxB,coeff*maxB,0)     ## lam_0, lam_1, ....
  ints<-c(locknot_1[1],locknot_1[2],locknot_1[4],locknot_1[5])   #t_o t_1, ...
  int_ends<-ints[-1]
  ## there are 3 segments, 4 knots involved
  arrs<-c()
  a<-c(  (lam[1]-lam[2])/(ints[1]-ints[2]) , 0 , (lam[3]-lam[4])/(ints[3]-ints[4])  )
  b<-c(lam[1] - ints[1]*a[1],lam[2], lam[2] - ints[3]*a[3])

  u_s<-c(1 - exp( -(a[1]/2)*(ints[2]^2-ints[1]^2)-b[1]*k_diff ),1-exp( -(a[2]/2)*(ints[3]^2-ints[2]^2)-b[2]*k_diff*2 ),1-exp( -(a[3]/2)*(ints[4]^2-ints[3]^2)-b[3]*k_diff ))

  u<-runif(1,0,1)
  k<-max(which(ints<=t_curr))   #which interval we fall in
  u_k<-  1-exp( -(a[k]/2)*(int_ends[k]^2-t_curr^2)-b[k]*(int_ends[k]-t_curr) )

  if(u<=u_k){
    if(a[k]==0){t_next<-t_curr - log(1-u)/b[k]}
    else{
      t_next<- (-b[k] + sqrt( b[k]^2 + a[k]^2*t_curr^2 + 2*a[k]*b[k]*t_curr - 2*a[k]*log(1-u)))/a[k]
    }
    return(t_next)
  }
  else{
    while(u>u_k){
      u<-(u-u_k)/(1-u_k)
      t_curr<-int_ends[k]
      if(t_curr==int_ends[length(int_ends)]){break}
      k<-k+1
      u_k<-u_s[k]
    }
    if(t_curr==int_ends[length(int_ends)])
    {
      return(NULL)
    }else if(a[k]==0){
      t_next<-t_curr - log(1-u)/b[k]
      return(t_next)
    }else{
      t_next<- (-b[k] + sqrt( b[k]^2 + a[k]^2*t_curr^2 + 2*a[k]*b[k]*t_curr - 2*a[k]*log(1-u)))/a[k]
      return(t_next)
    }
  }
}

# generates multiple arrivals t initialised at t=0
PL_maj_inversion<-function(k_diff,maxB,locknot_1,coeff){
  t<-locknot_1[1]
  arrs<-c()
  while( t <  max(locknot_1)){
    a<-PL_maj_single_arr(k_diff,maxB,locknot_1,coeff,t)
    if(is.null(a)==TRUE){break}
    arrs<-c(arrs,a)
    t<-a
  }
  return(arrs)
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

