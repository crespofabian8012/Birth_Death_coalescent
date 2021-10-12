#R code accompanying the paper "Coalescent Models Derived from a
#Birth-Death Process" by Crespo, Posada and Wiuf

# Author of R code: Fausto Fabian Crespo Fernandez
################################################################################################################
rm(list = ls())
################################################################################################################
################################################################################################################
#' relative
#' relative population size at time t given time origin(TOrigin) for model M
#' @param TOrigin: time of origin
#' @param Delta: scaled growth rate (scaled by N/birth_rate )
#' @param t: time in mode time(scaled in units of N/birth_rate )
################################################################################################################
relative <- function(Torigin, Delta, t)
  exp(-Delta * t) * (1 - exp(-Delta * (Torigin - t))) ^ 2 / (1 - exp(-Delta *
                                                                       Torigin)) ^ 2

################################################################################################################
#' lambdaB
#'
#' The  relative population size at time u given time origin(t) for model M*:
#' @param u The time in model time M*
#' @param d The parameter Delta
#' @param g The parameter Gamma
#' @param t The time of origin
#' @return the relative population size for mode M*
################################################################################################################
lambdaB <- function(u, d, g, t) {
  x <- exp(d * u) * (1 - exp(-d * t)) ^ 2 / (1 - exp(-d * (t - u))) ^ 2
  x <-
    (2.0 * d / g) * x * (1 - exp(-g * (1 / (1 - exp(
      -d * u
    )) - 1 / (1 - exp(
      -d * t
    )))))
  return(x)
}
################################################################################################################
#' LambdaB
#'
#' The cumulative relative population size for model M*:
#' @param s The time in model time M*
#' @param Delta The parameter Delta
#' @param Gamma The parameter Gamma
#' @param Torigin The time of origin
#' @return the cummulative relative population size for mode M*
################################################################################################################
LambdaB <- function(s, Delta, Gamma, Torigin) {
  if (s > 0 & s < Torigin) {
    upper.limit = min(s, Torigin)
    res <-
      integrate(
        Vectorize(lambdaB, c("u")),
        d = Delta,
        g = Gamma,
        t = 0.99999999999 * Torigin,
        lower = 0,
        upper = s,
        stop.on.error = FALSE
      )$value
    return(res)
  }
  else{
    if (s > Torigin)
      print(paste("Error!, s=", s, " is greater than Torigin=", Torigin))
    return(0)
  }
}
################################################################################################################
#' functionToZeroB
#'
#' The auxiliar function to solve for the inverse of  cumulative relative population size for model M*:
#' Lambda(x, Delta, Gamma, TOrigin) -y = 0 or Lambda(x, Delta, Gamma, TOrigin) =y
#' @param s The time in model time M*
#' @param Delta The parameter Delta
#' @param Gamma The parameter Gamma
#' @param Torigin The time of origin
#' @return return the function to f(x)= Lambda(x, Delta, Gamma, TOrigin) -y
################################################################################################################
functionToZeroB <- function(y, Delta, Gamma, Torigin)
{
  f <- function(x)
    LambdaB(x, Delta, Gamma, Torigin) - y
  return(f)
}
################################################################################################################
#' InverseLambdaB
#'
#' Inverse of the cumulative relative population size Lambda^{-1}(x, Delta, Gamma, TOrigin) for model M*
#' @param s The time in model time M*
#' @param Delta The parameter Delta
#' @param Gamma The parameter Gamma
#' @param Torigin The time of origin
#' @return return the Lambda^{-1}(x, Delta, Gamma, TOrigin)
################################################################################################################
InverseLambdaB <- function(s, Delta, Gamma, Torigin)
{
  result <-
    cmna::bisection(functionToZeroB(s, Delta, Gamma, Torigin),
                    0,
                    Torigin,
                    tol = 10 ^ (-6))
  return(result)
}
################################################################################################################
#' standard2model
#'
#' Kingsman coalescent  time to  M_0 model time
#' @param Torigin The time of origin
#' @param Delta The parameter Delta
#' @param Gamma The parameter Gamma
#' @param u The time variable u in model time
#' @return return the model time in model M_k with K=0 equivalent to Kingsman coalescent time u
########################################################################################
standard2model <- function(u,Torigin, Delta, Gamma)
  (2*Delta/Gamma)*(log(1 + Delta * u - exp(-Delta * Torigin)) - log(1 + (Delta * u - 1) *
                                                      exp(-Delta * Torigin))) / Delta


################################################################################################################
#' standard2modelB_K
#'
#' from Kingsman coalescent  time to M_K model time
#' @param u The time variable u is in Kingsman coalescent time
#' @param Torigin The time of origin
#' @param Delta The parameter Delta
#' @param Gamma The parameter Gamma
#' @param K The parameter K of the model M_K
#' @return return the model time in model M_K equivalent to Kingsman coalescent time u
################################################################################################################
standard2modelB_K <- function(u, Torigin, Delta, Gamma, K) {
  minus.exp.Delta.Torigin = exp(-Delta * Torigin)
  a = (K / Gamma) * (1 - minus.exp.Delta.Torigin) - minus.exp.Delta.Torigin
  b = 1 - (K / Gamma) * (1 - minus.exp.Delta.Torigin)
  x = exp(u * K / 2.0)
  y = (x - b) / (x * exp(-Delta * Torigin) + a)
  result = (1.0 / Delta) * log(y)
  return(result)
}
################################################################################################################
#' modelB_K2standard
#'
#' from  M_K model time to Kingsman coal time
#' @param u The time variable u is in M_K model time 
#' @param Torigin The time of origin
#' @param Delta The parameter Delta
#' @param Gamma The parameter Gamma
#' @param K The parameter K of the model M_K
#' @return return the Kingsman coal time  equivalent to M_K model time u
################################################################################################################
modelB_K2standard <- function(u, Torigin, Delta, Gamma, K) {
  minus.exp.Delta.Torigin = exp(-Delta * Torigin)
  a = (K / Gamma) * (1 - minus.exp.Delta.Torigin) - minus.exp.Delta.Torigin
  b = 1 - (K / Gamma) * (1 - minus.exp.Delta.Torigin)

  y = (a*exp(Delta*u) + b) / (1-minus.exp.Delta.Torigin *exp(Delta*u))

  result = (2.0 / K) * log(y)
  return(result)
}
#check standard2modelB_K(modelB_K2standard(0.025, 0.039, 100, 100, 0.8), 0.039, 100, 100, 0.8) == 0.025
################################################################################################################
#' modelB2standard
#'
#' from  M_0 model time to Kingsman coal time
#' @param u The time variable u is in M_0 model time 
#' @param Torigin The time of origin
#' @param Delta The parameter Delta
#' @param Gamma The parameter Gamma
#' @return return the Kingsman coal time  equivalent to M_0 model time u
################################################################################################################
modelB02standard <- function(u, Torigin, Delta, Gamma) {
 
  a = exp(Delta*u) -1.0
  b = 1 - exp(-1.0*Delta*Torigin)
  c = 1-  exp(-1.0*Delta*(Torigin -u))
  y = (a* b) / (Gamma *c)
  
  result = y
  return(result)
}
#check standard2modelB0(modelB02standard(0.025, 0.039, 100, 100), 0.039, 100, 100) == 0.025
################################################################################################################
#' standard2modelB0
#'
#' from    Kingsman coal time to M_0 model time
#' @param u The time variable u is in Kingsman coal  time 
#' @param Torigin The time of origin
#' @param Delta The parameter Delta
#' @param Gamma The parameter Gamma
#' @return return the  M_0 model time u
################################################################################################################
standard2modelB0 <- function(u, Torigin, Delta, Gamma) {
  
  a = exp(-1.0 * Delta * Torigin)
  y = log(1 + Delta * u - a) - log(1 + (Delta * u - 1) * a);
  
  result = (1.0/Gamma)*y
  return(result)
}
#check standard2modelB0(modelB02standard(0.025, 0.039, 100, 100), 0.039, 100, 100) == 0.025
################################################################################################################
#' standard2modelB
#'
#' Kingsman coalescent  time to time in model  M_K
#' @param u The time variable u in Kingsman coalescent time
#' @param Torigin The time of origin
#' @param Delta The parameter Delta
#' @param Gamma The parameter Gamma
#' @param K The parameter K of the model M*
#' @return return the model time in model M* equivalent to Kingsman coalescent time u
################################################################################################################
standard2modelB <- function(u, Torigin, Delta, Gamma, K) {
  minus.exp.Delta.Torigin = exp(-Delta * Torigin)
  a = (K / Gamma) * (1 - minus.exp.Delta.Torigin) - minus.exp.Delta.Torigin
  b = 1 - (K / Gamma) * (1 - minus.exp.Delta.Torigin)
  x = (exp(u * K / 2.0) - b) / (exp(u * K / 2.0 - Delta * Torigin) + a)
  result = (1.0 / Delta) * log(x)
  return(result)
}
################################################################################################################
#' standard2modelB_Mstar
#'
#' Kignman coalescent  time to time in model  M*
#' @param u The time variable u in Kingsman coalescent time
#' @param Torigin The time of origin
#' @param Delta The parameter Delta
#' @param Gamma The parameter Gamma
#' @return return the model time in model M* equivalent to Kingsman coalescent time  u
################################################################################################################
standard2modelB_Mstar <- function(u, Torigin, Delta, Gamma) {
  result = InverseLambdaB(u, Delta, Gamma, Torigin)
  return(result)
}
################################################################################################################
#' standardCoal
#'
#' Kignman coalescent(with constant population size)  inter event times for sample size n
#' @param n The sample size at time 0(present time)
#' @return return the list of coalescent times inter event times
################################################################################################################
standardCoal <-
  function(n)
    unlist(lapply(2:n, function(i)
      rexp(1, rate = choose(i, 2))))
################################################################################################################
#' standardCoalFrom
#'
#' Kignman coalescent inter event times for sample size n to sample size p(p<n)
#' @param p The sample size at some time earlier than present time(p<n)
#' @param n The sample size at time 0(present time)
#' @return return the list of coalescent times inter event times
################################################################################################################
standardCoalFrom <-
  function(p, n)
    unlist(lapply(p:n, function(i)
      rexp(1, rate = choose(i, 2))))

################################################################################################################
#' modelCoal
#'
#' Coalescent  times from sample size n to sample size 1 for model M_0
#' @param n The sample size at time 0(present time)
#' @param Torigin The time of origin
#' @param Delta The parameter Delta
#' @param Gamma The parameter Gamma
#' @return return the list of coalescent times in model M_0
################################################################################################################
modelCoal <- function(n, Torigin, Delta, Gamma) {
  w <- standardCoal(n)
  # time until i=n-1,...,1 ancestors
  cumulativeSTDCoaltimes <-
    unlist(lapply(1:(n - 1), function(i)
      sum(w[i:(n - 1)])))
  u <-
    unlist(lapply(
      cumulativeSTDCoaltimes,
      FUN = function(x)
        standard2model(x, Torigin, Delta, Gamma)
    ))
  return(u)
}
################################################################################################################
#' modelCoal_race
#'
#' Coalescent  times from sample size n to sample size 1 for model M_0
#' @param n The sample size at time 0(present time)
#' @param Torigin The time of origin
#' @param Delta The parameter Delta
#' @param Gamma The parameter Gamma
#' @return return the list of coalescent times in model M_0
################################################################################################################
modelCoal_race <- function(n, Torigin, Delta, Gamma) {
  #example modelCoal_race(5, 0.039, 100, 100)
  
  Kingman.coal.times<- unlist(lapply(1:choose(n, 2), function(i)
    rexp(1, rate =1)))
  idx<- 1:n
  pairs<-combn(idx, 2, paste, collapse=",")
  
  coal.model.times <- unlist(lapply(
    Kingman.coal.times,
    FUN = function(x)
      standard2modelB0(x, Torigin, Delta, Gamma)
  ))
  idx.next = n+1
  active.ind<- unlist(lapply(1:n, FUN=function(x)  toString(x)))
  result= list()
  pos =1
  curr.model.time = 0
  curr.Kingman.model.time  = 0
  iteration = 1
  while(iteration < n)
  {
    numpairs = choose(n-iteration+1, 2)
    stopifnot(length(active.ind)==(n-iteration+1))
    stopifnot(length(coal.model.times)==numpairs)
    stopifnot(length(pairs)==numpairs)
    
    winner <- min(coal.model.times)
    pos.winner <-which(coal.model.times == winner)[1]
    winner.pair <-pairs[pos.winner]
    subSplits=  unlist(strsplit(winner.pair, ","))
    
    left = subSplits[1]
    right = subSplits[2]
    active.ind<-active.ind[-which(active.ind==left)]
    active.ind<-active.ind[-which(active.ind==right)]
    pairs<- pairs[-pos.winner]; 
    result[[pos]] <- coal.model.times[pos.winner]
    pos = pos +1
    
    
    curr.model.time = coal.model.times[pos.winner]
    curr.Kingman.model.time = modelB02standard(curr.model.time, Torigin, Delta, Gamma)
    
    coal.model.times <- coal.model.times[-pos.winner]
    
    dropout.idxs<-unlist(lapply( pairs, FUN = function(x, left, right){
      list.left.right = unlist(strsplit(x, ","))
      #grepl( left, x, fixed = TRUE)|| grepl( right, x, fixed = TRUE)
      (left %in% list.left.right) || (right %in% list.left.right)
    },left, right))
    
    
    stopifnot(length(pairs)==length(dropout.idxs))
    if (length(pairs)==0){
      
      break
    }
    pairs = pairs[!dropout.idxs]
    coal.model.times <- coal.model.times[!dropout.idxs]
    str.new.ind = toString(idx.next)
    new.pairs = unlist(lapply(active.ind, FUN= function(x, str.new.ind) paste(x, str.new.ind, sep=","), str.new.ind ))
    active.ind= append(active.ind, toString(idx.next))
    pairs = append(pairs, new.pairs)
    new.Kingman.model.times<-curr.Kingman.model.time+unlist(lapply(1:length(new.pairs), function(i)
      rexp(1, rate =1)))
    new.coal.model.times<- unlist(lapply(
      new.Kingman.model.times,
      FUN = function(x)
        standard2modelB0(x, Torigin, Delta, Gamma)
    ))
    coal.model.times<-append(coal.model.times, new.coal.model.times)
    idx.next= idx.next+1 
    iteration = iteration +1  
  }
  result1 = unlist(result)
  if (length(result1)!=(n-1)){
    
    print(paste("length ",length(result1)))
    print(paste("result ",result1))
    print(paste("Delta ",Delta))
    print(paste("Gamma ",Gamma))
    print(paste("Torigin ",Torigin))
    print(paste("n ",n))
    print(paste("K ",K))
  }
  stopifnot(length(result1)==(n-1))
  return(result1)
}
################################################################################################################
#' modelCoalB_K
#'
#' Coalescent  times from sample size n to sample size 1
#' @param n The sample size at time 0(present time)
#' @param Torigin The time of origin
#' @param Delta The parameter Delta
#' @param Gamma The parameter Gamma
#' @param K The parameter K of the model M_K
#' @return return the list of coalescent times in model M_K
################################################################################################################
modelCoalB_K <- function(n, Torigin, Delta, Gamma, K) {
  w <- standardCoal(n)
  # time until i=n-1,...,1 ancestors
  cumulativeSTDCoaltimes <-
    unlist(lapply(1:(n - 1), function(i)
      sum(w[i:(n - 1)])))
  u <-
    unlist(lapply(
      cumulativeSTDCoaltimes,
      FUN = function(x)
        standard2modelB_K(x, Torigin, Delta, Gamma, K)
    ))
  return(u)
}
################################################################################################################
#' modelCoalB_K_race
#'
#' Coalescent  times from sample size n to sample size 1
#' @param n The sample size at time 0(present time)
#' @param Torigin The time of origin
#' @param Delta The parameter Delta
#' @param Gamma The parameter Gamma
#' @param K The parameter K of the model M_K
#' @return return the list of coalescent times in model M_K
################################################################################################################
modelCoalB_K_race <- function(n, Torigin, Delta, Gamma, K) {
  #example modelCoalB_K_race(5, 0.039, 100, 100, 0.8)
  # modelCoalB_K_race(10, 0.2882, 10, 10, 0.8) 
  
  Kingman.coal.times<- unlist(lapply(1:choose(n, 2), function(i)
    rexp(1, rate =1)))
  idx<- 1:n
  pairs<-combn(idx, 2, paste, collapse=",")
       
  coal.model.times <- unlist(lapply(
    Kingman.coal.times,
    FUN = function(x)
      standard2modelB_K(x, Torigin, Delta, Gamma, K)
  ))
  idx.next = n+1
  active.ind<- unlist(lapply(1:n, FUN=function(x)  toString(x)))
  result= list()
  pos =1
  curr.model.time = 0
  curr.Kingman.model.time  = 0
  iteration = 1
  while(iteration < n)
    {
      numpairs = choose(n-iteration+1, 2)
      stopifnot(length(active.ind)==(n-iteration+1))
      stopifnot(length(coal.model.times)==numpairs)
      stopifnot(length(pairs)==numpairs)
      
      winner <- min(coal.model.times)
      pos.winner <-which(coal.model.times == winner)[1]
      winner.pair <-pairs[pos.winner]
      subSplits=  unlist(strsplit(winner.pair, ","))
      
      left = subSplits[1]
      right = subSplits[2]
      active.ind<-active.ind[-which(active.ind==left)]
      active.ind<-active.ind[-which(active.ind==right)]
      pairs<- pairs[-pos.winner]; 
      result[[pos]] <- coal.model.times[pos.winner]
      pos = pos +1
    
         
      curr.model.time = coal.model.times[pos.winner]
      curr.Kingman.model.time = modelB_K2standard(curr.model.time, Torigin, Delta, Gamma, K)
      
      coal.model.times <- coal.model.times[-pos.winner]

      dropout.idxs<-unlist(lapply( pairs, FUN = function(x, left, right){
        list.left.right = unlist(strsplit(x, ","))
        #grepl( left, x, fixed = TRUE)|| grepl( right, x, fixed = TRUE)
        (left %in% list.left.right) || (right %in% list.left.right)
      },left, right))
      
      
      stopifnot(length(pairs)==length(dropout.idxs))
      if (length(pairs)==0){
        
        break
      }
      pairs = pairs[!dropout.idxs]
      coal.model.times <- coal.model.times[!dropout.idxs]
      str.new.ind = toString(idx.next)
      new.pairs = unlist(lapply(active.ind, FUN= function(x, str.new.ind) paste(x, str.new.ind, sep=","), str.new.ind ))
      active.ind= append(active.ind, toString(idx.next))
      pairs = append(pairs, new.pairs)
      new.Kingman.model.times<-curr.Kingman.model.time+unlist(lapply(1:length(new.pairs), function(i)
        rexp(1, rate =1)))
      new.coal.model.times<- unlist(lapply(
        new.Kingman.model.times,
        FUN = function(x)
          standard2modelB_K(x, Torigin, Delta, Gamma, K)
      ))
      coal.model.times<-append(coal.model.times, new.coal.model.times)
      idx.next= idx.next+1 
      iteration = iteration +1  
      }
  result1 = unlist(result)
  if (length(result1)!=(n-1)){

     print(paste("length ",length(result1)))
     print(paste("result ",result1))
     print(paste("Delta ",Delta))
     print(paste("Gamma ",Gamma))
     print(paste("Torigin ",Torigin))
     print(paste("n ",n))
     print(paste("K ",K))
  }
  stopifnot(length(result1)==(n-1))
  return(result1)
}
################################################################################################################
#' modelCoalB_MAster
#'
#' Coalescent  times from sample size n to sample size 1
#' @param n The sample size at time 0(present time)
#' @param Torigin The time of origin
#' @param Delta The parameter Delta
#' @param Gamma The parameter Gamma
#' @return return the list of coalescent times in model M*
################################################################################################################
modelCoalB_MAster <- function(n, Torigin, Delta, Gamma) {
  w <- standardCoal(n)
  # time until i=n-1,...,1 ancestors
  cumulativeSTDCoaltimes <-
    unlist(lapply(1:(n - 1), function(i)
      sum(w[i:(n - 1)])))
  u <-
    unlist(lapply(
      cumulativeSTDCoaltimes,
      FUN = function(x)
        standard2modelB_Mstar(x, Torigin, Delta, Gamma)
    ))
  return(u)
}
################################################################################################################
#' modelCoalHybrid
#'
#' Coalescent  times from sample size n to sample size nB for Kingsman coalescent  model
#' @param n The sample size at time 0(present time)
#' @param Torigin The time of origin
#' @param Delta The parameter Delta
#' @param Gamma The parameter Gamma
#' @param nB The sample size when the Hybrid model switchs from M_K to BD coalescent
#' @return return the list of coalescent times from sample size n to sample size nB for Kingsman coalescent model
################################################################################################################
modelCoalHybrid <- function(n, Torigin, Delta, Gamma, nB) {
  w <- standardCoalFrom(nB + 1, n)
  u <-  unlist(lapply(1:(n - nB), function(i)
    sum(w[i:(n - nB)])))
  return(u)
}
################################################################################################################
#' expCoal
#'
#' Coalescent waiting times from sample size n to sample size 1
#' @param n The sample size at time 0(present time)
#' @param Delta The exponent parameter in N(t)=N(0) exp(-Delta*t) (time backwards)
#' @param Gamma The exponent parameter in N(t)=N(0) exp(-Delta*t) (time backwards)
#' @return return the list of coalescent times when the population is growing exponentially
################################################################################################################
expCoalTimes <- function(n, Delta, Gamma) {
  w <- standardCoal(n)
  # time until i=n-1,...,1 ancestors
  u <- unlist(lapply(1:(n - 1), function(i)
    sum(w[i:(n - 1)])))
  u <- log(1 + 0.5*Gamma * u) / Delta
  return(unlist(u))
}
################################################################################################################
#' expCoalWaitingTimes
#'
#' Coalescent waiting  times from sample size n to sample size 1
#' @param n The sample size at time 0(present time)
#' @param beta The exponent parameter in N(t)=N(0) exp(-beta*t) (time backwards)
#' @return return the list of coalescent waiting times when the population is growing exponentially
################################################################################################################
expCoalWaitingTimes <- function(n, beta) {
  w <- standardCoal(n)
  # time until i=n-1,...,1 ancestors
  u <- unlist(lapply(1:(n - 1), function(i)
    sum(w[i:(n - 1)])))
  u <- log(1 + beta * u) / beta
  u <- c(u, 0)
  unlist(lapply(1:(n - 1), function(i)
    u[i] - u[i + 1]))
}
################################################################################################################
#' log.prod.between
#'
#' The logarithm of the product of the positive  numbers between from and to(including them)
#' i.e log(k*(k+1)*...*n) where k is the from parameter and n is the to parameter
#' i.e the sum of the logarithm of the factors in the product
#' @param from The first number in the product(if it is less than or equal to 0 then the initial number is
#' taken as 1)
#' @param to The last number of the product
#' @return return the sum of the logarithm of the numbers(including the from and to)
################################################################################################################
log.prod.between = function(from, to) {
  result = 0
  stopifnot(from >= 0)
  stopifnot(to >= 0)
  
  
  if (from > to)
    return(result)
  
  if (min(from, to) > 0)
  {
    if (to >= from) {
      return(sum(log(seq(
        max(from, 1), to, by = 1
      ))))
    }
  }
  else
    return(result)
}
################################################################################################################
#' log.prod.from.To
#'
#' The logarithm of  (n!/k!) where k is the  parameter from and n is the to parameter to
#' i.e the sum of the logarithm of the factors from k+1 to n
#' @param from The first number in the product(if it is less than or equal to 0 then the initial
#' number is taken as 1)
#' @param to The last number of the product
#' @return return the sum of the logarithm of the numbers(including  the last one but not the first one)
################################################################################################################
log.prod.from.To = function(from, to) {
  result = 0
  stopifnot(from >= 0)
  stopifnot(to >= 0)
  
  
  if (min(from, to) >= 0)
  {
    if (to > from) {
      return(sum(log(seq(
        max(from + 1, 1), to, by = 1
      ))))
    }
    else if (from > to)
    {
      return(-1.0 * sum(log(seq(
        max(to + 1, 1), from, by = 1
      ))))
    }
    else
    {
      result = 0
      return(result)
    }
  }
  else
    return(result)
}
################################################################################################################
#' sample.conditional.coalescent.time.distribution
#'
#' Sample the coalescent time from the Beta distribution
#' @param Delta The parameter Delta
#' @param alpha The last number of the product
#' @param k The current sample size
#' @param m The current population size
#' @return return the coalescent time when
################################################################################################################
sample.conditional.coalescent.time.distribution = function(Delta, alpha, k, m) {
  require(stats)
  stopifnot(k > 0)
  stopifnot(m - k + 1 > 0)
  U = rbeta(1, k, m - k + 1)
  stopifnot(U / (alpha + (1 - alpha) * U) > 0)
  T = (-1.0 / Delta) * log(U / (alpha + (1 - alpha) * U))
  return(T)
}
################################################################################################################
#' sample.first.coalescent.time.distribution
#'
#' Sample the coalescent time from the Gamma distribution
#' @param Delta The parameter Delta
#' @param Gamma The last number of the product
#' @param k The current sample size
#' @return return the coalescent time
################################################################################################################
sample.first.coalescent.time.distribution = function(Delta, Gamma, k) {
  require(stats)
  U = stats::rgamma(1, shape = k, rate = 1)
  T = (-1.0 / Delta) * log(1.0 - (Gamma / (U + Gamma)))
  return(T)
}
################################################################################################################
#' conditional.density.first.time.k.ancestors
#'
#' the conditional density of the coalescent time when there are k  population ancestors
#' @param delta The parameter Delta
#' @param gamma The parameter Gamma
#' @param n The sample size
#' @param time.Origin The time of origin
#' @param times
#' @return return the density of the coalescent time
################################################################################################################
conditional.density.first.time.k.ancestors = function(delta, gamma, n, time.Origin, times) {
  exp.delta.t = exp(-1 * delta * time.Origin)
  output = ((delta * gamma) ** (n - 1)) * (1 + (gamma - 1) * exp.delta.t) **
    (n - 1) / (1 - exp.delta.t) ** (n - 1)
  
  for (i in  2:n)
  {
    if (i < n) {
      output = i * output * exp(-1 * delta * time.Origin) /  (1 + (gamma - 1) *
                                                                (exp.delta.t) ** 2)
    }
    else{
      output = output * exp(-1 * delta * time.Origin) /  (1 + (gamma - 1) * (exp.delta.t) **
                                                            2)
    }
  }
  return(output)
}
################################################################################################################
#' density.first.time.k.ancestors
#'
#' Sample the coalescent time from the Beta distribution
#' @param Delta The parameter Delta
#' @param alpha The last number of the product
#' @param k The current sample size
#' @param m The current population size
#' @return return the coalescent time when
################################################################################################################
density.first.time.k.ancestors = function(delta, gamma, n, time.Origin, times) {
  exp.delta.t = exp(-1 * delta * time.Origin)
  output = ((delta * gamma) ** (n))
  for (i in  1:n)
  {
    exp.delta.t = exp(-1 * delta * times[i])
    output = i * output * exp(-1 * delta * time.Origin) /  (1 + (gamma -
                                                                   1) * ((exp.delta.t) ** 2))
    
  }
  return(output)
}
solve.cumulative.density.number.ancestors.population.present.time = function(u, sample.size)
{
  fun <-
    function (x)
      exp(log.prod.between(x - sample.size + 2, x) - log.prod.between(x + 2, x +
                                                                        sample.size)) - u
  fun2 <-
    function (x)
      (u - exp(
        log.prod.between(x - sample.size + 2, x) - log.prod.between(x + 2, x + sample.size)
      )) ^ 2
  
  result <- search_infinite_from(0, (sample.size - 2), fun)
  
  result <- ceiling(result)
  return(result)
  
}
solve.cumulative.density.number.ancestors.population = function(u, k, m)
{
  fun <- function (x) {
    result = exp(
      log.prod.between(x - k + 2, x) - log.prod.between(x + 2, x + k) + log.prod.from.To(m -
                                                                                           1, m + k - 1) +
        log.prod.from.To(m, m - k)
    ) - u
    return(result)
  }
  fun2 <- function (x) {
    result = (
      u - exp(
        log.prod.between(x - k + 2, x) - log.prod.between(x + 2, x + k) + log.prod.from.To(m -
                                                                                             1, m + k - 1) +
          log.prod.from.To(m, m - k)
      )
    ) ^ 2
    return(result)
  }
  
  interval <- seq((k - 2), m, by = 1)
  
  result <- nearest_value_search(0, interval, fun)
  
  result <- ceiling(result)
  
  return(result)
}
################################################################################################################
#' get_number_of_ancestors_population_when_sample_size_minus_1_sample
#' sample the  number of ancestors in the coalescent approximation
#' when there are sample_size-1 ancestors in the coalescent appproximation
#' in the sample. The number ancestors in the population is assumed to be infinity
#' @param sample_size: current number ancestors in the sample(present time)
#' @return: the number of ancestors when there are k-1 ancestors in the sample  the output >=sample_size-1
################################################################################################################
get_number_of_ancestors_population_when_sample_size_minus_1_sample = function(sample_size)
{
  accepted = FALSE
  while (!accepted) {
    u = runif(1, 0, 1)
    mprime = solve.cumulative.density.number.ancestors.population.present.time(u, sample_size)
    if (mprime >= (sample_size - 1)) {
      accepted = TRUE
    }
  }
  return(mprime)
}
################################################################################################################
#' get_number_of_ancestors_population_when_k_minus_1_ancestors_sample
#' sample the number of ancestors in the BD coalescent model approximation
#' when there are k-1 ancestors in the coalescent appproximation
#' in the sample
#' @param k: current number ancestors in the sample
#' @param m: current number ancestors in the population
#' @return: the number of ancestors when there are k-1 ancestors in the sample the output >=k-1 and < m
################################################################################################################
get_number_of_ancestors_population_when_k_minus_1_ancestors_sample = function(k, m)
{
  accepted = FALSE
  while (!accepted)
  {
    u = runif(1, 0, 1)
    mprime = solve.cumulative.density.number.ancestors.population(u, k, m)
    
    if (mprime < m & mprime >= (k - 1))
    {
      accepted = TRUE
    }
  }
  return(mprime)
}
################################################################################################################
#' simulate.list.number.ancestors.population
#' simulate the number of ancestors in the BD coalescent model approximation when there are  n-1, n-2,...,1
#' in the sample
#' where n is the sample size
#' @param sample.size: the sample size
#' @return: the list of population sizes
################################################################################################################
simulate.list.number.ancestors.population = function(sample.size) {
  population.present.time = get_number_of_ancestors_population_when_sample_size_minus_1_sample(sample.size)
  current.population.size.when.sample.minus1 = population.present.time
  list.populations.sizes = list(current.population.size.when.sample.minus1)
  m = current.population.size.when.sample.minus1
  rest.population.sizes = list(current.population.size.when.sample.minus1)
  p = 2
  for (k in (sample.size - 1):1)
  {
    m = get_number_of_ancestors_population_when_k_minus_1_ancestors_sample(k, m)
    list.populations.sizes[[p]] = m
    p = p + 1
  }
  
  return(unlist(list.populations.sizes))
}
################################################################################################################
#' simulate.list.number.ancestors.population.from
#' simulate the number of ancestors in the BD coalescent model approximation when there are  current.sample,
#' current.sample-1,...,1 in the sample
#' @param sample.size: the sample size at t=0
#' @param current.population.size: the sample size
#' @param current.sample: the current sample size
#' @return: the list of population sizes
################################################################################################################
simulate.list.number.ancestors.population.from = function(sample.size,
                                                          current.population.size,
                                                          current.sample) {
  stopifnot(current.sample <= sample.size)
  
  list.populations.sizes = list(current.population.size)
  if (sample.size == current.sample) {
    current.population.size.when.sample.minus1 =
      get_number_of_ancestors_population_when_sample_size_minus_1_sample(sample.size)
    m = current.population.size.when.sample.minus1
  }
  else{
    m = current.population.size
    
  }
  list.populations.sizes = list(m)
  
  p = 2
  for (k in (current.sample):2)
  {
    m = get_number_of_ancestors_population_when_k_minus_1_ancestors_sample(k, m)
    list.populations.sizes[[p]] = m
    p = p + 1
  }
  
  return(unlist(list.populations.sizes))
}
################################################################################################################
#' simulate.coalescent.times.A
#' simulate the coalescent times of the sample from sample.size to 1  using BD coalescent approximation
#' @param Gamma: Gamma parameter
#' @param Delta: Delta parameter
#' @param sample.size: the  sample size  at the t=0
#' @param list.number.ancestors.population: the  list of number ancestor population when there are n, n-1,..,1
#' in the sample
#' @return: the list of coalescent times
################################################################################################################
simulate.coalescent.times.A = function(Gamma,
                                       Delta,
                                       sample.size,
                                       list.number.ancestors.population)
{
  list.coal.times = list()
  m = list.number.ancestors.population[1]
  u = unlist(sample.first.coalescent.time.distribution(Delta, Gamma, m +
                                                         1))
  u = u[[1]]
  current.time = u
  list.coal.times = list(u)
  p = 2
  for (i in  (sample.size - 2):0)
  {
    m = list.number.ancestors.population[sample.size - i]
    current.m = list.number.ancestors.population[sample.size - i - 1]
    alpha = 1 - exp(-1.0 * Delta * current.time)
    u = unlist(sample.conditional.coalescent.time.distribution(Delta, alpha, m +
                                                                 1, current.m))
    u = u[[1]]
    current.time = current.time + u
    list.coal.times[p] = current.time
    p = p + 1
  }
  return(unlist(list.coal.times))
}
################################################################################################################
#' simulate.coalescent.times.A.from
#' simulate the coalescent times of the sample from current.sample to 1 using BD coalescent approximation
#' @param Gamma: Gamma parameter
#' @param Delta: Delta parameter
#' @param sample.size: the  sample size at t=0
#' @param list.number.ancestors.population: the  list of number ancestor population when there are current.sample,
#' current.sample-1,..,1  in the sample
#' @param current.sample: the  current sample(less than sample.size)
#' @param from_time: the  current time
#' @return: the list of coalescent times
################################################################################################################
simulate.coalescent.times.A.from = function(Gamma,
                                            Delta,
                                            sample.size,
                                            list.number.ancestors.population,
                                            current.sample,
                                            from_time)
{
  stopifnot(current.sample <= sample.size)
  list.coal.times = list()
  m = list.number.ancestors.population[1]
  if (current.sample == sample.size) {
    u = unlist(sample.first.coalescent.time.distribution(Delta, Gamma, m + 1))
    u = u[[1]]
    list.coal.times = list(u)
    p = 2
    
  }
  else{
    u = from_time
    list.coal.times = list()
    p = 1
  }
  
  current.time = u
  
  for (i in  (current.sample - 2):0)
  {
    m = list.number.ancestors.population[current.sample - i]
    current.m = list.number.ancestors.population[current.sample - i - 1]
    alpha = 1 - exp(-1.0 * Delta * current.time)
    u = unlist(sample.conditional.coalescent.time.distribution(Delta, alpha, m +
                                                                 1, current.m))
    u = u[[1]]
    current.time = current.time + u
    list.coal.times[p] = current.time
    p = p + 1
  }
  return(unlist(list.coal.times))
}
################################################################################################################
#' simulate.coalescent.times.A.test
#' compute the theoretical expected coalescent times usint the simulated numbers of ancestors
#' @param Gamma: Gamma parameter
#' @param Delta: Delta parameter
#' @param sample.size: the  sample size at t=0
#' @param list.number.ancestors.population: the  list of number ancestor population when there are current.sample,
#' current.sample-1,..,1  in the sample
#' @return: the list of coalescent times
################################################################################################################
simulate.coalescent.times.A.test = function(Gamma,
                                            Delta,
                                            sample.size,
                                            list.number.ancestors.population)
{
  p = 1
  list.coal.times <- list()
  list.coal.times <-
    lapply((sample.size - 1):0,
           FUN = function(i,
                          sample.size,
                          list.number.ancestors.population)
           {
             2.0 * (1.0 / list.number.ancestors.population[sample.size - i])
             
           },
           sample.size,
           list.number.ancestors.population
    )
  
  return(unlist(list.coal.times))
}
################################################################################################################
#' expected_number_ancestors
#' compute the theoretical expected  numbers of ancestors
#' @param sample.size: the  sample size at t=0
#' @return: the list of theoretical number of ancestors
################################################################################################################
expected_number_ancestors = function(sample.size) {
  list.number.ancestors <-
    lapply((sample.size - 2):1,
           FUN = function(k, sample.size)
           {
             sample.size * k / (sample.size - k - 1)
           },
           sample.size
    )
  
  return(unlist(list.number.ancestors))
}
################################################################################################################
#' compute_stats_number_ancestors
#' compute stats for the numbers of ancestors
#' @param sample.size: the  sample size at t=0
#' @param sim: the  number of simulations used
#' @param DeltaList: the  list of Delta values used to simulate number.ancestors.simA
#' @param number.ancestors.simA: the  list of simulated number of ancestors
#' @return: the list of data frames with statistics
################################################################################################################
compute_stats_number_ancestors = function(sample.size,
                                          sim,
                                          DeltaList,
                                          number.ancestors.simA) {
  list.Data.Frames.Numbers.Ancestors <-
    parallel::mclapply(
      1:length(DeltaList),
      mc.cores = parallel::detectCores() - 1,
      FUN = function(k,
                     sample.size,
                     sim,
                     number.ancestors.simA) {
        quants <- c(0.025, 0.50, 0.975)
        
        matrixCurrentValue <-
          matrix(
            number.ancestors.simA[k, ],
            nrow = sim ,
            ncol = sample.size - 1,
            byrow = TRUE
          )
        
        quantiles <-
          apply(matrixCurrentValue , 2 , quantile , probs = quants)
        mean.Number.Ancestors <- colMeans(matrixCurrentValue)
        result <-
          data.frame(
            mean = mean.Number.Ancestors,
            LI = quantiles[1, ],
            median = quantiles[2, ],
            UI = quantiles[3, ]
          )
        result
      },
      sample.size,
      sim,
      number.ancestors.simA
    )
  
  return(list.Data.Frames.Numbers.Ancestors)
  
}
################################################################################################################
#' simulateKingmanCoal.parallel
#' simulate coalescent times for model Kingman parallel
#' @param GammaList: the  paramater Gamma
#' @param sim: the  number of simulations sim
#' @param sample.size: the  sample size at time t=0
#' @param coal.events.times.simKingman: a big matrix bigstatsr::FBM(length(DeltaList),sim*(sample.size-1))
#' of size to store the simulated coalescent times
#' @return: a list of data frame(sample.size -1 columns ) with the stats for the coalescent times
#' for Kingman coal(one data frame per each Gamma value)
################################################################################################################
simulateKingmanCoal.parallel = function(DeltaList,GammaList,
                                sim,
                                sample.size,
                                coal.events.times.simKingman)
{
  require(foreach)
  require(doParallel)
  require(doSNOW)
  require(parallel)
  require(doFuture)
  require(bigstatsr)
  require(stats)
  require(rgenoud)
  require(doRNG)
  require(doMC)
  
  RNGkind("L'Ecuyer-CMRG")
  a=as.numeric(Sys.time())
  set.seed(floor(runif(1)*a)) #set seed to something



  ncores = parallel::detectCores() - 1
  innerCluster <-
    parallel::makeCluster(ncores, type = "FORK", outfile = "")
  on.exit(parallel::stopCluster(innerCluster), add = TRUE)
  doParallel::registerDoParallel(innerCluster)
  
  
  rng <- RNGseq(length(DeltaList) * sim, floor(runif(1)*a))
  opts <- list(chunkSize = 2)
  
  tmp3 <- foreach::foreach(j = 1:sim, .combine = 'c') %:%
    foreach::foreach(i = 1:length(GammaList),
                     r = rng[(j - 1) * length(GammaList) + 1:length(GammaList)],
                     .combine = 'c') %dopar% {
                       rngtools::setRNG(r)

                       
                      coal.events.times = (2*GammaList[i] /DeltaList[i])*standardCoal(sample.size)
                       
                       print(
                         paste(
                           "finished Kingman coal in parallel for ",
                           "Gamma",
                           GammaList[i],
                           " sim ",
                           j,
                           sep = " "
                         )
                       )
                       positions <-
                         cbind(rep(i, (sample.size - 1)), ((j - 1) * (sample.size - 1) + 1):(j *
                                                                                               (sample.size - 1)))
                       coal.events.times.simKingman[positions] <- rev(coal.events.times)
                       NULL
                     }
  list.Data.Frames <-
    parallel::mclapply(
      1:length(GammaList),
      mc.set.seed = TRUE,
      mc.cores = parallel::detectCores() - 1,
      FUN = function(k,
                     sample.size,
                     sim,
                     coal.events.times.simKingman)
      {
        quants <- c(0.025, 0.50, 0.975)
        matrixCurrentValue <-
          matrix(
            coal.events.times.simKingman[k, ],
            nrow = sim ,
            ncol = sample.size - 1,
            byrow = TRUE
          )
        #quantiles<-apply( matrixCurrentValue , 2 , quantile , probs = quants , na.rm = TRUE )
        quantiles <-
          apply(matrixCurrentValue , 2 , quantile , probs = quants)
        meanCoalTimes <- colMeans(matrixCurrentValue)
        result <-
          data.frame(
            mean = meanCoalTimes,
            LI = quantiles[1, ],
            median = quantiles[2, ],
            UI = quantiles[3, ]
          )
        result
      },
      sample.size,
      sim,
      coal.events.times.simKingman
    )
  
  return(list.Data.Frames)
}
################################################################################################################
#' simulateExponentialCoal.parallel
#' simulate coalescent times for model Kingman parallel
#' @param GammaList: the  parameter Gamma
#' @param sim: the  number of simulations sim
#' @param sample.size: the  sample size at time t=0
#' @param coal.events.times.sim.Exp: a big matrix bigstatsr::FBM(length(DeltaList),sim*(sample.size-1))
#' of size to store the simulated coalescent times
#' @return: a list of data frame(sample.size -1 columns ) with the stats for the coalescent times
#' for Exponential coal(one data frame per each Gamma value)
################################################################################################################
simulateExponentialCoal.parallel = function(DeltaList, GammaList,
                                        sim,
                                        sample.size,
                                        coal.events.times.sim.Exp)
{
  require(foreach)
  require(doParallel)
  require(doSNOW)
  require(parallel)
  require(doFuture)
  require(bigstatsr)
  require(stats)
  require(rgenoud)
  require(doRNG)
  require(doMC)
  
  RNGkind("L'Ecuyer-CMRG")
  a=as.numeric(Sys.time())
  set.seed(floor(runif(1)*a)) #set seed to something
 
  ncores = parallel::detectCores() - 1
  innerCluster <-
    parallel::makeCluster(ncores, type = "FORK", outfile = "")
  on.exit(parallel::stopCluster(innerCluster), add = TRUE)
  doParallel::registerDoParallel(innerCluster)
  
  
  rng <- RNGseq(length(DeltaList) * sim, floor(runif(1)*a))
  opts <- list(chunkSize = 2)
  
  tmp3 <- foreach::foreach(j = 1:sim, .combine = 'c') %:%
    foreach::foreach(i = 1:length(GammaList),
                     r = rng[(j - 1) * length(GammaList) + 1:length(GammaList)],
                     .combine = 'c') %dopar% {
                       rngtools::setRNG(r)
                       
                       
                       coal.events.times = expCoalTimes(sample.size, DeltaList[i], GammaList[i])
                       
                       print(
                         paste(
                           "finished exponential coal in parallel for ",
                           "Gamma",
                           GammaList[i],
                           " sim ",
                           j,
                           sep = " "
                         )
                       )
                       positions <-
                         cbind(rep(i, (sample.size - 1)), ((j - 1) * (sample.size - 1) + 1):(j *
                                                                                               (sample.size - 1)))
                       coal.events.times.sim.Exp[positions] <- rev(coal.events.times)
                       NULL
                     }
  list.Data.Frames <-
    parallel::mclapply(
      1:length(GammaList),
      mc.set.seed = TRUE,
      mc.cores = parallel::detectCores() - 1,
      FUN = function(k,
                     sample.size,
                     sim,
                     coal.events.times.sim.Exp)
      {
        quants <- c(0.025, 0.50, 0.975)
        matrixCurrentValue <-
          matrix(
            coal.events.times.sim.Exp[k, ],
            nrow = sim ,
            ncol = sample.size - 1,
            byrow = TRUE
          )
        #quantiles<-apply( matrixCurrentValue , 2 , quantile , probs = quants , na.rm = TRUE )
        quantiles <-
          apply(matrixCurrentValue , 2 , quantile , probs = quants)
        meanCoalTimes <- colMeans(matrixCurrentValue)
        result <-
          data.frame(
            mean = meanCoalTimes,
            LI = quantiles[1, ],
            median = quantiles[2, ],
            UI = quantiles[3, ]
          )
        result
      },
      sample.size,
      sim,
      coal.events.times.sim.Exp
    )
  
  return(list.Data.Frames)
}
################################################################################################################
#' simulateB_K
#' simulate coalescent times for model M_K for a single combination of Delta and Gamma
#' @param i: the position of this value of Delta in the DeltaList
#' @param Delta: the  paratemer Delta
#' @param Gamma: the  paratemer Gamma
#' @param sim: the  number of simulations sim
#' @param sample.size: the  sample size at time t=0
#' @param Time.Origin.STD: the  time of origin
#' @param coal.events.times.simB: a big matrix bigstatsr::FBM(length(DeltaList),sim*(sample.size-1))
#' of size to store the simulated coalescent times
#' @param K: the  parameter K of model M_K
#' @return: a data frame(sample.size -1 columns ) with the stats for the coalescent times for model M_K
################################################################################################################
simulateB_K = function(i,
                       Delta,
                       Gamma,
                       sim,
                       sample.size,
                       Time.Origin.STD,
                       coal.events.times.simB,
                       K)
{


  
  lapply(
    1:sim,
    FUN = function(j,
                   sample.size,
                   Time.Origin.STD,
                   Delta,
                   Gamma,
                   coal.events.times.simB,
                   i)
    {
      TimeOrigin = Time.Origin.STD[i, j]

      if (K>0){
                         
          coal.events.times = modelCoalB_K(sample.size, TimeOrigin, Delta, Gamma, K)
      }
      else{#K=0
          coal.events.times = modelCoal(sample.size, TimeOrigin, Delta, Gamma)
        }

     
     positions <-cbind(rep(i, (sample.size - 1)), ((j - 1) * (sample.size - 1) + 1):(j *
                  (sample.size - 1)))
      coal.events.times.simB[positions] <- rev(coal.events.times)
     
        print(paste("finished B M_k  for Delta",
                           Delta,
                           "Gamma",
                           Gamma,
                           "K",
                           K,
                           " sim ",
                           j,
                           sep = " "
                         )
                       )
    },
    sample.size,
    Time.Origin.STD,
    Delta,
    Gamma,
    coal.events.times.simB,
    i
  )
  
  
  quants <- c(0.025, 0.50, 0.975)
  matrixCurrentValue <-
          matrix(
            coal.events.times.simB[i, ],
            nrow = sim ,
            ncol = sample.size - 1,
            byrow = TRUE
          )
  quantiles <-
          apply(matrixCurrentValue , 2 , quantile , probs = quants)
  
  meanCoalTimes <-
          colMeans(matrixCurrentValue)
  
  result <-
          data.frame(
            mean = meanCoalTimes,
            LI = quantiles[1, ],
            median = quantiles[2, ],
            UI = quantiles[3, ]
          )
     
  return(result)
}
################################################################################################################
#' simulateB_K.parallel
#' simulate coalescent times for model M_K in parallel
#' @param DeltaList: the  paratemer Delta
#' @param GammaList: the  paratemer Gamma
#' @param sim: the  number of simulations sim
#' @param sample.size: the  sample size at time t=0
#' @param Time.Origin.STD: the  time of origin
#' @param coal.events.times.simB: a big matrix bigstatsr::FBM(length(DeltaList),sim*(sample.size-1))
#' of size to store the simulated coalescent times
#' @param K: the  parameter K>=0 of model M_K
#' @return: a list of data frame(sample.size -1 columns ) with the stats for the coalescent times
#' for model M_K(one data frame per each Delta value)
################################################################################################################
simulateB_K.parallel = function(DeltaList,
                                GammaList,
                                sim,
                                sample.size,
                                Time.Origin.STD,
                                coal.events.times.simB,
                                K)
{
  require(foreach)
  require(doParallel)
  require(doSNOW)
  require(parallel)
  require(doFuture)
  require(bigstatsr)
  require(stats)
  require(rgenoud)
  require(doRNG)
  require(doMC)
  
  RNGkind("L'Ecuyer-CMRG")
  a=as.numeric(Sys.time())
  set.seed(floor(runif(1)*a)) #set seed to something

  ncores = parallel::detectCores() - 1
  innerCluster <-
    parallel::makeCluster(ncores, type = "FORK", outfile = "")
  on.exit(parallel::stopCluster(innerCluster), add = TRUE)
  doParallel::registerDoParallel(innerCluster)
  
  
  rng <- RNGseq(length(DeltaList) * sim, floor(runif(1)*a))
  opts <- list(chunkSize = 2)
  
  tmp3 <- foreach::foreach(j = 1:sim, .combine = 'c') %:%
    foreach::foreach(i = 1:length(DeltaList),
                     r = rng[(j - 1) * length(DeltaList) + 1:length(DeltaList)],
                     .combine = 'c') %dopar% {
                       rngtools::setRNG(r)
                       
                       TimeOrigin = Time.Origin.STD[i, j]
                       
                       if (K>0){
                         
                         coal.events.times = modelCoalB_K(sample.size, TimeOrigin, DeltaList[i], GammaList[i], K)
                       }
                       else{#K=0
                         coal.events.times = modelCoal(sample.size, TimeOrigin, DeltaList[i], GammaList[i])
                       }
                       
                       print(
                         paste(
                           "finished B M_k in parallel for Delta",
                           DeltaList[i],
                           "Gamma",
                           GammaList[i],
                           "K",
                           K,
                           " sim ",
                           j,
                           sep = " "
                         )
                       )
                       positions <-
                         cbind(rep(i, (sample.size - 1)), ((j - 1) * (sample.size - 1) + 1):(j *
                                                                                               (sample.size - 1)))
                       coal.events.times.simB[positions] <- rev(coal.events.times)
                       NULL
                     }
  list.Data.FramesB <-
    parallel::mclapply(
      1:length(DeltaList),
      mc.set.seed = TRUE,
      mc.cores = parallel::detectCores() - 1,
      FUN = function(k,
                     sample.size,
                     sim,
                     coal.events.times.simB)
      {
        quants <- c(0.025, 0.50, 0.975)
        matrixCurrentValue <-
          matrix(
            coal.events.times.simB[k, ],
            nrow = sim ,
            ncol = sample.size - 1,
            byrow = TRUE
          )
        #quantiles<-apply( matrixCurrentValue , 2 , quantile , probs = quants , na.rm = TRUE )
        quantiles <-
          apply(matrixCurrentValue , 2 , quantile , probs = quants)
        meanCoalTimes <- colMeans(matrixCurrentValue)
        result <-
          data.frame(
            mean = meanCoalTimes,
            LI = quantiles[1, ],
            median = quantiles[2, ],
            UI = quantiles[3, ]
          )
        result
      },
      sample.size,
      sim,
      coal.events.times.simB
    )
  
  return(list.Data.FramesB)
}
################################################################################################################
#' simulateB_K_Race.parallel
#' simulate coalescent times for model M_K in parallel
#' @param DeltaList: the  paratemer Delta
#' @param GammaList: the  paratemer Gamma
#' @param sim: the  number of simulations sim
#' @param sample.size: the  sample size at time t=0
#' @param Time.Origin.STD: the  time of origin
#' @param coal.events.times.simB: a big matrix bigstatsr::FBM(length(DeltaList),sim*(sample.size-1))
#' of size to store the simulated coalescent times
#' @param K: the  parameter K>=0 of model M_K
#' @return: a list of data frame(sample.size -1 columns ) with the stats for the coalescent times
#' for model M_K(one data frame per each Delta value)
################################################################################################################
simulateB_K_Race.parallel = function(DeltaList,
                                GammaList,
                                sim,
                                sample.size,
                                Time.Origin.STD,
                                coal.events.times.simB.race,
                                K)
{
  require(foreach)
  require(doParallel)
  require(doSNOW)
  require(parallel)
  require(doFuture)
  require(bigstatsr)
  require(stats)
  require(rgenoud)
  require(doRNG)
  require(doMC)
  
  RNGkind("L'Ecuyer-CMRG")
  a=as.numeric(Sys.time())
  set.seed(floor(runif(1)*a)) #set seed to something
  
  ncores = parallel::detectCores() - 1
  innerCluster <-
    parallel::makeCluster(ncores, type = "FORK", outfile = "")
  on.exit(parallel::stopCluster(innerCluster), add = TRUE)
  doParallel::registerDoParallel(innerCluster)
  
  
  rng <- RNGseq(length(DeltaList) * sim, floor(runif(1)*a))
  opts <- list(chunkSize = 2)
  
  tmp3 <- foreach::foreach(j = 1:sim, .combine = 'c') %:%
    foreach::foreach(i = 1:length(DeltaList),
                     r = rng[(j - 1) * length(DeltaList) + 1:length(DeltaList)],
                     .combine = 'c') %dopar% {
                       rngtools::setRNG(r)
                       
                       TimeOrigin = Time.Origin.STD[i, j]
                       
                       if (K>0){
                         
                         coal.events.times = modelCoalB_K_race(sample.size, TimeOrigin, DeltaList[i], GammaList[i], K)
                       }
                       else{#K=0
                         coal.events.times = modelCoal_race(sample.size, TimeOrigin, DeltaList[i], GammaList[i])
                       }
                       
                       print(
                         paste(
                           "finished B M_k as race in parallel for Delta",
                           DeltaList[i],
                           "Gamma",
                           GammaList[i],
                           "K",
                           K,
                           " sim ",
                           j,
                           sep = " "
                         )
                       )
                       positions <-
                         cbind(rep(i, (sample.size - 1)), ((j - 1) * (sample.size - 1) + 1):(j *
                                                                                               (sample.size - 1)))
                       coal.events.times.simB.race[positions] <- coal.events.times
                       NULL
                     }
  list.Data.FramesB <-
    parallel::mclapply(
      1:length(DeltaList),
      mc.set.seed = TRUE,
      mc.cores = parallel::detectCores() - 1,
      FUN = function(k,
                     sample.size,
                     sim,
                     coal.events.times.simB.race)
      {
        quants <- c(0.025, 0.50, 0.975)
        matrixCurrentValue <-
          matrix(
            coal.events.times.simB.race[k, ],
            nrow = sim ,
            ncol = sample.size - 1,
            byrow = TRUE
          )
        #quantiles<-apply( matrixCurrentValue , 2 , quantile , probs = quants , na.rm = TRUE )
        quantiles <-
          apply(matrixCurrentValue , 2 , quantile , probs = quants)
        meanCoalTimes <- colMeans(matrixCurrentValue)
        result <-
          data.frame(
            mean = meanCoalTimes,
            LI = quantiles[1, ],
            median = quantiles[2, ],
            UI = quantiles[3, ]
          )
        result
      },
      sample.size,
      sim,
      coal.events.times.simB.race
    )
  
  return(list.Data.FramesB)
}
################################################################################################################
#' simulateB_MAster
#' simulate coalescent times for model M* for a combination of Gamma and Delta
#' @param i: the  position of Delta in the list of Delta values
#' @param Delta: the  paratemer Delta
#' @param Gamma: the  paratemer Gamma
#' @param sim: the  number of simulations sim
#' @param sample.size: the  sample size at time t=0
#' @param Time.Origin.STD: the  time of origin
#' @param coal.events.times.simB: a big matrix bigstatsr::FBM(length(DeltaList),sim*(sample.size-1))
#' of size to store the simulated coalescent times
#' @return: a  data frame(sample.size -1 columns ) with the stats for the coalescent times
#' for model M_K(one data frame per each Delta value)
################################################################################################################
simulateB_MAster = function(i,      Delta,
                                     Gamma,
                                     sim,
                                     sample.size,
                                     Time.Origin.STD,
                                     coal.events.times.simB)
{
 
  
  RNGkind("L'Ecuyer-CMRG")
   a=as.numeric(Sys.time())
  set.seed(floor(runif(1)*a))


  
  lapply(1:sim, FUN=function(j, Delta,
                                     Gamma, sample.size,
                                     sim, Time.Origin.STD, coal.events.times.simB, i

    ){
          TimeOrigin = Time.Origin.STD[i, j]
                       
          coal.events.times = modelCoalB_MAster(sample.size, TimeOrigin, Delta, Gamma)
                       
                print(
                   paste(
                           "finished B M Aster  for Delta",
                           Delta,
                           "Gamma",
                           Gamma,
                           " sim ",
                           j,
                           sep = " "
                         )
                       )


            positions <-
                         cbind(rep(i, (sample.size - 1)), ((j - 1) * (sample.size - 1) + 1):(j *
                                                                                               (sample.size - 1)))
                       #print(positions)
            coal.events.times.simB[positions] <- rev(coal.events.times)
          

  },  Delta,
      Gamma, sample.size,
  sim, Time.Origin.STD, coal.events.times.simB, i)
  
  

  quants <- c(0.025, 0.50, 0.975)
  matrixCurrentValue <-
          matrix(
            coal.events.times.simB[i, ],
            nrow = sim ,
            ncol = sample.size - 1,
            byrow = TRUE
          )
  quantiles <-
          apply(matrixCurrentValue , 2 , quantile , probs = quants)
  
  meanCoalTimes <-
          colMeans(matrixCurrentValue)
  
  result <-
          data.frame(
            mean = meanCoalTimes,
            LI = quantiles[1, ],
            median = quantiles[2, ],
            UI = quantiles[3, ]
          )
       
  
  return(result)
}
################################################################################################################
#' simulateB_MAster.parallel
#' simulate coalescent times for model M* in parallel
#' @param DeltaList: the  paratemer Delta
#' @param GammaList: the  paratemer Gamma
#' @param sim: the  number of simulations sim
#' @param sample.size: the  sample size at time t=0
#' @param Time.Origin.STD: the  time of origin
#' @param coal.events.times.simB: a big matrix bigstatsr::FBM(length(DeltaList),sim*(sample.size-1))
#' of size to store the simulated coalescent times
#' @return: a list of data frame(sample.size -1 columns ) with the stats for the coalescent times
#' for model M_K(one data frame per each Delta value)
################################################################################################################
simulateB_MAster.parallel = function(DeltaList,
                                     GammaList,
                                     sim,
                                     sample.size,
                                     Time.Origin.STD,
                                     coal.events.times.simB)
{
  require(foreach)
  require(doParallel)
  require(doSNOW)
  require(parallel)
  require(doFuture)
  require(bigstatsr)
  require(stats)
  require(rgenoud)
  require(doRNG)
  require(doMC)
  
  RNGkind("L'Ecuyer-CMRG")
   a=as.numeric(Sys.time())
  set.seed(floor(runif(1)*a))

  ncores = parallel::detectCores() - 1
  innerCluster <-
    parallel::makeCluster(ncores, type = "FORK", outfile = "")
  on.exit(parallel::stopCluster(innerCluster), add = TRUE)
  doParallel::registerDoParallel(innerCluster)
  
  
  rng <- RNGseq(length(DeltaList) * sim, floor(runif(1)*a))
  opts <- list(chunkSize = 2)
  
  tmp3 <- foreach::foreach(j = 1:sim, .combine = 'c') %:%
    foreach::foreach(i = 1:length(DeltaList),
                     r = rng[(j - 1) * length(DeltaList) + 1:length(DeltaList)],
                     .combine = 'c') %dopar% {
                       rngtools::setRNG(r)
                       
                       TimeOrigin = Time.Origin.STD[i, j]
                       
                       coal.events.times = modelCoalB_MAster(sample.size, TimeOrigin, DeltaList[i], GammaList[i])
                       
                       
                       print(
                         paste(
                           "finished B M Aster in parallel for Delta",
                           DeltaList[i],
                           "Gamma",
                           GammaList[i],
                           " sim ",
                           j,
                           sep = " "
                         )
                       )
                       positions <-
                         cbind(rep(i, (sample.size - 1)), ((j - 1) * (sample.size - 1) + 1):(j *
                                                                                               (sample.size - 1)))
                       #print(positions)
                       coal.events.times.simB[positions] <- rev(coal.events.times)
                       NULL
                     }
  
  
  
  list.Data.FramesB <-
    parallel::mclapply(
      1:length(DeltaList),
      mc.set.seed = TRUE,
      mc.cores = parallel::detectCores() - 1,
      FUN = function(k,
                     sample.size,
                     sim,
                     coal.events.times.simB)
      {
        quants <- c(0.025, 0.50, 0.975)
        matrixCurrentValue <-
          matrix(
            coal.events.times.simB[k, ],
            nrow = sim ,
            ncol = sample.size - 1,
            byrow = TRUE
          )
        quantiles <-
          apply(matrixCurrentValue , 2 , quantile , probs = quants)
        meanCoalTimes <-
          colMeans(matrixCurrentValue)
        result <-
          data.frame(
            mean = meanCoalTimes,
            LI = quantiles[1, ],
            median = quantiles[2, ],
            UI = quantiles[3, ]
          )
        result
      },
      sample.size,
      sim,
      coal.events.times.simB
    )
  
  return(list.Data.FramesB)
}
simulateC = function(Delta, Gamma,  sim, sample.size)
{
  require(parallel)
  coal.events.times.simC <- array(0, dim = c(sim, sample.size - 1))
  
  lapply(
    1:sim,
    FUN = function(i,
                   sample.size,
                   Delta,
                   Gamma,
                   coal.events.times.simC)
    {
      coal.events.times = modelCoalC(sample.size, Delta, Gamma)
      coal.events.times.simC[i, ] <<- rev(coal.events.times)
      print(paste0("finished C sim", i, sep = " "))
    },
    sample.size,
    Delta,
    Gamma,
    coal.events.times.simC
  )
  
  quants <- c(0.025, 0.50, 0.975)
  #quantiles<-apply( coal.events.times.simC , 2 , quantile , probs = quants , na.rm = TRUE )
  quantiles <-
    apply(coal.events.times.simC , 2 , quantile , probs = quants)
  meanCoalTimes <- colMeans(coal.events.times.simC)
  result <-
    data.frame(
      mean = meanCoalTimes,
      LI = quantiles[1, ],
      median = quantiles[2, ],
      UI = quantiles[3, ]
    )
  return(result)
}
################################################################################################################
#' simulateC.parallel
#' simulate coalescent times for model C in parallel
#' @param DeltaList: the  paratemer Delta
#' @param GammaList: the  paratemer Gamma
#' @param sim: the  number of simulations sim
#' @param sample.size: the  sample size at time t=0
#' @param Time.Origin.STD: the  time of origin
#' @param coal.events.times.simC: a big matrix bigstatsr::FBM(length(DeltaList),sim*(sample.size-1))
#' of size to store the simulated coalescent times
#' @return: a list of data frame(sample.size -1 columns ) with the stats for the coalescent times
#' for model M_K(one data frame per each Delta value)
################################################################################################################
simulateC.parallel = function(DeltaList,
                              GammaList,
                              sim,
                              sample.size,
                              coal.events.times.simC)
{
  require(foreach)
  require(doParallel)
  require(doSNOW)
  require(parallel)
  require(doFuture)
  require(bigstatsr)
  require(stats)
  require(rgenoud)
  require(doRNG)
  require(doMC)
  require(pracma)
  require(NLRoot)
  require(cmna)
  RNGkind("L'Ecuyer-CMRG")
  a=as.numeric(Sys.time())
  set.seed(floor(runif(1)*a))
  ncores = parallel::detectCores() - 1
  innerCluster <-
    parallel::makeCluster(ncores, type = "FORK", outfile = "")
  on.exit(parallel::stopCluster(innerCluster), add = TRUE)
  doParallel::registerDoParallel(innerCluster)
  rng <- RNGseq(length(DeltaList) * sim, floor(runif(1)*a))
  
  
  opts <- list(chunkSize = 2)
  
  tmp3 <- foreach::foreach(j = 1:sim, .combine = 'c') %:%
    foreach::foreach(i = 1:length(DeltaList),
                     r = rng[(j - 1) * length(DeltaList) + 1:length(DeltaList)],
                     .combine = 'c') %dopar% {
                       rngtools::setRNG(r)
                       
                       coal.events.times = modelCoalC(sample.size, DeltaList[i], GammaList[i])
                       
                       positions <-
                         cbind(rep(i, (sample.size - 1)), ((j - 1) * (sample.size - 1) + 1):(j *
                                                                                               (sample.size - 1)))
                       
                       coal.events.times.simC[positions] <- rev(coal.events.times)
                       print(
                         paste(
                           "finished C in parallel for Delta",
                           DeltaList[i],
                           "Gamma",
                           GammaList[i],
                           " sim ",
                           j,
                           sep = " "
                         )
                       )
                       NULL
                     }
  
  list.Data.FramesC <-
    parallel::mclapply(
      1:length(DeltaList),
      mc.cores = parallel::detectCores() - 1,
      FUN = function(k,
                     sample.size,
                     sim,
                     coal.events.times.simC) {
        quants <- c(0.025, 0.50, 0.975)
        matrixCurrentValue <-
          matrix(
            coal.events.times.simC[k, ],
            nrow = sim ,
            ncol = sample.size - 1,
            byrow = TRUE
          )
        #quantiles<-apply( matrixCurrentValue , 2 , quantile , probs = quants , na.rm = TRUE )
        quantiles <-
          apply(matrixCurrentValue , 2 , quantile , probs = quants)
        meanCoalTimes <- colMeans(matrixCurrentValue)
        result <-
          data.frame(
            mean = meanCoalTimes,
            LI = quantiles[1, ],
            median = quantiles[2, ],
            UI = quantiles[3, ]
          )
        result
      },
      sample.size,
      sim,
      coal.events.times.simC
    )
  
  
  return(list.Data.FramesC)
  
}
################################################################################################################
#' simulateA
#' simulate coalescent times for model BD for a given Delta and Gamma value
#' @param Delta: the  parameter Delta
#' @param Gamma: the  parameter Gamma
#' @param sim: the  number of simulations sim
#' @param sample.size: the  sample size at time t=0
#' @param Time.Origin.STD: the  time of origin
#' @param i: index of the Delta and Gamma value in the GammaList and Deltalist
#' @param coal.events.times.simA: a big matrix bigstatsr::FBM(length(DeltaList),sim*(sample.size-1))
#' of size to store the simulated coalescent times
#' @param number.ancestors.simA: a big matrix bigstatsr::FBM(length(DeltaList),sim*(sample.size-1))
#' of size to store the simulated number of ancestors in the population
#' @param number.ancestors.Transition: a big matrix bigstatsr::FBM(length(DeltaList),sim)
#' of size to store the simulated number of ancestors in the population when sample size is floor(0.5*sample.size)
#' @return: a list of data frame(sample.size -1 columns ) with the stats for the coalescent times
#' for model M_K(one data frame per each Delta value)
#' @return: a  data frame(sample.size -1 columns ) with the stats for the coalescent times
################################################################################################################
simulateA = function(Delta,
                     Gamma,
                     sim,
                     sample.size,
                     Time.Origin.STD,
                     i,
                     coal.events.times.simA,
                     number.ancestors.simA,
                     number.ancestors.Transition)
{
  
   RNGkind("L'Ecuyer-CMRG")
   a=as.numeric(Sys.time())
   set.seed(floor(runif(1)*a))


  lapply(1:sim, FUN=function(j, Delta,
                     Gamma,
                     sim,
                     sample.size,
                     Time.Origin.STD,
                     coal.events.times.simA,
                     number.ancestors.simA,  number.ancestors.Transition, i )
                 {
               
                
                       list.number.ancestors.population = simulate.list.number.ancestors.population(sample.size)
                       
                       
                       list.coal.times = simulate.coalescent.times.A(Gamma,
                                                                     Delta,
                                                                     sample.size,
                                                                     list.number.ancestors.population)
                       
                       coal.events.times = unlist(list.coal.times)
                       
                       Time.Origin = coal.events.times[length(coal.events.times)]
                       
                       Time.Origin.STD[i, j] <- Time.Origin
                       
                       positions <-
                         cbind(rep(i, (sample.size - 1)), ((j - 1) * (sample.size - 1) + 1):(j *
                                                                                               (sample.size - 1)))
                       
                       number.ancestors.simA[positions] <-
                         list.number.ancestors.population[2:(length(coal.events.times))]
                       
                       coal.events.times.simA[positions] <-
                         coal.events.times[1:(length(coal.events.times) - 1)]
                       m = list.number.ancestors.population[floor(length(list.number.ancestors.population) /
                                                                    2)]
                       number.ancestors.Transition[cbind(i, j)] <- m
                       
                       
                       print(
                         paste(
                           "finished A  for Delta",
                           Delta,
                           "Gamma",
                           Gamma,
                           " sim ",
                           j,
                           sep = " "
                         )
                       )
                    
            }, 
            Delta,
            Gamma,
            sim,
            sample.size,
            Time.Origin.STD,
            coal.events.times.simA,
            number.ancestors.simA,  number.ancestors.Transition, i )


  
  quants <- c(0.025, 0.50, 0.975)
  matrixCurrentValue <-
          matrix(
            coal.events.times.simB[i, ],
            nrow = sim ,
            ncol = sample.size - 1,
            byrow = TRUE
          )
  quantiles <-
          apply(matrixCurrentValue , 2 , quantile , probs = quants)
  
  meanCoalTimes <-
          colMeans(matrixCurrentValue)
  
  result <-
          data.frame(
            mean = meanCoalTimes,
            LI = quantiles[1, ],
            median = quantiles[2, ],
            UI = quantiles[3, ]
          )
     
  return(result)
}
################################################################################################################
#' simulateA.parallel
#' simulate coalescent times for model BD in parallel
#' @param DeltaList: the list of  Delta parameters
#' @param GammaList: the list of   Gamma parameters
#' @param sim: the  number of simulations sim
#' @param sample.size: the  sample size at time t=0
#' @param Time.Origin.STD: the  time of origin
#' @param coal.events.times.simA: a big matrix bigstatsr::FBM(length(DeltaList),sim*(sample.size-1))
#' of size to store the simulated coalescent times
#' @param number.ancestors.simA: a big matrix bigstatsr::FBM(length(DeltaList),sim*(sample.size-1))
#' of size to store the simulated number of ancestors in the population
#' @param number.ancestors.Transition: a big matrix bigstatsr::FBM(length(DeltaList),sim)
#' of size to store the simulated number of ancestors in the population when sample size is floor(0.5*sample.size)
#' @return: a list of data frame(sample.size -1 columns ) with the stats for the coalescent times
#' for model M_K(one data frame per each Delta value)
################################################################################################################
simulateA.parallel = function(DeltaList,
                              GammaList,
                              sim,
                              sample.size,
                              Time.Origin.STD,
                              coal.events.times.simA,
                              number.ancestors.simA,
                              number.ancestors.Transition)
{
  require(foreach)
  require(doParallel)
  require(doSNOW)
  require(parallel)
  require(doFuture)
  require(bigstatsr)
  require(stats)
  require(rgenoud)
  require(doRNG)
  require(doMC)
  
  RNGkind("L'Ecuyer-CMRG")
   a=as.numeric(Sys.time())
  set.seed(floor(runif(1)*a))
  ncores = parallel::detectCores() - 1
  innerCluster <-
    parallel::makeCluster(ncores, type = "FORK", outfile = "")
  on.exit(parallel::stopCluster(innerCluster))
  doParallel::registerDoParallel(innerCluster)
  #rng<- RNGseq( length(DeltaList)* sim, 6452913)
  
  rng <- RNGseq(length(DeltaList) * sim, floor(runif(1)*a))
  
  opts <- list(chunkSize = 2)
  
  tmp3 <- foreach::foreach(j = 1:sim, .combine = 'c') %:%
    foreach::foreach(i = 1:length(DeltaList),
                     r = rng[(j - 1) * length(DeltaList) + 1:length(DeltaList)],
                     .combine = 'c') %dopar% {
                       rngtools::setRNG(r)
                       
                       list.number.ancestors.population = simulate.list.number.ancestors.population(sample.size)
                       
                       
                       list.coal.times = simulate.coalescent.times.A(GammaList[i],
                                                                     DeltaList[i],
                                                                     sample.size,
                                                                     list.number.ancestors.population)
                       
                       coal.events.times = unlist(list.coal.times)
                       
                       Time.Origin = coal.events.times[length(coal.events.times)]
                       
                       Time.Origin.STD[i, j] <- Time.Origin
                       
                       positions <-
                         cbind(rep(i, (sample.size - 1)), ((j - 1) * (sample.size - 1) + 1):(j *
                                                                                               (sample.size - 1)))
                       
                       number.ancestors.simA[positions] <-
                         list.number.ancestors.population[2:(length(coal.events.times))]
                       
                       coal.events.times.simA[positions] <-
                         coal.events.times[1:(length(coal.events.times) - 1)]
                       m = list.number.ancestors.population[floor(length(list.number.ancestors.population) /
                                                                    2)]
                       number.ancestors.Transition[cbind(i, j)] <- m
                       
                       
                       print(
                         paste(
                           "finished A in parallel for Delta",
                           DeltaList[i],
                           "Gamma",
                           GammaList[i],
                           " sim ",
                           j,
                           sep = " "
                         )
                       )
                       NULL
                     }
  
  list.Data.FramesA <-
    parallel::mclapply(
      1:length(DeltaList),
      mc.cores = parallel::detectCores() - 1,
      FUN = function(k,
                     sample.size,
                     sim,
                     coal.events.times.simA) {
        quants <- c(0.025, 0.50, 0.975)
        
        matrixCurrentValue <-
          matrix(
            coal.events.times.simA[k, ],
            nrow = sim ,
            ncol = sample.size - 1,
            byrow = TRUE
          )
        #quantiles<-apply( matrixCurrentValue , 2 , quantile , probs = quants , na.rm = TRUE )
        quantiles <-
          apply(matrixCurrentValue , 2 , quantile , probs = quants)
        meanCoalTimes <- colMeans(matrixCurrentValue)
        result <-
          data.frame(
            mean = meanCoalTimes,
            LI = quantiles[1, ],
            median = quantiles[2, ],
            UI = quantiles[3, ]
          )
        result
      },
      sample.size,
      sim,
      coal.events.times.simA
    )
  
  
  return(list.Data.FramesA)
  
}
################################################################################################################
#' simulateHybrid
#' simulate coalescent times for  Hybrid model 
#' @param i: the  position of Delta in Delta list
#' @param Delta: the  parameter Delta
#' @param Gamma: the  parameter Gamma
#' @param sim: the  number of simulations sim
#' @param sample.size: the  sample size at time t=0
#' @param Time.Origin.STD: the  time of origin
#' @param coal.events.times.simHybrid: a big matrix bigstatsr::FBM(length(DeltaList),sim*(sample.size-1))
#' of size to store the simulated coalescent times
#' @param nB: the sample size when the Hybrid model pass from M_k model  to BD model
#' @param number.Fail.hybrid: a big matrix bigstatsr::FBM(length(mprimeList),length(GammaList)) where mprimeList
#' is the list of nB values,
#' to store the number of times we have to resample the number of ancestors in the population m_nB when
#' sample is nB because m_nB < nB
#' @param posMprime position of nB in the list of mprimeList
#' @param K parameter K of the model M_K(applied backwards in time from sample.size  to nB)
#' @param number.ancestors.simHybrid: a big matrix bigstatsr::FBM(length(DeltaList),sim*(sample.size-1))
#' of size to store the simulated number of ancestors in the population
#' @param number.ancestors.Transition : a big matrix bigstatsr::FBM(length(DeltaList),sim)
#' of size to store the simulated number of ancestors in the population when the Hybrid model pass from M_k model
#' to BD model
#' @return: a  data frame(sample.size -1 columns ) with the stats for the coalescent times
#' for model M_K(one data frame per each Delta value)
################################################################################################################
simulateHybrid = function(i, Delta,
                                   Gamma,
                                   sim,
                                   sample.size,
                                   Time.Origin.STD,
                                   coal.events.times.simHybrid,
                                   nB,
                                   number.Fail.hybrid,
                                   posMprime,
                                   K,
                                   number.ancestors.simHybrid,
                                   number.ancestors.Transition)
{

  
  RNGkind("L'Ecuyer-CMRG")
   a=as.numeric(Sys.time())
  set.seed(floor(runif(1)*a))
 
  
  
  tmp3 <- lapply(1:sim, FUN=function(j, sample.size,
    Time.Origin.STD,
    Delta,
    Gamma,
    coal.events.times.simHybrid,
    nB,
    number.Fail.hybrid,
    posMprime,
    K,
    number.ancestors.simHybrid,
    number.ancestors.Transition,
    i
    ){

                     TimeOrigin <- Time.Origin.STD[i, j]
                       
                       m = -1
                       while (m < nB) {
                         if (m != -1) {
                           number.Fail.hybrid[cbind(posMprime, i)] = number.Fail.hybrid[cbind(posMprime, i)] +
                             1
                         }
                         
                         std.coal.events.times.until.B = rev(modelCoalHybrid(sample.size, TimeOrigin, Delta, Gamma, nB))
                         
                         coal.events.times.until.nB = unlist(lapply(
                           std.coal.events.times.until.B,
                           FUN = function(x)
                             standard2modelB_K(x, TimeOrigin, Delta, Gamma, K)
                         ))
                         
                         t = std.coal.events.times.until.B[length(std.coal.events.times.until.B)]
                         
                         m = ceiling(2.0 / t)
                         
                         
                       }
                       number.ancestors.Transition[cbind(i, j)] <- m
                       
                       
                       list.number.ancestors.population = simulate.list.number.ancestors.population.from(sample.size, m, nB)
                       
                       s <-
                         standard2modelB_K(t, TimeOrigin, Delta, Gamma, K)
                       
                       
                       list.coal.times = simulate.coalescent.times.A.from(Gamma,
                                                                          Delta,
                                                                          sample.size,
                                                                          list.number.ancestors.population,
                                                                          nB,
                                                                          s)
                       
                       coal.events.times = c(unlist(coal.events.times.until.nB),
                                             unlist(list.coal.times))
                       
                    

                      positions <-
                         cbind(rep(i, (sample.size - 1)), ((j - 1) * (sample.size - 1) + 1):(j *
                                                                                               (sample.size - 1)))
                       
                       coal.events.times.simHybrid[positions] <-
                         unlist(coal.events.times)
                       
                       number.ancestors.simHybrid[positions] <-
                         c(rep(0, sample.size - nB - 2),
                           unlist(list.number.ancestors.population),
                           0)
                       
                       
                       print(
                         paste(
                           "finished hybrid scenario for Delta",
                           Delta,
                           "Gamma",
                           Gamma,
                           " k0 ",
                           nB,
                           " sim ",
                           j,

                           sep = " "
                         )
                       )
                       NULL




  },
  sample.size,
    Time.Origin.STD,
    Delta,
    Gamma,
    coal.events.times.simHybrid,
    nB,
    number.Fail.hybrid,
    posMprime,
    K,
    number.ancestors.simHybrid,
    number.ancestors.Transition,
    i)
    

   
  

  quants <- c(0.025, 0.50, 0.975)
  matrixCurrentValue <-
          matrix(
            coal.events.times.simB[i, ],
            nrow = sim ,
            ncol = sample.size - 1,
            byrow = TRUE
          )
  quantiles <-
          apply(matrixCurrentValue , 2 , quantile , probs = quants)
  
  meanCoalTimes <-
          colMeans(matrixCurrentValue)
  
  result <-
          data.frame(
            mean = meanCoalTimes,
            LI = quantiles[1, ],
            median = quantiles[2, ],
            UI = quantiles[3, ]
          )
  
  return(result)
  
}
################################################################################################################
#' simulateHybrid.parallel
#' simulate coalescent times for  Hybrid model in parallel
#' @param DeltaList: the  parameter Delta
#' @param GammaList: the  parameter Gamma
#' @param sim: the  number of simulations sim
#' @param sample.size: the  sample size at time t=0
#' @param Time.Origin.STD: the  time of origin
#' @param coal.events.times.simHybrid: a big matrix bigstatsr::FBM(length(DeltaList),sim*(sample.size-1))
#' of size to store the simulated coalescent times
#' @param nB: the sample size when the Hybrid model pass from M_k model  to BD model
#' @param number.Fail.hybrid: a big matrix bigstatsr::FBM(length(mprimeList),length(GammaList)) where mprimeList
#' is the list of nB values,
#' to store the number of times we have to resample the number of ancestors in the population m_nB when
#' sample is nB because m_nB < nB
#' @param posMprime position of nB in the list of mprimeList
#' @param K parameter K of the model M_K(applied backwards in time from sample.size  to nB)
#' @param number.ancestors.simHybrid: a big matrix bigstatsr::FBM(length(DeltaList),sim*(sample.size-1))
#' of size to store the simulated number of ancestors in the population
#' @param number.ancestors.Transition : a big matrix bigstatsr::FBM(length(DeltaList),sim)
#' of size to store the simulated number of ancestors in the population when the Hybrid model pass from M_k model
#' to BD model
#' @return: a list of data frame(sample.size -1 columns ) with the stats for the coalescent times
#' for model M_K(one data frame per each Delta value)
################################################################################################################
simulateHybrid.parallel = function(DeltaList,
                                   GammaList,
                                   sim,
                                   sample.size,
                                   Time.Origin.STD,
                                   coal.events.times.simHybrid,
                                   nB,
                                   number.Fail.hybrid,
                                   posMprime,
                                   K,
                                   number.ancestors.simHybrid,
                                   number.ancestors.Transition)
{
  require(foreach)
  require(doParallel)
  require(doSNOW)
  require(parallel)
  require(doFuture)
  require(bigstatsr)
  require(stats)
  require(rgenoud)
  require(doRNG)
  require(doMC)
  
  RNGkind("L'Ecuyer-CMRG")
   a=as.numeric(Sys.time())
  set.seed(floor(runif(1)*a))
  ncores = parallel::detectCores() - 1
  innerCluster2 <-
    parallel::makeCluster(ncores, type = "FORK", outfile = "")
  
  #on.exit(parallel::stopCluster(innerCluster), add = TRUE)
  on.exit(parallel::stopCluster(innerCluster2))
  doParallel::registerDoParallel(innerCluster2)
  rng <- RNGseq(length(DeltaList) * sim, floor(runif(1)*a))
  
  opts <- list(chunkSize = 2)
  
  
  tmp3 <- foreach::foreach(j = 1:sim, .combine = 'c') %:%
    foreach::foreach(i = 1:length(DeltaList),
                     r = rng[(j - 1) * length(DeltaList) + 1:length(DeltaList)],
                     .combine = 'c') %dopar% {
                       rngtools::setRNG(r)
                       
                       TimeOrigin <- Time.Origin.STD[i, j]
                       
                       m = -1
                       while (m < nB) {
                         if (m != -1) {
                           number.Fail.hybrid[cbind(posMprime, i)] = number.Fail.hybrid[cbind(posMprime, i)] +
                             1
                         }
                         
                         std.coal.events.times.until.B = rev(modelCoalHybrid(sample.size, TimeOrigin, DeltaList[i], GammaList[i], nB))
                         
                         coal.events.times.until.nB = unlist(lapply(
                           std.coal.events.times.until.B,
                           FUN = function(x)
                             standard2modelB_K(x, TimeOrigin, DeltaList[i], GammaList[i], K)
                         ))
                         
                         t = std.coal.events.times.until.B[length(std.coal.events.times.until.B)]
                         
                         m = ceiling(2.0 / t)
                         
                         
                       }
                       number.ancestors.Transition[cbind(i, j)] <- m
                       
                       
                       list.number.ancestors.population = simulate.list.number.ancestors.population.from(sample.size, m, nB)
                       
                       s <-
                         standard2modelB_K(t, TimeOrigin, DeltaList[i], GammaList[i], K)
                       
                       
                       list.coal.times = simulate.coalescent.times.A.from(GammaList[i],
                                                                          DeltaList[i],
                                                                          sample.size,
                                                                          list.number.ancestors.population,
                                                                          nB,
                                                                          s)
                       
                       coal.events.times = c(unlist(coal.events.times.until.nB),
                                             unlist(list.coal.times))
                       
                       
                       positions <-
                         cbind(rep(i, (sample.size - 1)), ((j - 1) * (sample.size - 1) + 1):(j *
                                                                                               (sample.size - 1)))
                       
                       coal.events.times.simHybrid[positions] <-
                         unlist(coal.events.times)
                       
                       number.ancestors.simHybrid[positions] <-
                         c(rep(0, sample.size - nB - 2),
                           unlist(list.number.ancestors.population),
                           0)
                       
                       
                       print(
                         paste(
                           "finished hybrid scenario in parallel for Delta",
                           DeltaList[i],
                           "Gamma",
                           GammaList[i],
                           " sim ",
                           j,
                           sep = " "
                         )
                       )
                       NULL
                     }
  
  list.Data.FramesHybrid <-
    parallel::mclapply(
      1:length(DeltaList),
      mc.cores = parallel::detectCores() - 1,
      FUN = function(k,
                     sample.size,
                     sim,
                     coal.events.times.simHybrid) {
        quants <- c(0.025, 0.50, 0.975)
        
        matrixCurrentValue <-
          matrix(
            coal.events.times.simHybrid[k, ],
            nrow = sim ,
            ncol = sample.size - 1,
            byrow = TRUE
          )
        #quantiles<-apply( matrixCurrentValue , 2 , quantile , probs = quants , na.rm = TRUE )
        quantiles <-
          apply(matrixCurrentValue , 2 , quantile , probs = quants)
        meanCoalTimes <- colMeans(matrixCurrentValue)
        result <-
          data.frame(
            mean = meanCoalTimes,
            LI = quantiles[1, ],
            median = quantiles[2, ],
            UI = quantiles[3, ]
          )
        result
      },
      sample.size,
      sim,
      coal.events.times.simHybrid
    )
  
  
  return(list.Data.FramesHybrid)
  
}
simulateA.test = function(Delta,
                          Gamma,
                          sim,
                          sample.size,
                          Time.Origin.STD,
                          i,
                          coal.events.times.simA)
{
  list.number.ancestors.population.sim <<-
    array(0, dim = c(sim, sample.size))
  Time.Origin.STD.sim <- array(0, dim = c(sim))
  lapply(
    1:sim,
    FUN = function(j,
                   sample.size,
                   Gamma,
                   Delta,
                   coal.events.times.simA,
                   Time.Origin.STD.sim)
    {
      list.number.ancestors.population = simulate.list.number.ancestors.population(sample.size)
      
      list.number.ancestors.population[list.number.ancestors.population ==
                                         0] <- 1
      list.number.ancestors.population.sim[j, ] <<-
        list.number.ancestors.population
      list.coal.times = simulate.coalescent.times.A.test(Gamma,
                                                         Delta,
                                                         sample.size,
                                                         list.number.ancestors.population)
      
      coal.events.times = unlist(list.coal.times)
      
      coal.events.times.simA[j, ] <<- coal.events.times
      print(paste0("finished A test sim ", j, sep = " "))
      
    },
    sample.size,
    Gamma,
    Delta,
    coal.events.times.simA,
    Time.Origin.STD.sim = Time.Origin.STD.sim
  )
  
  
  quants <- c(0.025, 0.50, 0.975)
  quantiles <-
    apply(coal.events.times.simA , 2 , quantile , probs = quants)
  meanCoalTimes <- colSums(coal.events.times.simA) / sim
  result <-
    data.frame(
      mean = meanCoalTimes,
      LI = quantiles[1, ],
      median = quantiles[2, ],
      UI = quantiles[3, ]
    )
  return(result)
}
simulateA.test.parallel = function(DeltaList,
                                   GammaList,
                                   sim,
                                   sample.size,
                                   Time.Origin.STD.test,
                                   coal.events.times.sim.testA)
{
  require(foreach)
  require(doParallel)
  require(doSNOW)
  require(parallel)
  require(doFuture)
  require(bigstatsr)
  require(stats)
  require(rgenoud)
  require(doMC)
  
  RNGkind("L'Ecuyer-CMRG")
   a=as.numeric(Sys.time())
  set.seed(floor(runif(1)*a))
 
  
  clust <-
    parallel::makeCluster(parallel::detectCores() - 1,
                          type = "FORK",
                          outfile = "")
  
  on.exit(parallel::stopCluster(clust), add = TRUE)
  doParallel::registerDoParallel(clust)
  rng <- RNGseq(length(DeltaList) * sim, floor(runif(1)*a))
  
  opts <- list(chunkSize = 2)
  
  
  results <- foreach::foreach(j = 1:sim, .combine = 'c') %:%
    foreach::foreach(i = 1:length(DeltaList),
                     r = rng[(j - 1) * length(DeltaList) + 1:length(DeltaList)],
                     .combine = 'c') %dopar% {
                       rngtools::setRNG(r)
                       
                       list.number.ancestors.population = simulate.list.number.ancestors.population(sample.size)
                       list.coal.times = simulate.coalescent.times.A.test(GammaList[i],
                                                                          DeltaList[i],
                                                                          sample.size,
                                                                          list.number.ancestors.population)
                       coal.events.times = unlist(list.coal.times)
                       Time.Origin = coal.events.times[length(coal.events.times)]
                       Time.Origin.STD.test[i, j] <- Time.Origin
                       positions <-
                         cbind(rep(i, (sample.size - 1)), ((j - 1) * (sample.size - 1) + 1):(j *
                                                                                               (sample.size - 1)))
                       coal.events.times.sim.testA[positions] <-
                         coal.events.times[1:(length(coal.events.times) - 1)]
                       print(
                         paste(
                           "finished A test for Delta",
                           DeltaList[i],
                           "Gamma",
                           GammaList[i],
                           " sim ",
                           j,
                           sep = " "
                         )
                       )
                       NULL
                     }
  
  list.Data.FramesA <-
    parallel::mclapply(
      1:length(DeltaList),
      mc.set.seed = TRUE,
      mc.cores = parallel::detectCores() - 1,
      FUN = function(k,
                     sample.size,
                     sim,
                     coal.events.times.sim.testA) {
        quants <- c(0.025, 0.50, 0.975)
        matrixCurrentValue <-
          matrix(
            coal.events.times.sim.testA[k, ],
            nrow = sim ,
            ncol = sample.size - 1,
            byrow = TRUE
          )
        quantiles <-
          apply(matrixCurrentValue , 2 , quantile , probs = quants)
        meanCoalTimes <- colMeans(matrixCurrentValue)
        result <-
          data.frame(
            mean = meanCoalTimes,
            LI = quantiles[1, ],
            median = quantiles[2, ],
            UI = quantiles[3, ]
          )
        result
      },
      sample.size,
      sim,
      coal.events.times.sim.testA
    )
  
  return(list.Data.FramesA)
}
################################################################################################################
#C scenario functions
Integrand <- function(u, Delta, Gamma) {
  integrand <- function(x) {
    term <- exp(-1.0 * Delta * x)
    
    relative(x, Delta, u) * Gamma * Delta * term * exp(-1.0 * Gamma * term /
                                                         (1 - term)) / (1 - term) ^ 2
  }
  return(integrand)
}
h_integral <- function(u, Delta, Gamma) {
  integrand <- function(x, u, Delta, Gamma)
  {
    term <- exp(-1.0 * Delta * x)
    
    return(relative(x, Delta, u) * Gamma * Delta * term * exp(-1.0 * Gamma * term /
                                                                (1 - term)) / (1 - term) ^ 2)
  }
  return(
    integrate(
      integrand,
      u = u,
      Delta = Delta,
      Gamma = Gamma,
      lower = u,
      upper = Inf,
      stop.on.error = FALSE
    )$value
  )
}
vh_integral <- function(u, Delta, Gamma) {
  integrand <- function(x, u, Delta, Gamma)
  {
    term <- exp(-1.0 * Delta * x)
    
    return(relative(x, Delta, u) * Gamma * Delta * term * exp(-1.0 * Gamma * term /
                                                                (1 - term)) / (1 - term) ^ 2)
  }
  return(
    integrate(
      integrand,
      u = u,
      Delta = Delta,
      Gamma = Gamma,
      lower = u,
      upper = Inf,
      stop.on.error = FALSE
    )$value
  )
}
hFinal_function <- function(Delta, Gamma) {
  res <- function(x)
  {
    return(h_integral(x, Delta, Gamma))
  }
  return(res)
}
h_analytic <- function(Delta, Gamma) {
  res <- function(x) {
    exp_minus_delta_u <- exp(-1.0 * Delta * x)
    
    v <- exp_minus_delta_u / (1 - exp_minus_delta_u)
    exp_minus_gamma_v <- exp(-1.0 * Gamma * v)
    
    result <-
      -1.0 * (1.0 / exp_minus_delta_u) * ((1.0 - exp_minus_delta_u) ^ 2) * (exp_minus_gamma_v *
                                                                              (v ^ 2 + 2.0 * v / Gamma + 2.0 / (Gamma ^ 2)) - 2.0 / (Gamma ^ 2))
    result <-
      result + 2.0 * (1 - exp_minus_delta_u) * (exp_minus_gamma_v * (v + 1.0 /
                                                                       Gamma) - 1.0 / Gamma)
    result <- result + exp_minus_delta_u * (1 - exp_minus_gamma_v)
    return(result)
  }
  return(res)
}
lambdaFinal <- function(u, Delta, Gamma) {
  factor = 2.0 * Delta / Gamma
  return(factor * (1.0 / h_integral(u, Delta, Gamma)))
}
lambda_analytic <- function(u, Delta, Gamma) {
  factor = 2.0 * Delta / Gamma
  funcX <- h_analytic(Delta, Gamma)
  return(factor * (1.0 / funcX(u)))
}
LambdaFinal <- function(s, Delta, Gamma) {
  if (s > 0) {
    return(
      integrate(
        Vectorize(lambda_analytic, c("u")),
        Delta = Delta,
        Gamma = Gamma,
        lower = 0,
        upper = s,
        stop.on.error = FALSE
      )$value
    )
    
  }
  else{
    return(0)
  }
}
TransformedLambda <- function(z, Delta, Gamma) {
  if (Gamma < 100)
    return(
      integrate(
        Vectorize(lambdaFinal, c("u")),
        Delta = Delta,
        Gamma = Gamma,
        lower = 0,
        upper = exp(z),
        stop.on.error = FALSE
      )$value
    )
  else{
    return(
      integrate(
        Vectorize(lambda_analytic, c("u")),
        Delta = Delta,
        Gamma = Gamma,
        lower = 0,
        upper =  exp(z),
        stop.on.error = FALSE
      )$value
    )
  }
  
}
Lambda_integral <- function(s, Delta, Gamma) {
  f <- function(x)
  {
    return(
      integrate(
        Vectorize(lambdaFinal, c("u")),
        Delta = Delta,
        Gamma = Gamma,
        lower = 0,
        upper = s,
        stop.on.error = FALSE
      )$value
    )
  }
  return(f)
}
log_functionToZero <- function(y, Delta, Gamma)
{
  f <- function(x) {
    return(log(LambdaFinal(x, Delta, Gamma)) - log(y))
  }
  return(f)
}
derivativeLogLambda <- function(Delta, Gamma)
{
  f <- function(x)
  {
    return(lambdaFinal(x, Delta, Gamma) / LambdaFinal(x, Delta, Gamma))
  }
  return(f)
}
functionToZero <- function(y, Delta, Gamma)
{
  f <- function(x)
    LambdaFinal(x, Delta, Gamma) - y
  return(f)
}
TransformedfunctionToZero <- function(y, Delta, Gamma)
{
  f <- function(x) {
    w <- exp(x)
    return(LambdaFinal(exp(x), Delta, Gamma) - y)
  }
  return(f)
}
derivativeTransformedLambda <- function(Delta, Gamma)
{
  f <- function(x)
  {
    return(lambda_analytic(exp(x), Delta, Gamma) * exp(x))
  }
  return(f)
}
################################################################################################################
#' InverseTransformedLambda
#' Computes Lambda^(-1)(s) for scenario scenario C by transforming the
#' variables Lambda^(-1)(exp(s)) to allow the argument of Lambda^(-1)
#' to be non negative
#' @param s
#' @param Delta Delta parameter
#' @param Gamma Gamma parameter
################################################################################################################
InverseTransformedLambda <- function(s, Delta, Gamma)
{
  upper <- search_infinite_from(0, 0, functionToZero(s, Delta, Gamma))
  result <-
    cmna::bisection(functionToZero(s, Delta, Gamma),
                    min(0, upper - 1),
                    upper,
                    tol = 10 ^ (-6))
  return(result)
}
################################################################################################################
#' standard2modelC
#'
#' Kingsnman coalescent  time to model time C
#' @param u The time variable u in model time
#' @param Delta The parameter Delta
#' @param Gamma The parameter Gamma
#' @return return the model time in model C equivalent to Kingsman coalescent time u
################################################################################################################
standard2modelC <- function(u, Delta, Gamma)
{
  result = InverseTransformedLambda(u, Delta, Gamma)
  return(result)
  
}
################################################################################################################
#' modelCoalC
#' computes the list of coalescent times in model C
#' @param sample.size The time variable u in model time
#' @param Delta The parameter Delta
#' @param Gamma The parameter Gamma
#' @return return the list of coalescent times in model C
################################################################################################################
modelCoalC <- function(sample.size, Delta, Gamma) {
  w <- standardCoal(sample.size)
  cumulativeSTDCoaltimes <-
    unlist(lapply(1:(sample.size - 1), function(i)
      sum(w[i:(sample.size - 1)])))
  u <-
    unlist(lapply(
      cumulativeSTDCoaltimes,
      FUN = function(x)
        standard2modelC(x, Delta, Gamma)
    ))
  return(u)
}
################################################################################################################
#' nearest_value_search
#' finds the first element p in the increasing ordered list increasing_list such that
#' monotone_fun(p)>x
#' @param x value to compare the monotone_fun(tipically 0)
#' @param increasing_list: increasing ordered list
#' @param monotone_function: monotone function(increasing or decreasing)
################################################################################################################
nearest_value_search = function(x, increasing_list, monotone_function) {
  left = 1
  right = length(increasing_list)
  while (right - left > 1) {
    middle = floor((left + right) / 2)
    if (x < monotone_function(increasing_list[middle])) {
      right = middle
    }
    else{
      left = middle
    }
  }
  return(increasing_list[right])
  
}
################################################################################################################
#' search_infinite_from
#' finds the first interval [2^{p-1}*start, 2^{p}*start] with p=1,2,.., in [start, infinity]
#' such that monotone_fun(2^{p}*start)>x and  monotone_fun(2^{p-1}*start)<=x
#' @param x value to compare the monotone_fun(tipically 0)
#' @param start: start of the interval to search
#' @param monotone_function: monotone function(increasing or decreasing)
#' @output  the first positive integer 2^{p}*start such that monotone_fun(2^{p}*start)>x
################################################################################################################
search_infinite_from <- function(x, start,  monotone_function) {
  if (start > 0) {
    from = start
    to = start
    
  }
  else{
    from = 0
    to = 1
  }
  while (monotone_function(to) < x)
  {
    from = to
    to = 2 * to
  }
  increasing_list <- seq(from, to, by = 1)
  result <-
    nearest_value_search(x, increasing_list, monotone_function)
  return(result)
}
################################################################################################################
#' get_simulated_kth_coal_times
#' gets the simulated  kth coalescence times
#' @param coal.events.times.sim is a big matrix bigstatsr::FBM(length(DeltaList),sim*(sample.size-1))
#' @param k index of coalescence time  from present time backwards
#' @param sim: number of simulations
#' @param sample.size: sample size
#' @output kth_coal_times: kth coalescent times
################################################################################################################
get_simulated_kth_coal_times <-
  function(coal.events.times.sim,
           pos,
           sim,
           sample.size)
  {
    num <- length(rows_along(coal.events.times.sim))
    kth_coal_times <-
      lapply(
        1:num,
        FUN = function(k,
                       sample.size,
                       sim,
                       coal.events.times.sim,
                       pos) {
          matrixCurrentValue <-
            matrix(
              coal.events.times.sim[k, ],
              nrow = sim ,
              ncol = sample.size - 1,
              byrow = TRUE
            )
          kth_times <- matrixCurrentValue[, pos]
          
          result <- data.frame(coal_times = kth_times)
          colnames(result) <- paste("coal_times_", pos, sep = "")
          return(result)
          
        },
        sample.size,
        sim,
        coal.events.times.sim,
        pos
      )
    
    return(kth_coal_times)
  }
################################################################################################################
#' get_list_simulated_trees
#' gets a simulate tree with some  set of parameter value
#' @param coal.events.times.sim is a big matrix bigstatsr::FBM(length(DeltaList),sim*(sample.size-1))
#' @param parameter_index index of parameter value
#' @param number.sim.trees: number of simulations
#' @param sample.size: the sample size used
#' @output list_trees: the list of simulated trees
################################################################################################################
get_list_simulated_trees <-
  function(coal.events.times.sim,
           parameter_index,
           number.sim.trees,
           sample.size)
  {
    require(dplyr)
    require(ape)
    
    matrixCurrentValue <-
      matrix(
        coal.events.times.sim[parameter_index, ],
        nrow = number.sim.trees ,
        ncol = sample.size - 1,
        byrow = TRUE
      )
    
    list_trees <-
      lapply(1:number.sim.trees, function(k, sample.size, matrixCurrentValue) {
        random_simulation <- unlist(matrixCurrentValue[k, ])
        coalescent_tree = rcoal(sample.size, br = "coalescent")
        x <- dplyr::as_tibble(coalescent_tree)
        
        df <- as.data.frame(x)
        
        node_times <- c(rep(0, sample.size), rev(random_simulation))
        
        x[, "branch.length"] <-
          unlist(lapply(1:nrow(x), function(idx, node_times, x) {
            if (as.numeric(x[idx, "parent"]) != as.numeric(x[idx, "node"]))
              return(node_times[as.numeric(x[idx, "parent"])] - node_times[as.numeric(x[idx, "node"])])
            else
              return(NA)
            
          }, node_times, x))
        
        tree <- ape::as.phylo(x)
        
        return(tree)
        
      }, sample.size, matrixCurrentValue)
    
    
    return(list_trees)
  }
################################################################################################################
#' install_required_packages
#' install required packages
#' @param list.of.packages list of packages to install
#' @output None
################################################################################################################
install_required_packages = function(list.of.packages) {
  new.packages <-
    list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
  if (length(new.packages) > 0)
    install.packages(new.packages, dependencies = TRUE, repos = "http://cran.us.r-project.org")
  
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  if ("ggtree"  %in%  new.packages)
    BiocManager::install("ggtree")
  
}
################################################################################################################
#' load_required_packages
#' load required packages
#' @param list_packages list of packages to load
#' @output None
################################################################################################################
load_required_packages = function(list_packages) {
  lapply(list_packages,
         library,
         character.only = TRUE)
  
  
}
################################################################################################################
#' getCurrentFileLocation
#' get the path of the current R script
#' @output the path of the current R script
################################################################################################################
getCurrentFileLocation <-  function()
{
  require(tibble)
  this_file <- commandArgs() %>%
    tibble::enframe(name = NULL) %>%
    tidyr::separate(
      col = value,
      into = c("key", "value"),
      sep = "=",
      fill = 'right'
    ) %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
  if (length(this_file) == 0)
  {
    this_file <- rstudioapi::getSourceEditorContext()$path
  }
  return(dirname(this_file))
}
################################################################################################################
#' computeCorrelationMatrixModel
#' computes correlation matriz between coalescent times for a model
#' @param coal.events.times.sim the simulated coalescent times bigstatsr::FBM(length(DeltaList),sim*(sample.size-1))
#' @param DeltaList list of Delta values
#' @param sample.size the sample size at time t=0
#' @param sim the  number of simulations used to generate coal.events.times.sim
#' @output the list of correlation matrices(one for each Delta value)
################################################################################################################
computeCorrelationMatrixModel <-
  function(coal.events.times.sim,
           DeltaList,
           sample.size,
           sim)
  {
    require(lineup)
    list_corr_mat <-
      parallel::mclapply(
        1:length(DeltaList),
        mc.cores = parallel::detectCores() - 1,
        FUN = function(k,
                       sample.size,
                       sim,
                       coal.events.times.sim)
        {
          matrixCurrentValue <-
            matrix(
              coal.events.times.sim[k, ],
              nrow = sim ,
              ncol = sample.size - 1,
              byrow = TRUE
            )
          
          result <-
            cor(matrixCurrentValue)
          return(as.matrix(result))
          
        },
        sample.size,
        sim,
        coal.events.times.sim
      )
    
    return(list_corr_mat)
    
  }
