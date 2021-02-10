#R code accompanying the paper "Coalescent Models Derived from a
#Birth-Death Process" by Crespo and Wiuf

# Author of R code: Fausto Crespo
#####################################################
rm(list=ls())
####################################################
#####################################################################
#relative
#relative population size at time t given time origin(TOrigin) for scenario B
#TOrigin: time of origin 
#Delta: scaled growth rate (scaled by N/birth_rate )
#t: time in mode time(scaled in units of N/birth_rate )
#####################################################################
relative <- function(Torigin,Delta,t)
  exp(-Delta*t)*(1-exp(-Delta*(Torigin-t)))^2/(1-exp(-Delta*Torigin))^2

#####################################################################
#reciprocRel
#relative population size at time t given time origin(TOrigin) for scenario B(New version)
#TOrigin: time of origin 
#Delta: scaled growth rate (scaled by N/birth_rate )
#t: time in mode time(scaled in units of N/birth_rate )
#####################################################################
lambdaB <- function(u,d, g, t) {
  x <- exp(d*u)*(1-exp(-d*t))^2/( 1-exp(-d*(t-u)) )^2
  x <- (2.0*d/g) *x*(1 - exp(-g*( 1/(1-exp(-d*u))-1/(1-exp(-d*t)) )))
 # x <- exp(Delta*u)*(1-exp(-Delta*Torigin))^2/( 1-exp(-Delta*(Torigin-u)) )^2
  #x <- (2 * Delta / Gamma)* x *(1 - exp(-Gamma*( 1.0/(1.0-exp(-Delta*u))-1.0/(1.0-exp(-Delta*Torigin)) )))
  return(x)
}

LambdaB<-function(s,Delta, Gamma, Torigin){
  if (s>0 & s<Torigin){
    upper.limit= min(s, Torigin)
    res <- integrate(Vectorize(lambdaB, c("u")),d=Delta, g=Gamma, t=0.99999999999*Torigin, lower = 0, upper = s,stop.on.error = FALSE)$value
    return(res)
  }
  else{
    if (s > Torigin)
      print(paste("s=", s, " is greater than Torigin=", Torigin))
    return(0)
  }
}
functionToZeroB<-function( y, Delta,Gamma, Torigin)
{
  f<-function(x) LambdaB(x, Delta,Gamma, Torigin)- y
  return(f)
}

InverseLambdaB<-function(s, Delta,Gamma, Torigin)
{
  #upper<-search_infinite_from(0, 0, functionToZeroB( s, Delta,Gamma, Torigin))
  #result<-cmna::bisection(functionToZeroB( s, Delta, Gamma, Torigin), min(0,upper-1),upper, tol=10^(-6))
  result <- cmna::bisection(functionToZeroB( s, Delta, Gamma, Torigin), 0, Torigin, tol=10^(-6))
  return(result)
}
# Standard time to model time B
standard2model <- function(Torigin,Delta,Gamma,u) 
  ( log(1+Delta*u-exp(-Delta*Torigin)) - log(1+(Delta*u-1)*exp(-Delta*Torigin)) )/Delta

# Standard time to model time B version M_K
standard2modelB_K <- function(u,Torigin,Delta,Gamma, K) {
  minus.exp.Delta.Torigin = exp(-Delta*Torigin)
  a = (K / Gamma) * (1 - minus.exp.Delta.Torigin )- minus.exp.Delta.Torigin
  #c= 1- (2.0/ Gamma)*(1 - minus.exp.Delta.Torigin)
  b= 1- (K/ Gamma)*(1 - minus.exp.Delta.Torigin)
  #d= a*(1-minus.exp.Delta.Torigin)*(1-minus.exp.Delta.Torigin)
  #K= exp(Gamma* u * (a*b+c) / (2.0 * d) + log((b+c)/(a-1)))
  #x=(K*a-c) /(K+b)
  #x=(exp(u*K/2.0)- b ) /( exp(u*K /2.0 -Delta*Torigin) + a)#version of December 31, 2020
  #result=(1.0/Delta)*log(x)#version of December 31, 2020
  x=exp(u*K/2.0)
  y=(x- b ) /( x*exp(-Delta*Torigin) + a)
  result=(1.0/Delta)*log(y)
  return(result)
}
# Standard time to model B version time M*
standard2modelB <- function(u,Torigin,Delta,Gamma, K) {
  minus.exp.Delta.Torigin = exp(-Delta*Torigin)
  a = (K / Gamma) * (1 - minus.exp.Delta.Torigin )- minus.exp.Delta.Torigin
  #c= 1- (2.0/ Gamma)*(1 - minus.exp.Delta.Torigin)
  b= 1- (K/ Gamma)*(1 - minus.exp.Delta.Torigin)
  #d= a*(1-minus.exp.Delta.Torigin)*(1-minus.exp.Delta.Torigin)
  #K= exp(Gamma* u * (a*b+c) / (2.0 * d) + log((b+c)/(a-1)))
  #x=(K*a-c) /(K+b)
  x=(exp(u*K/2.0)- b ) /( exp(u*K /2.0 -Delta*Torigin) + a)
  result=(1.0/Delta)*log(x)
  return(result)
}
standard2modelB_MAster <- function(u, Torigin,Delta,Gamma) {
  
  result=InverseLambdaB(u, Delta, Gamma, Torigin)
  return(result)
}

# Standard time to model time B version 2
standard2modelBOld <- function(Torigin,Delta,Gamma, u) {
  a = exp(-1* Delta*Torigin)
  b = Gamma *u/ (2*(1-a))
  x=(b+1)/(a*b+1)
  result=(1.0/Delta)*log(x)
  return(result)
}

# Coalescence times with constant population size
standardCoal <- function(n) unlist(lapply(2:n,function(i) rexp(1,rate=choose(i,2))))

# Coalescence times with constant population size from some sample size p to n
standardCoalFrom <- function(p,n) unlist(lapply(p:n,function(i) rexp(1,rate=choose(i,2))))

# Coalescence times in model time(inter event times)
modelCoal <- function(n,Torigin,Delta, Gamma) {
  w <- standardCoal(n)
  # time until i=n-1,...,1 ancestors
  u <- standard2model( Torigin,Delta,Gamma, unlist( lapply(1:(n-1),function(i) sum(w[i:(n-1)])) ) )
  u <- c(u,0)
  unlist(lapply(1:(n-1),function(i) u[i]-u[i+1]))
}
# Coalescence event  times in model time
modelCoalB_K <- function(n,Torigin,Delta, Gamma, K) {
  w <- standardCoal(n)
  # time until i=n-1,...,1 ancestors
  cumulativeSTDCoaltimes <- unlist( lapply(1:(n-1),function(i) sum(w[i:(n-1)])) )
  u <- unlist(lapply( cumulativeSTDCoaltimes, FUN=function(x) standard2modelB_K(x,Torigin,Delta,Gamma,K)))
  return(u)
}
# Coalescence event  times in model time
modelCoalB_MAster <- function(n,Torigin,Delta, Gamma) {
  w <- standardCoal(n)
  # time until i=n-1,...,1 ancestors
  cumulativeSTDCoaltimes <- unlist( lapply(1:(n-1),function(i) sum(w[i:(n-1)])) )
  u <- unlist(lapply( cumulativeSTDCoaltimes, FUN=function(x) standard2modelB_MAster(x,Torigin,Delta,Gamma)))
  return(u)
}
# Coalescence event  times in model time
modelCoalHybrid <- function(n,Torigin,Delta, Gamma, nB) {
  w <- standardCoalFrom(nB+1,n)
  u <-  unlist( lapply(1:(n-nB),function(i) sum(w[i:(n-nB)])) ) 
  return(u)
}
# Coalescence times with exponentially decreasing population size
expCoal <- function(n,beta) {
  w <- standardCoal(n)
  # time until i=n-1,...,1 ancestors
  u <- unlist( lapply(1:(n-1),function(i) sum(w[i:(n-1)])) )
  u <- log(1+beta*u)/beta
  u <- c(u,0)
  unlist(lapply(1:(n-1),function(i) u[i]-u[i+1]))
}

#######################################################################
compute.gamma=function(lambda, mu, rho)
{
  output= abs(lambda-mu)/max(rho*lambda, mu-(1-rho)*lambda)
  return(output)
}
log.prod.between=function(from,to){
  result=0
  stopifnot(from>=0)
  stopifnot(to>=0)
 

  if (from >to)
    return(result)
  
  if (min(from, to)>0)
  {
    if (to >=from ){
      return(sum(log(seq(max(from,1),to, by=1))))
    }
  }
 else 
   return(result)
}
log.prod.from.To=function(from,to){
  result=0
  stopifnot(from>=0)
  stopifnot(to>=0)

  
  if (min(from, to)>=0)
  {
    if (to >from ){
      return(sum(log(seq(max(from+1,1),to, by=1))))
    }
    else if (from >to )
    {
      return(-1.0*sum(log(seq(max(to+1,1),from, by=1))))
    }
    else
    {
      result=0
      return(result)
    }
  }
  else 
    return(result)
}
ToStandardTime<-function(t, Time1, Delta )
{
  a=exp(-1*Delta*Time1)
  c=exp(Delta*t)
  result<-(c-1)*(1-a)/(1-a*c)
  return(result)
}
ToModelTime<-function(u, Time1, Delta )
{
  a=exp(-1*Delta*Time1)
  b=1+u-a
  c=1+(u-1)*a
  if ((b== 0) || (c==0))
  {
    return((1.0/Delta)*(log(1+u)-(u*u*a)/(1+u)))
  } 
  else{
    return((1.0/Delta)*log(b)-log(c))
  }
}
sample.conditional.coalescent.time.distribution = function(Delta, alpha, k, m){
  require(stats)
  stopifnot(k>0)
  stopifnot(m-k+1>0)
  U= rbeta(1, k, m-k+1)
  stopifnot(U/(alpha+(1-alpha)*U)>0)
  T=(-1.0/Delta)* log(U/(alpha+(1-alpha)*U))
  return(T)
}
sample.first.coalescent.time.distribution = function(Delta, Gamma, k){
  require(stats)
  U= stats::rgamma(1, shape=k, rate=1)
  T=(-1.0 / Delta) * log(1.0 - (Gamma/(U+Gamma)))
  return(T)
}

conditional.density.first.time.k.ancestors = function(delta, gamma, n, time.Origin, times){
  exp.delta.t = exp(-1*delta* time.Origin)
  output = ((delta* gamma)**(n-1)) * (1+ (gamma-1)*exp.delta.t)**(n-1) / (1- exp.delta.t)**(n-1)
  
  for(i in  2:n)
  {
    if (i<n){
      
      output = i* output *exp(-1*delta* time.Origin) /  (1+ (gamma-1)*(exp.delta.t)**2)
    }
    else{
      
      output = output *exp(-1*delta* time.Origin) /  (1+ (gamma-1)*(exp.delta.t)**2)
    }
  }
  return(output)
}
density.first.time.k.ancestors = function(delta, gamma, n, time.Origin, times){
  exp.delta.t = exp(-1*delta* time.Origin)
  output = ((delta* gamma)**(n)) 
  for(i in  1:n)
  {
    exp.delta.t = exp(-1*delta* times[i])
    output = i* output *exp(-1*delta* time.Origin) /  (1+ (gamma-1)*((exp.delta.t)**2))
    
  }
  return(output)
}
solve.cumulative.density.number.ancestors.population.present.time=function(u, sample.size)
{
 
  fun <- function (x)  exp(log.prod.between(x-sample.size+2,x)-log.prod.between(x+2,x+sample.size))-u
  fun2 <- function (x)  (u-exp(log.prod.between(x-sample.size+2,x)-log.prod.between(x+2,x+sample.size)))^2
  
  result<-search_infinite_from( 0, (sample.size-2), fun ) 

  result<-ceiling(result)
  return(result)

}
solve.cumulative.density.number.ancestors.population=function(u, k,m)
{
   fun <- function (x) {
        result=exp(log.prod.between(x-k+2,x)-log.prod.between(x+2,x+k)+ log.prod.from.To(m-1,m+k-1) + log.prod.from.To(m,m-k))-u
      return(result)
   } 
   fun2 <- function (x) {
     result=(u-exp(log.prod.between(x-k+2,x)-log.prod.between(x+2,x+k)+ log.prod.from.To(m-1,m+k-1) + log.prod.from.To(m,m-k)))^2
     return(result)
   }
   
   interval<-seq((k-2),m, by=1)
   
   result<-nearest_value_search(0,interval,fun)

   result<-ceiling(result)
  
  return(result)
}
#####################################################################
#get_number_of_ancestors_population_when_sample_size_minus_1_sample
#sample the  number of ancestors in the coalescent approximation 
# when there are sample_size-1 ancestors in the coalescent appproximation 
# in the sample. The number ancestors in the population is assumed to be infinity
#sample_size: current number ancestors in the sample(present time) 
#output: the number of ancestors when there are k-1 ancestors in the sample
# the output >=sample_size-1 
#####################################################################
get_number_of_ancestors_population_when_sample_size_minus_1_sample=function(sample_size)
{
  
  accepted=FALSE
  while(!accepted){
    u=runif(1,0,1)
    mprime=solve.cumulative.density.number.ancestors.population.present.time(u, sample_size)
    if(mprime >=(sample_size - 1)){
     
      accepted=TRUE
    }
  }
  return(mprime)
}
#####################################################################
#get_number_of_ancestors_population_when_k_minus_1_ancestors_sample
#sample the number of ancestors in the coalescent approximation 
# when there are k-1 ancestors in the coalescent appproximation 
# in the sample
#k: current number ancestors in the sample 
#m: current number ancestors in the population
#output: the number of ancestors when there are k-1 ancestors in the sample
# the output >=k-1 and < m 
#####################################################################
get_number_of_ancestors_population_when_k_minus_1_ancestors_sample=function(k,m)
{
  accepted=FALSE
  while(!accepted)
  {
    u=runif(1,0,1)
    mprime=solve.cumulative.density.number.ancestors.population(u, k, m)
    
      if(mprime <m & mprime >= (k-1))
      {
        accepted=TRUE
      }
  }
  return(mprime)
}
simulate.list.number.ancestors.population=function(sample.size){

  population.present.time= get_number_of_ancestors_population_when_sample_size_minus_1_sample(sample.size)
  current.population.size.when.sample.minus1=population.present.time
  list.populations.sizes=list(current.population.size.when.sample.minus1)
  m=current.population.size.when.sample.minus1
  rest.population.sizes=list(current.population.size.when.sample.minus1)
  p=2
  for (k in (sample.size-1):1)
    {
    m=get_number_of_ancestors_population_when_k_minus_1_ancestors_sample(k,m)
    list.populations.sizes[[p]]=m
    p=p+1
  }

  return(unlist(list.populations.sizes))
}
simulate.list.number.ancestors.population.from=function(sample.size, current.population.size, current.sample){

  stopifnot(current.sample <=sample.size)
  
  list.populations.sizes=list(current.population.size)
  if (sample.size ==current.sample){
    
    current.population.size.when.sample.minus1= get_number_of_ancestors_population_when_sample_size_minus_1_sample(sample.size)
    m=current.population.size.when.sample.minus1
  }
  else{
    m=current.population.size
    
  }
  list.populations.sizes=list(m)
  
  p=2
  for (k in (current.sample):2)
  {
    m=get_number_of_ancestors_population_when_k_minus_1_ancestors_sample(k,m)
    list.populations.sizes[[p]]=m
    p=p+1
  }
  
  return(unlist(list.populations.sizes))
}

simulate.coalescent.times.A=function(Gamma, Delta, sample.size, list.number.ancestors.population)
{
  list.coal.times=list()
  m= list.number.ancestors.population[1]
  u=unlist(sample.first.coalescent.time.distribution(Delta, Gamma, m+1))
  u=u[[1]]
  current.time=u
  list.coal.times=list(u)
  p=2
  for(i in  (sample.size-2):0)
  {
    m= list.number.ancestors.population[sample.size-i]
    current.m=list.number.ancestors.population[sample.size-i-1]
    alpha=1 - exp(-1.0 * Delta * current.time)
    u=unlist(sample.conditional.coalescent.time.distribution(Delta, alpha, m+1, current.m))
    u=u[[1]]
    current.time= current.time+u
    list.coal.times[p]=current.time
    p=p+1
  }
  return(unlist(list.coal.times))
}
simulate.coalescent.times.A.from=function(Gamma, Delta, sample.size, list.number.ancestors.population, current.sample, from_time)
{
  stopifnot(current.sample<=sample.size)
  list.coal.times=list()
  m= list.number.ancestors.population[1]
   if (current.sample==sample.size){
     u=unlist(sample.first.coalescent.time.distribution(Delta, Gamma, m+1))
     u=u[[1]]
     list.coal.times=list(u)
     p=2
     
   }
  else{
     u=from_time
     list.coal.times=list()
     p=1
  }
  
  current.time=u
 
  for(i in  (current.sample-2):0)
  {
    m= list.number.ancestors.population[current.sample-i]
    current.m=list.number.ancestors.population[current.sample-i-1]
    alpha=1 - exp(-1.0 * Delta * current.time)
    u=unlist(sample.conditional.coalescent.time.distribution(Delta, alpha, m+1, current.m))
    u=u[[1]]
    current.time= current.time+u
    list.coal.times[p]=current.time
    p=p+1
  }
  return(unlist(list.coal.times))
}
simulate.coalescent.times.A.test=function(Gamma, Delta, sample.size, list.number.ancestors.population)
{
  p=1
  list.coal.times<-list()
  list.coal.times<-lapply((sample.size-1):0, FUN=function(i, sample.size,list.number.ancestors.population )
  {
   2.0 * ( 1.0 /list.number.ancestors.population[sample.size-i])
    
   }, sample.size, list.number.ancestors.population)

  return(unlist(list.coal.times))
}
expected_number_ancestors =function(sample.size){
  
  list.number.ancestors<-lapply((sample.size-2):1, FUN=function(k, sample.size)
    {
    
    sample.size*k/(sample.size-k-1)
   }, sample.size)
  
  return(unlist(list.number.ancestors))
}
compute_stats_number_ancestors=function(sample.size,sim,DeltaList, number.ancestors.simA){
  
  list.Data.Frames.Numbers.Ancestors<-parallel::mclapply(1:length(DeltaList),mc.cores=parallel::detectCores()-1, FUN=function(k, sample.size,sim, number.ancestors.simA ){
    
    quants <- c(0.025,0.50,0.975)
    
    matrixCurrentValue<-matrix(number.ancestors.simA[k,], nrow=sim , ncol=sample.size-1, byrow=TRUE)
    #quantiles<-apply( matrixCurrentValue , 2 , quantile , probs = quants , na.rm = TRUE )
    quantiles<-apply( matrixCurrentValue , 2 , quantile , probs = quants  )
    mean.Number.Ancestors<-colMeans(matrixCurrentValue)
    result <-data.frame(mean=mean.Number.Ancestors, LI= quantiles[1,], median=quantiles[2,], UI=quantiles[3,])
    result
  },sample.size,sim, number.ancestors.simA )
  
  return(list.Data.Frames.Numbers.Ancestors)

}
simulateB_K=function(i,Delta, Gamma,  sim, sample.size, Time.Origin.STD, K)
  {
  require(parallel)
  coal.events.times.simB <- array(0,dim=c(sim,sample.size-1))
  
  lapply(1:sim,
    FUN=function(j,sample.size, Time.Origin.STD, Delta, Gamma, coal.events.times.simB, i )
  {
    TimeOrigin=Time.Origin.STD[i,j]
    coal.events.times= modelCoalB_K(sample.size,TimeOrigin,Delta, Gamma, K)
    coal.events.times.simB[j,] <<- rev(coal.events.times)
    print(paste0("finished B sim",j, sep=" "))
  },  sample.size, Time.Origin.STD,Delta, Gamma, coal.events.times.simB, i)

  quants <- c(0.025,0.50,0.975)
  #quantiles<-apply( coal.events.times.simB , 2 , quantile , probs = quants , na.rm = TRUE )
  quantiles<-apply( coal.events.times.simB , 2 , quantile , probs = quants )
  meanCoalTimes<-colMeans(coal.events.times.simB)
  result <-data.frame(mean=meanCoalTimes, LI= quantiles[1,], median=quantiles[2,], UI=quantiles[3,])
  return(result)
}

simulateB_K.parallel=function(DeltaList,GammaList, sim, sample.size, Time.Origin.STD, coal.events.times.simB, K)
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
  set.seed(7596034) #set seed to something
  s <- .Random.seed
  ncores=parallel::detectCores()-1
  innerCluster <- parallel::makeCluster(ncores, type="FORK", outfile="")
   on.exit(parallel::stopCluster(innerCluster), add = TRUE)
  doParallel::registerDoParallel(innerCluster)
  
  #rng<- RNGseq( length(DeltaList)* sim, 5363131)
  #rng<- RNGseq( length(DeltaList)* sim, 32252)
  
  rng<- RNGseq( length(DeltaList)* sim, 53727)
  opts <- list(chunkSize=2)

   tmp3 <-foreach::foreach(j=1:sim, .combine = 'c') %:%
    foreach::foreach(i=1:length(DeltaList), r=rng[(j-1)*length(DeltaList) + 1:length(DeltaList)], .combine = 'c'
            ) %dopar% {
      rngtools::setRNG(r)

    TimeOrigin=Time.Origin.STD[i,j]
   
    coal.events.times= modelCoalB_K(sample.size,TimeOrigin,DeltaList[i], GammaList[i], K)
  
    #print(paste("time origin=",TimeOrigin,"coal event times=", coal.events.times ))
    #print(unlist(coal.events.times))
    
    print(paste("finished B M_k in parallel for Delta", DeltaList[i], "Gamma", GammaList[i], "K",K, " sim ",j, sep=" "))
    positions <-cbind(rep(i, (sample.size-1)),((j-1)*(sample.size-1) +1):(j*(sample.size-1)))
    #print(positions)
    coal.events.times.simB[positions] <- rev(coal.events.times)
    NULL
  }

   
 
  list.Data.FramesB<-parallel::mclapply(1:length(DeltaList),mc.set.seed = TRUE,mc.cores=parallel::detectCores()-1, 
    FUN=function(k, sample.size,sim, coal.events.times.simB )
    {

     quants <- c(0.025,0.50,0.975)
     matrixCurrentValue<-matrix(coal.events.times.simB[k,], nrow=sim , ncol=sample.size-1, byrow=TRUE)
     #quantiles<-apply( matrixCurrentValue , 2 , quantile , probs = quants , na.rm = TRUE )
     quantiles<-apply( matrixCurrentValue , 2 , quantile , probs = quants  )
     meanCoalTimes<-colMeans(matrixCurrentValue)
     result <-data.frame(mean=meanCoalTimes, LI= quantiles[1,], median=quantiles[2,], UI=quantiles[3,])
     result
     },sample.size,sim, coal.events.times.simB )

    return(list.Data.FramesB)
}

simulateB_MAster.parallel=function(DeltaList,GammaList, sim, sample.size, Time.Origin.STD, coal.events.times.simB)
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
  set.seed(7596034) #set seed to something
  s <- .Random.seed
  ncores=parallel::detectCores()-1
  innerCluster <- parallel::makeCluster(ncores, type="FORK", outfile="")
  on.exit(parallel::stopCluster(innerCluster), add = TRUE)
  doParallel::registerDoParallel(innerCluster)
  
  #rng<- RNGseq( length(DeltaList)* sim, 5363131)
  #rng<- RNGseq( length(DeltaList)* sim, 32252)
  
  rng<- RNGseq( length(DeltaList)* sim, 53727)
  opts <- list(chunkSize=2)
  
  tmp3 <-foreach::foreach(j=1:sim, .combine = 'c') %:%
    foreach::foreach(i=1:length(DeltaList), r=rng[(j-1)*length(DeltaList) + 1:length(DeltaList)], .combine = 'c'
    ) %dopar% {
      rngtools::setRNG(r)
      
      TimeOrigin=Time.Origin.STD[i,j]
      
      coal.events.times= modelCoalB_MAster(sample.size,TimeOrigin,DeltaList[i], GammaList[i])
      
      #print(paste("time origin=",TimeOrigin,"coal event times=", coal.events.times ))
      #print(unlist(coal.events.times))
      
      print(paste("finished B M Aster in parallel for Delta", DeltaList[i], "Gamma", GammaList[i],  " sim ",j, sep=" "))
      positions <-cbind(rep(i, (sample.size-1)),((j-1)*(sample.size-1) +1):(j*(sample.size-1)))
      #print(positions)
      coal.events.times.simB[positions] <- rev(coal.events.times)
      NULL
    }
  
  
  
  list.Data.FramesB<-parallel::mclapply(1:length(DeltaList),mc.set.seed = TRUE,mc.cores=parallel::detectCores()-1, 
                                        FUN=function(k, sample.size,sim, coal.events.times.simB )
                                        {
                                          
                                          quants <- c(0.025,0.50,0.975)
                                          matrixCurrentValue<-matrix(coal.events.times.simB[k,], nrow=sim , ncol=sample.size-1, byrow=TRUE)
                                          #quantiles<-apply( matrixCurrentValue , 2 , quantile , probs = quants , na.rm = TRUE )
                                          quantiles<-apply( matrixCurrentValue , 2 , quantile , probs = quants )
                                          meanCoalTimes<-colMeans(matrixCurrentValue)
                                          result <-data.frame(mean=meanCoalTimes, LI= quantiles[1,], median=quantiles[2,], UI=quantiles[3,])
                                          result
                                        },sample.size,sim, coal.events.times.simB )
  
  return(list.Data.FramesB)
}
simulateC=function(Delta, Gamma,  sim, sample.size)
{
  require(parallel)
  coal.events.times.simC <- array(0,dim=c(sim,sample.size-1))
  
  lapply(1:sim,FUN=function(i,sample.size,  Delta, Gamma, coal.events.times.simC )
  {
    coal.events.times= modelCoalC(sample.size, Delta, Gamma)
    coal.events.times.simC[i,] <<- rev(coal.events.times)
    print(paste0("finished C sim",i, sep=" "))
  },  sample.size, Delta, Gamma, coal.events.times.simC)
  
  quants <- c(0.025,0.50,0.975)
  #quantiles<-apply( coal.events.times.simC , 2 , quantile , probs = quants , na.rm = TRUE )
  quantiles<-apply( coal.events.times.simC , 2 , quantile , probs = quants  )
  meanCoalTimes<-colMeans(coal.events.times.simC)
  result <-data.frame(mean=meanCoalTimes, LI= quantiles[1,], median=quantiles[2,], UI=quantiles[3,])
  return(result)
}
simulateC.parallel=function(DeltaList,GammaList, sim, sample.size, coal.events.times.simC)
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
  set.seed(69012365) #set seed to something
  s <- .Random.seed
  ncores=parallel::detectCores()-1
  innerCluster <- parallel::makeCluster(ncores, type="FORK", outfile="")
   on.exit(parallel::stopCluster(innerCluster), add = TRUE)
  doParallel::registerDoParallel(innerCluster)
  rng<- RNGseq( length(DeltaList)* sim, 1234)
 
  
  opts <- list(chunkSize=2)

  tmp3 <-foreach::foreach(j=1:sim, .combine = 'c') %:%
    foreach::foreach(i=1:length(DeltaList), r=rng[(j-1)*length(DeltaList) + 1:length(DeltaList)], .combine = 'c'
            
            ) %dopar% {
      rngtools::setRNG(r)
     
     coal.events.times= modelCoalC(sample.size, DeltaList[i], GammaList[i])

     positions <-cbind(rep(i, (sample.size-1)),((j-1)*(sample.size-1) +1):(j*(sample.size-1)))
     
      coal.events.times.simC[positions] <- rev(coal.events.times)
      print(paste("finished C in parallel for Delta", DeltaList[i], "Gamma", GammaList[i], " sim ",j, sep=" "))
     NULL
  }
  
  list.Data.FramesC<-parallel::mclapply(1:length(DeltaList),mc.cores=parallel::detectCores()-1, FUN=function(k, sample.size,sim, coal.events.times.simC){

     quants <- c(0.025,0.50,0.975)
     matrixCurrentValue<-matrix(coal.events.times.simC[k,], nrow=sim , ncol=sample.size-1, byrow=TRUE)
     #quantiles<-apply( matrixCurrentValue , 2 , quantile , probs = quants , na.rm = TRUE )
     quantiles<-apply( matrixCurrentValue , 2 , quantile , probs = quants  )
     meanCoalTimes<-colMeans(matrixCurrentValue)
     result <-data.frame(mean=meanCoalTimes, LI= quantiles[1,], median=quantiles[2,], UI=quantiles[3,])
     result
     },sample.size,sim, coal.events.times.simC )

  
    return(list.Data.FramesC)
  
}
simulateA=function(Delta,Gamma, sim, sample.size, Time.Origin.STD,indexA, coal.events.times.simA)
  {
  list.number.ancestors.population.sim <<-array(0,dim=c(sim,sample.size))
  
  Time.Origin.STD.sim <- array(0,dim=c(sim))

  lapply(1:sim, FUN=function(j, sample.size, Gamma, Delta, coal.events.times.simA,  Time.Origin.STD.sim)
  {
    list.number.ancestors.population= simulate.list.number.ancestors.population(sample.size)

    list.number.ancestors.population.sim[j,] <<-list.number.ancestors.population
    list.coal.times= simulate.coalescent.times.A(Gamma, Delta,sample.size, list.number.ancestors.population)

    coal.events.times =unlist(list.coal.times)

    Time.Origin= coal.events.times[length(coal.events.times)]

    Time.Origin.STD.sim[j]<<- Time.Origin

    coal.events.times.simA[j,] <<- coal.events.times[1:(length(coal.events.times)-1)]
    print(paste0("finished A sim ",j, "for Delta=", Delta, "Gamma=", Gamma, sep=" "))
    
  },
  sample.size, Gamma, Delta, coal.events.times.simA,
  Time.Origin.STD.sim=Time.Origin.STD.sim)
  
  Time.Origin.STD[cbind(rep(indexA,sim), 1:sim)]<-Time.Origin.STD.sim
  quants <- c(0.025,0.50,0.975)
  #quantiles<-apply( coal.events.times.simA , 2 , quantile , probs = quants , na.rm = TRUE )
  quantiles<-apply( coal.events.times.simA , 2 , quantile , probs = quants  )
  meanCoalTimes<-colMeans(coal.events.times.simA)
  result <-data.frame(mean=meanCoalTimes, LI= quantiles[1,], median=quantiles[2,], UI=quantiles[3,])
  return(result)
}
simulateA.parallel=function(DeltaList,GammaList, sim, sample.size, Time.Origin.STD, coal.events.times.simA, number.ancestors.simA, number.ancestors.Transition)
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
  set.seed(728493) #set seed to something
  s <- .Random.seed
  ncores=parallel::detectCores()-1
  innerCluster <- parallel::makeCluster(ncores, type="FORK", outfile="")
  on.exit(parallel::stopCluster(innerCluster))
  doParallel::registerDoParallel(innerCluster)
  #rng<- RNGseq( length(DeltaList)* sim, 6452913)

  rng<- RNGseq( length(DeltaList)* sim, 53727)
  
  opts <- list(chunkSize=2)

  tmp3 <-foreach::foreach(j=1:sim, .combine = 'c') %:%
    foreach::foreach(i=1:length(DeltaList), r=rng[(j-1)*length(DeltaList) + 1:length(DeltaList)], .combine = 'c'
            
            ) %dopar% {
      rngtools::setRNG(r)
     
      list.number.ancestors.population= simulate.list.number.ancestors.population(sample.size)
      
  
      list.coal.times= simulate.coalescent.times.A(GammaList[i], DeltaList[i],sample.size, list.number.ancestors.population)

      coal.events.times =unlist(list.coal.times)
     
      Time.Origin= coal.events.times[length(coal.events.times)]
      
      Time.Origin.STD[i,j]<- Time.Origin

      positions <-cbind(rep(i, (sample.size-1)),((j-1)*(sample.size-1) +1):(j*(sample.size-1)))
     
      number.ancestors.simA[positions] <- list.number.ancestors.population[2:(length(coal.events.times))]
      
      coal.events.times.simA[positions] <- coal.events.times[1:(length(coal.events.times)-1)]
      m= list.number.ancestors.population[floor(length(list.number.ancestors.population)/2)]
      number.ancestors.Transition[cbind(i, j)] <- m
      
      
      print(paste("finished A in parallel for Delta", DeltaList[i], "Gamma", GammaList[i], " sim ",j, sep=" "))
      NULL
    }
  
   
   
   list.Data.FramesA<-parallel::mclapply(1:length(DeltaList),mc.cores=parallel::detectCores()-1, FUN=function(k, sample.size,sim, coal.events.times.simA ){

     quants <- c(0.025,0.50,0.975)

     matrixCurrentValue<-matrix(coal.events.times.simA[k,], nrow=sim , ncol=sample.size-1, byrow=TRUE)
     #quantiles<-apply( matrixCurrentValue , 2 , quantile , probs = quants , na.rm = TRUE )
     quantiles<-apply( matrixCurrentValue , 2 , quantile , probs = quants )
     meanCoalTimes<-colMeans(matrixCurrentValue)
     result <-data.frame(mean=meanCoalTimes, LI= quantiles[1,], median=quantiles[2,], UI=quantiles[3,])
     result
     },sample.size,sim, coal.events.times.simA )

  
    return(list.Data.FramesA)

}
simulateHybrid.parallel=function(DeltaList,GammaList, sim, sample.size, Time.Origin.STD, coal.events.times.simHybrid, nB, number.Fail.hybrid, posMprime, K, number.ancestors.simHybrid, number.ancestors.Transition )
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
  set.seed(69012365) #set seed to something
  s <- .Random.seed
  ncores=parallel::detectCores()-1
  innerCluster2 <- parallel::makeCluster(ncores, type="FORK", outfile="")

  #on.exit(parallel::stopCluster(innerCluster), add = TRUE)
  on.exit(parallel::stopCluster(innerCluster2))
  doParallel::registerDoParallel(innerCluster2)
  rng<- RNGseq( length(DeltaList)* sim, 234567)
  
  opts <- list(chunkSize=2)
  
  
  tmp3 <-foreach::foreach(j=1:sim, .combine = 'c') %:%
    foreach::foreach(i=1:length(DeltaList), r=rng[(j-1)*length(DeltaList) + 1:length(DeltaList)], .combine = 'c'
                     
    ) %dopar% {
      rngtools::setRNG(r)
      
      TimeOrigin<-Time.Origin.STD[i,j]
      
      m= -1
      while(m < nB){
        
        if (m!= -1){
          
          number.Fail.hybrid[cbind(posMprime, i)]=number.Fail.hybrid[cbind(posMprime, i)]+1
        }

        std.coal.events.times.until.B= rev(modelCoalHybrid(sample.size,TimeOrigin,DeltaList[i], GammaList[i], nB))
        
        #coal.events.times.until.nB=standard2modelB_K( unlist(std.coal.events.times.until.B), TimeOrigin,DeltaList[i],GammaList[i],K)

         coal.events.times.until.nB= unlist(lapply( std.coal.events.times.until.B, FUN=function(x) standard2modelB_K(x,TimeOrigin,DeltaList[i],GammaList[i],K)))

        t= std.coal.events.times.until.B[length(std.coal.events.times.until.B)]
        
        m= ceiling(2.0 / t)
        
        
       # m = sample.size*nB/(sample.size-nB-1)
        
      }
        number.ancestors.Transition[cbind(i, j)] <- m
      

        list.number.ancestors.population= simulate.list.number.ancestors.population.from(sample.size,m, nB)
        
        s<-standard2modelB_K(t, TimeOrigin,DeltaList[i],GammaList[i],K)

        
        list.coal.times=simulate.coalescent.times.A.from(GammaList[i], DeltaList[i], sample.size, list.number.ancestors.population, nB, s)
        
        coal.events.times =c(unlist(coal.events.times.until.nB), unlist(list.coal.times))
   
       
      positions <-cbind(rep(i, (sample.size-1)),((j-1)*(sample.size-1) +1):(j*(sample.size-1)))
      
      coal.events.times.simHybrid[positions] <- unlist(coal.events.times)
      
      number.ancestors.simHybrid[positions] <- c(rep(0, sample.size-nB-2), unlist(list.number.ancestors.population), 0)
      
      #list.number.ancestors.population[2:(length(coal.events.times))]
      
      print(paste("finished hybrid scenario in parallel for Delta", DeltaList[i], "Gamma", GammaList[i], " sim ",j, sep=" "))
      NULL
    }
  
  list.Data.FramesHybrid <- parallel::mclapply(1:length(DeltaList),mc.cores=parallel::detectCores()-1, FUN=function(k, sample.size,sim, coal.events.times.simHybrid ){
    
    quants <- c(0.025,0.50,0.975)

    matrixCurrentValue<-matrix(coal.events.times.simHybrid[k,], nrow=sim , ncol=sample.size-1, byrow=TRUE)
    #quantiles<-apply( matrixCurrentValue , 2 , quantile , probs = quants , na.rm = TRUE )
    quantiles<-apply( matrixCurrentValue , 2 , quantile , probs = quants )
    meanCoalTimes<-colMeans(matrixCurrentValue)
    result <-data.frame(mean=meanCoalTimes, LI= quantiles[1,], median=quantiles[2,], UI=quantiles[3,])
    result
  },sample.size,sim, coal.events.times.simHybrid )


  return(list.Data.FramesHybrid)
  
}
simulateA.test=function(Delta,Gamma, sim, sample.size, Time.Origin.STD,i, coal.events.times.simA)
{
  list.number.ancestors.population.sim <<-array(0,dim=c(sim,sample.size))
  Time.Origin.STD.sim <- array(0,dim=c(sim))
  lapply(1:sim, FUN=function(j, sample.size, Gamma, Delta, coal.events.times.simA,  Time.Origin.STD.sim)
  {
    list.number.ancestors.population= simulate.list.number.ancestors.population(sample.size)

    list.number.ancestors.population[list.number.ancestors.population==0]<-1
    list.number.ancestors.population.sim[j,] <<-list.number.ancestors.population
    list.coal.times= simulate.coalescent.times.A.test(Gamma, Delta,sample.size, list.number.ancestors.population)

    coal.events.times =unlist(list.coal.times)

    coal.events.times.simA[j,] <<- coal.events.times 
    print(paste0("finished A test sim ",j, sep=" "))
    
  },
  sample.size, Gamma, Delta, coal.events.times.simA,
  Time.Origin.STD.sim=Time.Origin.STD.sim)
  

  quants <- c(0.025,0.50,0.975)
  #quantiles<-apply( coal.events.times.simA , 2 , quantile , probs = quants , na.rm = TRUE )
  quantiles<-apply( coal.events.times.simA , 2 , quantile , probs = quants )
  meanCoalTimes<-colSums(coal.events.times.simA)/sim
  result <-data.frame(mean=meanCoalTimes, LI= quantiles[1,], median=quantiles[2,], UI=quantiles[3,])
  return(result)
}
simulateA.test.parallel=function(DeltaList,GammaList, sim, sample.size, Time.Origin.STD.test,coal.events.times.sim.testA)
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
  set.seed(69012365) #set seed to something
  s <- .Random.seed
  
  clust <- parallel::makeCluster(parallel::detectCores()-1, type="FORK", outfile="")

  on.exit(parallel::stopCluster(clust), add = TRUE)
  doParallel::registerDoParallel(clust)
  rng<- RNGseq( length(DeltaList)* sim, 234567)

  opts <- list(chunkSize=2)


  results <-foreach::foreach(j=1:sim, .combine = 'c') %:%
    foreach::foreach(i=1:length(DeltaList), r=rng[(j-1)*length(DeltaList) + 1:length(DeltaList)], .combine = 'c'
            
            ) %dopar% {

      rngtools::setRNG(r)

      list.number.ancestors.population= simulate.list.number.ancestors.population(sample.size)
      list.coal.times= simulate.coalescent.times.A.test(GammaList[i], DeltaList[i],sample.size, list.number.ancestors.population)
      coal.events.times =unlist(list.coal.times)
      Time.Origin= coal.events.times[length(coal.events.times)]
      Time.Origin.STD.test[i,j]<- Time.Origin
      positions <-cbind(rep(i, (sample.size-1)),((j-1)*(sample.size-1) +1):(j*(sample.size-1)))
      coal.events.times.sim.testA[positions] <- coal.events.times[1:(length(coal.events.times)-1)]
      print(paste("finished A test for Delta", DeltaList[i], "Gamma", GammaList[i], " sim ",j, sep=" "))
      NULL
    }

  list.Data.FramesA<-parallel::mclapply(1:length(DeltaList), mc.set.seed = TRUE, mc.cores=parallel::detectCores()-1,FUN=function(k, sample.size,sim, coal.events.times.sim.testA){
    
    quants <- c(0.025,0.50,0.975)
    matrixCurrentValue<-matrix(coal.events.times.sim.testA[k,], nrow=sim , ncol=sample.size-1, byrow=TRUE)
    #quantiles<-apply( matrixCurrentValue , 2 , quantile , probs = quants , na.rm = TRUE )
    quantiles<-apply( matrixCurrentValue , 2 , quantile , probs = quants  )
    meanCoalTimes<-colMeans(matrixCurrentValue)
    result <-data.frame(mean=meanCoalTimes, LI= quantiles[1,], median=quantiles[2,], UI=quantiles[3,])
    result
  },sample.size,sim, coal.events.times.sim.testA )
  
  return(list.Data.FramesA)
}
###################################################################################
#C scenario functions
Integrand <-function(u,Delta,Gamma){
  integrand<-function(x){
    term <- exp(-1.0 * Delta *x)

    relative(x,Delta,u)*Gamma * Delta * term * exp(-1.0 * Gamma * term /(1-term)) / (1-term)^2
  }
  return(integrand)
}
h_integral <-function(u,Delta,Gamma){
  integrand<-function(x, u,Delta,Gamma) 
    {
    term <- exp(-1.0 * Delta *x)

    return(relative(x,Delta,u)*Gamma * Delta * term * exp(-1.0 * Gamma * term /(1-term)) / (1-term)^2)
  }
    return(integrate(integrand, u=u, Delta=Delta, Gamma=Gamma, lower = u, upper = Inf, stop.on.error = FALSE)$value)
}
vh_integral<-function(u,Delta,Gamma){
  integrand<-function(x, u,Delta,Gamma) 
  {
    term <- exp(-1.0 * Delta *x)
    
    return(relative(x,Delta,u)*Gamma * Delta * term * exp(-1.0 * Gamma * term /(1-term)) / (1-term)^2)
  }
  return(integrate(integrand, u=u, Delta=Delta, Gamma=Gamma, lower = u, upper = Inf, stop.on.error = FALSE)$value)
}
hFinal_function <-function(Delta,Gamma){

  res<-function(x)
    {
    return(h_integral(x,Delta,Gamma))
  }
  return(res)
}
h_analytic<-function(Delta,Gamma){
  res<-function(x){
   exp_minus_delta_u<- exp(-1.0 * Delta *x)

   v<- exp_minus_delta_u / (1 - exp_minus_delta_u )
   exp_minus_gamma_v<- exp(-1.0*Gamma*v)
   
   result<- -1.0*(1.0/exp_minus_delta_u)*((1.0-exp_minus_delta_u)^2)*(exp_minus_gamma_v*(v^2 + 2.0*v/ Gamma + 2.0 / (Gamma^2)) -2.0/(Gamma^2))
   result<-result + 2.0*(1-exp_minus_delta_u)*(exp_minus_gamma_v*(v+ 1.0/Gamma )-1.0/Gamma)
   result<-result + exp_minus_delta_u*(1-exp_minus_gamma_v)
   return(result)
  }
  return(res)
}
lambdaFinal<-function(u,Delta,Gamma){
  factor=2.0 *Delta/Gamma
  return( factor * (1.0/h_integral(u, Delta, Gamma)))
}
lambda_analytic<-function(u,Delta,Gamma){
  factor=2.0 *Delta/Gamma
  funcX<-h_analytic(Delta, Gamma)
  return( factor * (1.0/funcX(u)))
}
LambdaFinal<-function(s,Delta,Gamma){
  if (s>0){
       return(integrate(Vectorize(lambda_analytic, c("u")),Delta=Delta, Gamma=Gamma,  lower = 0, upper = s,stop.on.error = FALSE)$value)

  }
  else{
    return(0)
  }
}
TransformedLambda<-function(z,Delta,Gamma){
  if (Gamma<100)
    return(integrate(Vectorize(lambdaFinal, c("u")),Delta=Delta, Gamma=Gamma,  lower = 0, upper = exp(z),stop.on.error = FALSE)$value)
  else{
    return(integrate(Vectorize(lambda_analytic, c("u")),Delta=Delta, Gamma=Gamma,  lower = 0, upper =  exp(z),stop.on.error = FALSE)$value)
  }
   
}
Lambda_integral<-function(s,Delta,Gamma){
  f<-function(x)
  {
    return(integrate(Vectorize(lambdaFinal, c("u")),Delta=Delta, Gamma=Gamma,  lower = 0, upper = s, stop.on.error = FALSE)$value)
  } 
  return(f)
}
log_functionToZero<-function( y, Delta,Gamma)
{
  f<-function(x) {
       return(log(LambdaFinal(x, Delta,Gamma))-log(y))
  }
  return(f)
}
derivativeLogLambda<-function(Delta, Gamma)
{
  f<-function(x)
  {
    return(lambdaFinal(x, Delta, Gamma) / LambdaFinal(x, Delta, Gamma))
  } 
  return(f)
}
functionToZero<-function( y, Delta,Gamma)
{
  f<-function(x) LambdaFinal(x, Delta,Gamma)- y
  return(f)
}
TransformedfunctionToZero<-function( y, Delta,Gamma)
{
  f<-function(x) {
    w<-exp(x)
    return(LambdaFinal(exp(x), Delta,Gamma)- y)
    }
  return(f)
}
derivativeTransformedLambda<-function(Delta, Gamma)
{
  f<-function(x)
  {
    return(lambda_analytic(exp(x), Delta, Gamma) *exp(x))
  } 
  return(f)
}
#####################################################################
# InverseTransformedLambda
#Computes Lambda^(-1)(s) for scenario scenario C by transforming the
#variables Lambda^(-1)(exp(s)) to allow the argument of Lambda^(-1)
#to be non negative
#####################################################################
InverseTransformedLambda<-function(s, Delta,Gamma)
{
  upper<-search_infinite_from(0, 0, functionToZero( s, Delta,Gamma))
  result<-cmna::bisection(functionToZero( s, Delta,Gamma), min(0,upper-1),upper, tol=10^(-6))
  return(result)
}
standard2modelC <- function(u, Delta, Gamma)
{
  result=InverseTransformedLambda(u, Delta, Gamma)
  return(result)

} 
modelCoalC <- function(sample.size, Delta, Gamma) {
  w <- standardCoal(sample.size)
  cumulativeSTDCoaltimes<-unlist(lapply(1:(sample.size-1),function(i) sum(w[i:(sample.size-1)])))
  u <- unlist(lapply( cumulativeSTDCoaltimes, FUN=function(x) standard2modelC(x,Delta, Gamma)))
  return(u)
}
#####################################################################
#sampleTOriginFormEVD8
#sample random time of origin for cancer model scenario B from Extreme Value Distribution(EVD)
#in particular Gumbel distribution
#####################################################################
sampleTOriginFormEVD8<-function(n, Delta, Gamma )
{
  require(evd)
  results<-list()
  sampled.gumbel <- evd::rgumbel(n, 0,1)
  lapply(1:n, FUN=function(i, Delta, Gamma, results){
    x=sampled.gumbel[i]
    y <- (1.0 / Delta) * log(1 + Gamma*exp(x) )
    results[[i]]<<-y
  }, Delta=Delta, Gamma=Gamma, results)
  return(unlist(results))
}
#####################################################################
#doParallel
#performs fun in parallel with progress bar 
#####################################################################
doParallel <- function(pop, fun, ncores = parallel::detectCores(),...) {
  require("doFuture")
  require(doSNOW)
  require(tcltk)
  registerDoFuture()
  ntasks <- 100
  pb <- tkProgressBar(max=ntasks)
  progress <- function(n) setTkProgressBar(pb, n)
  progress <- function(n) cat(sprintf("current Delta, Gamma %d is complete\n", n))
  opts <- list(progress=progress)
  cl <- parallel::makeCluster(ncores, outfile="")
  registerDoSNOW(cl)
  old_plan <- plan(cluster, workers = cl)
  on.exit({
    plan(old_plan)
    parallel::stopCluster(cl)
  })
  foreach(i = pop,
          .multicombine=TRUE) %dopar% fun(i,...)
}
#####################################################################
#nearest_value_search
#finds the first element p in the increasing ordered list increasing_list such that
# monotone_fun(p)>x 
#x value to compare the monotone_fun(tipically 0) 
#increasing_list: increasing ordered list
#monotone_function: monotone function(increasing or decreasing)
#####################################################################
nearest_value_search = function(x, increasing_list, monotone_function){
  left = 1
  right = length(increasing_list)
  while(right - left > 1){
    middle = floor((left + right) / 2)
    if(x < monotone_function(increasing_list[middle])){
      right = middle
    }
    else{
      left = middle
    }
  }
    return(increasing_list[right])
  
}
simulateB_all<-function(sample.size, sim, DeltaList, GammaList, Time.Origin.STD){

  listDataFramesB<-lapply(1:length(DeltaList), FUN=function(i, DeltaList,GammaList, sim, sample.size, Time.Origin.STD)
                       {
                        Delta = DeltaList[i]
                        Gamma = GammaList[i]
                        return(simulateB_K(i, Delta, Gamma,  sim, sample.size, Time.Origin.STD))
                        },DeltaList, GammaList, sim=sim, sample.size=sample.size, Time.Origin.STD, 2)
  return(listDataFramesB)
}
#####################################################################
#search_infinite_from
#finds the first interval [2^{p-1}*start, 2^{p}*start] with p=1,2,.., in [start, infinyty] 
# such that monotone_fun(2^{p}*start)>x and  monotone_fun(2^{p-1}*start)<=x
#x value to compare the monotone_fun(tipically 0) 
#start: start of the interval to search
#monotone_function: monotone function(increasing or decreasing)
#####################################################################
search_infinite_from<-function( x, start,  monotone_function ){ 
    if (start >0){
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
       to = 2*to   
    }
   increasing_list<-seq(from,to, by=1)
   result<-nearest_value_search(x, increasing_list, monotone_function)
   return(result)
}
#####################################################################
#get_simulated_kth_coal_times
#gets the simulated  kth coalescence times
#coal.events.times.simA is a big matrix bigstatsr::FBM(length(DeltaList),sim*(sample.size-1))
#k index of coalescence time  from present time backwards
#sim: number of simulations
#####################################################################
get_simulated_kth_coal_times<-function(coal.events.times.sim, pos, sim, sample.size)
{
  num<-length(rows_along(coal.events.times.sim))
  kth_coal_times<-lapply(1:num, FUN=function(k, sample.size,sim, coal.events.times.sim, pos ){
    
    matrixCurrentValue<-matrix(coal.events.times.sim[k,], nrow=sim , ncol=sample.size-1, byrow=TRUE)
    kth_times<- matrixCurrentValue[,pos]
    
    result <-data.frame(coal_times=kth_times)
    colnames(result)<-paste("coal_times_",pos, sep="")
    return(result)
    
  },sample.size,sim, coal.events.times.sim, pos )
  
  return(kth_coal_times)
}
#####################################################################
#get_list_simulated_trees
#gets a simulate tree with some  set of paramater value
#coal.events.times.sim is a big matrix bigstatsr::FBM(length(DeltaList),sim*(sample.size-1))
#parameter_index index of parameter value
#number.sim.trees: number of simulations
#sample.size: the sample size used
#####################################################################
get_list_simulated_trees<-function(coal.events.times.sim, parameter_index, number.sim.trees, sample.size)
{
   require(dplyr)
   require(ape)
  
    matrixCurrentValue<-matrix(coal.events.times.sim[parameter_index,], nrow=number.sim.trees , ncol=sample.size-1, byrow=TRUE)
  
    list_trees<-lapply(1:number.sim.trees, function(k,sample.size, matrixCurrentValue ){
      
      random_simulation<- unlist(matrixCurrentValue[k,])
      coalescent_tree = rcoal(sample.size, br="coalescent")
      x <- tibble::as_tibble(coalescent_tree)
      
      df<-as.data.frame(x)
      
      node_times<- c(rep(0,sample.size),rev(random_simulation))
      
      x[,"branch.length"]<-unlist(lapply(1:nrow(x), function(idx,node_times,x){
        if (as.numeric(x[idx, "parent"]) != as.numeric(x[idx, "node"]))
          return(node_times[as.numeric(x[idx, "parent"])]-node_times[as.numeric(x[idx, "node"])])
        else
          return(NA)
        
      }, node_times,x))
      
      tree<-ape::as.phylo(x)
      
      return(tree)
      
    }, sample.size, matrixCurrentValue)


  return(list_trees)
}
install_required_packages=function(list.of.packages){
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0) install.packages(new.packages, dependencies=TRUE, repos = "http://cran.us.r-project.org")
  
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  if ("ggtree"  %in%  new.packages )
     BiocManager::install("ggtree")
  
}
load_required_packages=function(list_packages){

  lapply(list_packages,
                   library,
                   character.only = TRUE)

  
}
getCurrentFileLocation <-  function()
{
  require(tibble)
  this_file <- commandArgs() %>% 
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
  if (length(this_file)==0)
  {
    this_file <- rstudioapi::getSourceEditorContext()$path
  }
  return(dirname(this_file))
}
computeCorrelationMatrix  <-  function(coal.events.times.simA, coal.events.times.simB, DeltaList, sample.size, sim, path_to_save)
{

  require(lineup)
  list.corr.vect<-parallel::mclapply(1:length(DeltaList),mc.cores=parallel::detectCores()-1, 
                                        FUN=function(k, sample.size,sim, coal.events.times.simB, coal.events.times.simA )
                                        {
                                          
                                          matrixCurrentValueA<-matrix(coal.events.times.simA[k,], nrow=sim , ncol=sample.size-1, byrow=TRUE)
                                          matrixCurrentValueB<-matrix(coal.events.times.simB[k,], nrow=sim , ncol=sample.size-1, byrow=TRUE)
                                          
                                          result <- lineup::corbetw2mat(matrixCurrentValueA, matrixCurrentValueB, what = "paired")
                                          return(list(as.vector(result)))
                                          
                                        },sample.size,sim, coal.events.times.simB, coal.events.times.simA )
  
  return(list.corr.vect)

}
computeCorrelationMatrixModel <- function(coal.events.times.sim, DeltaList, sample.size, sim, path_to_save)
{
  
  require(lineup)
  # list_corr_mat<-parallel::mclapply(1:length(DeltaList),mc.cores=parallel::detectCores()-1, 
  #                                    FUN=function(k, sample.size,sim, coal.events.times.sim )
  #                                    {
  #                                      
  #                                      matrixCurrentValue<-matrix(coal.events.times.sim[k,], nrow=sim , ncol=sample.size-1, byrow=TRUE)
  #                       
  #                                      result <- cor(matrixCurrentValue)
  #                                      return(as.matrix(result))
  #                                      
  #                                    },sample.size,sim, coal.events.times.sim )
  list_corr_mat<-lapply(1:length(DeltaList), 
                                    FUN=function(k, sample.size,sim, coal.events.times.sim )
                                    {
                                      
                                      matrixCurrentValue<-matrix(coal.events.times.sim[k,], nrow=sim , ncol=sample.size-1, byrow=TRUE)
                                      
                                      result <- cor(matrixCurrentValue)
                                      return(as.matrix(result))
                                      
                                    },sample.size,sim, coal.events.times.sim )
  return(list_corr_mat)
  
}
#############################################################################################
sim = 10000
sample.size = 10
runScenarioA = TRUE
runScenarioC = FALSE
runScenarioB = TRUE
doTestA = FALSE
varyDelta = FALSE
doJustBenchMark = FALSE
runHybridScenario = TRUE
number.sim.trees=5
simulate.Trees = FALSE
K=1

list_packages<-c("parallel", "doFuture", "bigstatsr", "doRNG", "microbenchmark", "matrixStats", "cmna", "NLRoot", "ggplot2", "future.apply",
                 "ape",  "tibble",  "ggtree", "geiger", "dplyr", "pracma", "foreach", "doParallel",  "stats",
                 "doSNOW", "tcltk", "ggplot2", "tidyverse", "pryr", "lineup")


install_required_packages(list_packages)
lapply(list_packages,
       library,
       character.only = TRUE)



path_to_save =getCurrentFileLocation()

GammaList=c(0.001,  0.1,  10)


DeltaList= rep(1, length(GammaList))

Time.Origin.STD <-bigstatsr::FBM(length(DeltaList), sim , type="double", init=0)
coal.events.times.simA=bigstatsr::FBM(length(DeltaList),sim*(sample.size-1))

number.ancestors.simA=bigstatsr::FBM(length(DeltaList),sim*(sample.size-1))

Time.Origin.STD.trees <-bigstatsr::FBM(length(DeltaList), number.sim.trees , type="double", init=0)
coal.events.times.simA.trees=bigstatsr::FBM(length(DeltaList),number.sim.trees*(sample.size-1))

coal.events.times.simB=bigstatsr::FBM(length(DeltaList),sim*(sample.size-1))

coal.events.times.simHybrid=bigstatsr::FBM(length(DeltaList),sim*(sample.size-1))

nB= floor(sample.size / 2)


mprimeList<-c(99,95,90,80,70,50)

number.Fail.hybrid=bigstatsr::FBM(length(mprimeList),length(GammaList))


if (doTestA)
{
  Time.Origin.STD.test <-bigstatsr::FBM(length(DeltaList), sim , type="double", init=0)
  coal.events.times.sim.testA=bigstatsr::FBM(length(DeltaList),sim*(sample.size-1))
 
  listDataFramesA_test<-simulateA.test.parallel(DeltaList,GammaList, sim, sample.size, Time.Origin.STD.test,coal.events.times.sim.testA)
  saveRDS(listDataFramesA_test, paste(path_to_save,"/","listDataFramesA_",sample.size,"_test.rds", sep=""))
  meansTorigin.test <- bigstatsr::big_apply(Time.Origin.STD.test, a.FUN = function(Time.Origin.STD.test, ind) rowMeans(Time.Origin.STD.test[ind, ]),
                            ind = rows_along(Time.Origin.STD.test), a.combine = 'c',
                            block.size = 100)
  
  saveRDS(meansTorigin.test, file=paste(path_to_save,"/", "meansToriginA_",sample.size,"_test.rds", sep=""))

  
  pryr::mem_change(rm(Time.Origin.STD.test))
  pryr::mem_change(rm(coal.events.times.sim.testA))

}
if (!doJustBenchMark){
######################################################################################
# Scenario A with stochastic population size

if (runScenarioA){
  number.ancestors.Transition=bigstatsr::FBM(length(DeltaList),sim)
  
  listDataFramesA<-simulateA.parallel(DeltaList, GammaList, sim, sample.size, Time.Origin.STD, coal.events.times.simA, number.ancestors.simA, number.ancestors.Transition)
  saveRDS(listDataFramesA, file=paste(path_to_save,"/", "listDataFramesA_",sample.size,"_", sim, ".rds", sep=""))
  
  number.ancestors.Transition.Matrix<-as.matrix(number.ancestors.Transition[])
  saveRDS(number.ancestors.Transition.Matrix, file=paste(path_to_save,"/","number.ancestors.Transition_", sample.size,"_", sim,".rds", sep=""))
  
  
  
  listDataFramesAncestors <- compute_stats_number_ancestors(sample.size,sim,DeltaList, number.ancestors.simA)
  saveRDS(listDataFramesAncestors, file=paste(path_to_save,"/", "listDataFramesA_NumberAncestors",sample.size,"_", sim,".rds", sep=""))
  
  saveRDS(Time.Origin.STD, file=paste(path_to_save,"/", "Time.Origin.STD_",sample.size,"_", sim,".rds", sep=""))
  
  meansTorigin <- bigstatsr::big_apply(Time.Origin.STD, a.FUN = function(Time.Origin.STD, ind) rowMeans(Time.Origin.STD[ind, ]),
                                       ind = rows_along(Time.Origin.STD), a.combine = 'c',
                                       block.size = 500)
  
  saveRDS(meansTorigin, file=paste(path_to_save,"/", "meansToriginA_",sample.size,"_", sim,".rds", sep=""))
  
  listToriginCI <- bigstatsr::big_apply(Time.Origin.STD, a.FUN = function(Time.Origin.STD, ind){
    quants <- c(0.025,0.50,0.975)
    matrixStats::rowQuantiles(Time.Origin.STD[ind, ],  probs = quants)
  } ,
  ind = rows_along(Time.Origin.STD),a.combine = "rbind",
  block.size = 500)
  saveRDS(listToriginCI, file=paste(path_to_save,"/", "listToriginCI_A_",sample.size,"_", sim,".rds", sep=""))
  
  pos=sample.size-1
  simul_coal_times_A<-get_simulated_kth_coal_times(coal.events.times.simA, pos, sim, sample.size)
  
  saveRDS(simul_coal_times_A, file=paste(path_to_save,"/", "simul_coal_times_A_",pos,"_",sample.size,"_", sim,".rds", sep=""))
  
  list_corr_mat <- computeCorrelationMatrixModel(coal.events.times.simA, DeltaList, sample.size, sim, path_to_save)
  
  pryr::mem_change(rm(coal.events.times.simA))
  pryr::mem_change(rm(number.ancestors.simA))
  pryr::mem_change(rm(number.ancestors.Transition))
  
  saveRDS(list_corr_mat, file=paste(path_to_save,"/", "list_corr_mat_A","_",sample.size,"_", sim,".rds", sep=""))
}
else{
  print(paste(path_to_save,"/", "listDataFramesA_",sample.size,".rds", sep=""))
  listDataFramesA<-  readRDS(paste(path_to_save,"/", "listDataFramesA_",sample.size,"_", sim,".rds", sep=""))
  
  print(paste(path_to_save,"/", "Time.Origin.STD_",sample.size,"_",nB,".rds", sep=""))
  Time.Origin.STD <- readRDS(paste(path_to_save,"/", "Time.Origin.STD_",sample.size,"_", sim,".rds", sep=""))
  
  print(paste(path_to_save,"/", "meansToriginA_",sample.size,".rds", sep=""))
  meansTorigin <- readRDS(paste(path_to_save,"/", "meansToriginA_",sample.size,"_", sim,".rds", sep=""))
}


#######################################################################################
#hybrid scenario

if (runHybridScenario){
  
  coal.events.times.simHybrid=bigstatsr::FBM(length(DeltaList),sim*(sample.size-1))
  number.ancestors.simHybrid=bigstatsr::FBM(length(DeltaList),sim*(sample.size-1))
  number.ancestors.Transition=bigstatsr::FBM(length(DeltaList),sim)
  
K=2

listDataFramesHybrid <- simulateHybrid.parallel(DeltaList,GammaList, sim, sample.size, Time.Origin.STD, coal.events.times.simHybrid, nB, number.Fail.hybrid, 1, K, number.ancestors.simHybrid, number.ancestors.Transition)

saveRDS(listDataFramesHybrid, file=paste(path_to_save,"/","listDataFramesHybrid_","K=", K, "_",sample.size,"_",nB,"_", sim,".rds", sep=""))

number.ancestors.Transition.Matrix<-as.matrix(number.ancestors.Transition[])
saveRDS(number.ancestors.Transition.Matrix, file=paste(path_to_save,"/","number.ancestors.Transition_","K=", K, "_",sample.size,"_",nB,"_", sim,".rds", sep=""))


print(number.Fail.hybrid[])

listDataFramesAncestors <- compute_stats_number_ancestors(sample.size,sim,DeltaList, number.ancestors.simHybrid)
saveRDS(listDataFramesAncestors, file=paste(path_to_save,"/", "listDataFramesHybrid_NumberAncestors",sample.size,"_", sim,".rds", sep=""))

list_corr_mat <- computeCorrelationMatrixModel(coal.events.times.simHybrid, DeltaList, sample.size, sim, path_to_save)
  
saveRDS(list_corr_mat, file=paste(path_to_save,"/", "list_corr_mat_Hybrid","_","K=", K,"_",sample.size,"_",nB,"_", sim,".rds", sep=""))

pryr::mem_change(rm(coal.events.times.simHybrid))
pryr::mem_change(rm(number.ancestors.simHybrid))
pryr::mem_change(rm(number.ancestors.Transition))
}

######################################################################################
#Scenario B(cancer model)
if (runScenarioB){
  
  coal.events.times.simB=bigstatsr::FBM(length(DeltaList),sim*(sample.size-1))
  print("starting simulation M*")
  #K.list= seq(0.5,2, by=0.1)
  #K.list= seq(2,2, by=1)
  K.list= c(0.8)
  
  listDataFramesB<-simulateB_MAster.parallel(DeltaList, GammaList, sim, sample.size, Time.Origin.STD, coal.events.times.simB)
  saveRDS(listDataFramesB, file=paste(path_to_save,"/","listDataFramesB_",sample.size,"_", sim,".rds", sep=""))
  
  print("starting simulation M_k")
  for(i in  1:length(K.list)){
   
     K=K.list[i]
    listDataFramesB<-simulateB_K.parallel(DeltaList, GammaList, sim, sample.size, Time.Origin.STD, coal.events.times.simB, K)
    saveRDS(listDataFramesB, file=paste(path_to_save,"/","listDataFramesB_",sample.size,"_K=",K, "_", sim,".rds", sep=""))
    
    
    list_corr_mat <- computeCorrelationMatrixModel(coal.events.times.simB, DeltaList, sample.size, sim, path_to_save)
    saveRDS(list_corr_mat, file=paste(path_to_save,"/", "list_corr_mat_B_","_",sample.size,"_K=",K,"_", sim, ".rds", sep=""))
    
  #  if (i < length(K.list)){
   #   listDataFramesA<-simulateA.parallel(DeltaList, GammaList, sim, sample.size, Time.Origin.STD, coal.events.times.simA)
   # }
    
  }
  
  pos=sample.size-1
  simul_coal_times_B<-get_simulated_kth_coal_times(coal.events.times.simB, pos, sim, sample.size)

  saveRDS(simul_coal_times_B, file=paste(path_to_save,"/", "simul_coal_times_B_",sample.size,".rds", sep=""))
  pryr::mem_change(rm(coal.events.times.simB))

  
}

######################################################################################
#Scenario C

if (runScenarioC)
{

  
  coal.events.times.simC=bigstatsr::FBM(length(DeltaList),sim*(sample.size-1))
  
  listDataFramesC<-simulateC.parallel(DeltaList, GammaList, sim, sample.size,  coal.events.times.simC)
  saveRDS(listDataFramesC, file=paste(path_to_save,"/listDataFramesC_",sample.size,"_",nB,".rds", sep=""))

  pryr::mem_change(rm(coal.events.times.simC))
 
}

##############################################################################################
#Draw some trees: number.sim.trees trees per Delta and Gamma combination

if (simulate.Trees){
    
    packageVersion("phangorn")
    packageVersion("phytools")
    
    print("Simulating some example trees...")
    
    listDataFramesA.trees<-simulateA.parallel(DeltaList, GammaList, number.sim.trees, sample.size, Time.Origin.STD.trees, coal.events.times.simA.trees)
    
    lapply(1:length(GammaList), function(param_index, coal.events.times.simA.trees, number.sim.trees, sample.size, path_to_save){
      
      random_simulated_trees<-get_list_simulated_trees(coal.events.times.simA.trees, param_index,  number.sim.trees, sample.size)
      
      lapply(1:number.sim.trees, function(idx, param_index, GammaList, path_to_save, random_simulated_trees ){
        
        tree=random_simulated_trees[[idx]]
        tree<- ape::as.phylo(tree)
        tiff(paste(path_to_save, "/",  "tree_n=",sample.size,"_Gamma=",GammaList[param_index],"_",idx,".tiff", sep=""), units="in", width=3.25, height=3.25, res=300)
        plot(tree, show.tip.label = FALSE)
        dev.off()
      }, param_index, GammaList, path_to_save, random_simulated_trees )
      
      
    }, coal.events.times.simA.trees, number.sim.trees, sample.size, path_to_save)
    
    pryr::mem_change(rm(coal.events.times.simA.trees))
    pryr::mem_change(rm(Time.Origin.STD.trees))
    
  }


}else{
  
#################################################################################################
#Benchmarking
#uncomment this if you want to benchmark the code
print("Starting simulation benchmark...")
  
listDataFramesA<-simulateA.parallel(DeltaList, GammaList, sim, sample.size, Time.Origin.STD, coal.events.times.simA)
K=1
mbm <- microbenchmark::microbenchmark("BD"=listDataFramesA<-simulateA.parallel(DeltaList, GammaList, sim, sample.size, Time.Origin.STD, coal.events.times.simA),
                         "M*"= listDataFramesB<-simulateB_MAster.parallel(DeltaList, GammaList, sim, sample.size, Time.Origin.STD, coal.events.times.simB),
                         "M0.8"= listDataFramesB<-simulateB_K.parallel(DeltaList, GammaList, sim, sample.size, Time.Origin.STD, coal.events.times.simB, 0.8),
                         "M2"=listDataFramesB<-simulateB_K.parallel(DeltaList, GammaList, sim, sample.size, Time.Origin.STD, coal.events.times.simB, 2),
                         "H50"= listDataFramesHybrid<-simulateHybrid.parallel(DeltaList,GammaList, sim, sample.size, Time.Origin.STD, coal.events.times.simHybrid, 50, number.Fail.hybrid, 6, 2, number.ancestors.simHybrid, number.ancestors.Transition),
                         times=100L)
  print(mbm)

  
  tiff(paste(path_to_save, "/",  "PlotBenchmark",".tiff", sep=""), units="in", width=3.25, height=3.25, res=300)
  plot_mbm<-ggplot2::autoplot(mbm, log=TRUE)+  theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  plot_mbm
  ggsave(paste(path_to_save, "/",  "plot_mbm_all.pdf", sep=""))
  ggsave(paste(path_to_save, "/", "plot_mbm_all.png", sep=""))
  ggsave(paste(path_to_save, "/", "plot_mbm_all.jpg", sep=""))
  
  dev.off()
  
  total.failed <- bigstatsr::big_apply(number.Fail.hybrid, a.FUN = function(X, ind) {
    rowSums(X[, ind, drop = FALSE])
  }, a.combine = 'plus', block.size = 10)
  
  print(total.failed)
  print(number.Fail.hybrid[])
  
  
}





