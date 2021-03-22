functions{
  real conditionalDensityTOrigin_rng(real delta){
    //this distribution is related with Extreme Value(EV) distribution
        real result =  (1.0 / delta) * log(1- delta / (log(uniform_rng(0,1))));
        return(result);
  }
  real conditionalDensityTOrigin_pdf( real y, real delta) {
      real result;
      real term1;
      real term2;
      real term3;
      real partial;
      term1 = exp(-1.0*delta*y);
      term2 = delta * term1;
      term3 = 1.0-term1;
      partial = delta * term2 /(term3 * term3);
      partial = partial * exp( -1 * term2/term3);
      result  = partial;
      return(result);
   }
  real conditionalDensityTOrigin_lpdf( real y, real delta) {
      real result;
      real term1;
      real term2;
      real term3;
      real partial;
      term1 = exp(-1.0*delta*y);
      term2 = delta * term1;
      term3 = 1.0-term1;
      partial = delta * term2 /(term3 * term3);
      partial = partial * exp( -1 * term2/term3);
      result  = log(partial);
      return(result);
   }
   real conditionalDensityTOrigin_cdf(real y, real delta) {
     real term1;
     real term2;
     real term3;
     real result;
     term1 = exp(-1.0 * delta * y);
     term2 = delta * term1;
     term3 = 1.0-term1;
     result = exp( -1 * term2/term3);
     return(result);
  }
  real conditionalDensityTOrigin_lcdf(real y, real delta) {
    real term1;
    real term2;
    real term3;
    real result;
    term1 = exp(-1.0 * delta * y);
    term2 = delta * term1;
    term3 = 1.0-term1;
    result = exp( -1 * term2/term3);
    return(log(result));
 }
  real log_h(real t,
         real torigin,
         real delta,
         real K){

    real a = 1.0 - exp(-1.0 * delta * (torigin - t));
    real first_term = 2.0 * log(a);
    real second_term = -1.0 * delta * t;
    real third_term = exp(delta * t);
    real above_term = first_term + second_term;
    real b = 1.0 - exp(-1.0 * delta * torigin);
    real  below_term = 2.0 * log(b);
    real extra_term = log(1+ (K /delta)*b*(third_term -1.0) / a );
    real logH;
    if (K==0)
        logH = above_term - below_term;
    else
        logH =  above_term - below_term + extra_term;
    return(logH);
  }
  real model_time_to_standard_time(real t,
         real torigin,
         real delta,
         real K)
         {
         //if (t==0)
        //    return (0.0);

          real model_time=0.0;
          real a = exp(delta * t) - 1.0;
          real b = 1.0 - exp(-1.0 * delta * torigin);
          real c = 1.0 - exp(-1.0 * delta * (torigin - t));
          real d;
          real numerator;
          real  denominator;
          if (K==0){
            model_time = a * b / (delta * c);
            return(2*model_time);
          }
         else{
           c = exp(-1.0 * delta * torigin);
           d = 1.0 -c;
           a = ((K / delta) * d) - c;
           b = 1.0 - (K/ delta)*d;

           numerator =(a*exp(delta*t)+b);
           denominator = 1-exp(-delta*(torigin-t));
           model_time = (2.0 / K )*log(numerator/denominator);
           return(model_time);
            }
         }
  real log_likelihood(int sample_size, real K, real delta,
       real torigin, vector sorted_coalescent_times_scaled_by_theta)
    {

    vector[sample_size] sortedCoalescentTimes;
    real current_time;
    real log_lik = 0.0;
    int alive_cells;
    real term_only_after_first_coal_event;

    //sortedCoalescentTimes = (1.0/theta) * sorted_coalescent_times_scaled_by_theta;
    sortedCoalescentTimes =  sorted_coalescent_times_scaled_by_theta;

    for(j in 1:(sample_size-1)){
      current_time = sortedCoalescentTimes[j];
      alive_cells =(sample_size-j+1);//alive_cells goes from N downto 2
      log_lik = log_lik + log(alive_cells*(alive_cells-1)/2.0);

      log_lik = log_lik + log(2)-log_h(current_time,torigin, delta, K);

      if (alive_cells==sample_size)
         term_only_after_first_coal_event=0.0;
      else
         term_only_after_first_coal_event= model_time_to_standard_time(sortedCoalescentTimes[j-1],torigin, delta, K);

      log_lik = log_lik - (alive_cells /2.0)*(alive_cells-1)*(model_time_to_standard_time(current_time,torigin, delta, K)-
                term_only_after_first_coal_event);
    }
    return(log_lik);
  }
}
data{
  real<lower=0.0, upper=2.0> K;
  int<lower=2> sample_size;
  vector<lower = 0.0>[sample_size] sorted_coalescent_times_scaled_by_theta;
}
parameters{
  //real<lower=0> theta;
  real<lower=0> delta;
  real<lower=0> torigin;
}
model{
  //theta~exponential(10);
  //theta=1
  delta~exponential(10);
  torigin~conditionalDensityTOrigin(delta);

  //target += log(parameter_exponential_torigin) - torigin * parameter_exponential_torigin;
  //target +=log_likelihood(sample_size, K, theta, delta, torigin, sorted_coalescent_times_scaled_by_theta);
  target +=log_likelihood(sample_size, K, delta, torigin, sorted_coalescent_times_scaled_by_theta);
}
