data {
          int<lower=0> m1; //number of active patients at interim from first stage
          int<lower=0> n1event; //number of patients with event out of n1-m1
          real mu; //prior mu for theta
          real<lower=0> sigma;//prior sd for theta
          real<lower=0> a; //prior a for kappa
          real<lower=0> b; //prior b for kappa
          vector[n1event] y_looktime; //time till event at looktime
          vector[m1] z_looktime; //censored b/c no event until looktime

          }
          parameters {
          real<lower=0> theta; //scale parameter to estimate
          real<lower=0> kappa; //shape parameter to estimate
          }

          model{
          target += weibull_lpdf(y_looktime | kappa, theta);
          target += weibull_lccdf(z_looktime | kappa, theta);

          theta ~ lognormal(mu, sigma); //prior for theta, Reference Ibrahim
          kappa ~ gamma(a, b); //prior for kappa, Reference Ibrahim
          }
