data {
          int<lower=0> nevent; // number of events at study end
          int<lower=0> ncens; // number of aC at study end
          real mu; //prior mu for theta
          real<lower=0> sigma;//prior sd for theta
          real<lower=0> a; //prior a for kappa
          real<lower=0> b; //prior b for kappa
          vector[nevent] y_studyend; //time till event at looktime
          vector[ncens] z_studyend; //aC time

          }
          parameters {
          real<lower=0> theta; //scale parameter to estimate
          real<lower=0> kappa; //shape parameter to estimate
          }

          model{
          target += weibull_lpdf(y_studyend | kappa, theta);
          target += weibull_lccdf(z_studyend | kappa, theta);

          theta ~ lognormal(mu, sigma); //prior for theta, Reference Ibrahim
          kappa ~ gamma(a, b); //prior for kappa, Reference Ibrahim
          }
