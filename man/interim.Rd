\name{interim}
\alias{interim}
\alias{print.interim}
\title{
Interim Analysis for a single-arm two-stage phase II clinical trial based on a Bayesian Weibull model.
}
\description{
Conducts the interim analysis for a Bayesian Weibull model based on the Bayesian Probability of Success.
Taking the survival data observed so far, the function utilizes the function \link[rstanPoS]{data_looktime} to determine all the important variables for the stan model. Then it saves the input variables from the data block of the stan model in a list, to then run the stan model to estimate the posterior values of the shape and scale parameter of the Weibull Model. The prior for the shape is given by a Gamma distribution G(a,b), and the prior for the scale by a lognormal distribution lognormal(mu,sigma). Next it extracts all the the posterior values of the Weibull parameters to then predict the survival data of the second stage patients and remaining first stage patients based on one random posterior value of shape and scale. With the help of the function \link[rstanPoS]{data_studyend} the important values for the predicted survival data are determined and furthermore it combines those with the already observed survival data at interim. Next, the observed and predicted survival data is utilized to once again estimate the posterior distribution. With the newly extracted posterior values the survival probabilities at the primary endpoint are calculated, which are then used to determine the posterior probability of a clinical relevant outcome and from this the success criterion is verified. This second step, from the prediction on, is repeated for d times, a feasible value for d is 100. This will then return the Bayesian Probabilit of Success, as the relative frequency of successes out of d samples.
}
\usage{
interim(n.cores,rddata,discrete,data,n,n1,p0,p2_0,k,tau,omega1,omega2,kappa,theta,d,eta_star,lambda,mu,sigma,a,b,till,to,to2)
\method{print}{interim}(x,...)
}
\arguments{
  \item{n.cores}{number of cores on computer minus 1 - to run the code parallel on multiple cores.}
  \item{rddata}{a indicator variable - 1 if a random data set is simulated, 0 if an external data set is utilized.}
  \item{discrete}{set to 1 if discrete times are determined, 0 if continous times are determined}
  \item{data}{if rddata = 0, data is the external data set, a dataframe containing the observed accrual time variable e, the time-to-event variable y, the event time variable y_e
  and the event indicator variable delta at interim for the first stage patients*. Set this to any number if a simulated random data set is used.}
  \item{n}{the total number of patients in the study.}
  \item{n1}{the number of first stage patients.}
  \item{p0}{the survival probability of the corresponding primary endpoint.}
  \item{p2_0}{the survival probability of the corresponding intermediate endpoint.}
  \item{k}{the intermediate endpoint.}
  \item{tau}{the primary endpoint.}
  \item{omega1}{the accrual period for the first stage patients.}
  \item{omega2}{the accrual period for the second stage patients.}
  \item{kappa}{the assumed shape parameter of the Weibull distribution.}
  \item{theta}{the assumed scale parameter of the Weibull distribution.}
  \item{d}{the number of repetitions/simulations, use at least 1000.}
  \item{eta_star}{the threshold for determining success.}
  \item{lambda}{the threshold for making an interim decision.}
  \item{mu}{the prior mu for the lognormal distribution of the scale parameter.}
  \item{sigma}{the prior sigma for the lognormal distribution of the scale parameter.}
  \item{a}{the prior shape parameter for the gamma distribution of the shape parameter.}
  \item{b}{the prior scale parameter for the gamma distribution of the shape parameter.}
   \item{till}{if the observation is taken every 2 months until the 6th month, till = 6, if even after 6 months, the last observation is taken after 2 months, set till to any number.}
  \item{to}{the time steps at which the observations are taken e.g. 2,4,6 months so to = 2}
  \item{to2}{if the first observations are taken every two months but the last one, three months after the one before the last, set to2 = 3, otherwise to2 is equal to "to".}

}
\details{
*The time-to-event variable of the observed survival dataset must not include any NA or missing values for the patients without an event. Instead a high number, as presened in the example section shall be utilized.
}
\value{
Prints the following items with class "interim":
  \item{Results}{a data frame consisting of the assumed survival probability for the primary and intermediate endpoint, the assumed shape and scale parameters of the Weibull model, the Bayesian Probability of Success and the interim decision.}
  \item{Summary of Posterior Values at looktime}{a data frame displaying the median, 5 percent and 95 percent quantile of the posterior values at looktime for both shape and scale.}
   \item{Summary of Posterior Values at study end}{a data frame displaying the median, 5 percent and 95 percent quantile of the posterior values at studyend for both shape and scale. These is the summary of the posterior values based on the observed and predicted survival data.}
}

\author{
Simrit Rattan
}

\examples{
for a simulated random data set:
interim(n.cores = 7,rddata = 1,discrete = 1,data = 0,n = 43,n1 = 24,p0 = 0.127,p2_0 = 0.4,k = 4,tau = 9,omega1 = 24,omega2 = 19,kappa = 1.001144,theta = 4.36499,d = 1000,eta_star = 0.95,lambda = 0.10,mu = 0,sigma = 10,a=1,b = 0.001, till = 6, to = 2, to2 = 3)

for an external data set:
data = data.frame("e" = c(4,8,10,3.5,12.5,21), "y" = c(2,4,1000,1000,1000,1000), "y_e" = c(6,12,1010,1003.5,1012.5,1021),"delta" = c(1,1,0,0,0,0))
#y = 1000 indicates that no TTE data is available as no event has occured so far for this patient.
interim(n.cores = 7, rddata = 0, discrete = 1, data = data ,n = 43,n1 = 12,p0 = 0.127,p2_0 = 0.4,k = 4,tau = 9,omega1 = 24,omega2 = 19,kappa = 1.001144,theta = 4.36499,d = 1000,eta_star = 0.95,lambda = 0.10,mu = 0,sigma = 10,a=1,b = 0.001, till = 6, to = 2, to2 = 3)
}
