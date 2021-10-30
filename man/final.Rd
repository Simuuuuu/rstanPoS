\name{final}
\alias{final}
\alias{print.final}
\title{
Final Analysis for a single-arm two-stage phase II clinical trial based on a Bayesian Weibull model.
}
\description{
Conducts the final analysis of a Bayesian Weibull model, based on the posterior probability of a clinical relevant outcome.
As a first step, the full data set of all n patients will be utilized to get all the important variables for the Stan model via the function \link[rstanPoS]{data_studyend2}. Next, all the input factors stated in the data block of the Stan model will be saved into a list to estimate the posterior values of the Weibull parameters. The prior for the shape is given by a Gamma distribution G(a,b), and the prior for the scale by a lognormal distribution lognormal(mu,sigma). Then the survival probabilities at the primary endpoint for all posterior values of shape and scale are calculated. These are utilized to determine the posterior probability of a clinical relevant outcome. Last, it will be checked whether the success criterion is met to make a final decision.
}
\usage{
final(rddata,interim,data,n,n1,p0,p2_0,k,tau,omega2,kappa,theta,eta_star,mu,sigma,a,b,till,to,to2)
\method{print}{final}(x,...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rddata}{a indicator variable - 1 if a random data set is simulated, 0 if a external data set is utilized.}
  \item{interim}{if rdata = 1, interim is the output of the function \link[rstanPoS]{interim}, from which the random interim data set is needed. Set this to any number if an external data set is used.}
  \item{data}{if rddata = 0, data is the external data set, a dataframe containing the observed accrual time variable e, the              time-to-event variable y,and the event indicator variable delta at studyend for all n patients*.Set this to any number if a simulated random data set is used.}
  \item{n}{the total number of patients in the study.}
  \item{n1}{the number of first stage patients.}
  \item{p0}{the survival probability of the corresponding primary endpoint.}
  \item{p2_0}{the survival probability of the corresponding intermediate endpoint.}
  \item{k}{the intermediate endpoint.}
  \item{tau}{the primary endpoint.}
  \item{omega2}{the accrual period for the second stage patients.}
  \item{kappa}{the assumed shape parameter of the Weibull distribution.}
  \item{theta}{the assumed scale parameter of the Weibull distribution.}
  \item{eta_star}{the threshold for determining success.}
  \item{mu}{the prior mu for the lognormal distribution of the scale parameter.}
  \item{sigma}{the prior sigma for the lognormal distribution of the scale parameter.}
  \item{a}{the prior shape parameter for the gamma distribution of the shape parameter.}
  \item{b}{the prior scale parameter for the gamma distribution of the shape parameter.}
   \item{till}{if the observation is taken every 2 months until the 6th month, till = 6, if even after 6 months, the last observation is taken after 2 months, set till to any number.}
  \item{to}{the time steps at which the observations are taken e.g. 2,4,6 months so to = 2}
  \item{to2}{if the first observations are taken every two months but the last one, three months after the one before the last, set to2 = 3, otherwise to2 is equal to "to".}

}

\details{
*The time-to-event variable of the observed survival dataset at interim must not include any NA or missing values for the patients without an event. Instead a high number, as presened in the example section shall be utilized.
}
\value{
 \item{Results}{a dataframe consisting of the assumed survival probability for the primary and intermediate endpoint, the assume shape and scale parameters of the Weibull model, the posterior probability of a clinical relevant outcome and the final decision.}
 \item{RStan Sampling Output}{the ouput of the RStan function \link[rstan]{sampling} }
  \item{Summary of Posterior Values at studyend}{a dataframe displaying the median, 5 percent and 95 percent quantile of the posterior values at study end for both shape and scale.}

}

\author{
Simrit Rattan}

\examples{
for a simulated random data set:
interim = interim(n.cores = 7, rddata = 1, discrete = 1,data = 0,n = 43,n1 = 24,p0 = 0.127,p2_0 = 0.4,k = 4,tau = 9,omega1 = 24,omega2 = 19,kappa = 1.001144,theta = 4.36499,d = 1000,eta_star = 0.95,lambda = 0.10,mu = 0,sigma = 10,a=1,b = 0.001, till = 6, to = 2, to2 = 3)
final(rddata = 1, discrete = 1, interim = interim, data = 0,n = 43,n1 = 24,p0 = 0.127,p2_0 = 0.4,k = 4,tau = 9,omega2 = 19,kappa = 1.001144,theta = 4.36499,eta_star = 0.95,mu = 0,sigma = 10,a = 1,b = 0.001, till = 6, to = 2, to2 = 3)

for an external data set:
data = data.frame("e" = c(4,8,10,3.5,12.5,21), "y" = c(2,4,1000,1000,1000,1000), "y_e" = c(6,12,1010,1003.5,1012.5,1021),"delta" = c(1,1,0,0,0,0))
#y = 1000 indicates that no TTE data is available as no event has occured so far for this patient.

final(rddata = 0, discrete = 1, interim = 0, data = data , n = 43,n1 = 12,p0 = 0.127,p2_0 = 0.4,k = 4,tau = 9,omega2 = 19,kappa = 1.001144,theta = 4.36499,eta_star = 0.95,mu = 0,sigma = 10,a = 1,b = 0.001, till = 6, to = 2, to2 = 3)

}
