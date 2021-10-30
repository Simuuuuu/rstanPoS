\name{prediction_samples}
\alias{prediction_samples}

\title{
Generates the predicted survival data based on the posterior parameter values.
}
\description{
The function generates survival data based on the posterior values obtained from
the sampling done in RStan. It is a prediction of the survival data in the trial
where the outcome is not known so far. This includes the survival data for the m1
patients without an event from the first stage, as well as the survival data for
the n-n1 second stage patients.
}
\usage{
prediction_samples(data,n,n1,m1,tau,omega2,kappa,theta,k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a dataframe containing the observed accrual time variable e, the time-to-event variable y, the event time variable y_e
  and the event indicator variable delta at interim for the first stage patients.*}
  \item{n}{the total number of patients in the study.}
  \item{n1}{the number of first stage patients.}
  \item{m1}{the number of censored patients out of n1 patients.}
  \item{tau}{the primary endpoint.}
  \item{omega2}{the accrual period in the second stage.}
  \item{kappa}{posterior value of the shape parameter of the Weibull model.}
  \item{theta}{posterior value of the scale parameter of the Weibull model.}
  \item{k}{the intermediate endpoint.}
}
\details{
The following five variables are generated: e,y,y_e,z,delta.
The accrual time e consists of the accrual time for the m1 first stage patients which is taken from the dataframe data and of the accrual times for the n-n1 second stage patients, for which a Uniform distribution U(max(e1) + looktime, max(e1) + looktime + omega2) is assumed. This indicates that the second stage patients are recruited once the intermediate endpoint of the last first stage patient has been observed. The time-to-event data for the remaining m1 first stage patients is determined by generating random numbers from a Weibull distribution W(kappa,theta) and adding to it the time spent in the study without an event.
Furthermore the time-to-event variable y for the second stage patients, follows a Weibull distribution W(kappa,theta). For all the predicted time-to-event data kappa and theta are taken from the estimated posterior distributions.

The function is used to determine the Bayesian probability of success, which is based on the observed data until looktime and the predicted data of the remaining patients. To estimate the Bayesian probability of success, data will be predicted for a larger number of times, therefore the function will be used multiple of times in the function `interim()`, serving as a supportive function.

*The time-to-event variable of the observed survival dataset at interim must not include any NA or missing values for the patients without an event. Instead a high number, as presened in the example section shall be utilized.
}
\value{
A dataframe containing the following variables:
\item{e}{the predicted accrual times.}
  \item{y}{the predicted time-to-event data.}
  \item{y_e}{the predicted actual time of the event.}
  \item{z}{the predicted time of administrative censoring, which is equal to tau.}
  \item{delta}{the event indicator.}
}

\author{
Simrit Rattan
}

\examples{
data = data.frame("e" = c(4,8,10,3.5,12.5,21), "y" = c(2,4,1000,1000,1000,1000), "y_e" = c(6,12,1010,1003.5,1012.5,1021),"delta" = c(1,1,0,0,0,0))
#y = 1000 indicates that no TTE data is available as no event has occured so far for this patient.
prediction_samples(data = data,n = 25, n1 = 6, m1 = 4,tau = 9,omega2 = 19,kappa = 1.2,theta = 5.9,k = 4)

}

