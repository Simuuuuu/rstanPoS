\name{secondstage_samples}
\alias{secondstage_samples}
\title{
The function generates Weibull survival data fort the second stage of a phase II clinical trial.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
secondstage_samples(data,n,n1,tau,omega2,kappa,theta,k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a dataframe containing the observed accrual time variable e, the time-to-event variable y, the event time variable y_e
  and the event indicator variable delta at interim for the first stage patients.*}
  \item{n}{the total number of patients in the study.}
  \item{n1}{the number of first stage patients.}
  \item{tau}{the primary endpoint.}
  \item{omega2}{the accrual period for the second stage patients.}
  \item{kappa}{the assumed shape parameter of the Weibull distribution.}
  \item{theta}{the assumed scale parameter of the Weibull distribution.}
  \item{k}{the intermediate endpoint.}
}
\details{
The following five variables are generated: e,y,y_e,z,delta.
The accrual time e follows a uniform distribution U(max(e1) + looktime,max(e1) + looktime + omega2). The max(e1) + looktime indicates that the second stage patients are recruited once the intermediate endpoint of the last first stage patient has been observed.
Furthermore the survival time y, which is the time-to-event data, follows a Weibull distribution W(kappa,theta).

*The time-to-event variable of the observed survival dataset at interim must not include any NA or missing values for the patients without an event. Instead a high number, as presened in the example section shall be utilized.
}


\value{
A dataframe containing the following variables:
  \item{e}{the accrual times.}
  \item{y}{the time-to-event data.}
  \item{y_e}{the actual time of the event.}
  \item{z}{the time of administrative censoring, which is the time spent in the study without experiencing an event.}
  \item{delta}{the event indicator.}
}

\author{
Simrit Rattan
}

\examples{
data = data.frame("e" = c(4,8,10,3.5,12.5,21), "y" = c(2,4,1000,1000,1000,1000), "y_e" = c(6,12,1010,1003.5,1012.5,1021),"delta" = c(1,1,0,0,0,0))
#y = 1000 indicates that no TTE data is available as no event has occured so far for this patient.
secondstage_samples(data = data,n = 43,n1 = 24,tau = 9,omega2 = 19,kappa = 1.15,theta = 6.25,k = 4)
}
