\name{data_studyend}
\alias{data_studyend}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Observed and predicted survival data.
}
\description{
The function has two purposes:
First it takes the predicted survival data and determines the number of events,
the number of censored patients, the time-to-event data rounded and unrounded,
and the censoring times rounded and unrounded for the specified primary endpoint at the time of the final analysis.
Second it adds the observed data at interim to the predicted data.
}
\usage{
data_studyend(data,prediction,tau,k,till,to,to2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a dataframe containing the observed accrual time variable e, the time-to-event variable y, the event time variable y_e
  and the event indicator variable delta at interim for the first stage patients.*}
  \item{prediction}{a dataframe containing the predicted accrual time variable e, the time-to-event variable y,
  and the event indicator variable delta at studyend for the remaining m1 first stage patients and the new n-n1 second stage patients.}
  \item{tau}{the primary endpoint.}
  \item{k}{the intermediate endpoint}
    \item{till}{if the observation is taken every 2 months until the 6th month, till = 6, if even after 6 months, the last observation is taken after 2 months, set till to any number.}
  \item{to}{the time steps at which the observations are taken e.g. 2,4,6 months so to = 2}
  \item{to2}{if the first observations are taken every two months but the last one, three months after the one before the last, set to2 = 3, otherwise to2 is equal to "to".}
}

\details{
*The time-to-event variable of the observed survival dataset at interim must not include any NA or missing values for the patients without an event. Instead a high number, as presened in the example section shall be utilized.
}

\value{
A list containing the following components:
  \item{mevent}{the number of events out of n - n1 + m1 patients.}
  \item{ncens}{the number of censored patients out of n patients (administrative censoring).}
  \item{nevent}{the number of events out of all n patients.}
  \item{y_studyend}{the time-to-event times for the number of events out of n patients.(observed + predicted)}
  \item{y_studyend_rounded}{the rounded time-to-event times for the number of events out of n patients.(observed + predicted)}
  \item{z_studyend}{the censoring times, therefore the time spent in the study without experiencing an event.}
  \item{z_studyend_rounded}{the rounded censoring times, therefore the time spent in the study without experiencing an event.}
}

\author{
Simrit Rattan
}

\examples{
data = data.frame("e" = c(4,8,10,3.5,12.5,21), "y" = c(2,4,1000,1000,1000,1000), "y_e" = c(6,12,1010,1003.5,1012.5,1021),"delta" = c(1,1,0,0,0,0))
#y = 1000 indicates that no TTE data is available as no event has occured so far for this patient.
prediction = data.frame("e" = c(1.5,3.5,12.5,21,25,25.5,40,42,26.5,33), "y" = c(5,9,6,5,3,12,14,15,2,1), "y_e" = c(6.5,12.5,18.5,26,28,37.5,54,47,28.5,34), "delta" = c(1,1,1,1,1,0,0,0,1,1))
data_studyend(data = data,prediction = prediction,tau = 9, k = 4, till = 6, to = 2, to2 = 3)
}
