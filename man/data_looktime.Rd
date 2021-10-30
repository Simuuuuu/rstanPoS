\name{data_looktime}
\alias{data_looktime}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Survival data up to looktime
}
\description{
The function takes the survival data and determines the number of events,
the number of censored patients, the time-to-event data and
the censoring times for a specified intermediate endpoint at the time of the interim analysis. Moreover, it also determines the rounded time-to-event data and censoring times.
}
\usage{
data_looktime(data,k,till,to,to2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a dataframe containing the observed accrual time variable e, the time-to-event variable y, the event time variable y_e
  and the event indicator variable delta at interim for the first stage patients.*}
  \item{k}{the intermediate endpoint.}
    \item{till}{if the observation is taken every 2 months until the 6th month, till = 6, if even after 6 months, the last observation is taken after 2 months, set till to any number.}
  \item{to}{the time steps at which the observations are taken e.g. 2,4,6 months so to = 2}
  \item{to2}{if the first observations are taken every two months but the last one, three months after the one before the last, set to2 = 3, otherwise to2 is equal to "to".}
}

\details{
*The time-to-event variable of the observed survival dataset at interim must not include any NA or missing values for the patients without an event. Instead a high number, as presented in the example section shall be utilized.
}

\value{
A list containing the following components:
    \item{m1}{the number of censored patients out of n1 patients at interim.}
    \item{n1event }{the number of events out of n1 patients at interim.}
    \item{y_looktime}{the time-to-event times for the number of events out of n1 patients at interim.}
    \item{y_looktime_round}{the rounded time-to-event times for the number of events out of n1 patients at interim.}
    \item{z_looktime}{the censoring times, therefore the time spent in the study without experiencing an event.}
    \item{z_looktime_round}{the rounded censoring times, therefore the time spent in the study without experiencing an event.}

}

\author{
Simrit Rattan
}

\examples{
data = data.frame("e" = c(4,8,10,3.5,12.5,21), "y" = c(2,4,1000,1000,1000,1000), "y_e" = c(6,12,1010,1003.5,1012.5,1021),"delta" = c(1,1,0,0,0,0))
#y = 1000 indicates that no TTE data is available as no event has occured so far for this patient.
data_looktime(data = data,k = 4, till = 6, to = 2, to2 = 3)

}
