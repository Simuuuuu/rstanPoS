\name{rdsamples_generator}
\alias{rdsamples_generator}
\title{Random samples generator for first stage survival data
}
\description{
The function generates Weibull survival data fort the first stage of a phase II clinical trial.
}
\usage{
rdsamples_generator(n1,tau,k,omega1,kappa,theta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n1}{the number of first stage patients.}
  \item{tau}{the primary endpoint.}
  \item{k}{the intermediate endpoint.}
  \item{omega1}{the accrual period for the first stage patients.}
  \item{kappa}{the assumed shape parameter of the Weibull distribution.}
  \item{theta}{the assumed scale parameter of the Weibull distribution.}

}
\details{
The following five variables are generated: e,y,y_e,z,delta.
The accrual time e follows a Uniform distribution U(0,omega1).
Furthermore the survival time y, which is the time-to-event data, follows a Weibull distribution W(kappa,theta).

}
\value{
A dataframe containing the following variables:
  \item{e}{the accrual times.}
  \item{y}{the time-to-event data.}
  \item{y_e}{the actual time of the event.}
  \item{z}{the administrative censoring times, which are equal to tau.}
  \item{delta}{the event indicator.}
}

\author{Simrit Rattan
}

\examples{
rdsamples_generator(n1 = 24,tau = 9,k = 2,omega1 = 24,kappa = 1.15,theta = 6.25)

}
