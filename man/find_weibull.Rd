\name{find_weibull}
\alias{find_weibull}
\alias{plot.find_weibull}


\title{Determination of Weibull Parameters given two quantiles.
%%  ~~function to do ... ~~
}
\description{For two different endpoints and their survival probabilities, the function
determines the shape and scale parameter of the Weibull distribution as well as the median survival time.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
find_weibull(PFS_1_time,PFS_1_probability,PFS_2_time,PFS_2_probability)
\method{plot}{find_weibull}(x, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PFS_1_time}{the primary endpoint.}
  \item{PFS_1_probability}{the survival probability of the corresponding primary endpoint.}
  \item{PFS_2_time}{the intermediate endpoint.}
  \item{PFS_2_probability}{the survival probability of the corresponding intermediate endpoint.}
  \item{...}{additional arguments affecting the plot method.}

}
\details{
With the plot method, the Weibull survival curve is displayed.
}
\value{
A list with class "find_weibull" containing the following
components:
\item{Shape}{the corresponding shape parameter of the Weibull Distribution.}
\item{Scale}{the corresponding scale parameter of the Weibull Distribution.}
\item{Median}{the corresponding median survival time.}
\item{PFS1}{the survival probability of the corresponding primary endpoint.}
\item{PFS2}{the survival probability of the corresponding intermediate endpoint}
\item{PFS1time}{the primary endpoint.}
\item{PFS2time}{the indermediate endpoint.}

}
\references{
Cook, J. D. (2010). Determining distribution parameters from quantiles.
}
\author{
Dario Zocholl
}

\examples{
find_weibull(PFS_1_time = 9, PFS_1_probability = 0.127, PFS_2_time = 4,
                         PFS_2_probability = 0.4)
}

