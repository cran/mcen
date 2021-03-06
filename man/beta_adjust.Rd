% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcen.r
\name{beta_adjust}
\alias{beta_adjust}
\title{Adjusts the value of the coefficients to account for the scaling of x and y.}
\usage{
beta_adjust(beta, sigma_x, sigma_y, mean_x, mean_y)
}
\arguments{
\item{beta}{The estiamte of beta with scaled data.}

\item{sigma_x}{Sample tandard deviations of the original predictors.}

\item{sigma_y}{Sample standard deviations of the orignal responses.}

\item{mean_x}{Sample means of the original predictors .}

\item{mean_y}{Sample means of the original responses.}
}
\value{
Returns the adjusted coefficients
}
\description{
Adjusts the value of the coefficients to account for the scaling of x and y.
}
\author{
Ben Sherwood <ben.sherwood@ku.edu>, Brad Price <brad.price@mail.wvu.edu>
}
