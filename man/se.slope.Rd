% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meta_se.R
\name{se.slope}
\alias{se.slope}
\title{Computes a slope and standard error}
\usage{
se.slope(cor, sdy, sdx, n)
}
\arguments{
\item{cor}{estimated Pearson correlation}

\item{sdy}{estimated standard deviation of the response variable}

\item{sdx}{estimated standard deviation of the predictor variable}

\item{n}{sample size}
}
\value{
Returns a one-row matrix:
\itemize{
\item Estimate - estimated slope
\item SE - standard error
}
}
\description{
This function computes a slope and its standard error
for a simple linear regression model (random-x model) using the estimated
Pearson correlation and the estimated standard deviations of the response
variable and predictor variable. This function is useful in a meta-analysis
of slopes of a simple linear regression model where some studies report
the Pearson correlation but not the slope.
}
\examples{
se.slope(.392, 4.54, 2.89, 60)

# Should return: 
#          Estimate        SE
# Slope:  0.6158062 0.1897647


}
\references{
\insertRef{Snedecor1980}{vcmeta}
}
