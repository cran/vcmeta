% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meta_ave.R
\name{meta.ave.mean.ps}
\alias{meta.ave.mean.ps}
\title{Confidence interval for an average mean difference from paired-samples studies}
\usage{
meta.ave.mean.ps(alpha, m1, m2, sd1, sd2, cor, n, bystudy = TRUE)
}
\arguments{
\item{alpha}{alpha level for 1-alpha confidence}

\item{m1}{vector of estimated means for measurement 1}

\item{m2}{vector of estimated means for measurement 2}

\item{sd1}{vector of estimated SDs for measurement 1}

\item{sd2}{vector of estimated SDs for measurement 2}

\item{cor}{vector of estimated correlations for paired measurements}

\item{n}{vector of sample sizes}

\item{bystudy}{logical to also return each study estimate (TRUE) or not}
}
\value{
Returns a matrix.  The first row is the average estimate across all studies.  If bystudy
is TRUE, there is 1 additional row for each study.  The matrix has the following columns:
\itemize{
\item Estimate - estimated effect size
\item SE - standard error
\item LL - lower limit of the confidence interval
\item UL - upper limit of the confidence interval
\item df - degrees of freedom
}
}
\description{
Computes the estimate, standard error, and confidence interval for an
average mean difference from two or more paired-samples studies.
A Satterthwaite adjustment to the degrees of freedom is used to improve
the accuracy of the confidence interval for the average effect size.
Equality of variances within or across studies is not assumed.
}
\examples{
m1 <- c(53, 60, 53, 57)
m2 <- c(55, 62, 58, 61)
sd1 <- c(4.1, 4.2, 4.5, 4.0)
sd2 <- c(4.2, 4.7, 4.9, 4.8)
cor <- c(.7, .7, .8, .85)
n <- c(30, 50, 30, 70)
meta.ave.mean.ps(.05, m1, m2, sd1, sd2, cor, n, bystudy = TRUE)

# Should return:
#        Estimate        SE        LL         UL      df
# Average   -3.25 0.2471557 -3.739691 -2.7603091 112.347
# Study 1   -2.00 0.5871400 -3.200836 -0.7991639  29.000
# Study 2   -2.00 0.4918130 -2.988335 -1.0116648  49.000
# Study 3   -5.00 0.5471136 -6.118973 -3.8810270  29.000
# Study 4   -4.00 0.3023716 -4.603215 -3.3967852  69.000


}
\references{
\insertRef{Bonett2009a}{vcmeta}
}
