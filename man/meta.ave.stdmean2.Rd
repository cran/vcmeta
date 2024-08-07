% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meta_ave.R
\name{meta.ave.stdmean2}
\alias{meta.ave.stdmean2}
\title{Confidence interval for an average standardized mean difference
from 2-group studies}
\usage{
meta.ave.stdmean2(alpha, m1, m2, sd1, sd2, n1, n2, stdzr, bystudy = TRUE)
}
\arguments{
\item{alpha}{alpha level for 1-alpha confidence}

\item{m1}{vector of estimated means for group 1}

\item{m2}{vector of estimated means for group 2}

\item{sd1}{vector of estimated SDs for group 1}

\item{sd2}{vector of estimated SDs for group 2}

\item{n1}{vector of group 1 sample sizes}

\item{n2}{vector of group 2 sample sizes}

\item{stdzr}{\itemize{
\item set to 0 for square root unweighted average variance standardizer
\item set to 1 for group 1 SD standardizer
\item set to 2 for group 2 SD standardizer
\item set to 3 for square root weighted average variance standardizer
}}

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
}
}
\description{
Computes the estimate, standard error, and confidence interval for an
average standardized mean difference from two or more 2-group studies.
Square root unweighted variances, square root weighted variances, and
single group standard deviation are options for the standardizer.
Equality of variances within or across studies is not assumed.
}
\examples{
m1 <- c(21.9, 23.1, 19.8)
m2 <- c(16.1, 17.4, 15.0)
sd1 <- c(3.82, 3.95, 3.67)
sd2 <- c(3.21, 3.30, 3.02)
n1 <- c(40, 30, 24)
n2 <- c(40, 28, 25)
meta.ave.stdmean2(.05, m1, m2, sd1, sd2, n1, n2, 0, bystudy = TRUE)

# Should return: 
#         Estimate        SE        LL       UL
# Average 1.526146 0.1734341 1.1862217 1.866071
# Study 1 1.643894 0.2629049 1.1286100 2.159178
# Study 2 1.566132 0.3056278 0.9671126 2.165152
# Study 3 1.428252 0.3289179 0.7835848 2.072919


}
\references{
\insertRef{Bonett2009a}{vcmeta}
}
