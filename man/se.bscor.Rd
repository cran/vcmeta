% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meta_se.R
\name{se.bscor}
\alias{se.bscor}
\title{Computes the standard error for a biserial correlation}
\usage{
se.bscor(m1, m2, sd1, sd2, n1, n2)
}
\arguments{
\item{m1}{estimated mean for level 1}

\item{m2}{estimated mean for level 2}

\item{sd1}{estimated standard deviation for level 1}

\item{sd2}{estimated standard deviation for level 2}

\item{n1}{sample size for level 1}

\item{n2}{sample size for level 2}
}
\value{
Returns a one-row matrix:
\itemize{
\item Estimate - estimated biserial correlation
\item SE - standard error
}
}
\description{
This function computes a biserial correlation and its standard error. A
biserial correlation can be used when one variable is quantitative and the
other variable has been artifically dichotmized. The biserial correlation
estimates the correlation between an observable quantitative variable and
an unobserved quantitative variable that is measured on a dichotomous
scale. This function requires the estimated mean, estimated standard
deviation, and samples size from each level of the dichotomized variable.
This function is useful in a meta-analysis of Pearson correlations where
some studies report a Pearson correlation and other studies report the
information needed to compute a biserial correlation. The biserial
correlation and standard error output from this function can be used as
input in the \link[vcmeta]{meta.ave.cor.gen} function.
}
\details{
This function computes a point-biserial correlation and its standard error
as a function of a standardized mean difference with a weighted variance
standardizer. Then the point-biserial estimate is transformed into a
biserial correlation using the traditional adjustment. The adjustment is
also applied to the point-biserial standard error to obtain the standard
error for the biserial correlation.

The biserial correlation assumes that the observed quantitative variable
and the unobserved quantitative variable have a bivariate normal
distribution. Bivariate normality is a crucial assumption underlying the
transformation of a point-biserial correlation to a biserial correlation.
Bivariate normality also implies equal variances of the observed
quantitative variable at each level of the dichotomized variable, and this
assumption is made in the computation of the standard error.
}
\examples{
se.bscor(21.9, 16.1, 3.82, 3.21, 40, 40)

#  Should return: 
#                          Estimate         SE
#  Biserial correlation:  0.8018318 0.07451665


}
\references{
\insertRef{Bonett2020b}{vcmeta}
}
