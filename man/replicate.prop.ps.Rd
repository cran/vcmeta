% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meta_rep.R
\name{replicate.prop.ps}
\alias{replicate.prop.ps}
\title{Compares and combines paired-samples proportion differences in original and
follow-up studies}
\usage{
replicate.prop.ps(alpha, f1, f2)
}
\arguments{
\item{alpha}{alpha level for 1-alpha confidence}

\item{f1}{vector of frequency counts for 2x2 table in original study}

\item{f2}{vector of frequency counts for 2x2 table in follow-up study}
}
\value{
A 4-row matrix. The rows are:
\itemize{
\item Row 1 summarizes the original study
\item Row 2 summarizes the follow-up study
\item Row 3 estimates the difference in proportion differences
\item Row 4 estimates the average proportion difference
}

The columns are:
\itemize{
\item Estimate - proportion difference estimate (single study, difference, average)
\item SE - standard error
\item z - z-value
\item p - p-value
\item LL - lower limit of the confidence interval
\item UL - upper limit of the confidence interval
}
}
\description{
This function computes confidence intervals from an original study and a
follow-up study where the effect size is a paired-samples proportion
difference. Confidence intervals for the difference and average of effect
sizes are also computed. The confidence level for the difference is
1 – 2*alpha, which is recommended for equivalence testing.
}
\examples{
f1 <- c(42, 2, 15, 61)
f2 <- c(69, 5, 31, 145)
replicate.prop.ps(.05, f1, f2)

# Should return:
#                           Estimate         SE           z            p
# Original:              0.106557377 0.03440159  3.09745539 1.951898e-03
# Follow-up:             0.103174603 0.02358274  4.37500562 1.214294e-05
# Original - Follow-up:  0.003852359 0.04097037  0.09402793 9.250870e-01
# Average:               0.105511837 0.02048519  5.15064083 2.595979e-07
#                                LL         UL
# Original:              0.03913151 0.17398325
# Follow-up:             0.05695329 0.14939592
# Original - Follow-up: -0.06353791 0.07124263
# Average:               0.06536161 0.14566206


}
\references{
\insertRef{Bonett2021}{vcmeta}
}
