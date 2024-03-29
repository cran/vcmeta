% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meta_model.R
\name{meta.lm.spear}
\alias{meta.lm.spear}
\title{Meta-regression analysis for Spearman correlations}
\usage{
meta.lm.spear(alpha, n, cor, X)
}
\arguments{
\item{alpha}{alpha level for 1-alpha confidence}

\item{n}{vector of sample sizes}

\item{cor}{vector of estimated Spearman correlations}

\item{X}{matrix of predictor values}
}
\value{
Returns a matrix.  The first row is for the intercept with one additional
row per predictor.  The matrix has the following columns:
\itemize{
\item Estimate - OLS estimate
\item SE - standard error
\item z - z-value
\item p - p-value
\item LL - lower limit of the confidence interval
\item UL - upper limit of the confidence interval
}
}
\description{
This function estimates the intercept and slope coefficients in a
meta-regression model where the dependent variable is a
Fisher-transformed Spearman correlation. The estimates are OLS estimates
with robust standard errors that accommodate residual heteroscedasticity.
The correlations are Fisher-transformed and hence the parameter
estimates do not have a simple interpretation. However, the hypothesis
test results can be used to decide if a population slope is either
positive or negative.
}
\examples{

n <- c(150, 200, 300, 200, 350)
cor <- c(.14, .29, .16, .21, .23)
x1 <- c(18, 25, 23, 19, 24)
X <- matrix(x1, 5, 1)
meta.lm.spear(.05, n, cor, X)

# Should return: 
#       Estimate         SE          z     p           LL         UL
# b0 -0.08920088 0.26686388 -0.3342561 0.738 -0.612244475 0.43384271
# b1  0.01370866 0.01190212  1.1517825 0.249 -0.009619077 0.03703639


}
