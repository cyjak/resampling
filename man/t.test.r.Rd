\name{t.test.r}
\alias{t.test.r}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
Using resampling methods, performs one and two sample t-test on vectors of data.
}
\usage{
t.test.r(x, ...)

t.test.r(x, y = NULL,
         alternative = c("two.sided", "less", "greater"),
         mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95,
         n.perm = 10000, seed = 0, boot.type = c("bca", "percentile"),
         pvalue = TRUE, confint = TRUE, boot.values = F, perm.values = F, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a (non-tempy) numeric vector of data values.}
  \item{y}{an optional (non-empty) numeric vector of data values.}
  \item{alternative}{a character string specifying the alternative hypothesis, must be one of "\code{two.sided}" (default), "\code{greater}" or "\code{less}".}
   \item{mu}{a number indicating the true value of the mean (or difference in means if you are performing a two sample test).}
   \item{paired}{a logical indicating whether you want a paired t-test.}
   \item{var.equal}{a logical variable indicating whether to treat the two variances as being equal. If TRUE then the pooled variance is used to estimate the variance otherwise the Welch (or Satterthwaite) approximation to the degrees of freedom is used.}
   \item{conf.level}{confidence level of the interval.}
   \item{n.perm}{an integer indicating the number of permutations.}
   \item{seed}{an integer indicating the number of the seed to ensure replicability of the random process. Set to \code{NULL} to remove the use of a seed.}
   \item{boot.type}{a character string specifying the bootstrap method to be used for computing the confidence intervals, must be one of "\code{bca}" (default) or "\code{percentile}".}
   \item{pvalue}{a logical indicating whether you want the p-value to be computed}
   \item{confint}{a logical indicating whether you want the confidence intervals to be computed}
   \item{boot.values}{a logical indicating whether the values of all bootstraps should be returned.}
   \item{perm.values}{a logical indicating whether the values of all permutations should be returned.}

}
}
\details{

If \code{paired} is \code{TRUE} then both x and y must be specified and they must be the same length. Missing values are silently removed (in pairs if \code{paired} is \code{TRUE}). If \code{var.equal} is \code{TRUE} then the pooled estimate of the variance is used. By default, if \code{var.equal} is \code{FALSE} then the variance is estimated separately for both groups and the Welch modification to the degrees of freedom is used.

If \code{boot.type} is \code{bca} (default), the bias-corrected and accelerated (BCa) bootstrap interval is used to determine the quantiles to be used for the corresponding \code{conf.level}. If \code{boot.type} is \code{percentile}, the percentile corresponding to the \code{conf.level} are directly used. The BCa method (default) is know as being more robust to extreme values.

\code{pvalue} and \code{confint} allow to specify whether the p-value and the confidence intervals respectively are computed. If one is not interested in the one of the two metrics, disabling it allows for a faster running of the function.

}
\value{
A list with class "htest" containing the following components:
  \item{delta}{difference between \code{mu} and the mean of \code{x} for one sample t-test, the means' difference between \code{x} and \code{y} for two independent samples t-test, or the mean of the difference between \code{x} and \code{y} for paired t-test.}
  \item{p.value}{the p-value for the test, calculated with a bootstrap method.}
  \item{conf.int}{a confidence interval for the mean appropriate to the specified alternative hypothesis. Only computed if \code{confint = TRUE}.}
  \item{estimate}{the estimated mean or difference in means depending on whether it was a one-sample test or a two-sample test.}
  \item{null.value}{the specified hypothesized value of the mean or mean difference depending on whether it was a one-sample test or a two-sample test.}
  \item{alternative}{a character string describing the alternative hypothesis.}
  \item{method}{a character string indicating what type of t-test was performed.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Generate data
set.seed(123)
x = rnorm(n = 100, mean = 1, sd = 10)
y = rnorm(n = 100, mean = 0, sd = 10)

## One sample t-test
t.test.r(x)

## Two independent sample t-test
t.test.r(x, y)

## Paired t-test
t.test.r(x, y, paired = TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
