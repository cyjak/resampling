\name{cor.test.r}
\alias{cor.test.r}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
Using resampling methods, test for association between paired samples, using one of Pearson's product moment correlation coefficient, Kendall's tau or Spearman rho.
}
\usage{
cor.test.r(x, y)

cor.test(x, y,
         alternative = c("two.sided", "less", "greater"),
         method = c("pearson", "kendall", "spearman"),
         conf.level = 0.95, continuity = FALSE,
         n.perm = 10000, seed = 0, boot.type = c("bca", "percentile"),
         pvalue = TRUE, confint = TRUE, boot.values = F, perm.values = F, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, y}{numeric vectors of data values. x and y must have the same length.}
  \item{alternative}{indicates the alternative hypothesis and must be one of "\code{two.sided}", "\code{greater}" or "\code{less}". "\code{greater}" corresponds to positive association, "\code{less}" to negative association.}
  \item{method}{a character string indicating which correlation coefficient is to be used for the test. One of "\code{pearson}", "\code{kendall}", or "\code{spearman}".}
  \item{conf.level}{confidence level for the returned confidence interval. Currently only used for the Pearson product moment correlation coefficient if there are at least 4 complete pairs of observations.}
  \item{n.perm}{an integer indicating the number of permutations.}
  \item{seed}{an integer indicating the number of the seed to ensure replicability of the random process. Set to \code{NULL} to remove the use of a seed.}
  \item{boot.type}{a character string specifying the bootstrap method to be used for computing the confidence intervals, must be one of "\code{bca}" (default) or "\code{percentile}".}
  \item{pvalue}{a logical indicating whether you want the p-value to be computed}
  \item{confint}{a logical indicating whether you want the confidence intervals to be computed}
  \item{boot.values}{a logical indicating whether the values of all bootstraps should be returned.}
  \item{perm.values}{a logical indicating whether the values of all permutations should be returned.}

}

\details{
If \code{boot.type} is \code{bca} (default), the bias-corrected and accelerated (BCa) bootstrap interval is used to determine the quantiles to be used for the corresponding \code{conf.level}. If \code{boot.type} is \code{percentile}, the percentile corresponding to the \code{conf.level} are directly used. The BCa method (default) is know as being more robust to extreme values.

\code{pvalue} and \code{confint} allow to specify whether the p-value and the confidence intervals respectively are computed. If one is not interested in the one of the two metrics, disabling it allows for a faster running of the function.
}
\value{
A list with class "htest" containing the following components:
  \item{statistic}{the value of the correlation coefficient.}
  \item{parameter}{the sample size.}
  \item{p.value}{the p-value of the test.}
  \item{estimate}{the estimated measure of association.}
  \item{null.value}{the value of the association measure under the null hypothesis, always \code{0}.}
  \item{alternative}{a character string describing the alternative hypothesis.}
  \item{method}{a character string indicating how the association was measured.}
  \item{conf.int}{a confidence interval for the measure of association. Only computed if \code{confint = TRUE}.}

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
x = rnorm(n = 100, mean = 0, sd = 1)
y = rnorm(n = 100, mean = 0, sd = 1) + x

## Pearson's correlation
cor.test.r(x, y)

## Spearman's correlation
cor.test.r(x, y, method = "spearman")

## Kendall's correlation
cor.test.r(x, y, method = "kendall")

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
