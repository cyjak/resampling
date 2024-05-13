\name{chisq.test.r}
\alias{chisq.test.r}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
chisq.test.r %%  ~~function to do ... ~~
}
\description{
Using resampling methods, performs chi-squared contingency table tests.}
\usage{
chisq.test.r(x, y = NULL, correct = TRUE,
             p = rep(1/length(x), length(x)), rescale.p = FALSE,
             n.perm = 10000, seed = 0, perm.values = F, pvalue = TRUE,
          confint = TRUE, conf.level = 0.95, boot.type = c("bca"), boot.values = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric vector or matrix. x and y can also both be factors.}
  \item{y}{a numeric vector; ignored if x is a matrix. If x is a factor, y should be a factor of the same length.}
  \item{correct}{a logical indicating whether to apply continuity correction when computing the test statistic for 2 by 2 tables: one half is subtracted from all |O - E|∣O−E∣ differences; however, the correction will not be bigger than the differences themselves. No correction is done if \code{simulate.p.value = TRUE}.}
  \item{p}{a vector of probabilities of the same length as \code{x}. An error is given if any entry of \code{p} is negative.}
  \item{rescale.p}{a logical scalar; if \code{TRUE} then \code{p} is rescaled (if necessary) to sum to 1. If \code{rescale.p} is \code{FALSE}, and \code{p} does not sum to 1, an error is given.}
  \item{n.perm}{an integer indicating the number of permutations.}
  \item{seed}{an integer indicating the number of the seed to ensure replicability of the random process. Set to \code{NULL} to remove the use of a seed.}
  \item{perm.values}{a logical indicating whether the values of all permutations should be returned.}
  \item{pvalue}{a logical indicating whether you want the p-value to be computed}
  \item{confint}{a logical indicating whether you want the confidence intervals to be computed. Only used for 2x2 tables.}
  \item{conf.level}{confidence level for the returned confidence interval.}
  \item{boot.type}{a character string specifying the bootstrap method to be used for computing the confidence intervals.}
  \item{boot.values}{a logical indicating whether the values of all bootstraps should be returned.}


}
\details{
When comparing two proportions (2x2 table), the difference in proportions and its confidence intervals can be calculated. If \code{boot.type} is \code{bca} (default), the bias-corrected and accelerated (BCa) bootstrap interval is used to determine the quantiles to be used for the corresponding \code{conf.level}. If \code{boot.type} is \code{percentile}, the percentile corresponding to the \code{conf.level} are directly used. The BCa method (default) is know as being more robust to extreme values.
}
\value{
A list with class "htest" containing the following components:
 \item{statistic}{the value of the test statistic.}
 \item{parameter}{sample size.}
 \item{p.value}{the p-value of the test.}
 \item{method}{ character string indicating the type of test performed.}
 \item{observed}{the observed counts.}
 \item{expected}{the expected counts under the null hypothesis.}
 \item{residuals}{the Pearson residuals, \code{(observed - expected) / sqrt(expected)}}
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.


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
