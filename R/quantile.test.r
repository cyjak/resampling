#' @export

quantile.test.r = function (x, y = NULL, quantile = 0.5,  alternative = c("two.sided", "less", "greater"),
                     mu = 0, conf.level = 0.95,
                     seed=0, n.perm=10000, confint = TRUE, pvalue = TRUE, boot.type = c("bca", "percentile"), boot.values = F, perm.values = F,
                     ...){

  boot.type <- match.arg(boot.type)
  alternative <- match.arg(alternative)
  if (!missing(mu) && (length(mu) != 1 || is.na(mu)))    stop("'mu' must be a single number")
  if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) || conf.level < 0 || conf.level > 1)){
    stop("'conf.level' must be a single number between 0 and 1")}

  if (!missing(quantile) && (length(quantile) != 1 || !is.finite(quantile) || quantile < 0 || quantile > 1)){
    stop("'quantile' must be a single number between 0 and 1")}

    xok <- !is.na(x)
    x <- x[xok]
    nx <- length(x)

  if (!is.null(y)) {
    dname <- paste(deparse1(substitute(x)), "and", deparse1(substitute(y)))
    yok <- !is.na(y)
    y <- y[yok]
    ny <- length(y)
    theta_hat <- quantile(y, quantile) - quantile(x, quantile) - mu
    estimate <- c(quantile(x, quantile), quantile(y, quantile))
    names(estimate) <- c(paste0("quantile ",quantile," of x"), paste0("quantile ",quantile," of y"))

  } else {
    dname <- deparse1(substitute(x))
    theta_hat <- quantile(x, quantile) - mu
    estimate = theta_hat
    names(estimate) <- c(paste0("quantile ",quantile," of x"))
  }

  if (is.null(y)) {
    if (nx < 2){
      stop("not enough 'x' observations")
    }
  } else {
    if (nx < 1)  stop("not enough 'x' observations")
    if (ny < 1)  stop("not enough 'y' observations")
  }

    method = paste0(quantile,' quantile test')


    ### permutation
    if (pvalue){
      if (!is.null(y)){
          all_values = c(x,y)
          set.seed(seed)
          perm_theta_hat = replicate(n.perm, {
            positions = sample(x=c(rep(1, length(x)), rep(2, length(y))), size=length(all_values), replace=FALSE)
            quantile(all_values[positions==2], quantile) - quantile(all_values[positions==1], quantile)
          })

      } else {
        set.seed(seed)
        perm_theta_hat = replicate(n.perm, { quantile(sign(runif(length(x), -1, 1)) * abs(x), quantile)  })
      }


        ### p-values

        if (alternative == 'two.sided'){
          temp_quantile = ecdf(perm_theta_hat) (theta_hat)
          temp_quantile = ifelse(temp_quantile>.5, 1-temp_quantile, temp_quantile)
          pval = temp_quantile*2

        } else if (alternative == "less") {
          temp_quantile = ecdf(perm_theta_hat) (theta_hat)
          temp_quantile = ifelse(temp_quantile<.5, 1-temp_quantile, temp_quantile)
          pval = ifelse(theta_hat<.0, 1-temp_quantile, temp_quantile)

        }  else if (alternative == "greater") {
          temp_quantile = ecdf(perm_theta_hat) (theta_hat)
          temp_quantile = ifelse(temp_quantile>.5, 1-temp_quantile, temp_quantile)
          pval = ifelse(theta_hat<.0, 1-temp_quantile, temp_quantile)

        }

    }






    ### bootstrap
    alpha <- 1 - conf.level
    if (alternative=='two.sided') sided=2 else sided=1
    set.seed(seed)
    if (confint){
      if (!is.null(y)){
        boot_theta_hat = replicate(n.perm, {
          quantile(sample(y, length(y), replace=TRUE), quantile) - quantile(sample(x, length(x), replace=TRUE), quantile) - mu })
        } else {
          boot_theta_hat = replicate(n.perm, {
            quantile(sample(x, length(x), replace=TRUE), quantile) - mu   })
        }


    quantiles = c(alpha/sided, 1-alpha/sided)

    if (boot.type == 'bca'){
      if (!is.null(y)){
        cint = boot.bca(x, y, boot.values=boot_theta_hat, theta_hat=theta_hat, paired=FALSE, quantiles=quantiles, fun=function(x,y) quantile(y, quantile) - quantile(x, quantile)-mu)
      } else {
        cint = boot.bca(x, boot.values=boot_theta_hat, theta_hat=theta_hat, paired=FALSE, quantiles=quantiles, fun=function(x) quantile(x, quantile)-mu)
      }
    }

    if (alternative=='less') cint[1] = -Inf
    if (alternative=='greater') cint[2] = Inf

    }




    method = paste0(method,' (permutation & bootstrap)')
    names(theta_hat) <- "delta"
    names(n) <- "n"


    RVAL = list(statistic = theta_hat, df = n,
                null.value = mu,
                alternative = alternative, method = method, estimate = estimate,
                boot.type = boot.type)

    if (pvalue) RVAL <- c(RVAL, list(p.value = as.numeric(pval)))
    if (boot.values) RVAL = c(RVAL, list(boot.values = boot_theta_hat))
    if (perm.values) RVAL = c(RVAL, list(perm.values = perm_theta_hat))
    if (confint){
      attr(cint, "conf.level") <- conf.level
      RVAL <- c(RVAL, list(conf.int = cint))
    }

    class(RVAL) <- "htest"
    RVAL
}


# set.seed(0) ; n=100
# x1 = rnorm(n,0,10)
# x2 = rnorm(n,5,10)
# median(x2)-median(x1)
# tt = quantile.test.r(x1,x2, mu=0) ; tt









































