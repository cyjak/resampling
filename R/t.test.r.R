#' @export

t.test.r = function (x, y = NULL, alternative = c("two.sided", "less", "greater"),
          mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95,
          seed=0, n.perm=10000, confint = TRUE, pvalue = TRUE, boot.type = c("bca", "percentile"), boot.values = F, perm.values = F,
          ...)
{
  boot.type <- match.arg(boot.type)
  alternative <- match.arg(alternative)
  if (!missing(mu) && (length(mu) != 1 || is.na(mu)))
    stop("'mu' must be a single number")
  if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) ||
                               conf.level < 0 || conf.level > 1))
    stop("'conf.level' must be a single number between 0 and 1")
  if (!is.null(y)) {
    dname <- paste(deparse1(substitute(x)), "and", deparse1(substitute(y)))
    if (paired)
      xok <- yok <- complete.cases(x, y)
    else {
      yok <- !is.na(y)
      xok <- !is.na(x)
    }
    y <- y[yok]
  }
  else {
    dname <- deparse1(substitute(x))
    if (paired)
      stop("'y' is missing for paired test")
    xok <- !is.na(x)
    yok <- NULL
  }
  x <- x[xok]
  if (paired) {
    x <- x - y
    y <- NULL
  }
  nx <- length(x)
  mx <- mean(x)
  vx <- var(x)

  if (is.null(y)) {
    if (nx < 2)
      stop("not enough 'x' observations")
    stderr <- sqrt(vx/nx)
    if (stderr < 10 * .Machine$double.eps * abs(mx))
      stop("data are essentially constant")
    method <- if (paired)
      "Paired t-test"
    else "One Sample t-test"
    estimate <- setNames(mx, if (paired)
      "mean difference"
      else "mean of x")
    n = c(nx)
    theta_hat <- (mx - mu)
  }
  else {
    ny <- length(y)
    if (nx < 1 || (!var.equal && nx < 2))
      stop("not enough 'x' observations")
    if (ny < 1 || (!var.equal && ny < 2))
      stop("not enough 'y' observations")
    if (var.equal && nx + ny < 3)
      stop("not enough observations")
    my <- mean(y)
    vy <- var(y)
    method <- paste(if (!var.equal)
      "Two Sample t-test")
    estimate <- c(mx, my)
    names(estimate) <- c("mean of x", "mean of y")
    if (var.equal) {
      df <- nx + ny - 2
      v <- 0
      if (nx > 1)
        v <- v + (nx - 1) * vx
      if (ny > 1)
        v <- v + (ny - 1) * vy
      v <- v/df
      stderr <- sqrt(v * (1/nx + 1/ny))
    }
    else {
      stderrx <- sqrt(vx/nx)
      stderry <- sqrt(vy/ny)
      stderr <- sqrt(stderrx^2 + stderry^2)
    }
    if (stderr < 10 * .Machine$double.eps * max(abs(mx),
                                                abs(my)))
      stop("data are essentially constant")
    theta_hat <- (mx - my - mu)
    n = c(nx, ny)
  }

  ### permutation
  if (pvalue){
    if (method == 'Two Sample t-test'){
      all_values = c(x,y)
      set.seed(seed)
      perm_theta_hat = replicate(n.perm, {
        positions = sample(x=c(rep(1, length(x)), rep(2, length(y))), size=length(all_values), replace=FALSE)
        mean(all_values[positions==1]) -   mean(all_values[positions==2])
      })


    } else {   # if one sample t-test or paired t-test
      set.seed(seed)
      perm_theta_hat = replicate(n.perm, { mean(sign(runif(length(x), -1, 1)) * abs(x)) })
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
  } else pval = '-'





  ### bootstrap
  alpha <- 1 - conf.level
  if (alternative=='two.sided') sided=2 else sided=1
  set.seed(seed)
  if (confint){
    if (method == 'Two Sample t-test'){
      boot_theta_hat = replicate(n.perm, {
        mean(sample(x, length(x), replace=TRUE)) - mean(sample(y, length(y), replace=TRUE)) - mu
      })} else {
        boot_theta_hat = replicate(n.perm, {
          mean(sample(x, length(x), replace=TRUE)) - mu   })
      }


    quantiles = c(alpha/sided, 1-alpha/sided)

    if (boot.type == 'bca'){
      # ### bca (accelerated attempt)
      # ### https://stats.stackexchange.com/questions/437477/calculate-accelerated-bootstrap-interval-in-r
      # # compute constants
      # z0 = qnorm(mean(boot_theta_hat <= theta_hat))
      # zq = qnorm(quantiles)
      #
      # # accelerated constant
      # if (method == 'Two Sample t-test'){
      #   I <- rep(NA, nx+ny)
      #   for(i in 1:nx){
      #     xnew <- x[-i]  #Remove ith data point
      #
      #     delta_jack <- mean(xnew) - mean(y) - mu #Estimate delta
      #     I[i] <- (nx+ny-1)*(theta_hat - delta_jack)
      #   }
      #   for (i in 1:ny){
      #     ynew <- y[-i]  #Remove ith data point
      #     delta_jack <- mean(x) - mean(ynew) - mu #Estimate delta
      #     I[i+nx] <- (nx+ny-1)*(theta_hat - delta_jack)
      #   }
      #
      # } else {
      #   I <- rep(NA, n)
      #   for(i in 1:n){
      #     xnew <- x[-i]  #Remove ith data point
      #     delta_jack <- mean(xnew) - mu #Estimate delta
      #     I[i] <- (n-1)*(theta_hat - delta_jack)
      #   }
      # }
      #
      # #Estimate a
      # a <- (sum(I^3) / sum(I^2)^1.5) / 6
      #
      # # adjusted quantiles
      # quantiles = pnorm(z0 + (z0+zq)/(1-a*(z0+zq)))

      if (method == 'Two Sample t-test'){
        cint = boot.bca(x, y, boot.values=boot_theta_hat, theta_hat=theta_hat, paired=FALSE, quantiles=quantiles, fun=function(x,y) mean(x)-mean(y)-mu)
      } else {
        cint = boot.bca(x, boot.values=boot_theta_hat, theta_hat=theta_hat, paired=FALSE, quantiles=quantiles, fun=function(x) mean(x)-mu)
      }

    }
    # cint = quantile(boot_theta_hat, quantiles)
    if (alternative=='less') cint[1] = -Inf
    if (alternative=='greater') cint[2] = Inf
  }







  method = paste0(method,' (permutation & bootstrap)')

  names(theta_hat) <- "delta"
  names(n) <- "n"
  names(mu) <- if (paired)
    "mean difference"
  else if (!is.null(y))
    "difference in means"
  else "mean"

  RVAL <- list(statistic = theta_hat, df = n,
               estimate = estimate, null.value = mu,
               stderr = stderr, alternative = alternative, method = method,
               data.name = dname, boot.type = boot.type)
  if (pvalue)  RVAL <- c(RVAL, list(p.value = pval))
  if (boot.values) RVAL = c(RVAL, list(boot.values = boot_theta_hat))
  if (perm.values) RVAL = c(RVAL, list(perm.values = perm_theta_hat))
  if (confint){
    attr(cint, "conf.level") <- conf.level
    RVAL <- c(RVAL, list(conf.int = cint))
  }
  class(RVAL) <- "htest"
  RVAL
}
