
cor.test.r = function (x, y, alternative = c("two.sided", "less", "greater"),
          method = c("pearson", "kendall", "spearman"), exact = NULL,
          conf.level = 0.95, continuity = FALSE, seed=0, n.perm=10000, confint=TRUE, pvalue=TRUE, boot.type = c("bca", "percentile"), boot.values = F, perm.values = F,
          ...)
{
  boot.type <- match.arg(boot.type)
  correlation = function(x,y) cov(x,y) / (sd(x)*sd(y))
  alternative <- match.arg(alternative)
  method <- match.arg(method)
  DNAME <- paste(deparse1(substitute(x)), "and", deparse1(substitute(y)))
  if (!is.numeric(x))
    stop("'x' must be a numeric vector")
  if (!is.numeric(y))
    stop("'y' must be a numeric vector")
  if (length(x) != length(y))
    stop("'x' and 'y' must have the same length")
  OK <- complete.cases(x, y)
  x <- x[OK]
  y <- y[OK]
  n <- length(x)
  NVAL <- 0
  conf.int <- FALSE

# Pearson -----------------------------------------------------------------

  if (method == "pearson") {
    if (n < 3L)
      stop("not enough finite observations")
    method <- "Pearson's product-moment correlation"
    names(NVAL) <- "correlation"
    theta_hat <- cor(x, y)
    df <- n - 2L
    ESTIMATE <- c(cor = theta_hat)
    PARAMETER <- c(df = df)
    STATISTIC <- c(t = sqrt(df) * theta_hat/sqrt(1 - theta_hat^2))
    if (n > 3) {
      if (!is.null(conf.level)){
        if (!missing(conf.level) && (length(conf.level) !=
                                     1 || !is.finite(conf.level) || conf.level < 0 ||
                                     conf.level > 1))
          stop("'conf.level' must be a single number between 0 and 1")
        conf.int <- TRUE


        ### bootstrap
        alpha <- 1 - conf.level
        if (confint){
          set.seed(seed)
          boot_theta_hat = replicate(n.perm, {
            positions = sample(1:n, n, replace=TRUE)
            correlation(x[positions], y[positions])
          })


          if (alternative=='two.sided') sided=2 else sided=1
          quantiles = c(alpha/sided, 1-alpha/sided)

          if (boot.type == 'bca'){
            # z0 = qnorm(mean(boot_theta_hat <= theta_hat))
            # zq = qnorm(quantiles)
            #
            # I <- rep(NA, n)
            # for(i in 1:n){
            #   xnew <- x[-i]
            #   ynew <- y[-i]
            #   r_jack <-  cor(xnew, ynew)
            #   I[i] <- (n-1)*(theta_hat - r_jack)
            # }
            # a <- (sum(I^3) / sum(I^2)^1.5) / 6
            #
            # quantiles = pnorm(z0 + (z0+zq)/(1-a*(z0+zq)))
            cint = boot.bca(x, y, boot.values=boot_theta_hat, theta_hat=theta_hat, paired=TRUE, quantiles=quantiles, fun=function(x,y) cor(x,y, method="pearson"))
          }

          #cint = quantile(boot_theta_hat, quantiles)
          if (alternative=='less') cint[1] = -1
          if (alternative=='greater') cint[2] = 1
          attr(cint, "conf.level") <- conf.level
        }
      }
    }


    ### permutation
    if (pvalue){
      set.seed(seed)
      perm.cor = replicate(n.perm, {
        positions_x = sample(1:n, n, replace=FALSE)
        positions_y = sample(1:n, n, replace=FALSE)
        correlation(x[positions_x], y[positions_y])
      })

      # calculates p-value

      if (alternative == 'two.sided'){
        temp_quantile = ecdf(perm.cor) (theta_hat)
        temp_quantile = ifelse(temp_quantile>.5, 1-temp_quantile, temp_quantile)
        PVAL = temp_quantile*2

      } else if (alternative == "less") {
        temp_quantile = ecdf(perm.cor) (theta_hat)
        temp_quantile = ifelse(temp_quantile<.5, 1-temp_quantile, temp_quantile)
        PVAL = ifelse(theta_hat<.0, 1-temp_quantile, temp_quantile)

      }  else if (alternative == "greater") {
        temp_quantile = ecdf(perm.cor) (theta_hat)
        temp_quantile = ifelse(temp_quantile>.5, 1-temp_quantile, temp_quantile)
        PVAL = ifelse(theta_hat<.0, 1-temp_quantile, temp_quantile)

      }
    }







  }   # end Pearson
  else {
    if (n < 2)
      stop("not enough finite observations")
    PARAMETER <- NULL

# Kendall -----------------------------------------------------------------

    if (method == "kendall") {
      method <- "Kendall's rank correlation tau"
      names(NVAL) <- "tau"
      theta_hat <- cor(x, y, method = "kendall")
      ESTIMATE <- c(tau = theta_hat)
      if (!is.finite(ESTIMATE)) {
        ESTIMATE[] <- NA
        STATISTIC <- c(T = NA)
        PVAL <- NA
      }
      else {

        if (!is.null(conf.level)){
          ### bootstrap
          conf.int = TRUE
          alpha <- 1 - conf.level
          set.seed(seed)
          boot_theta_hat = replicate(n.perm, {
            positions = sample(1:n, n, replace=TRUE)
            cor(x[positions], y[positions], method='kendall', use='complete.obs')
          })

          if (alternative=='two.sided') sided=2 else sided=1
          quantiles = c(alpha/sided, 1-alpha/sided)

          if (boot.type == 'bca'){
            # z0 = qnorm(mean(boot_theta_hat <= theta_hat))
            # zq = qnorm(quantiles)
            #
            # I <- rep(NA, n)
            # for(i in 1:n){
            #   xnew <- x[-i]
            #   ynew <- y[-i]
            #   r_jack <-  cor(xnew, ynew, method = "kendall")
            #   I[i] <- (n-1)*(theta_hat - r_jack)
            # }
            # a <- (sum(I^3) / sum(I^2)^1.5) / 6
            #
            # quantiles = pnorm(z0 + (z0+zq)/(1-a*(z0+zq)))

            cint = boot.bca(x, y, boot.values=boot_theta_hat, theta_hat=theta_hat, paired=TRUE, quantiles=quantiles, fun=function(x,y) cor(x, y, method="kendall"))
          }

          #cint = quantile(boot_theta_hat, quantiles)
          if (alternative=='less') cint[1] = -1
          if (alternative=='greater') cint[2] = 1
          attr(cint, "conf.level") <- conf.level
        }



      ### permutation
      set.seed(seed)
      perm.cor = replicate(n.perm, {
        positions_x = sample(1:n, n, replace=FALSE)
        positions_y = sample(1:n, n, replace=FALSE)
        cor(x[positions_x], y[positions_y], method='kendall')
      })

      # calculates p-value

      if (alternative == 'two.sided'){
        temp_quantile = ecdf(perm.cor) (theta_hat)
        temp_quantile = ifelse(temp_quantile>.5, 1-temp_quantile, temp_quantile)
        PVAL = temp_quantile*2

      } else if (alternative == "less") {
        temp_quantile = ecdf(perm.cor) (theta_hat)
        temp_quantile = ifelse(temp_quantile<.5, 1-temp_quantile, temp_quantile)
        PVAL = ifelse(theta_hat<.0, 1-temp_quantile, temp_quantile)

      }  else if (alternative == "greater") {
        temp_quantile = ecdf(perm.cor) (theta_hat)
        temp_quantile = ifelse(temp_quantile>.5, 1-temp_quantile, temp_quantile)
        PVAL = ifelse(theta_hat<.0, 1-temp_quantile, temp_quantile)

      }



      }
    }
    else {

# Spearman ----------------------------------------------------------------

      method <- "Spearman's rank correlation rho"
      if (is.null(exact))
        exact <- TRUE
      names(NVAL) <- "rho"
      theta_hat <- cor(rank(x), rank(y))
      ESTIMATE <- c(rho = theta_hat)
      if (!is.finite(ESTIMATE)) {
        ESTIMATE[] <- NA
        STATISTIC <- c(S = NA)
        PVAL <- NA
      }
      else {


        if (!is.null(conf.level)){
          ### bootstrap
          conf.int = TRUE
          alpha <- 1 - conf.level
          set.seed(seed)
          boot_theta_hat = replicate(n.perm, {
            positions = sample(1:n, n, replace=TRUE)
            correlation(rank(x[positions]), rank(y[positions]))
          })

          if (alternative=='two.sided') sided=2 else sided=1
          quantiles = c(alpha/sided, 1-alpha/sided)

          if (boot.type == 'bca'){
            # z0 = qnorm(mean(boot_theta_hat <= theta_hat))
            # zq = qnorm(quantiles)
            #
            # I <- rep(NA, n)
            # for(i in 1:n){
            #   xnew <- x[-i]
            #   ynew <- y[-i]
            #   r_jack <-  cor(rank(xnew), rank(ynew))
            #   I[i] <- (n-1)*(theta_hat - r_jack)
            # }
            # a <- (sum(I^3) / sum(I^2)^1.5) / 6
            #
            # quantiles = pnorm(z0 + (z0+zq)/(1-a*(z0+zq)))
            cint = boot.bca(x, y, boot.values=boot_theta_hat, theta_hat=theta_hat, paired=TRUE, quantiles=quantiles, fun=function(x,y) cor(rank(x), rank(y), method="pearson"))
          }

          #cint = quantile(boot_theta_hat, quantiles)
          if (alternative=='less') cint[1] = -1
          if (alternative=='greater') cint[2] = 1
          attr(cint, "conf.level") <- conf.level
        }


      ### permutation
      set.seed(seed)
      perm.cor = replicate(n.perm, {
        positions_x = sample(1:n, n, replace=FALSE)
        positions_y = sample(1:n, n, replace=FALSE)
        correlation(rank(x[positions_x]), rank(y[positions_y]))
      })

      # calculates p-value

      if (alternative == 'two.sided'){
        temp_quantile = ecdf(perm.cor) (theta_hat)
        temp_quantile = ifelse(temp_quantile>.5, 1-temp_quantile, temp_quantile)
        PVAL = temp_quantile*2

      } else if (alternative == "less") {
        temp_quantile = ecdf(perm.cor) (theta_hat)
        temp_quantile = ifelse(temp_quantile<.5, 1-temp_quantile, temp_quantile)
        PVAL = ifelse(theta_hat<.0, 1-temp_quantile, temp_quantile)

      }  else if (alternative == "greater") {
        temp_quantile = ecdf(perm.cor) (theta_hat)
        temp_quantile = ifelse(temp_quantile>.5, 1-temp_quantile, temp_quantile)
        PVAL = ifelse(theta_hat<.0, 1-temp_quantile, temp_quantile)

      }




      }
    }  # end Spearman
  }



  statistic = c(theta_hat=theta_hat)
  PARAMETER = c(n=n) ; names(PARAMETER) = 'n'
  method = paste0(method,' (permutation & bootstrap)')

  RVAL <- list(statistic = statistic, parameter = PARAMETER,
               estimate = ESTIMATE, null.value = NVAL,
               alternative = alternative, method = method, data.name = DNAME, boot.type = boot.type)
  if (conf.int & confint) RVAL <- c(RVAL, list(conf.int = cint))
  if (boot.values) RVAL = c(RVAL, list(boot.values = boot_theta_hat))
  if (perm.values) RVAL = c(RVAL, list(perm.values = perm_theta_hat))
  if (pvalue) RVAL <- c(RVAL, list(p.value = as.numeric(PVAL)))
  class(RVAL) <- "htest"
  RVAL
}
