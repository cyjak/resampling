#' @export

cor_test_r = function (x, y, alternative = c("two.sided", "less", "greater"),
          method = c("pearson", "kendall", "spearman"),
          conf.level = 0.95, seed=0, n.perm=10000, confint=TRUE, pvalue=TRUE, boot.type = c("bca", "percentile"), pvalue.type = c("CI.inversion", "permutation"), boot.values = F, perm.values = F,
          ...)
{
  boot.type <- match.arg(boot.type)
  pvalue.type <- match.arg(pvalue.type)
  alternative <- match.arg(alternative)
  method <- match.arg(method)
  if (pvalue && pvalue.type=='CI.inversion' && alternative!='two.sided') stop("The CI.inversion method for computing the p-value is only availble for two-sided tests. Choose instead the 'permutation' method.")

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







  if (n < 3L)      stop("not enough finite observations")


  names(NVAL) <- "correlation"
  theta_hat <- cor(x, y, method=method)
  df <- n - 2L
  ESTIMATE <- c(cor = theta_hat)
  PARAMETER <- c(df = df)
  STATISTIC <- c(t = sqrt(df) * theta_hat/sqrt(1 - theta_hat^2))
  if (n > 3) {
    if (!is.null(conf.level)){
      if (!missing(conf.level) && (length(conf.level) !=1 || !is.finite(conf.level) || conf.level < 0 || conf.level > 1))
        stop("'conf.level' must be a single number between 0 and 1")
      conf.int <- TRUE


      ### bootstrap
      alpha <- 1 - conf.level
      if (confint){
        set.seed(seed)
        boot_theta_hat = replicate(n.perm, {
          positions = sample(1:n, n, replace=TRUE)
          cor(x[positions], y[positions], method=method)
        })


        if (alternative=='two.sided') sided=2 else sided=1
        quantiles = c(alpha/sided, 1-alpha/sided)


        if (boot.type == 'bca'){
          pval_precision = 1/n.perm
          alpha_seq = seq(1e-16, 1 - 1e-16, pval_precision)
          quantiles_seq = c(alpha_seq/sided, 1-alpha_seq/sided)

          cint = boot.bca(x, y, boot.values=boot_theta_hat, theta_hat=theta_hat, paired=TRUE, quantiles=quantiles, fun=function(x,y) cor(x,y, method=method))

          if (pvalue && pvalue.type=='CI.inversion'){
            cint_seq = boot.bca(x, y, boot.values=boot_theta_hat, theta_hat=theta_hat, paired=TRUE, quantiles=quantiles_seq, fun=function(x,y) cor(x,y, method=method))
          }


          if (pvalue && pvalue.type=='CI.inversion'){
            ll_seq = cint_seq[1:length(alpha_seq)]
            ul_seq = cint_seq[(length(alpha_seq)+1): (2*length(alpha_seq))]

            pval = NULL
            pval = suppressWarnings( alpha_seq[which(ll_seq==min(na.omit(ll_seq[na.omit(ll_seq)>=0])) & !is.na(ll_seq))][1] )
            if(length(pval)==0 | is.na(pval)) pval = alpha_seq[which(ul_seq==max(na.omit(ul_seq[na.omit(ul_seq)<=0])) & !is.na(ul_seq))][1]
          }

        }

        #cint = quantile(boot_theta_hat, quantiles)
        if (alternative=='less') cint[1] = -1
        if (alternative=='greater') cint[2] = 1
        attr(cint, "conf.level") <- conf.level
      }
    }
  }



  ### permutation
  if (pvalue && pvalue.type=='permutation'){
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
      pval = temp_quantile*2

    } else if (alternative == "less") {
      temp_quantile = ecdf(perm.cor) (theta_hat)
      temp_quantile = ifelse(temp_quantile<.5, 1-temp_quantile, temp_quantile)
      pval = ifelse(theta_hat<.0, 1-temp_quantile, temp_quantile)

    }  else if (alternative == "greater") {
      temp_quantile = ecdf(perm.cor) (theta_hat)
      temp_quantile = ifelse(temp_quantile>.5, 1-temp_quantile, temp_quantile)
      pval = ifelse(theta_hat<.0, 1-temp_quantile, temp_quantile)

    }
  }









  if (method=='pearson'){         method_label <- "Pearson's product-moment correlation"
  } else if (method=='spearman'){  method_label <- "Spearman's rank correlation rho"
  } else if (method=='kendall')   method_label <- "Kendall's rank correlation tau"

  if (!pvalue) { method_label = paste0(method_label,' (bootstrap)')
  } else if (pvalue && pvalue.type=='permutation'){ method_label = paste0(method_label,' (permutation & bootstrap)')
  } else if (pvalue && pvalue.type=='CI.inversion'){ method_label = paste0(method_label,' (bootstrap & CI inversion)')}

  statistic = c(theta_hat=theta_hat)
  PARAMETER = c(n=n) ; names(PARAMETER) = 'n'

  RVAL <- list(statistic = statistic, parameter = PARAMETER,
               estimate = ESTIMATE, null.value = NVAL,
               alternative = alternative, method = method_label, data.name = DNAME, boot.type = boot.type)
  if (conf.int & confint) RVAL <- c(RVAL, list(conf.int = cint))
  if (boot.values) RVAL = c(RVAL, list(boot.values = boot_theta_hat))
  if (perm.values) RVAL = c(RVAL, list(perm.values = perm_theta_hat))
  if (pvalue) RVAL <- c(RVAL, list(p.value = as.numeric(pval)))
  class(RVAL) <- "htest"
  RVAL
}



# set.seed(0);n=100
# y1 = rnorm(n, 0, 1)
# y2 = rnorm(n, 0, 1) + .1*y1
#
# cor.test.r(y1,y2, pvalue=T, confint=T, method='pearson')
# cor.test(y1,y2, method='kendall')
