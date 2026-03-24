#' @export

lm_r = function(formula, data, conf.level = 0.95,
                seed=0, n.perm=10000, confint = TRUE, pvalue = TRUE, boot.type = c("bca", "percentile"), pvalue.type = c("CI.inversion", "permutation"), boot.values = F, perm.values = F,
                ...){

  boot.type <- match.arg(boot.type)
  pvalue.type <- match.arg(pvalue.type)
  alpha <- 1 - conf.level

  all_variables = as.character( attr(terms(formula), 'variables')[-1])
  IVs = all_variables[-1]
  namy = c('(Intercept)', IVs)
  dname = paste(deparse1(substitute(data)))
  data = data[complete.cases(data[, all_variables]), all_variables]
  mod = lm(formula, data)
  # theta_hat = data.frame(summary(mod)$coefficient[, 1])
  theta_hat = c(my.lm(formula, data)$coeff)
  formula = as.formula(formula)
  n = nrow(data)


  if (n < length(all_variables)) stop("The number of complete case observation must be greater or equal to the number of parameters being estimated.")



  if (confint){
    sided=2
    quantiles = c(alpha/sided, 1-alpha/sided)

    set.seed(seed)

    boot.values = t(as.data.frame(replicate(n.perm, {my.lm(formula, data[sample(n, n, replace=TRUE), ])$coeff})))
    output = boot.bca.lm(data = data, n = n, theta_hat = theta_hat, boot.values = boot.values, quantiles = quantiles, alpha = alpha, IVs = IVs)
    coefficients = data.frame(row.names = namy,
                              Estimate = theta_hat,
                              lower = as.numeric(output[output$quantiles==quantiles[1], 3:ncol(output)]),
                              upper = as.numeric(output[output$quantiles==quantiles[2], 3:ncol(output)]) )

  }


  if (pvalue && pvalue.type=='CI.inversion'){

    pval_precision = 1/n.perm
    alpha_seq = seq(1e-16, 1 - 1e-16, pval_precision)
    quantiles_seq = c(alpha_seq/sided, 1-alpha_seq/sided)

    cint_seq = boot.bca.lm(data = data, n = n, theta_hat = theta_hat, boot.values = boot.values, quantiles = quantiles_seq, alpha = alpha_seq, IVs = IVs)

    ll_seq = cint_seq[1:length(alpha_seq), ]
    ul_seq = cint_seq[(length(alpha_seq)+1): (2*length(alpha_seq)), ]

    pval = c()
    for (i in 1:length(all_variables)){
      pval[i] = suppressWarnings( alpha_seq[which(ll_seq[i+2]==min(na.omit(ll_seq[i+2][na.omit(ll_seq[i+2])>=0])) & !is.na(ll_seq[i+2]))][1] )
      if(length(pval[i])==0 | is.na(pval[i])) pval[i] = alpha_seq[which(ul_seq[i+2]==max(na.omit(ul_seq[i+2][na.omit(ul_seq[i+2])<=0])) & !is.na(ul_seq[i+2]))][1]
    }

    coefficients$p.value = pval

  }





  if (pvalue && pvalue.type=='permutation'){
    # see https://cran.r-project.org/web//packages//permuco/vignettes/permuco_tutorial.pdf
    # actually, for a given predictor, one should predict the outcome with the other predictors, and use those residuals for the p-value of the predictor of interest and repeat the process for the other predictors...
    temp_df = data
    set.seed(seed)
    perm_theta_hat = data.frame(matrix(ncol=length(all.vars(formula)), nrow=n.perm))
    for (i in 1:n.perm){
      temp_df[[formula[[2]]]] = data[[formula[[2]]]] + sample(mod$residuals, replace=FALSE)
      perm_theta_hat[i, ] = my.lm(formula, temp_df)$coeff[, 1]
    }

    pval = rep(NA, length(all_variables))
    for (i in 1:length(all_variables)){
      temp_quantile = ecdf(perm_theta_hat[, i]) (theta_hat[i])
      temp_quantile = ifelse(temp_quantile>.5, 1-temp_quantile, temp_quantile)
      pval[i] = temp_quantile*2
    }

    coefficients$p.value = pval

  } # end pvalue






  RVAL <- list(coefficients = apply(as.matrix(coefficients), 1:2, as.numeric), n = n,
               stderr = stderr,
               data.name = dname, boot.type = boot.type,
               call = match.call())
  # if (pvalue)  RVAL <- c(RVAL, list(p.value = pval))
  # if (boot.values) RVAL = c(RVAL, list(boot.values = boot_theta_hat))
  # if (perm.values) RVAL = c(RVAL, list(perm.values = perm_theta_hat))
  # if (confint){
  #   attr(cint, "conf.level") <- conf.level
  #   RVAL <- c(RVAL, list(conf.int = cint))
  # }

  class(RVAL) <- 'lm'
  RVAL

}






# n=100
# df = data.frame(x1=rnorm(n,0,1), x5=rnorm(n,0,1)); df$y = rnorm(n,10,20)+-1*df$x1+3*df$x5
# formula = as.formula('y~x1+x5')
# data = df; conf.level = 0.95; seed=0; n.perm=10000; confint = TRUE; pvalue = TRUE; boot.type = "bca"; boot.values = F; perm.values = F; pvalue.type = "CI.inversion"
#
# tt=lm.r(y~x1+x5, df, confint = T, n.perm=10000, pvalue = T)$coefficients; tt
# summary(lm(y~x1+x5, df))$coefficient
# confint(lm(y~x1+x5, df))
#
# tictoc::tic()
# tt = lm.r(y~x1+x5, df, confint = T, n.perm=10000, pvalue = T); tt
# tt = lm.r(y~x1+x5, df, confint = T, n.perm=10000, pvalue = T, pvalue.type='permutation'); tt
# tictoc::toc()
#
# # mody = lm(y~x1+x5, df)
# # summary(mody)$coefficients ; confint(mody)
#
# n=100
# df = data.frame(y=rnorm(n,0,5), x=rep(0:1, each=(n/2)))
# df$y = df$y + df$x
#
# tt=t.test.r(df$y[df$x==1], df$y[df$x==0])
# tt
# lm.r(y~x, df, confint = T, n.perm=10000, pvalue = F)$coefficients
# lm.r(y~x, df, confint = T, n.perm=10000, pvalue = F)





