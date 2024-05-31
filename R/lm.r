#' @export



lm.r = function(formula, data, conf.level = 0.95,
                seed=0, n.perm=10000, confint = TRUE, pvalue = TRUE, boot.type = c("bca", "percentile"), boot.values = F, perm.values = F,
                ...){

  boot.type <- match.arg(boot.type)

  alpha <- 1 - conf.level

  mod = lm(formula, data)

  coefficients = data.frame(summary(mod)$coefficient[,1:2])



  if (confint){

    sided=2
    quantiles = c(alpha/sided, 1-alpha/sided)

    set.seed(seed)
    func = function(data, indices){
      lm.chol(as.formula(formula), data[indices,])[, 1]
    }
    foo = boot::boot(data, R=n.perm, statistic=func)

    if (boot.type == 'bca'){
      cint = boot.bca.lm(x=data[, all.vars(formula)[2:length(all.vars(formula))]], y=data[[formula[[2]]]], theta_hat=coefficients[, 'Estimate'], boot.values=foo$t, quantiles=quantiles)
    }

    coefficients$lower = cint[, 2]
    coefficients$upper = cint[, 3]


  } # end confint




  if (pvalue){
    # see https://cran.r-project.org/web//packages//permuco/vignettes/permuco_tutorial.pdf
    # actually, for a given predictor, one should predict the outcome with the other predictors, and use those residuals for the p-value of the predictor of interest and repeat the process for the other predictors...
    temp_df = data
    set.seed(seed)
    perm_theta_hat = data.frame(matrix(ncol=length(all.vars(formula)), nrow=n.perm))
    for (i in 1:n.perm){
      temp_df[[formula[[2]]]] = data[[formula[[2]]]] + sample(mod$residuals, replace=FALSE)
      perm_theta_hat[i, ] = lm.chol(as.formula(formula), temp_df)[,'Estimate']
    }

    ecdf(as.matrix(perm_theta_hat+theta_hat)) (theta_hat)

    temp_quantile = ecdf(as.matrix(perm_theta_hat)) (theta_hat)
    temp_quantile = ifelse(temp_quantile>.5, 1-temp_quantile, temp_quantile)
    pval = temp_quantile*2

    coefficients$p.value = pval

  } # end pvalue


  coefficients$'Std..Error' = NULL
  coefficients

}


df = data.frame(x1=rnorm(100,0,1), x5=rnorm(100,0,1)); df$y = rnorm(100,10,20)+df$x1+10*df$x5
lm.r(y~x1+x5, df, confint = T, n.perm=10000, pvalue = F)

mody = lm(y~x1+x5, df)
summary(mody)$coefficients ; confint(mody)





