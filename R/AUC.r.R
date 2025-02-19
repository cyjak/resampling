#' @export



AUC.r = function(formula, data, alternative = c("two.sided", "less", "greater"), conf.level = 0.95,
                 seed=0, n.perm=10000, confint = TRUE, pvalue = TRUE, boot.type = c("bca", "percentile"), pvalue.type = c("CI.inversion", "permutation"), boot.values = F, perm.values = F,
                 ...){

  boot.type <- match.arg(boot.type)
  alternative <- match.arg(alternative)
  pvalue.type <- match.arg(pvalue.type)
  if (pvalue && pvalue.type=='CI.inversion' && alternative!='two.sided') stop("The CI.inversion method for computing the p-value is only availble for two-sided tests. Choose instead the 'permutation' method.")

  quick.auc = function(continuous, binary){

    tempy = data.frame(continuous, binary)
    tempy = tempy[order(tempy$continuous, decreasing = TRUE), ]
    continuous = tempy$continuous
    binary = tempy$binary

    modas = unique(na.omit(binary))
    TN = cumsum(binary == modas[1])
    TP = cumsum(binary == modas[2])

    FP = sum(binary == modas[1]) - TN
    FN = sum(binary == modas[2]) - TP

    sensi = TP / sum(binary == modas[1])
    speci = TN / sum(binary == modas[2])

    auc = sum(diff(speci) * (sensi[-1] + sensi[-length(sensi)])) / 2
    return(auc)
  }


  continuous = data[[formula[[2]]]]
  binary = data[[formula[[3]]]]
  temp = data.frame(continuous, binary)
  temp = temp[complete.cases(temp), ]
  continuous = temp$continuous
  binary = temp$binary
  theta_hat = quick.auc(temp$continuous, temp$binary)
  n = nrow(temp)


  ### permutation
  if (pvalue && pvalue.type=='permutation'){
    set.seed(seed)
    perm_theta_hat = replicate(n.perm, {
      positions_x = sample(1:n, n, replace=FALSE)
      positions_y = sample(1:n, n, replace=FALSE)
      quick.auc(continuous = continuous[positions_x], binary = binary[positions_y])
    })
    temp_quantile = ecdf(perm_theta_hat) (theta_hat)
    temp_quantile = ifelse(temp_quantile>.5, 1-temp_quantile, temp_quantile)
    pval = temp_quantile*2
  }







  ### bootstrap
  alpha <- 1 - conf.level
  if (confint){
    set.seed(seed)
    boot_theta_hat = replicate(n.perm, {
      positions = sample(1:n, n, replace=TRUE)
      quick.auc(continuous = continuous[positions], binary = binary[positions])
    })


    if (alternative=='two.sided') sided=2 else sided=1
    quantiles = c(alpha/sided, 1-alpha/sided)


    if (boot.type == 'bca'){

      pval_precision = 1/n.perm
      alpha_seq = seq(1e-16, 1 - 1e-16, pval_precision)
      quantiles_seq = c(alpha_seq/sided, 1-alpha_seq/sided)

      x = continuous ; y = binary
      cint = boot.bca(x, y, boot.values=boot_theta_hat, theta_hat=theta_hat, paired=TRUE, quantiles=quantiles, fun=function(x,y) quick.auc(continuous = x, binary = y))

      if (pvalue && pvalue.type=='CI.inversion'){
        cint_seq = boot.bca(x, y, boot.values=boot_theta_hat, theta_hat=theta_hat, paired=TRUE, quantiles=quantiles_seq, fun=function(x,y) quick.auc(continuous = x, binary = y))

        ll_seq = cint_seq[1:length(alpha_seq)]
        ul_seq = cint_seq[(length(alpha_seq)+1): (2*length(alpha_seq))]

        pval = NULL
        pval = suppressWarnings( alpha_seq[which(ll_seq==min(na.omit(ll_seq[na.omit(ll_seq)>=.5])) & !is.na(ll_seq))][1] )
        if(length(pval)==0 | is.na(pval)) pval = alpha_seq[which(ul_seq==max(na.omit(ul_seq[na.omit(ul_seq)<=.5])) & !is.na(ul_seq))][1]
      }


    }


    cint = quantile(boot_theta_hat, quantiles)
    if (alternative=='less') cint[1] = -1
    if (alternative=='greater') cint[2] = 1
    attr(cint, "conf.level") <- conf.level
  }



  if (theta_hat < .5){
    theta_hat = 1 - theta_hat
    if (confint){
      cint = c(1-cint[2], 1-cint[1])
    }
  }



  method_label = 'AUC' ; null.value = 0.5
  if (!pvalue) { method_label = paste0(method_label,' (bootstrap)')
  } else if (pvalue && pvalue.type=='permutation'){ method_label = paste0(method_label,' (permutation & bootstrap)')
  } else if (pvalue && pvalue.type=='CI.inversion'){ method_label = paste0(method_label,' (bootstrap & CI inversion)')}

  data.name = paste(deparse1(substitute(data)))

  RVAL <- list(
               estimate = theta_hat, null.value = null.value, method = method_label, data.name = data.name
               )

  if (pvalue) RVAL = c(RVAL, list(p.value = pval))
  if (confint) RVAL = c(RVAL, list(boot.type = boot.type,
                                   conf.int = cint))
  if (boot.values) RVAL = c(RVAL, list(boot.values = boot_theta_hat))
  if (perm.values) RVAL = c(RVAL, list(perm.values = perm_theta_hat))


  class(RVAL) = "htest"
  RVAL


}



# set.seed(1)
# df = data.frame(cbind(
#   outy = c(rnorm(50, 0, 1), rnorm(50, 1, 10)),
#   groupy = rep(0:1, each=50)
# ))
# AUC.r(outy ~ groupy, df, boot.type = 'bca', pvalue=T, confint=T)
#
#
#
# pROC::auc(pROC::roc(df$groupy, df$outy))
# brunnermunzel::brunnermunzel.test(df$outy[df$groupy==0], df$outy[df$groupy==1])

















