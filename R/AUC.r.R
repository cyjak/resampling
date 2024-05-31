#' @export
# set.seed(1)
# df = data.frame(cbind(
#   outy = c(rnorm(50, 0, 1), rnorm(50, 1, 10)),
#   groupy = rep(0:1, each=50)
# ))
#sdfasefsfy

AUC.r = function(formula, data, alternative = c("two.sided", "less", "greater"), pval = TRUE, conf.level = 0.95,
                 seed=0, n.perm=10000, confint = TRUE, pvalue = TRUE, boot.type = c("bca", "percentile"), boot.values = F, perm.values = F,
                 ...){

  boot.type <- match.arg(boot.type)
  alternative <- match.arg(alternative)

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
  if (pvalue){
    set.seed(seed)
    perm_theta_hat = replicate(n.perm, {
      positions_x = sample(1:n, n, replace=FALSE)
      positions_y = sample(1:n, n, replace=FALSE)
      quick.auc(continuous = continuous[positions_x], binary = binary[positions_y])
    })
    temp_quantile = ecdf(perm_theta_hat) (theta_hat)
    temp_quantile = ifelse(temp_quantile>.5, 1-temp_quantile, temp_quantile)
    p.value = temp_quantile*2
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
      # z0 = qnorm(mean(boot_theta_hat < theta_hat))
      # zq = qnorm(quantiles)
      #
      # I <- rep(NA, n)
      # for(i in 1:n){
      #   contin_new <- continuous[-i]
      #   binary_new <- binary[-i]
      #   auc_jack <-  quick.auc(continuous = contin_new, binary = binary_new)
      #   I[i] <- (n-1)*(theta_hat - auc_jack)
      # }
      # a <- (sum(I^3) / sum(I^2)^1.5) / 6
      #
      # quantiles = pnorm(z0 + (z0+zq)/(1-a*(z0+zq)))

      x = continuous ; y = binary
      cint = boot.bca(x, y, boot.values=boot_theta_hat, theta_hat=theta_hat, paired=TRUE, quantiles=quantiles, fun=function(x,y) quick.auc(continuous = x, binary = y))

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



  method = 'AUC (permutation & bootstrap)' ; null.value = 0.5
  data.name = paste(deparse1(substitute(data)))

  RVAL <- list(
               estimate = theta_hat, null.value = null.value, method = method, data.name = data.name
               )

  if (pval) RVAL = c(RVAL, list(p.value = p.value))
  if (confint) RVAL = c(RVAL, list(boot.type = boot.type,
                                   conf.int = cint))
  if (boot.values) RVAL = c(RVAL, list(boot.values = boot_theta_hat))
  if (perm.values) RVAL = c(RVAL, list(perm.values = perm_theta_hat))


  class(RVAL) = "htest"
  RVAL


}


# AUC.r(outy ~ groupy, df, boot.type = 'bca', pval=T, confint=T)
#
#
#
# pROC::auc(pROC::roc(df$groupy, df$outy))
# brunnermunzel::brunnermunzel.test(df$outy[df$groupy==0], df$outy[df$groupy==1])

















