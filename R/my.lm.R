
my.lm = function(formula, data, estimator=T, SE=F, R2=F){

  IV =  as.character(attr(terms(formula), 'variables')[-(1:2)])
  DV = as.character(attr(terms(formula), 'variables')[2])

  X = as.matrix(cbind('intercept' = rep(1, times=nrow(data)), data[,IV]), ncol=length(IV))


  XtX_inv = solve(t(X) %*% X)

  # coefficients
  coeff = NA
  if (estimator==T || SE==T){
    coeff = XtX_inv %*% t(X) %*% as.matrix(data[[DV]], ncol=1)  # estimators
  }


  # Standard error
  seOLS=NA
  if (SE==T){
    residus = data[[DV]] - X%*%coeff
    varcovOLS = sum(residus^2) / (nrow(data)-2) * XtX_inv
    seOLS = sqrt(diag(varcovOLS))
  }



  # R2
  if (R2==T){
    H = X %*% XtX_inv %*% t(X)
    idn = diag(rep(1, nrow(data)))
    RSS = t(as.matrix(data[[DV]])) %*% (idn - H) %*% as.matrix(data[[DV]])   # RSS

    H0 = matrix(X[,1]) %*% solve(t(X[,1]) %*% X[,1]) %*% t(matrix(X[,1]))
    TSS= t(as.matrix(data[[DV]])) %*% (idn - H0) %*% as.matrix(data[[DV]])   # TSS

    R2 = 1 - RSS/TSS
  }


  return(list(
    coeff = coeff,
    SE = seOLS,
    R2 = R2
  ))

}



