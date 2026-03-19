
my.lm = function(formula, data){

  IV =  as.character(attr(terms(formula), 'variables')[-(1:2)])
  DV = as.character(attr(terms(formula), 'variables')[2])
  X = as.matrix(cbind('intercept' = rep(1, times=nrow(data)), data[,IV]), ncol=length(IV))

  XtX_inv = solve(t(X) %*% X)

  coeff = XtX_inv %*% t(X) %*% as.matrix(data[[DV]], ncol=1)

  return(list(
    coeff = coeff
  ))

}



