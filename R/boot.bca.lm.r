

boot.bca.lm = function(formula, data, n, theta_hat, boot.values, quantiles, alpha, IVs){

  theta_jack = data.frame(matrix(nrow=n, ncol=ncol(data)))
  for(i in 1:n)      theta_jack[i, ] <- c(my.lm(formula, data[-i, ])$coeff)

  I <- sweep(theta_jack, 2, colMeans(theta_jack), FUN='-')
  a <- colSums(I^3) / (6 * colSums(I^2)^1.5)

  z0 = qnorm(colMeans(sweep(boot.values, 2, theta_hat, FUN='<=')))
  zq = qnorm(quantiles)

  # output = data.frame(matrix(nrow=ncol(data), ncol=5))
  output = data.frame(matrix(nrow=length(quantiles), ncol=ncol(data)+2))
  namy = c('(Intercept)', IVs)
  output[, 1] = quantiles
  output[, 2] = alpha
  for (i in 1:(ncol(data))){
    temp_quant = pnorm(z0[i] + (z0[i]+zq)/(1-a[i]*(z0[i]+zq)))
    temp_ci = quantile(boot.values[,i], temp_quant)

    output[, i+2] = temp_ci

    # output[i, ] = c(namy[i], temp_ci, temp_quant)
  }
  names(output) = c('quantiles','alpha',namy)

  # names(output) = c('variable','ll_value','ul_value','ll_quantile','ul_quantile')
  return(output)
}








