

boot.bca.lm = function(x, y, theta_hat, boot.values, quantiles){

  # quantiles = c(.025, .975)
  z0 = qnorm(colMeans(boot.values <= theta_hat))
  zq = qnorm(quantiles)

  n = length(y)
  #I <- data.frame(matrix(nrow=n, ncol=ncol(x)+1))

  theta_jack = data.frame(matrix(nrow=n, ncol=ncol(x)+1))

  for(i in 1:n){
    xnew <- x[-i, ]
    ynew <- y[-i]
    temp_df = data.frame(cbind(ynew, xnew))
    theta_jack[i, ] <- lm.chol(ynew~., temp_df)[, 1]
  }
  I <- colMeans(theta_jack) - theta_jack
  a <- colSums(I^3) / (6 * colSums(I^2)^1.5)

  output = data.frame(matrix(nrow=ncol(x)+1, ncol=5))
  namy = c('(Intercept)', names(x))

  for (i in 1:(ncol(x)+1)){
    temp_quant = pnorm(z0[i] + (z0[i]+zq)/(1-a[i]*(z0[i]+zq)))
    temp_ci = quantile(boot.values[,i], temp_quant)
    output[i, ] = c(namy[i], temp_ci, temp_quant)
  }

  output

}








