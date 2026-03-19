
### https://stats.stackexchange.com/questions/437477/calculate-accelerated-bootstrap-interval-in-r

boot.bca = function(x, y=NULL, theta_hat, boot.values, fun, paired, quantiles){


  if (is.null(y)){
    n = length(x)
    theta_jack = rep(NA, n)
    for(i in 1:n){
      xnew <- x[-i]
      theta_jack[i] <- fun(xnew)
    }


  } else if (!paired) { # if y is not null
    nx = length(x) ; ny = length(y)
    theta_jack = rep(NA, nx+ny)
    for(i in 1:nx){
      xnew <- x[-i]
      theta_jack[i] <- fun(xnew, y)
    }

    for (i in 1:ny){
      ynew <- y[-i]
      theta_jack[nx+i] <- fun(x, ynew)
    }


  } else {  # if paired
    n = length(x)
    theta_jack = rep(NA, n)
    for(i in 1:n){
      xnew <- x[-i]
      ynew <- y[-i]
      theta_jack[i] <- fun(xnew, ynew)
    }


  }


  I <- mean(theta_jack) - theta_jack
  a <- sum(I^3) / (6 * sum(I^2)^1.5)

  z0 = qnorm(mean(boot.values <= theta_hat))  # bias-correction
  zq = qnorm(quantiles)

  quantiles_corrected = pnorm(z0 + (z0+zq)/(1-a*(z0+zq)))
  quantile(boot.values, quantiles_corrected)
}















