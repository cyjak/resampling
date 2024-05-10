
### https://stats.stackexchange.com/questions/437477/calculate-accelerated-bootstrap-interval-in-r

boot.bca = function(x, y=NULL, theta_hat, boot.values, fun, paired, quantiles){

  #quantiles = c(.025, .975)
  z0 = qnorm(mean(boot.values <= theta_hat))
  zq = qnorm(quantiles)


  if (is.null(y)){

    n = length(x)
    I <- rep(NA, n)

    theta_jack = rep(NA, n)
    for(i in 1:n){
      xnew <- x[-i]
      theta_jack[i] <- fun(xnew)
    }
    I <- mean(theta_jack) - theta_jack
    a <- sum(I^3) / (6 * sum(I^2)^1.5)
    quantiles = pnorm(z0 + (z0+zq)/(1-a*(z0+zq)))
    quantile(boot.values, quantiles)


  } else if (!paired) { # if y is not null

    nx = length(x) ; ny = length(y)

    I <- rep(NA, nx+ny)

    theta_jack = rep(NA, nx+ny)
    for(i in 1:nx){
      xnew <- x[-i]
      theta_jack[i] <- mean(xnew) - mean(y)
    }

    for (i in 1:ny){
      ynew <- y[-i]
      theta_jack[nx+i] <- mean(x) - mean(ynew)
    }

    I <- mean(theta_jack) - theta_jack
    a <- sum(I^3) / (6 * sum(I^2)^1.5)
    quantiles = pnorm(z0 + (z0+zq)/(1-a*(z0+zq)))
    quantile(boot.values, quantiles)

  } else {  # if paired

    n = length(x)
    I <- rep(NA, n)

    theta_jack = rep(NA, n)
    for(i in 1:n){
      xnew <- x[-i]
      ynew <- y[-i]
      theta_jack[i] <- fun(xnew, ynew)
    }
    I <- mean(theta_jack) - theta_jack
    a <- sum(I^3) / (6 * sum(I^2)^1.5)
    quantiles = pnorm(z0 + (z0+zq)/(1-a*(z0+zq)))
    quantile(boot.values, quantiles)

  }
}















