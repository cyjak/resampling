# https://stackoverflow.com/questions/38601794/fit-many-formulae-at-once-faster-options-than-lapply/38604244#38604244

# performs fast linear regression

## linear model estimation based on pivoted Cholesky factorization with Jacobi preconditioner
lm.chol <- function(formula, data) {
  ## stage0: get response vector and model matrix
  ## we did not follow the normal route: match.call, model.frame, model.response, model matrix, etc
  y <- data[[as.character(formula[[2]])]]
  X <- model.matrix(formula, data)
  n <- nrow(X); p <- ncol(X)
  ## stage 1: XtX and Jacobi diagonal preconditioner
  XtX <- crossprod(X)
  D <- 1 / sqrt(diag(XtX))
  ## stage 2: pivoted Cholesky factorization
  R <- suppressWarnings(chol(t(D * t(D * XtX)), pivot = TRUE))
  piv <- attr(R, "pivot")
  r <- attr(R, "rank")
  if (r < p) {
    warning("Model is rank-deficient!")
    piv <- piv[1:r]
    R <- R[1:r, 1:r]
  }
  ## stage 3: solve linear system for coefficients
  D <- D[piv]
  b <- D * crossprod(X, y)[piv]
  z <- forwardsolve(t(R), b)
  RSS <- sum(y * y) - sum(z * z)
  sigma <- sqrt(RSS / (n - r))
  para <- D * backsolve(R, z)
  beta.hat <- rep(NA, p)
  beta.hat[piv] <- para
  ## stage 4: get standard error
  Rinv <- backsolve(R, diag(r))
  se <- rep(NA, p)
  se[piv] <- D * sqrt(rowSums(Rinv * Rinv)) * sigma
  ## stage 5: t-statistic and p-value
  t.statistic <- beta.hat / se
  p.value <- 2 * pt(-abs(t.statistic), df = n - r)
  ## stage 6: construct coefficient summary matrix
  coefficients <- matrix(c(beta.hat, se, t.statistic, p.value), ncol = 4L)
  # colnames(coefficients) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  rownames(coefficients) <- colnames(X)


  # ## stage 7: compute adjusted R.squared
  # adj.R2 <- 1 - sigma * sigma / var(y)
  # ## return model fitting results
  # attr(coefficients, "sigma") <- sigma
  # attr(coefficients, "adj.R2") <- adj.R2


  coefficients = data.frame(coefficients)
  names(coefficients) = c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  coefficients
}


