
my.glm <- function(formula, data,
                   dist = "gaussian", link = "identity",
                   max.iter = 100, tol = 1e-8) {
  # -------------------------------------------------------------------------
  # 1. Build model frame & design matrix (handles any number of predictors)
  # -------------------------------------------------------------------------
  mf   <- model.frame(formula, data, na.action = na.pass)
  y    <- model.response(mf)                     # response vector
  X    <- model.matrix(attr(mf, "terms"), mf)    # includes intercept automatically

  # -------------------------------------------------------------------------
  # 2. Create the family object with the requested link
  # -------------------------------------------------------------------------
  fam <- switch(dist,
                gaussian = gaussian(),
                binomial = binomial(),
                poisson  = poisson(),
                gamma    = Gamma(),
                stop("Unsupported distribution"))
  fam$link <- make.link(link)

  linkinv   <- fam$link$linkinv      # μ = g⁻¹(η)
  mu.eta    <- fam$link$mu.eta       # dμ/dη
  variance  <- fam$variance          # V(μ)

  # -------------------------------------------------------------------------
  # 3. Initialise
  # -------------------------------------------------------------------------
  n <- length(y)
  p <- ncol(X)
  beta <- rep(0, p)               # start at zero (or use lm fit for faster start)
  eta  <- X %*% beta
  mu   <- linkinv(eta)

  # -------------------------------------------------------------------------
  # 4. IR‑IRLS loop
  # -------------------------------------------------------------------------
  for (iter in seq_len(max.iter)) {
    mu_eta   <- mu.eta(eta)                     # dμ/dη
    var_mu   <- variance(mu)                   # V(μ)
    w        <- (mu_eta^2) / var_mu             # weights
    z        <- eta + (y - mu) / mu_eta         # working response

    WX <- sweep(X, 1L, sqrt(w), FUN = "*")   # multiply each row of X by sqrt(w)
    Wz <- sqrt(w) * z                        # z is already a vector
    beta_new <- solve(t(WX) %*% WX, t(WX) %*% Wz)

    if (max(abs(beta_new - beta)) < tol) break
    beta <- beta_new
    eta  <- X %*% beta
    mu   <- linkinv(eta)
  }

  # -------------------------------------------------------------------------
  # 5. Return results
  # -------------------------------------------------------------------------
  list(coefficients = beta,
       iterations   = iter,
       converged    = iter < max.iter,
       formula      = formula)
}
