#' @export

chisq.test.r = function (x, y = NULL, correct = TRUE, p = rep(1/length(x), length(x)),
          rescale.p = FALSE, seed=0, n.perm=10000, perm.values = F, pvalue = TRUE,
          confint = TRUE, conf.level = 0.95, boot.type = c("bca"), boot.values = F)
{

  boot.type <- match.arg(boot.type)
  get.statistic = function(x, YATES){
    row_sums  = rowSums(x)
    col_sums  = colSums(x)
    total_obs = sum(x)
    E = outer(row_sums, col_sums) / total_obs
    sum((abs(x - E) - YATES)^2/E)
  }


  DNAME <- deparse(substitute(x))
  if (is.data.frame(x))
    x <- as.matrix(x)
  if (is.matrix(x)) {
    if (min(dim(x)) == 1L)
      x <- as.vector(x)
  }

  if (!is.matrix(x) && !is.null(y)) {
    if (length(x) != length(y))
      stop("'x' and 'y' must have the same length")
    DNAME2 <- deparse(substitute(y))
    xname <- if (length(DNAME) > 1L || nchar(DNAME, "w") >
                 30){
      ""
    } else DNAME
    yname <- if (length(DNAME2) > 1L || nchar(DNAME2, "w") >
                 30){
      ""
    } else DNAME2
    OK <- complete.cases(x, y)
    xok <- factor(x[OK])
    yok <- factor(y[OK])
    if ((nlevels(xok) < 2L) || (nlevels(yok) < 2L))
      stop("'x' and 'y' must have at least 2 levels")
    x <- table(xok, yok)
    names(dimnames(x)) <- c(xname, yname)
    DNAME <- paste(paste(DNAME,  collapse = "\n"), "and",
                   paste(DNAME2, collapse = "\n"))
  }

  if (any(x < 0) || anyNA(x))
    stop("all entries of 'x' must be nonnegative and finite")
  if ((n <- sum(x)) == 0)
    stop("at least one entry of 'x' must be positive")


  if (is.matrix(x)) {
    METHOD <- "Pearson's Chi-squared test"
    nr <- as.integer(nrow(x))
    nc <- as.integer(ncol(x))
    if (is.na(nr) || is.na(nc) || is.na(nr * nc)){
      stop("invalid nrow(x) or ncol(x)", domain = NA)}




  } else {
    if (length(dim(x)) > 2L)
      stop("invalid 'x'")
    if (length(x) == 1L)
      stop("'x' must at least have 2 elements")
    if (length(x) != length(p))
      stop("'x' and 'p' must have the same number of elements")
    if (any(p < 0))
      stop("probabilities must be non-negative.")
    if (abs(sum(p) - 1) > sqrt(.Machine$double.eps)) {
      if (rescale.p){
        p <- p/sum(p)
      } else stop("probabilities must sum to 1.")
    }
    x = as.matrix(table(x,y))
  }


  row_sums  = rowSums(x)
  col_sums  = colSums(x)
  total_obs = sum(x)
  E = outer(row_sums, col_sums) / total_obs

  if (is.null(colnames(x))) colnames(x) = 1:ncol(x)
  if (is.null(rownames(x))) rownames(x) = 1:nrow(x)

  PARAMETER <- NA
  PVAL <- 1



  if (correct && nrow(x) >= 2L && ncol(x) >= 2L) {
    YATES <- min(0.5, abs(x - E))
    if (YATES > 0)
      METHOD <- paste(METHOD, "with Yates' continuity correction")

    }  else YATES <- 0
  STATISTIC <- sum((abs(x - E) - YATES)^2/E)
  PARAMETER <- (nr - 1L) * (nc - 1L)






  ### permutation
  if (pvalue){
    set.seed(seed)
    x1 = c() ; y1 = c()
    for (col in 1:ncol(x))    x1 = c(x1, rep(colnames(x)[col], col_sums[col]))
    for (row in 1:nrow(x))    y1 = c(y1, rep(rownames(x)[row], row_sums[row]))

    perm_theta_hat = replicate(n.perm, {
      positions_x = sample(1:total_obs, total_obs, replace=FALSE)
      positions_y = sample(1:total_obs, total_obs, replace=FALSE)
      get.statistic(as.matrix(table(x1[positions_x], y1[positions_y])), YATES=YATES)
    })

    temp_quantile = ecdf(-abs(perm_theta_hat)) (-abs(STATISTIC))
    PVAL = temp_quantile
  }




  p1 = x[1,1] / sum(x[1,])
  p2 = x[2,1] / sum(x[2,])
  theta_hat <- p1 - p2
  estimate <- c(p2, p1, theta_hat)
  names(estimate) <- c("proportion A", "proportion B", "delta")




  ### bootstrap
  set.seed(seed)
  confint_present = FALSE
  if (sum(dim(x))==4){ # if 2x2 table
    if (confint){

      confint_present = TRUE
      if(!is.null(y)){
        temp_df = data.frame(x1=as.numeric(as.factor(xok))-1, x2=as.numeric(as.factor(yok))-1)

      } else { # if input is a table, we recreate the full dataset

        temp_df = data.frame(matrix(ncol=2) )

        for (position_x in 1:2){
          for (position_y in 1:2){
            nb = x[position_x, position_y]
            if (nb > 0){
              for (i in 1:nb){
                temp_df = rbind(temp_df, c(position_x-1, position_y-1))
              }

            }
          }
        }
         temp_df = temp_df[2:(sum(x)+1) , ]
      }



       alpha <- 1 - conf.level
       #if (alternative=='two.sided') sided=2 else sided=1
       sided = 2
       quantiles = c(alpha/sided, 1-alpha/sided)


       nb_row = nrow(temp_df)

       set.seed(seed)
       boot_theta_hat = replicate(n.perm, {
         tempy <- temp_df[sample(1:nb_row, replace=TRUE), ]
         sum(tempy[tempy[,1]==1, 2]) / length(tempy[tempy[,1]==1, 2]) - sum(tempy[tempy[,1]==0, 2]) / length(tempy[tempy[,1]==0, 2])
       })

       x_table = x
       if (boot.type == 'bca'){
         x = temp_df[, 1]
         y = temp_df[, 2]
          cint = boot.bca(x, y, boot.values=boot_theta_hat, theta_hat=theta_hat, paired=TRUE, quantiles=quantiles, fun=function(x,y) sum(x[y==1]) / length(x[y==1]) - sum(x[y==0]) / length(x[y==0]) )
        }
        x = x_table





    }
  }





  PARAMETER = total_obs
  METHOD = paste0(METHOD,' permutation & bootstrap')
  names(STATISTIC) <- "X-squared"
  names(PARAMETER) <- "n"
  if (any(E < 5) && is.finite(PARAMETER))
    warning("Chi-squared approximation may be incorrect")

  RVAL = list(parameter = PARAMETER, estimate = estimate,
                 method = METHOD, data.name = DNAME, observed = x,
                 expected = E, residuals = (x - E)/sqrt(E))

  if (pvalue) RVAL = c(RVAL, list(statistic = STATISTIC, p.value = PVAL))

  if (perm.values) RVAL = c(RVAL, list(perm.values = perm_theta_hat))

  if (confint_present) RVAL = c(RVAL, list(conf.int = cint))
  if (confint_present & boot.values) RVAL = c(RVAL, list(boot.values = boot_theta_hat))

  class(RVAL) <- "htest"
  RVAL
}




#
# x = data.frame(x1=c(20, 80), x2=c(50,50))
#
# #chisq.test.r(x)
# set.seed(6)
# x=c(rbinom(50,1,.2), rbinom(50,1,.1))
# y=rep(0:1, each=50)
# chisq.test.r(x, y, pvalue=T, confint=T, correct = F)
# chisq.test(x, y, correct=F)
#
#
# tt = chisq.test.r(x=data.frame(x1=c(20, 80), x2=c(50,50)), pvalue=T, boot.values = T)
# tt
#
#
# tt$parameter
#
# x = as.matrix(table(x,y))





