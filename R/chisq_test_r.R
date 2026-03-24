#' @export

chisq_test_r = function (x, y = NULL, correct = TRUE, p = rep(1/length(x), length(x)),
          rescale.p = FALSE, seed=0, n.perm=10000, perm.values = F, pvalue = TRUE,
          confint = TRUE, conf.level = 0.95, boot.type = c("bca"), pvalue.type = c("CI.inversion", "permutation"), boot.values = F)
{

  boot.type <- match.arg(boot.type)
  pvalue.type <- match.arg(pvalue.type)
  get.statistic = function(x, YATES){
    row_sums  = rowSums(x)
    col_sums  = colSums(x)
    total_obs = sum(x)
    E = outer(row_sums, col_sums) / total_obs
    sum((abs(x - E) - YATES)^2/E)
  }
  if (pvalue && pvalue.type=='CI.inversion' && alternative!='two.sided') stop("The CI.inversion method for computing the p-value is only availble for two-sided tests. Choose instead the 'permutation' method.")


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
    method_label <- "Pearson's Chi-squared test"
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
  pval <- 1



  if (correct && nrow(x) >= 2L && ncol(x) >= 2L) {
    YATES <- min(0.5, abs(x - E))
    if (YATES > 0)
      method_label <- paste(method_label, "with Yates' continuity correction")

    }  else YATES <- 0
  STATISTIC <- sum((abs(x - E) - YATES)^2/E)
  PARAMETER <- (nr - 1L) * (nc - 1L)






  ### permutation
  if (pvalue && pvalue.type=='permutation'){
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
    pval = temp_quantile
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

         x1 = xok[yok==0]
         x2 = xok[yok==1]


         boot_theta_hat_valid = boot_theta_hat[!is.nan(boot_theta_hat)]
         z0 = qnorm(mean(boot_theta_hat_valid <= theta_hat))
         zq = qnorm(quantiles)
         nx = length(x1) ; ny = length(x2)

         I <- rep(NA, nx+ny)

         theta_jack = rep(NA, nx+ny)
         for(i in 1:nx){
           xnew <- x1[-i]
           theta_jack[i] <- sum(xnew==1)/length(xnew) - sum(x2==1)/length(x2)
         }

         for (i in 1:ny){
           ynew <- x2[-i]
           theta_jack[nx+i] <- sum(ynew==1)/length(ynew) - sum(x1==1)/length(x1)
         }

         I <- mean(theta_jack) - theta_jack
         a <- sum(I^3) / (6 * sum(I^2)^1.5)
         quantiles = pnorm(z0 + (z0+zq)/(1-a*(z0+zq)))
         cint = quantile(boot_theta_hat_valid, quantiles)

         if (pvalue && pvalue.type=='CI.inversion'){
           pval_precision = 1/n.perm
           alpha_seq = seq(1e-16, 1 - 1e-16, pval_precision)
           quantiles_seq = c(alpha_seq/sided, 1-alpha_seq/sided)
           zq_seq = qnorm(quantiles_seq)
           quantiles_seq = pnorm(z0 + (z0+zq_seq)/(1-a*(z0+zq_seq)))

           cint_seq = quantile(boot_theta_hat_valid, quantiles_seq)


           ll_seq = cint_seq[1:length(alpha_seq)]
           ul_seq = cint_seq[(length(alpha_seq)+1): (2*length(alpha_seq))]

           pval = NULL
           pval = suppressWarnings( alpha_seq[which(ll_seq==min(na.omit(ll_seq[na.omit(ll_seq)>=0])) & !is.na(ll_seq))][1] )
           if(length(pval)==0 | is.na(pval)) pval = alpha_seq[which(ul_seq==max(na.omit(ul_seq[na.omit(ul_seq)<=0])) & !is.na(ul_seq))][1]
         }


      }
    }
  }




  PARAMETER = total_obs


  if (!pvalue) { method_label = paste0(method_label,' (bootstrap)')
  } else if (pvalue && pvalue.type=='permutation'){ method_label = paste0(method_label,' (permutation & bootstrap)')
  } else if (pvalue && pvalue.type=='CI.inversion'){ method_label = paste0(method_label,' (bootstrap & CI inversion)')}

  names(STATISTIC) <- "X-squared"
  names(PARAMETER) <- "n"
  if (any(E < 5) && is.finite(PARAMETER))
    warning("Chi-squared approximation may be incorrect")

  RVAL = list(parameter = PARAMETER, estimate = estimate,
                 method = method_label, data.name = DNAME, observed = x,
                 expected = E, residuals = (x - E)/sqrt(E))

  if (pvalue) RVAL = c(RVAL, list(statistic = STATISTIC, p.value = pval))

  if (perm.values) RVAL = c(RVAL, list(perm.values = perm_theta_hat))

  if (confint_present) RVAL = c(RVAL, list(conf.int = cint))
  if (confint_present & boot.values) RVAL = c(RVAL, list(boot.values = boot_theta_hat))

  class(RVAL) <- "htest"
  RVAL
}





# x = data.frame(x1=c(20, 80), x2=c(50,50))
#
# #chisq.test.r(x)
# set.seed(6)
# x=c(rbinom(50,1,.2), rbinom(50,1,.1))
# y=rep(0:1, each=50)
# chisq.test.r(x, y, pvalue=T, confint=T, correct = F)
# chisq.test(x, y, correct=F)
# chisq.test.r(x, y)
#
# tt = chisq.test.r(x=data.frame(x1=c(20, 80), x2=c(50,50)), pvalue=T, boot.values = T)
# tt
#
#
# tt$parameter
#
# x = as.matrix(table(x,y))





