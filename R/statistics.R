#------------------------------------------------------------------------------
# R file             --Tianfu Yang
# Type:              Functions
# Subtype/Project:   Statistical analysis functions
# Descripetions:     
# Last Update:       Aug 27, 2015 12:14 PM
# Contents:
# 1
# Function:     Stat.qvalue(p.value)
# Description:  calculate the q-value based on a set of p-value in MultiTest
#------------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# UPDATED Aug 27, 2015 12:14 PM
# FUNCTION:     Stat.qvalue(p.value)
#' @title       calculate q-values based on a set of p-values in Multiple Tests
#' @param       p.value The p-values obtained from multiple testing
#' @return      q.value with the same length of p.value
# -----------------------------------------------------------------------------
#' @export      
#' @note        This is a revised version of the qvalue() function in 
#'              R/rrBLUP::GWAS(). 
#' @references  Estimating the False Discovery Rate using SAS (Paper 190-31)
#' @examples    
#' p = runif(1000)
#' q2 = Stat.qvalue(p)
# -----------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Unit test 
# 
# # The comparison with qvalue() in rrBLUP::GWAS()
# qvalue <- function(p) {           ## copy from R/rrBLUP::GWAS
#   smooth.df = 3
#   if (min(p) < 0 || max(p) > 1) {
#     print("ERROR: p-values not in valid range.")
#     return(0)
#   }
#   lambda = seq(0, 0.9, 0.05)
#   m <- length(p)
#   pi0 <- rep(0, length(lambda))
#   for (i in 1:length(lambda)) {
#     pi0[i] <- mean(p >= lambda[i])/(1 - lambda[i])
#   }
#   spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
#   pi0 <- predict(spi0, x = max(lambda))$y
#   pi0 <- min(pi0, 1)
#   if (pi0 <= 0) {
#     print("ERROR: The estimated pi0 <= 0. Check that you have valid p-values.")
#     return(0)
#   }
#   u <- order(p)
#   qvalue.rank <- function(x) {
#     idx <- sort.list(x)
#     fc <- factor(x)
#     nl <- length(levels(fc))
#     bin <- as.integer(fc)
#     tbl <- tabulate(bin)
#     cs <- cumsum(tbl)
#     tbl <- rep(cs, tbl)
#     tbl[idx] <- tbl
#     return(tbl)
#   }
#   v <- qvalue.rank(p)
#   qvalue <- pi0 * m * p/v
#   qvalue[u[m]] <- min(qvalue[u[m]], 1)
#   for (i in (m - 1):1) {
#     qvalue[u[i]] <- min(qvalue[u[i]], qvalue[u[i + 1]], 
#                         1)
#   }
#   return(qvalue)
# }
# 
# p = runif(1000)
# q1 = qvalue(p)
# q2 = Stat.qvalue(p)
# any(q1 != q2)
# 
# [1] FALSE
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#   qvalue.rank <- function(x) {    ## input p
#     idx <- sort.list(x)           ## the same as order(p)
#     fc <- factor(x)               ## trun to factors
#     bin <- as.integer(fc)         ## integer for factors (possible tie)
#     tbl <- tabulate(bin)          ## count the number
#     cs <- cumsum(tbl)             ## cumulative sum
#     tbl <- rep(cs, tbl)           ## 
#     tbl[idx] <- tbl
#     return(tbl)
#   } # it is a function to give the rank allowing tie
#   
#------------------------------------------------------------------------------
# Mar 12, 2015 12:24
# Is it similar to match(p,sort(p))?
#   In this expression the rank would be "higher" than qvalue.rank()
#   rank(p, ties.method = "max")  may be a better choice
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# Unit test 
# # The test here did not consider the situation when two P-val is equal
# # All equal.
#   for (i in 1: 10000){
#     p = rbinom(100,50,0.45) / 50
#     r1 = rank(p, ties.method = "max")
#     r2 = qvalue.rank(p)
#     if (any(r1 != r2)) {
#       stop()
#     }
#   }
#   for (i in 1:5000){
#     p = rnorm(1000)
#     r1 = rank(p,ties.method = "max")
#     r2 = qvalue.rank(p)
#     if (any(r1 != r2)) {
#       stop()
#     }
#   }
#------------------------------------------------------------------------------   

Stat.qvalue = function(p.value){
  ## check p data
  if (min(p.value,na.rm = T) < 0 || max(p.value,na.rm = T) > 1) {
    stop("p-values should be in the range of [0,1]")
  }
  m = length(p.value)
  ## calculate pi0
  lambda = seq(0, 0.9, 0.05) 
  pi0    = rep(0, length(lambda))
  for (i in 1:length(lambda)) {               ## negative value
    pi0[i] = mean(p.value >= lambda[i]) / (1 - lambda[i])
  }
  
  spi0 = smooth.spline(lambda, pi0, df = 3)
  pi0  = predict(spi0, x = max(lambda))$y
  pi0  = min(pi0, 1)
  if (pi0 <= 0) {
    stop("pi0 should be larger than 0")
  }
  
  u = order(p.value)      ## the index of p-value, greater for larger value
  v = rank(p.value, ties.method = "max")    ## ranking allowing ties
  qvalue = pi0 * m * p.value / v            ##  
  qvalue[u[m]] = min(qvalue[u[m]], 1) ## no larger than 1
  for (i in (m - 1):1) {
    qvalue[u[i]] <- min(qvalue[u[i]], qvalue[u[i + 1]], 1)
  }
  return(qvalue)
}
