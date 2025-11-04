#' Rejection Region for Two-Arm Trials with Binary Endpoints
#'
#' @description
#' This function calculates the rejection region for two-arm superiority trials
#' with binary endpoints. It supports five types of one-sided statistical tests:
#' the Pearson chi-squared test, the Fisher exact test, the Fisher mid-P test,
#' the Z-pooled exact unconditional test, and the Boschloo exact unconditional test.
#' The rejection region is determined based on the specified sample sizes and
#' significance level.
#'
#' @details
#' The function implements the methodology described in Homma and Yoshida (20XX)
#' for calculating rejection regions in the context of 2-in-1 adaptive designs
#' with binary endpoints. For the Pearson chi-squared test and Z-pooled test,
#' the rejection region is based on the standardized test statistic. For Fisher-based
#' tests, the rejection region is determined using the hypergeometric distribution.
#' The Boschloo test maximizes the type I error rate over the nuisance parameter
#' under the null hypothesis to obtain an exact unconditional test.
#'
#' @param N1 Sample size for group 1 (treatment group).
#' @param N2 Sample size for group 2 (control group).
#' @param alpha One-sided level of significance (e.g., 0.025).
#' @param Test Type of statistical test. Must be one of:
#'   \itemize{
#'     \item "Chisq": One-sided Pearson chi-squared test
#'     \item "Fisher": Fisher exact test
#'     \item "Fisher-midP": Fisher mid-P test
#'     \item "Z-pool": Z-pooled exact unconditional test
#'     \item "Boschloo": Boschloo exact unconditional test
#'   }
#'
#' @return A matrix of logical values (TRUE/FALSE) of dimension (N1 + 1) x (N2 + 1),
#'   where TRUE indicates that the combination of observed responders (x1, x2)
#'   falls in the rejection region. Rows correspond to the number of responders
#'   in group 1 (0 to N1) and columns correspond to the number of responders
#'   in group 2 (0 to N2).
#'
#' @importFrom fpCompare %<<%  %<=%  %>=%  %>>%  %!=%
#' @importFrom stats dbinom pnorm phyper dhyper
#'
#' @export
#'
#' @references
#' Homma, G. and Yoshida, T. (20XX). A 2-in-1 adaptive design for binary endpoints.
#'
#' @examples
#' # Example 1: Boschloo exact unconditional test
#' N1 <- 20
#' N2 <- 10
#' alpha <- 0.025
#' RR <- RR.binary(N1, N2, alpha, Test = "Boschloo")
#' print(RR)
#'
#' # Example 2: Fisher exact test
#' RR_fisher <- RR.binary(N1 = 30, N2 = 30, alpha = 0.025, Test = "Fisher")
#'
#' # Example 3: Pearson chi-squared test
#' RR_chisq <- RR.binary(N1 = 40, N2 = 20, alpha = 0.025, Test = "Chisq")
RR.binary <- function(N1, N2, alpha, Test) {
  if((Test == 'Chisq') | (Test == 'Z-pool')) {
    # Test statistics for the Chisq over all combinations of i(=0,...,N1) and j(=0,...,N2)
    Z.ij <- '/'(
      outer(N2 * (0:N1), N1 * (0:N2), '-') / (N1 * N2),
      sqrt(outer(0:N1, 0:N2, '+') / (N1 * N2) * (1 - outer(0:N1, 0:N2, '+') / (N1 + N2)))
    )
    Z.ij[is.na(Z.ij)] <- 0
    if(Test == 'Chisq') {
      # p-values for the Chisq
      p.val <- 1 - pnorm(Z.ij)
    } else {
      # Since zero and negative values of test statistics must not be statistically significant, they are omitted
      Z.ij.posi <- Z.ij[Z.ij %>>% 0]
      order.Z.ij.posi <- order(Z.ij.posi, decreasing = TRUE)
      i <- (row(Z.ij)[Z.ij %>>% 0] - 1)[order.Z.ij.posi]
      j <- (col(Z.ij)[Z.ij %>>% 0] - 1)[order.Z.ij.posi]
      # Calculate P_H0(X1 = i, X2 = j | theta)
      uniq.i <- sort(unique(i))
      dbinom.i <- sapply(seq(0, 1, l = 100), function(theta) dbinom(uniq.i, N1, theta))
      uniq.j <- sort(unique(j))
      dbinom.j <- sapply(seq(0, 1, l = 100), function(theta) dbinom(uniq.j, N2, theta))
      P_H0 <- dbinom.i[i, ] * dbinom.j[j + 1, ]
      # p-values for all possible values of theta
      p.ij <- apply(apply(P_H0, 2, cumsum), 1, max)
      # p-values for the Z-pool
      p.val <- 0 * Z.ij + 1
      p.val[cbind(i + 1, j + 1)] <- p.ij
    }
  } else if(Test == 'Fisher-midP') {
    # p-values for the Fisher-midP
    p.val <- outer(0:N1, 0:N2, function(i, j) {
      phyper(i, N1, N2, i + j, lower.tail = FALSE) + 0.5 * dhyper(i, N1, N2, i + j)
    })
  } else if((Test == 'Fisher') | (Test == 'Boschloo')) {
    # p-values for the Fisher over all combinations of i(=0,...,N1) and j(=0,...,N2)
    p.fisher <- outer(0:N1, 0:N2, function(i, j) phyper(i - 1, N1, N2, i + j, lower.tail = FALSE))
    if(Test == 'Fisher') {
      # p-values for the Fisher
      p.val <- p.fisher
    } else {
      # Since p-values satisfying hat{p}_{1} - hat{p}_{2} <= 0 must not be statistically significant, they are omitted
      p.max.boschloo <- min(p.fisher[(outer(N2 * (0:N1), N1 * (0:N2), '-') / (N1 * N2)) %<=% 0])
      p.fisher.posi <- p.fisher[p.fisher %<<% p.max.boschloo]
      order.p.fisher.posi <- order(p.fisher.posi, decreasing = FALSE)
      i <- (row(p.fisher)[p.fisher %<<% p.max.boschloo] - 1)[order.p.fisher.posi]
      j <- (col(p.fisher)[p.fisher %<<% p.max.boschloo] - 1)[order.p.fisher.posi]
      # Calculate P_H0(X1 = i, X2 = j | theta)
      uniq.i <- sort(unique(i))
      dbinom.i <- sapply(seq(0, 1, l = 100), function(theta) dbinom(uniq.i, N1, theta))
      uniq.j <- sort(unique(j))
      dbinom.j <- sapply(seq(0, 1, l = 100), function(theta) dbinom(uniq.j, N2, theta))
      P_H0 <- dbinom.i[i - min(uniq.i) + 1, ] * dbinom.j[j + 1, ]
      # p-values for all possible values of theta
      p.ij <- apply(apply(P_H0, 2, cumsum), 1, max)
      # p-values for the Boschloo
      p.val <- 0 * p.fisher + 1
      p.val[cbind(i + 1, j + 1)] <- p.ij
    }
  }
  # Rejection region
  RR <- (p.val %<<% alpha)
  # Return RR
  return(RR)
}
