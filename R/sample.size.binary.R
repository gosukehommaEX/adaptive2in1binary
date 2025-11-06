#' Sample Size Calculation for Two-Arm Trials with Binary Endpoints
#'
#' @description
#' This function calculates the required sample size for two-arm superiority
#' trials with binary endpoints to achieve a target power. It supports five
#' types of one-sided statistical tests: the Pearson chi-squared test, the
#' Fisher exact test, the Fisher mid-P test, the Z-pooled exact unconditional
#' test, and the Boschloo exact unconditional test. The function uses a grid
#' search algorithm to find the minimum sample size.
#'
#' @details
#' The function implements the sample size calculation methodology described in
#' Homma and Yoshida (20XX). The algorithm starts with an initial sample size
#' estimate based on the normal approximation for the Pearson chi-squared test,
#' then uses a grid search to find the exact minimum sample size that achieves
#' the target power. The exact power is calculated using the \code{power.binary}
#' function without approximation.
#'
#' @param p1 True probability of responders for group 1 (treatment group).
#' @param p2 True probability of responders for group 2 (control group).
#' @param r Allocation ratio to group 1 (i.e., group 1:group 2 = r:1, r > 0).
#' @param alpha One-sided level of significance (e.g., 0.025).
#' @param tar.power Target power (e.g., 0.8 or 0.9).
#' @param Test Type of statistical test. Must be one of:
#'   \itemize{
#'     \item "Chisq": One-sided Pearson chi-squared test
#'     \item "Fisher": Fisher exact test
#'     \item "Fisher-midP": Fisher mid-P test
#'     \item "Z-pool": Z-pooled exact unconditional test
#'     \item "Boschloo": Boschloo exact unconditional test
#'   }
#'
#' @return A data frame with the following columns:
#'   \itemize{
#'     \item \code{p1}: True probability of responders for group 1
#'     \item \code{p2}: True probability of responders for group 2
#'     \item \code{r}: Allocation ratio
#'     \item \code{alpha}: One-sided level of significance
#'     \item \code{tar.power}: Target power
#'     \item \code{Test}: Type of statistical test
#'     \item \code{Power}: Achieved power
#'     \item \code{N1}: Required sample size for group 1
#'     \item \code{N2}: Required sample size for group 2
#'     \item \code{N}: Required total sample size
#'   }
#'
#' @importFrom fpCompare %>=%  %<<%
#' @importFrom stats qnorm
#'
#' @export
#'
#' @examples
#' # Example 1: Pearson chi-squared test
#' sample.size.binary(p1 = 0.4, p2 = 0.2, r = 2, alpha = 0.025,
#'                    tar.power = 0.8, Test = "Chisq")
#'
#' # Example 2: Fisher exact test
#' sample.size.binary(p1 = 0.5, p2 = 0.2, r = 3, alpha = 0.025,
#'                    tar.power = 0.9, Test = "Fisher")
#'
#' # Example 3: Boschloo exact unconditional test
#' sample.size.binary(p1 = 0.7, p2 = 0.2, r = 4, alpha = 0.025,
#'                    tar.power = 0.8, Test = "Boschloo")
#'
#' # Example 4: Compare sample sizes across different tests
#' tests <- c("Chisq", "Fisher", "Fisher-midP", "Z-pool", "Boschloo")
#' results <- lapply(tests, function(test) {
#'   sample.size.binary(p1 = 0.6, p2 = 0.3, r = 1, alpha = 0.025,
#'                      tar.power = 0.8, Test = test)
#' })
#' do.call(rbind, results)
sample.size.binary <- function(p1, p2, r, alpha, tar.power, Test) {
  # Step 0 (calculate the required sample size for the one-sided Pearson chi-squared test)
  p <- (r * p1 + p2) / (1 + r)
  init_N2 <- '*'(
    (1 + 1 / r) / ((p1 - p2) ^ 2),
    (qnorm(alpha) * sqrt(p * (1 - p)) + qnorm(1 - tar.power) * sqrt((p1 * (1 - p1) / r + p2 * (1 - p2)) / (1 + 1 / r))) ^ 2
  )
  # Step 1 (power calculation given initial sample size)
  N2 <- ceiling(init_N2)
  N1 <- ceiling(r * N2)
  Power <- power.binary(p1, p2, N1, N2, alpha, Test)
  # Step 2 (sample size calculation via a grid search algorithm)
  if(Power %>=% tar.power) {
    while(Power %>=% tar.power) {
      N2 <- N2 - 1
      N1 <- ceiling(r * N2)
      Power <- power.binary(p1, p2, N1, N2, alpha, Test)
    }
    N2 <- N2 + 1
  } else {
    while(Power %<<% tar.power) {
      N2 <- N2 + 1
      N1 <- ceiling(r * N2)
      Power <- power.binary(p1, p2, N1, N2, alpha, Test)
    }
  }
  # Step 3 (determine the final sample size)
  N1 <- ceiling(r * N2)
  N <- N1 + N2
  Power <- power.binary(p1, p2, N1, N2, alpha, Test)
  # Return result
  result <- data.frame(p1, p2, r, alpha, tar.power, Test, Power, N1, N2, N)
  return(result)
}
