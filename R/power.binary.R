#' Power Calculation for Two-Arm Trials with Binary Endpoints
#'
#' @description
#' This function calculates the statistical power for two-arm superiority trials
#' with binary endpoints. It supports five types of one-sided statistical tests:
#' the Pearson chi-squared test, the Fisher exact test, the Fisher mid-P test,
#' the Z-pooled exact unconditional test, and the Boschloo exact unconditional test.
#' The power is computed exactly using the binomial distribution without
#' approximation.
#'
#' @details
#' The function implements the exact power calculation methodology described in
#' Homma and Yoshida (20XX). Power is calculated by summing the probabilities
#' of all outcomes in the rejection region under the alternative hypothesis.
#' The rejection region is obtained from the \code{RR.binary} function.
#' If vectors of \code{p1} and \code{p2} are provided, the function returns
#' a vector of powers corresponding to each (p1, p2) combination.
#'
#' @param p1 True probability of responders for group 1 (treatment group).
#'   Can be a scalar or a vector.
#' @param p2 True probability of responders for group 2 (control group).
#'   Must have the same length as \code{p1}.
#' @param N1 Sample size for group 1.
#' @param N2 Sample size for group 2.
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
#' @return A numeric vector of power values. If \code{p1} and \code{p2} are
#'   scalars, returns a single power value. If vectors are provided,
#'   returns a vector of powers with the same length.
#'
#' @importFrom fpCompare %!=%
#' @importFrom stats dbinom pbinom
#'
#' @export
#'
#' @references
#' Homma, G. and Yoshida, T. (20XX). A 2-in-1 adaptive design for binary endpoints.
#'
#' @examples
#' # Example 1: Single power calculation
#' p1 <- 0.5
#' p2 <- 0.2
#' N1 <- 10
#' N2 <- 40
#' alpha <- 0.025
#'
#' # Pearson chi-squared test
#' power.binary(p1, p2, N1, N2, alpha, Test = "Chisq")
#'
#' # Fisher exact test
#' power.binary(p1, p2, N1, N2, alpha, Test = "Fisher")
#'
#' # Example 2: Multiple power calculations
#' p1_vec <- c(0.5, 0.6, 0.7, 0.8)
#' p2_vec <- c(0.2, 0.2, 0.2, 0.2)
#'
#' # Calculate powers for all combinations
#' powers <- power.binary(p1_vec, p2_vec, N1, N2, alpha, Test = "Boschloo")
#' print(powers)
#'
#' # Example 3: Compare different tests
#' tests <- c("Chisq", "Fisher", "Fisher-midP", "Z-pool", "Boschloo")
#' for (test in tests) {
#'   pwr <- power.binary(0.6, 0.3, 30, 30, 0.025, Test = test)
#'   cat(sprintf("%s: Power = %.4f\n", test, pwr))
#' }
power.binary <- function(p1, p2, N1, N2, alpha, Test) {
  # Check that p1 and p2 are the same length.
  if(length(p1) %!=% length(p2)) stop('p1 and p2 should be the same length')
  # Set rejection region
  RR <- RR.binary(N1, N2, alpha, Test)
  # Return power
  sapply(seq(length(p1)), function(i) sum(dbinom(0:N1, N1, p1[i]) * pbinom(rowSums(RR) - 1, N2, p2[i])))
}
