#' Test Size Calculation for 2-in-1 Adaptive Design with Binary Endpoints
#'
#' @description
#' This function calculates the test size (maximum type I error rate) for a
#' 2-in-1 adaptive phase 2/3 design with binary endpoints. The test size is
#' defined as the maximum type I error rate over the entire range of the
#' nuisance parameter theta under the composite null hypothesis H0: p1 = p2 = theta.
#' The function provides results for five types of one-sided statistical tests.
#'
#' @details
#' The function implements the test size calculation methodology described in
#' Homma and Yoshida (20XX) for 2-in-1 adaptive designs with binary endpoints.
#' While the point null hypothesis assumes H0: p1 = p2 = p for a specific value
#' of p, the true response probability might differ from p. Therefore, the exact
#' type I error rate must be evaluated under the composite null hypothesis
#' H0: p1 = p2 by searching over all possible values of theta in [0, 1].
#'
#' The test size is calculated by:
#' \enumerate{
#'   \item Computing power under H0 (i.e., p1 = p2) for a grid of theta values
#'   \item Taking the maximum power across all theta values
#' }
#'
#' This provides the maximum type I error rate, which can be used to verify
#' whether the trial controls the type I error rate at the nominal significance level.
#'
#' @param N11 Sample size for group 1 at stage 1 (interim analysis).
#' @param N21 Sample size for group 2 at stage 1 (interim analysis).
#' @param N12 Sample size for group 1 at stage 2 (Phase 2 part).
#' @param N22 Sample size for group 2 at stage 2 (Phase 2 part).
#' @param N13 Sample size for group 1 at stage 3 (Phase 3 part).
#' @param N23 Sample size for group 2 at stage 3 (Phase 3 part).
#' @param cutpoint A cutpoint of the risk difference for the interim decision.
#'   Can be NA (to return results for all possible cutpoints) or a single value
#'   within -1 < cutpoint <= 1.
#' @param alpha2 One-sided level of significance at satge 2 (e.g., 0.025).
#' @param alpha3 One-sided level of significance at satge 3 (e.g., 0.025).
#'
#' @return A tibble with the following columns:
#'   \itemize{
#'     \item \code{cutpoint}: Cutpoint(s) of the risk difference over x11 (0,...,N11) and x21 (0,...,N21)
#'     \item \code{Test}: Names of the statistical tests (Chisq, Fisher, Fisher-midP, Z-pool, Boschloo)
#'     \item \code{Test.size}: Test size (maximum type I error rate over the entire range of theta)
#'   }
#'
#' @importFrom dplyr select group_by reframe
#'
#' @export
#'
#' @references
#' Homma, G. and Yoshida, T. (20XX). A 2-in-1 adaptive design for binary endpoints.
#'
#' @examples
#' # Example 1: All possible cutpoints
#' test.size.2in1.binary(
#'   N11 = 40, N21 = 20, N12 = 50, N22 = 25, N13 = 100, N23 = 50,
#'   cutpoint = NA, alpha2 = 0.025, alpha3 = 0.025
#' )
#'
#' # Example 2: Specific cutpoint
#' test.size.2in1.binary(
#'   N11 = 40, N21 = 20, N12 = 50, N22 = 25, N13 = 100, N23 = 50,
#'   cutpoint = 0.2, alpha2 = 0.025, alpha3 = 0.025
#' )
#'
#' # Example 3: Verify type I error rate control
#' # For a balanced design with equal sample sizes
#' test_size_result <- test.size.2in1.binary(
#'   N11 = 30, N21 = 30, N12 = 40, N22 = 40, N13 = 80, N23 = 80,
#'   cutpoint = 0.15, alpha2 = 0.025, alpha3 = 0.025
#' )
#' print(test_size_result)
#'
#' # Check if test size is controlled at alpha = 0.025
#' max_test_size <- max(test_size_result$Test.size)
#' cat(sprintf("Maximum test size: %.4f\n", max_test_size))
#' if (max_test_size <= 0.025) {
#'   cat("Type I error rate is controlled at alpha = 0.025\n")
#' } else {
#'   cat("Type I error rate exceeds alpha = 0.025\n")
#' }
test.size.2in1.binary <- function(N11, N21, N12, N22, N13, N23, cutpoint, alpha2, alpha3) {

  # Calculate test size
  test.size <- power.2in1.binary(
    seq(0, 1, l = 100), seq(0, 1, l = 100), N11, N21, N12, N22, N13, N23, cutpoint, alpha2, alpha3
  ) %>%
    select(cutpoint, Test, Power.Total) %>%
    group_by(cutpoint, Test) %>%
    reframe(
      Test.size = max(Power.Total)
    )

  # Return output
  return(test.size)
}
