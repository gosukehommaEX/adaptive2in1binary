#' Sample Size Calculation for 2-in-1 Adaptive Design with Binary Endpoints
#'
#' @description
#' This function calculates the required sample size for a 2-in-1 adaptive
#' phase 2/3 design with binary endpoints to achieve a target power. The 2-in-1
#' adaptive design allows a trial to maintain a small trial or to expand to a
#' large trial adaptively based on decisions made at the interim analysis.
#' The function provides results for five types of one-sided statistical tests
#' and compares the required sample size with a traditional design.
#'
#' @details
#' The function implements the sample size calculation methodology described in
#' Homma and Yoshida (20XX) for 2-in-1 adaptive designs with binary endpoints.
#' The algorithm proceeds as follows:
#' \enumerate{
#'   \item Calculate the required sample size for a traditional design (Phase 3 only)
#'     using the \code{sample.size.binary} function
#'   \item Set initial sample sizes for the 2-in-1 design based on the allocation
#'     ratios pi1 and pi2
#'   \item Iteratively increase the sample size until the calculated power meets
#'     or exceeds the target power
#' }
#'
#' The function returns comprehensive information including:
#' \itemize{
#'   \item Required sample sizes for both traditional and 2-in-1 designs
#'   \item Expected sample size (ESS)
#'   \item Go probability
#'   \item Power for Phase 2 part, Phase 3 part, and total power
#' }
#'
#' @param p1 True probability of responders for group 1 (treatment group).
#' @param p2 True probability of responders for group 2 (control group).
#' @param r Allocation ratio to group 1 (i.e., group 1:group 2 = r:1, r > 0).
#' @param pi1 Allocation ratio of total sample size to stage 1.
#'   N21 = ceiling(pi1 * N2) and N11 = ceiling(r * N21), where N2 is the
#'   total sample size for group 2.
#' @param pi2 Allocation ratio of total sample size to stage 2.
#'   N22 = ceiling(pi2 * N2) and N12 = ceiling(r * N22).
#'   Note: To control type I error rate, pi2 should be less than or equal to 1 - pi1.
#' @param cutpoint A cutpoint of the risk difference for the interim decision.
#'   The value should be within -1 < cutpoint <= 1.
#' @param alpha One-sided level of significance (e.g., 0.025).
#' @param tar.power Target power (e.g., 0.8 or 0.9).
#'
#' @return A tibble with the following columns:
#'   \itemize{
#'     \item \code{p1}: True probability of responders for group 1
#'     \item \code{p2}: True probability of responders for group 2
#'     \item \code{r}: Allocation ratio
#'     \item \code{pi1}: Allocation ratio to stage 1
#'     \item \code{pi2}: Allocation ratio to stage 2
#'     \item \code{alpha}: One-sided level of significance
#'     \item \code{tar.power}: Target power
#'     \item \code{Test}: Statistical testing approach
#'     \item \code{cutpoint}: Cutpoint of the risk difference
#'     \item \code{N1.trad}: Required sample size of group 1 for a traditional design
#'     \item \code{N2.trad}: Required sample size of group 2 for a traditional design
#'     \item \code{N.trad}: Required total sample size for a traditional design
#'     \item \code{N1.2in1}: Required sample size of group 1 for a 2-in-1 design
#'     \item \code{N2.2in1}: Required sample size of group 2 for a 2-in-1 design
#'     \item \code{N.2in1}: Required total sample size for a 2-in-1 design
#'     \item \code{ESS}: Expected sample size
#'     \item \code{Go.prob}: Go probability of expanding to a Phase 3 trial
#'     \item \code{Power.2}: Power of tests for Phase 2 part
#'     \item \code{Power.3}: Power of tests for Phase 3 part
#'     \item \code{Power.Total}: Total power of tests (sum of Power.2 and Power.3)
#'   }
#'
#' @importFrom dplyr tibble group_by reframe mutate ungroup select filter rename left_join bind_rows pull
#' @importFrom tidyr unnest
#' @importFrom fpCompare %<<%
#'
#' @export
#'
#' @references
#' Homma, G. and Yoshida, T. (20XX). A 2-in-1 adaptive design for binary endpoints.
#'
#' @examples
#' # Example 1: Basic usage
#' sample.size.2in1.binary(
#'   p1 = 0.4, p2 = 0.2, r = 2, pi1 = 0.2, pi2 = 0.5,
#'   cutpoint = 0.2, alpha = 0.025, tar.power = 0.8
#' )
#'
#' # Example 2: Balanced allocation
#' sample.size.2in1.binary(
#'   p1 = 0.6, p2 = 0.3, r = 1, pi1 = 0.3, pi2 = 0.3,
#'   cutpoint = 0.3, alpha = 0.025, tar.power = 0.9
#' )
#'
#' # Example 3: Compare different cutpoints
#' cutpoints <- c(0.15, 0.20, 0.25)
#' results <- lapply(cutpoints, function(c) {
#'   sample.size.2in1.binary(
#'     p1 = 0.5, p2 = 0.3, r = 1, pi1 = 0.2, pi2 = 0.3,
#'     cutpoint = c, alpha = 0.025, tar.power = 0.8
#'   )
#' })
#' do.call(rbind, results)
#'
#' # Example 4: High power requirement
#' sample.size.2in1.binary(
#'   p1 = 0.7, p2 = 0.3, r = 1, pi1 = 0.3, pi2 = 0.3,
#'   cutpoint = 0.5, alpha = 0.025, tar.power = 0.9
#' )
sample.size.2in1.binary <- function(p1, p2, r, pi1, pi2, cutpoint, alpha, tar.power) {
  # Step 0: sample size calculation for a traditional design of a phase 3 trial ensuring power = tar.power
  info.trad <- tibble(
    p1 = p1,
    p2 = p2
  ) %>%
    group_by_all() %>%
    reframe(
      Test = factor(
        c('Chisq', 'Fisher', 'Fisher-midP', 'Z-pool', 'Boschloo'),
        levels = c('Chisq', 'Fisher', 'Fisher-midP', 'Z-pool', 'Boschloo')
      )
    ) %>%
    group_by(p1, p2, Test) %>%
    mutate(
      sample.size = list(
        sample.size.binary(p1, p2, r, alpha, tar.power, Test)
      )
    ) %>%
    ungroup() %>%
    select(sample.size) %>%
    unnest(sample.size)
  N.min <- info.trad %>%
    filter(N == min(N)) %>%
    slice(1) %>%
    select(N1, N2, N)
  N2 <- N.min[['N2']]
  N1 <- N.min[['N1']]
  N21 <- ceiling(pi1 * N2)
  N11 <- ceiling(r * N21)
  N22 <- ceiling(pi2 * N2)
  N12 <- ceiling(r * N22)
  N23 <- N2 - N21
  N13 <- N1 - N11
  # Step 1: power calculation for a 2-in-1 adaptive phase 2/3 design given initial sample sizes
  info.2in1 <- power.2in1.binary(p1, p2, N11, N21, N12, N22, N13, N23, cutpoint, alpha) %>%
    mutate(
      N1 = N1,
      N2 = N2,
      N = N1 + N2
    )
  Power <- info.2in1 %>%
    pull(Power.Total) %>%
    min()
  # Step 2: sample size will be increased until calculated power satisfying the target power
  while(Power %<<% tar.power) {
    N2 <- N2 + 1
    N1 <- ceiling(r * N2)
    N <- N1 + N2
    N21 <- ceiling(pi1 * N2)
    N11 <- ceiling(r * N21)
    N22 <- ceiling(pi2 * N2)
    N12 <- ceiling(r * N22)
    N23 <- N2 - N21
    N13 <- N1 - N11
    info.2in1 <- info.2in1 %>%
      bind_rows(
        power.2in1.binary(p1, p2, N11, N21, N12, N22, N13, N23, cutpoint, alpha) %>%
          mutate(
            N1 = N1,
            N2 = N2,
            N = N1 + N2
          )
      )
    Power <- info.2in1 %>%
      filter(N == .env$N) %>%
      pull(Power.Total) %>%
      min()
  }
  # Step 3: print summary result
  result <- info.2in1 %>%
    group_by(p1, p2, Test) %>%
    filter(
      Power.Total >= tar.power
    ) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(
      pi1 = pi1,
      pi2 = pi2,
      .after = p2
    ) %>%
    rename(N1.2in1 = N1, N2.2in1 = N2, N.2in1 = N) %>%
    left_join(
      info.trad %>%
        rename(N1.trad = N1, N2.trad = N2, N.trad = N) %>%
        select(-Power),
      by = c('p1', 'p2', 'Test')
    ) %>%
    select(
      p1, p2, r, pi1, pi2, alpha, tar.power, Test, cutpoint, contains('trad'), contains('2in1'), ESS, Go.prob, contains('Power.')
    )

  # Return result
  return(result)
}
