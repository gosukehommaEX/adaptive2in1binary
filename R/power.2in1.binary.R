#' Power Calculation for 2-in-1 Adaptive Design with Binary Endpoints
#'
#' @description
#' This function calculates the statistical power for a 2-in-1 adaptive phase 2/3
#' design with binary endpoints using closed-form formulas. The 2-in-1 adaptive
#' design allows a trial to maintain a small trial or to expand to a large trial
#' adaptively based on decisions made at the interim analysis. The function
#' simultaneously provides results for five types of one-sided statistical tests.
#'
#' @details
#' The function implements the exact power calculation methodology described in
#' Homma and Yoshida (20XX) for 2-in-1 adaptive designs with binary endpoints.
#' The design proceeds in three stages: Stage 1 (interim analysis), Stage 2
#' (small trial part), and Stage 3 (large trial part). At the end of Stage 1,
#' if the observed risk difference is greater than or equal to the cutpoint C,
#' the trial is expanded to Stage 3. Otherwise, it proceeds to Stage 2.
#'
#' The power calculation accounts for:
#' \itemize{
#'   \item Go probability: Probability of expanding to a large trial
#'   \item Power.2: Power for the Phase 2 part (Stage 1 + Stage 2)
#'   \item Power.3: Power for the Phase 3 part (Stage 1 + Stage 3)
#'   \item Power.Total: Total power (sum of Power.2 and Power.3)
#'   \item ESS: Expected sample size
#' }
#'
#' The exact power is obtained without Monte Carlo simulations by evaluating
#' all possible combinations of responders at each stage.
#'
#' @param p1 True probability of responders for group 1 (treatment group).
#'   Can be a scalar or a vector.
#' @param p2 True probability of responders for group 2 (control group).
#'   Must have the same length as \code{p1}.
#' @param N11 Sample size for group 1 at stage 1 (interim analysis).
#' @param N21 Sample size for group 2 at stage 1 (interim analysis).
#' @param N12 Sample size for group 1 at stage 2 (Phase 2 part).
#'   Total sample size for Phase 2 is N11 + N12.
#' @param N22 Sample size for group 2 at stage 2 (Phase 2 part).
#'   Total sample size for Phase 2 is N21 + N22.
#' @param N13 Sample size for group 1 at stage 3 (Phase 3 part).
#'   Total sample size for Phase 3 is N11 + N13.
#' @param N23 Sample size for group 2 at stage 3 (Phase 3 part).
#'   Total sample size for Phase 3 is N21 + N23.
#' @param cutpoint A cutpoint of the risk difference for the interim decision.
#'   Can be NA (to return results for all possible cutpoints) or a single value
#'   within -1 < cutpoint <= 1.
#' @param alpha2 One-sided level of significance at satge 2 (e.g., 0.025).
#' @param alpha3 One-sided level of significance at satge 3 (e.g., 0.025).
#'
#' @return A tibble with the following columns:
#'   \itemize{
#'     \item \code{p1}: True probability of responders for group 1
#'     \item \code{p2}: True probability of responders for group 2
#'     \item \code{cutpoint}: Cutpoint(s) of the risk difference
#'     \item \code{Test}: Names of the statistical tests (Chisq, Fisher, Fisher-midP, Z-pool, Boschloo)
#'     \item \code{Go.prob}: Go probability of expanding to a Phase 3 trial
#'     \item \code{Power.2}: Power of tests for Phase 2 part
#'     \item \code{Power.3}: Power of tests for Phase 3 part
#'     \item \code{Power.Total}: Total power of tests (sum of Power.2 and Power.3)
#'     \item \code{ESS}: Expected sample size
#'   }
#'
#' @importFrom dplyr tibble group_by mutate reframe ungroup filter select bind_cols rename across arrange if_any everything if_else cross_join join_by left_join inner_join n
#' @importFrom tidyr unnest pivot_longer pivot_wider
#' @importFrom purrr map imap
#' @importFrom fpCompare %!=%  %<=% %==% %<<%
#' @importFrom stats dbinom
#' @importFrom magrittr %>%
#' @importFrom tidyselect contains
#'
#' @export
#'
#' @examples
#' # Example 1: Single cutpoint
#' power.2in1.binary(
#'   p1 = 0.4, p2 = 0.2, N11 = 40, N21 = 20, N12 = 50, N22 = 25,
#'   N13 = 100, N23 = 50, cutpoint = 0.2, alpha2 = 0.025, alpha3 = 0.025
#' )
#'
#' # Example 2: All possible cutpoints
#' power.2in1.binary(
#'   p1 = 0.4, p2 = 0.2, N11 = 40, N21 = 20, N12 = 50, N22 = 25,
#'   N13 = 100, N23 = 50, cutpoint = NA, alpha2 = 0.025, alpha3 = 0.025
#' )
#'
#' # Example 3: Multiple (p1, p2) combinations
#' p1_vec <- c(0.4, 0.5, 0.6)
#' p2_vec <- c(0.2, 0.2, 0.2)
#' power.2in1.binary(
#'   p1 = p1_vec, p2 = p2_vec, N11 = 30, N21 = 30, N12 = 40, N22 = 40,
#'   N13 = 80, N23 = 80, cutpoint = 0.2, alpha2 = 0.025, alpha3 = 0.025
#' )
power.2in1.binary <- function(p1, p2, N11, N21, N12, N22, N13, N23, cutpoint, alpha2, alpha3) {

  # Check that p1 and p2 are the same length
  if(length(p1) %!=% length(p2)) stop('p1 and p2 should be the same length')

  # Probability mass functions of the binomial distribution for group j(=1,2) at stage k(=1,2,3)
  dbinom.js <- lapply(seq(3), function(k) {
    lapply(seq(2), function(j) {
      sapply(seq(length(p1)), function(i) {
        dbinom(
          seq(0, list(c(N11, N12, N13), c(N21, N22, N23))[[j]][[k]]),
          list(c(N11, N12, N13), c(N21, N22, N23))[[j]][[k]],
          list(p1, p2)[[j]][[i]]
        )
      })
    })
  })

  ### Run numerical investigations
  run.num.inv <- tibble(
    ## Set sample sizes for each treatment group at each stage
    # Stage s(s=2,3)
    s = 2:3,
    # Sample sizes at the interim analysis
    N11 = N11,
    N21 = N21,
    # Sample sizes at stage s(=2,3)
    N1s = c(N12, N13),
    N2s = c(N22, N23),
    # Total sample sizes for stage s(=2,3) analyses
    N1s.ast = N11 + N1s,
    N2s.ast = N21 + N2s,
    alpha = c(alpha2, alpha3)
  ) %>%
    group_by(s) %>%
    mutate(
      # Set rejection regions for all tests
      RR = map(list(1), ~ {
        test.name = c('Chisq', 'Fisher', 'Fisher-midP', 'Z-pool', 'Boschloo')
        do.call(
          bind_cols,
          lapply(test.name, function(m) {
            RR = RR.binary(N1s.ast, N2s.ast, alpha, m)
            tibble(
              x1s.ast = c((row(RR) - 1) * NA ^ (RR %<=% 0)),
              x2s.ast = c((col(RR) - 1) * NA ^ (RR %<=% 0))
            ) %>%
              setNames(paste0(m, '.', names(.))) %>%
              mutate(
                test.result = as.double(Reduce(`*`, .) >= 0)
              )
          }) %>%
            setNames(test.name) %>%
            imap(~ rename(.x, '{.y}' := test.result))
        ) %>%
          filter(if_any(everything(), ~ !is.na(.))) %>%
          mutate(
            x1s.ast = do.call(pmax, c(select(., contains('x1')), na.rm = TRUE)),
            x2s.ast = do.call(pmax, c(select(., contains('x2')), na.rm = TRUE))
          ) %>%
          select(
            x1s.ast, x2s.ast, matches(paste0('^', test.name, '$'))
          )
      })
    ) %>%
    mutate(
      # Dataset including Go probabilities and test results
      data.all = map(list(1), ~ {
        # Dataset of all possible combinations of x11 and x21 and corresponding hat.delta.1
        data.hat.delta.1 = tibble(
          x11 = rep(0:N11, each = N21 + 1),
          x21 = rep(0:N21, times = N11 + 1),
          hat.delta.1 = (N21 * x11 - N11 * x21) / (N11 * N21),
          Go.s = if_else(hat.delta.1 %<<% cutpoint, 2, 3)
        ) %>%
          filter((Go.s %==% s) | (is.na(Go.s))) %>%
          select(-Go.s) %>%
          mutate(
            hat.delta.1 = '+'(
              is.na(cutpoint) * hat.delta.1,
              if_else(is.na(cutpoint), 0, cutpoint)
            )
          )
        # Dataset of go probability
        data.Go.prob = data.hat.delta.1 %>%
          group_by(hat.delta.1) %>%
          reframe(
            p1 = p1,
            p2 = p2,
            Go.prob = colSums(
              '*'(
                dbinom.js[[1]][[1]][x11 + 1, , drop = FALSE],
                dbinom.js[[1]][[2]][x21 + 1, , drop = FALSE]
              )
            )
          ) %>%
          group_by(p1, p2) %>%
          reframe(
            hat.delta.1 = hat.delta.1,
            # Go probability calculation
            Go.prob = (s == 2) * 0 + (s == 3) * rev(cumsum(rev(Go.prob)))
          )
        # Dataset of test results
        data.hat.delta.1 %>%
          cross_join(
            tibble(
              min.x1s = 0,
              max.x1s = N1s,
              min.x2s = 0,
              max.x2s = N2s
            )
          ) %>%
          mutate(
            min.x1s.ast = x11 + min.x1s,
            max.x1s.ast = x11 + max.x1s,
            min.x2s.ast = x21 + min.x2s,
            max.x2s.ast = x21 + max.x2s
          ) %>%
          inner_join(
            RR[[1]],
            by = join_by(
              min.x1s.ast <= x1s.ast,
              max.x1s.ast >= x1s.ast,
              min.x2s.ast <= x2s.ast,
              max.x2s.ast >= x2s.ast
            )
          ) %>%
          mutate(
            x1s = x1s.ast - x11,
            x2s = x2s.ast - x21
          ) %>%
          select(
            x11, x21, hat.delta.1, x1s, x2s, c('Chisq', 'Fisher', 'Fisher-midP', 'Z-pool', 'Boschloo')
          ) %>%
          group_by(hat.delta.1) %>%
          reframe(
            p1 = p1,
            p2 = p2,
            across(
              contains(
                c('Chisq', 'Fisher', 'Fisher-midP', 'Z-pool', 'Boschloo')),
              ~ colSums(. * '*'(
                dbinom.js[[1]][[1]][x11 + 1, , drop = FALSE] * dbinom.js[[1]][[2]][x21 + 1, , drop = FALSE],
                dbinom.js[[s]][[1]][x1s + 1, , drop = FALSE] * dbinom.js[[s]][[2]][x2s + 1, , drop = FALSE]
              ), na.rm = TRUE)
            )
          ) %>%
          left_join(
            data.Go.prob, by = c('hat.delta.1', 'p1', 'p2')
          ) %>%
          rename(cutpoint = hat.delta.1)
      })
    ) %>%
    select(s, data.all) %>%
    unnest(data.all) %>%
    group_by(s, p1, p2) %>%
    reframe(
      cutpoint = cutpoint,
      Go.prob = Go.prob,
      across(
        # Power calculation
        contains(c('Chisq', 'Fisher', 'Fisher-midP', 'Z-pool', 'Boschloo')),
        ~ '+'(
          (s %==% 2) * ((n() > 1) * cumsum(c(0, rev(rev(.)[-1]))) + (n() == 1) * .),
          (s %==% 3) * rev(cumsum(rev(.)))
        )
      )
    ) %>%
    group_by(p1, p2, cutpoint) %>%
    pivot_longer(
      cols = -c(s, p1, p2, cutpoint, Go.prob), names_to = 'Test', values_to = 'Power'
    ) %>%
    # Result of power calculation
    mutate(
      Test = factor(Test, levels = c('Chisq', 'Fisher', 'Fisher-midP', 'Z-pool', 'Boschloo')),
      Go.prob = max(Go.prob)
    ) %>%
    arrange(p1, p2, cutpoint, Test) %>%
    ungroup() %>%
    pivot_wider(
      names_from = s, values_from = Power, names_prefix = 'Power.'
    ) %>%
    mutate(
      'Power.2' = {if('Power.2' %in% names(.)) pmax(0, Power.2, na.rm = TRUE) else 0},
      'Power.3' = {if('Power.3' %in% names(.)) pmax(0, Power.3, na.rm = TRUE) else 0},
      Power.Total = Power.2 + Power.3,
      # Calculate an expected sample size
      ESS = (N11 + N21) + (1 - Go.prob) * (N12 + N22) + Go.prob * (N13 + N23)
    ) %>%
    select(
      p1, p2, cutpoint, Test, Go.prob, Power.2, Power.3, Power.Total, ESS
    )

  # Return result
  return(run.num.inv)
}
