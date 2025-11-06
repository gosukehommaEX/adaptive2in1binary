########################################################################################################################
## Script for Table 1: Total maximum sample size comparison across different scenarios
##
## Description:
##   This script generates Table 1 from Homma and Yoshida (20XX), which shows
##   the total maximum sample size and its corresponding ESS for each test
##   across three cases (A, B, C) with different allocation ratios (r=1, 2).
##
## Design parameters:
##   - Case A: p1=0.6, p2=0.2
##   - Case B: p1=0.6, p2=0.3
##   - Case C: p1=0.6, p2=0.4
##   - Allocation ratios: r = 1 (balanced), r = 2 (imbalanced)
##   - pi1 = pi2 = 0.3
##   - alpha2 = alpha3 = 0.025
##   - target power = 0.8
##   - rho = 0.5, 1, 1.5
##
## Output:
##   - table1.tex: LaTeX table file
##
## Usage:
##   source("table1.R")
##
## Note:
##   This script requires the following source file in the same directory:
##   - sample.size.2in1.binary.R
##
##   Required packages: kableExtra
########################################################################################################################
library(dplyr)
library(tidyr)
library(kableExtra)
library(adaptive2in1binary)
options(knitr.table.format = 'latex')
p1 = c(0.6, 0.6, 0.6)
p2 = c(0.2, 0.3, 0.4)
r = c(1, 2)
pi1 = 0.3
pi2 = 0.3
alpha2 = 0.025
alpha3 = 0.025
tar.power = 0.8
rho = c(0.5, 1, 1.5)
# Results of the maximum sample sizes
results = tibble(
  p1 = p1,
  p2 = p2,
  Case = factor(
    LETTERS[seq(length(p1))],
    levels = LETTERS[seq(length(p1))]
  )
) %>%
  group_by_all() %>%
  reframe(
    r = r
  ) %>%
  group_by_all() %>%
  reframe(
    rho = rho,
    cutpoint = rho * (p1 - p2)
  ) %>%
  group_by_all() %>%
  reframe(
    sample.size = list(
      sample.size.2in1.binary(
        p1 = p1, p2 = p2, r, pi1, pi2, cutpoint = cutpoint, alpha2, alpha3, tar.power
      )
    )
  ) %>%
  select(Case, rho, sample.size) %>%
  unnest(
    sample.size
  )

# Sample size table
latex.table = results %>%
  mutate(
    N.and.ESS = paste(N.2in1, if_else(nchar(N.2in1) == 3, ' (', '   ('), sprintf('%.1f', ESS), ')', sep = '')
  ) %>%
  select(Case, r, Test, cutpoint, Go.prob, N.and.ESS) %>%
  group_by(Case, r, cutpoint) %>%
  mutate(
    Go.prob = paste0('(', sprintf('%.2f', min(Go.prob)), '-', sprintf('%.2f', max(Go.prob)), ')')
  ) %>%
  arrange(Case, r, cutpoint, Test) %>%
  pivot_wider(
    names_from = Test, values_from = N.and.ESS
  ) %>%
  ungroup() %>%
  rename(
    '$r$' = r,
    '$C$' = cutpoint,
    'Range of $p_{{\\rm Go}}$' = Go.prob
  ) %>%
  kbl(format = 'latex', booktabs = TRUE, linesep = '', escape = F, align = 'l') %>%
  landscape(margin = NULL) %>%
  add_header_above(c(' ' = 4, 'The total maximum sample size (ESS)' = 5)) %>%
  collapse_rows(columns = 1:2, latex_hline = 'major', row_group_label_position = 'first')

# Save latex table file
latex.table %>%
  save_kable('table1.tex')
