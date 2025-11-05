########################################################################################################################
## Script for Table 2: Power and test size for pyruvate kinase deficiency trial
##
## Description:
##   This script generates Table 2 from Homma and Yoshida (20XX), which shows
##   power (test size) and expected sample size for each test under different
##   combinations of pi1 and C for a trial in pyruvate kinase deficiency.
##
## Design parameters:
##   - p1 = 0.35 (expected response rate in treatment group)
##   - p2 = 0.05 (expected response rate in control group)
##   - N1 = N2 = 63 (fixed total sample size per group)
##   - r = 1 (balanced allocation ratio)
##   - pi1 = 0.1, 0.2, 0.3 (allocation ratios to stage 1)
##   - pi2 = 0.3 (allocation ratio to stage 2)
##   - alpha2 = alpha3 = 0.025
##   - C = 0.15, 0.30, 0.45 (cutpoints)
##
## Output:
##   - table2.tex: LaTeX table file
##
## Reference:
##   Homma G and Yoshida T (20XX). A 2-in-1 adaptive design for binary endpoints.
##   Statistics in Medicine (under review).
##
## Usage:
##   source("table2.R")
##
## Note:
##   This script requires the following source files in the same directory:
##   - test.size.2in1.binary.R
##   - power.2in1.binary.R
##
##   Required packages: kableExtra
########################################################################################################################
library(kableExtra)
options(knitr.table.format = 'latex')
source('test.size.2in1.binary.R')
source('power.2in1.binary.R')
p1 = 0.35
p2 = 0.05
r = 1
N1 = N2 = 63
pi1 = c(0.1, 0.2, 0.3)
pi2 = 0.3
alpha2 = 0.025
alpha3 = 0.025
N21 = ceiling(pi1 * N2)
N11 = ceiling(r * N21)
N22 = ceiling(pi2 * N2)
N12 = ceiling(r * N22)
N23 = N2 - N21
N13 = N1 - N11
cutpoint = c(0.15, 0.3, 0.45)
# Results for the trial of pyruvate kinase deficiency
results.power = tibble(
  pi1 = pi1,
  N11 = N11,
  N21 = N21,
  N12 = N12,
  N22 = N22,
  N13 = N13,
  N23 = N23
) %>%
  group_by_all() %>%
  reframe(
    cutpoint = cutpoint
  ) %>%
  group_by_all() %>%
  mutate(
    power = list(
      power.2in1.binary(p1, p2, N11, N21, N12, N22, N13, N23, cutpoint, alpha2, alpha3)
    )
  ) %>%
  ungroup() %>%
  select(-cutpoint) %>%
  unnest(power) %>%
  select(
    pi1, cutpoint, Test, Go.prob, Power.Total, ESS
  )
results.test.size = tibble(
  pi1 = pi1,
  N11 = N11,
  N21 = N21,
  N12 = N12,
  N22 = N22,
  N13 = N13,
  N23 = N23
) %>%
  group_by_all() %>%
  reframe(
    cutpoint = cutpoint
  ) %>%
  group_by_all() %>%
  mutate(
    test.size = list(
      test.size.2in1.binary(N11, N21, N12, N22, N13, N23, cutpoint, alpha2, alpha3)
    )
  ) %>%
  ungroup() %>%
  select(pi1, test.size) %>%
  unnest(test.size)
results = results.power %>%
  left_join(
    results.test.size, by = c('pi1', 'cutpoint', 'Test')
  )

# Table2
latex.table = results %>%
  mutate(
    Power.Testsize = paste(sprintf('%.3f', Power.Total), ' (', sprintf('%.4f', Test.size), ')', sep = ''),
    Go.prob = sprintf('%.2f', Go.prob),
    ESS = sprintf('%.1f', ESS)
  ) %>%
  select(pi1, cutpoint, Test, Go.prob, Power.Testsize, ESS) %>%
  pivot_wider(
    names_from = Test, values_from = Power.Testsize
  ) %>%
  rename(
    '$\\pi_{1}$' = pi1,
    '$C$' = cutpoint,
    '$p_{{\\rm Go}}$' = Go.prob
  ) %>%
  kbl(format = 'latex', booktabs = TRUE, linesep = '', escape = F, align = 'l') %>%
  landscape(margin = NULL) %>%
  add_header_above(c(' ' = 4, 'Power (Test size)' = 5)) %>%
  collapse_rows(columns = 1, latex_hline = 'major', row_group_label_position = 'first')

# Save latex table file
latex.table %>%
  save_kable('table2.tex')
