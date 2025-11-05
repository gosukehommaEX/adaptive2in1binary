########################################################################################################################
## Script for Figure 5: Test size and power with different pi1, pi2, C combinations (r=2)
##
## Description:
##   This script generates Figure 5 from Homma and Yoshida (20XX), which shows
##   test size and power for each test with different combinations of pi1 (=0.2, 0.3, 0.4),
##   pi2 (=0.2, 0.3, 0.4), and C (=0.5*delta, delta, 1.5*delta) for the imbalanced
##   allocation ratio scenario (r=2) under alpha_2 = alpha_3 = 0.025.
##
## Design parameters:
##   - Case A: p1=0.6, p2=0.2, N3*=60
##   - Case B: p1=0.6, p2=0.3, N3*=120
##   - Case C: p1=0.6, p2=0.4, N3*=300
##   - Imbalanced allocation ratio (r=2)
##   - pi1, pi2 = 0.2, 0.3, 0.4 (all combinations)
##   - alpha2 = alpha3 = 0.025
##   - rho = 0.5, 1, 1.5
##
## Output:
##   - fig5.eps: Publication-quality heatmap figure (800 dpi)
##
## Reference:
##   Homma G and Yoshida T (20XX). A 2-in-1 adaptive design for binary endpoints.
##   Statistics in Medicine (under review).
##
## Usage:
##   source("fig5.R")
##
## Note:
##   This script requires the following source files in the same directory:
##   - power.2in1.binary.R
##   - test.size.2in1.binary.R
##
##   Required packages: ggplot2, ggh4x, scales, patchwork
########################################################################################################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggh4x)
library(scales)
library(patchwork)
library(adaptive2in1binary)
p1 = c(0.6, 0.6, 0.6)
p2 = c(0.2, 0.3, 0.4)
r = 2
rho = c(0.5, 1, 1.5)
N2 = c(20, 40, 100)
pi1 = c(0.2, 0.3, 0.4)
pi2 = c(0.2, 0.3, 0.4)
alpha2 = 0.025
alpha3 = 0.025
# Results of power and test size
results.fig5 = tibble(
  p1 = p1,
  p2 = p2,
  r = r,
  N2 = N2,
  N1 = r * N2,
  N = N1 + N2,
  Case = factor(
    paste0('Case ', LETTERS[seq(length(N))]),
    levels = paste0('Case ', LETTERS[seq(length(N))])
  )
) %>%
  group_by_all() %>%
  reframe(
    rho = rho,
    cutpoint = rho * (p1 - p2)
  ) %>%
  group_by_all() %>%
  reframe(
    pi1 = pi1,
    N21 = ceiling(pi1 * N2),
    N11 = ceiling(r * N21)
  ) %>%
  group_by_all() %>%
  reframe(
    pi2 = pi2,
    N22 = ceiling(pi2 * N2),
    N12 = ceiling(r * N22),
    N23 = N2 - N21,
    N13 = N1 - N11
  ) %>%
  group_by_all() %>%
  mutate(
    power = list(
      power.2in1.binary(p1, p2, N11, N21, N12, N22, N13, N23, cutpoint, alpha2, alpha3)
    ),
    test.size = list(
      test.size.2in1.binary(N11, N21, N12, N22, N13, N23, cutpoint, alpha2, alpha3)
    )
  ) %>%
  ungroup()

# Plot test size
fig5.test.size = results.fig5 %>%
  select(Case, rho, pi1, pi2, test.size) %>%
  unnest(test.size) %>%
  select(Case, rho, pi1, pi2, cutpoint, Test, Test.size) %>%
  mutate(
    cutpoint = as.factor(paste0('C = ', cutpoint)),
    y.axis.lab = factor(paste0(Case, '.', rho), levels = unique(paste0(Case, '.', rho))),
    Label = 'Test size'
  ) %>%
  ggplot(aes(pi1, pi2, fill = Test.size)) +
  facet_nested(
    y.axis.lab ~ Label + Test,
    nest_line = element_line(colour = 'black'),
    labeller = labeller(
      rho  = as_labeller(rho, label_parsed)
    )
  ) +
  theme_bw() +
  labs(x = expression(pi[1]), y = expression(pi[2])) +
  geom_tile(aes(fill = Test.size), color = 'gray50') +
  labs(fill = 'Test size') +
  scale_fill_stepsn(
    colors = c('#fff5f0', '#fee0d2', '#fcbba1', '#a50f15'),
    breaks = c(0.015, 0.02, 0.025),
    limits = c(0, 1),
    labels = label_number(accuracy = 0.005),
    values = rescale(c(0, 0.015, 0.02, 0.025, 1))
  ) +
  guides(
    color = 'none'
  ) +
  scale_x_continuous(
    breaks = seq(0.2, 0.4, by = 0.1)
  ) +
  scale_y_continuous(
    breaks = seq(0.2, 0.4, by = 0.1)
  ) +
  theme(
    plot.margin = margin(0, 0, 0, 0, 'pt'),
    strip.text.x = element_text(size = 16),
    strip.text.y = element_blank(),
    text = element_text(size = 20),
    panel.spacing = unit(-0.1, 'lines'),
    legend.key.width = unit(2, 'cm'),
    legend.text = element_text(size = 20),
    legend.title = element_text(hjust = 0.5),
    legend.title.position = 'top',
    legend.position = 'bottom'
  )

# Plot power
fig5.power = results.fig5 %>%
  select(Case, rho, pi1, pi2, power) %>%
  unnest(power) %>%
  select(Case, rho, pi1, pi2, cutpoint, Test, Power.Total) %>%
  group_by(Case, rho) %>%
  mutate(
    rho = factor(
      deparse(bquote(italic(C)==.(unique(if_else(rho == 1, '', as.character(rho))))~ delta)),
      levels = deparse(bquote(italic(C)==.(unique(if_else(rho == 1, '', as.character(rho))))~ delta))
    ),
    Label = 'Power'
  ) %>%
  ggplot(aes(pi1, pi2, fill = Power.Total)) +
  facet_nested(
    Case + rho ~ Label + Test,
    nest_line = element_line(colour = 'black'),
    labeller = labeller(
      rho  = as_labeller(rho, label_parsed)
    )
  ) +
  theme_bw() +
  labs(x = expression(pi[1]), y = expression(pi[2])) +
  geom_tile(aes(fill = Power.Total), color = 'gray50') +
  labs(fill = 'Power') +
  scale_fill_stepsn(
    colors = c('#ffffe5', '#f7fcb9', '#d9f0a3', '#006837'),
    breaks = c(0.6, 0.7, 0.8),
    limits = c(0, 1),
    labels = label_number(accuracy = 0.1),
    values = rescale(c(0, 0.6, 0.7, 0.8, 1))
  ) +
  guides(
    color = 'none'
  ) +
  scale_x_continuous(
    breaks = seq(0.2, 0.4, by = 0.1)
  ) +
  scale_y_continuous(
    breaks = seq(0.2, 0.4, by = 0.1)
  ) +
  theme(
    plot.margin = margin(0, 0, 0, 0, 'pt'),
    strip.text.x = element_text(size = 16),
    strip.text.y = element_text(size = 16),
    text = element_text(size = 20),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.spacing = unit(-0.1, 'lines'),
    legend.key.width = unit(2, 'cm'),
    legend.text = element_text(size = 20),
    legend.title = element_text(hjust = 0.5),
    legend.title.position = 'top',
    legend.position = 'bottom'
  )

# Save figures
ggsave(
  file = paste('fig5', 'eps', sep = '.'),
  plot = fig5.test.size + fig5.power +
    plot_layout(
      nrow = 1,
      widths = c(1, 1),
      guides = 'collect',
      axis_titles = 'collect'
    ) & theme(legend.position = 'bottom'),
  dpi = 800,
  width = 16,
  height = 16
)
