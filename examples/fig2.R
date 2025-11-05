########################################################################################################################
## Script for Figure 2: Test size and power under balanced allocation (r=1)
##
## Description:
##   This script generates Figure 2 from Homma and Yoshida (20XX), which shows
##   test size and power for each test, as well as Go probability and ESS under H1,
##   evaluated using the proposed formulas under alpha_2 = alpha_3 = 0.025,
##   pi_1 = pi_2 = 0.3, and balanced allocation ratio (r=1) over -1 <= C <= 1.
##
## Design parameters:
##   - Case A: p1=0.6, p2=0.2, N3*=60
##   - Case B: p1=0.6, p2=0.3, N3*=120
##   - Case C: p1=0.6, p2=0.4, N3*=300
##   - Balanced allocation ratio (r=1)
##   - pi1 = pi2 = 0.3
##   - alpha2 = alpha3 = 0.025
##   - rho ranges from 1 to 1.5
##
## Output:
##   - fig2.eps: Publication-quality figure (800 dpi)
##
## Reference:
##   Homma G and Yoshida T (20XX). A 2-in-1 adaptive design for binary endpoints.
##   Statistics in Medicine (under review).
##
## Usage:
##   source("fig2.R")
##
## Note:
##   This script requires the following source files in the same directory:
##   - test.size.2in1.binary.R
##   - power.2in1.binary.R
##
##   Required packages: ggplot2, patchwork
########################################################################################################################
library(ggplot2)
library(patchwork)
library(adaptive2in1binary)
p1 = c(0.6, 0.6, 0.6)
p2 = c(0.2, 0.3, 0.4)
r = 1
N2 = c(30, 60, 150)
N1 = ceiling(r * N2)
N = N1 + N2
pi1 = 0.3
pi2 = 0.3
N21 = ceiling(pi1 * N2)
N11 = ceiling(r * N21)
N22 = ceiling(pi2 * N2)
N12 = ceiling(r * N22)
N23 = N2 - N21
N13 = N1 - N11
range.rho = c(1, 1.5)
alpha2 = 0.025
alpha3 = 0.025
targe.alpha = 0.025
# Results of power and test size
results = tibble(
  p1 = p1,
  p2 = p2,
  N = N,
  N11 = N11,
  N21 = N21,
  N12 = N12,
  N22 = N22,
  N13 = N13,
  N23 = N23,
  Case = factor(
    paste0('Case ', LETTERS[seq(length(N))]),
    levels = paste0('Case ', LETTERS[seq(length(N))])
  )
) %>%
  group_by(Case) %>%
  mutate(
    power = list(
      power.2in1.binary(p1, p2, N11, N21, N12, N22, N13, N23, NA, alpha2, alpha3)
    ),
    test.size = list(
      test.size.2in1.binary(N11, N21, N12, N22, N13, N23, NA, alpha2, alpha3)
    )
  ) %>%
  ungroup() %>%
  select(-p1, -p2)

# Plot Go probability and ESS
scale = max(N1 + N2)
data.go.prob.ESS = results %>%
  select(-test.size) %>%
  unnest(power) %>%
  group_by(Case, cutpoint) %>%
  reframe(
    min.cutpoint = range.rho[1] * (p1 - p2)[n()],
    max.cutpoint = range.rho[2] * (p1 - p2)[n()],
    Go.prob = Go.prob[n()],
    ESS = ESS[n()]
  )
fig.go.prob.ESS = data.go.prob.ESS %>%
  mutate(
    Go.prob = if_else(Go.prob > 1, 1, Go.prob),
    Label = 'Go.Prob & ESS'
  ) %>%
  ggplot(aes(x = cutpoint, y = Go.prob)) +
  geom_rect(
    aes(xmin = min.cutpoint, xmax = max.cutpoint),
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.1,
    fill = 'gray'
  ) +
  geom_line(aes(linetype = factor('Go probability', levels = 'Go probability')), linewidth = 1.2, colour = 'grey50') +
  geom_line(aes(y = ESS/scale, linetype = 'ESS'), linewidth = 1.2, colour = 'grey50') +
  facet_grid(
    Case ~ Label
  ) +
  theme_bw() +
  scale_x_continuous(
    breaks = seq(-1, 1, by = 0.5)
  ) +
  scale_y_continuous(
    sec.axis = sec_axis(~ . * scale, name = 'ESS'),
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.1)
  ) +
  labs(
    x = expression(italic(C)),
    y = 'Go probability'
  ) +
  theme(
    strip.text.x = element_text(size = 20),
    strip.text.y = element_blank(),
    text = element_text(size = 20),
    panel.spacing = unit(-0.1, 'lines'),
    legend.key.width = unit(2, 'cm'),
    legend.text = element_text(size = 20),
    legend.title = element_blank(),
    legend.background = element_rect(fill = 'white', color = 'black'),
    legend.position = 'bottom'
  )

# Plot test size
fig.test.size = results %>%
  select(-power) %>%
  unnest(test.size) %>%
  left_join(
    data.go.prob.ESS %>% select(Case, cutpoint, min.cutpoint, max.cutpoint),
    by = c('Case', 'cutpoint')
  ) %>%
  mutate(
    Label = as.factor('Test size')
  ) %>%
  ggplot(aes(x = cutpoint, y = Test.size)) +
  geom_rect(
    aes(xmin = min.cutpoint, xmax = max.cutpoint),
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.1,
    fill = 'gray'
  ) +
  geom_line(aes(color = Test, linetype = Test), linewidth = 1.2) +
  facet_grid(
    Case ~ Label
  ) +
  theme_bw() +
  geom_hline(
    yintercept = targe.alpha,
    color = 'gray',
    linetype = 'longdash',
    linewidth = 1.2
  ) +
  scale_linetype_manual(
    values = c(
      'solid',
      'dotted',
      'twodash',
      'dashed',
      'dotdash'
    )
  ) +
  scale_color_manual(
    values = c(
      'Chisq' = 'black',
      'Fisher' = 'darkorchid1',
      'Fisher-midP' = 'blue',
      'Z-pool' = 'deepskyblue',
      'Boschloo' = 'chartreuse2'
    ),
    labels =  c('Chisq', 'Fisher', 'Fisher-midP', 'Z-pool', 'Boschloo')
  ) +
  scale_x_continuous(
    breaks = seq(-1, 1, by = 0.5)
  ) +
  scale_y_continuous(
    limits = c(0.005, 0.045),
    breaks = seq(0.005, 0.045, length = 5)
  ) +
  labs(
    x = expression(italic(C)),
    y = 'Test size'
  ) +
  theme(
    strip.text.x = element_text(size = 20),
    strip.text.y = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 20),
    panel.spacing = unit(-0.1, 'lines'),
    legend.key.width = unit(2, 'cm'),
    legend.text = element_text(size = 20),
    legend.background = element_rect(fill = 'white', color = 'black'),
    legend.title = element_text(hjust = 0.5),
    legend.title.position = 'top',
    legend.position = 'bottom'
  )

# Plot power
fig.power = results %>%
  select(-test.size) %>%
  unnest(power) %>%
  left_join(
    data.go.prob.ESS %>% select(Case, cutpoint, min.cutpoint, max.cutpoint),
    by = c('Case', 'cutpoint')
  ) %>%
  mutate(
    Label = 'Power'
  ) %>%
  ggplot(aes(x = cutpoint, y = Power.Total)) +
  geom_rect(
    aes(xmin = min.cutpoint, xmax = max.cutpoint),
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.1,
    fill = 'gray'
  ) +
  geom_line(aes(color = Test, linetype = Test), linewidth = 1.2) +
  facet_grid(
    Case ~ Label
  ) +
  theme_bw() +
  scale_linetype_manual(
    values = c(
      'solid',
      'dotted',
      'twodash',
      'dashed',
      'dotdash'
    )
  ) +
  scale_color_manual(
    values = c(
      'Chisq' = 'black',
      'Fisher' = 'darkorchid1',
      'Fisher-midP' = 'blue',
      'Z-pool' = 'deepskyblue',
      'Boschloo' = 'chartreuse2'
    ),
    labels =  c('Chisq', 'Fisher', 'Fisher-midP', 'Z-pool', 'Boschloo')
  ) +
  scale_x_continuous(
    breaks = seq(-1, 1, by = 0.5)
  ) +
  scale_y_continuous(
    limits = c(0.5, 1),
    breaks = seq(0.5, 1, length = 6)
  ) +
  labs(
    x = expression(italic(C)),
    y = 'Power'
  ) +
  theme(
    strip.text.x = element_text(size = 20),
    strip.text.y = element_text(size = 20),
    axis.title.y = element_blank(),
    text = element_text(size = 20),
    panel.spacing = unit(-0.1, 'lines'),
    legend.key.width = unit(2, 'cm'),
    legend.text = element_text(size = 20),
    legend.background = element_rect(fill = 'white', color = 'black'),
    legend.title = element_text(hjust = 0.5),
    legend.title.position = 'top',
    legend.position = 'bottom'
  )

# Save figures
ggsave(
  file = paste('fig2', 'eps', sep = '.'),
  plot = fig.go.prob.ESS +
    plot_spacer() +
    fig.test.size +
    fig.power +
    plot_layout(nrow = 1, widths = c(2, 0.1, 3, 3), guides = 'collect') & theme(legend.position = 'bottom'),
  device = cairo_pdf,
  dpi = 800,
  width = 16,
  height = 16
)
