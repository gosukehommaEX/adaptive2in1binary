########################################################################### README ###########################################################################
## The code was written by Gosuke Homma (gosuke.homma@boehringer-ingelheim.com), and was tested using R version 4.4.2.
#
## The "sample.size.2in1.binary" function aims to calculate the required sample size
## for a 2-in-1 adaptive phase 2/3 design with binary endpoints.
#  The following five tests can be specified:
#  (1) The one-sided Pearson chi-squared test (Chisq)
#  (2) The Fisher exact test (Fisher)
#  (3) The Fisher mid-p test (Fisher-midP)
#  (4) The Z-pooled exact unconditional test (Z-pool)
#  (5) The Boschloo exact unconditional test (Boschloo)
#  Note: The function only covers one-sided tests.
#
## sample.size.2in1.binary has the following arguments.
#  p1:         true probability of responders for group 1
#  p2:         true probability of responders for group 2
#  r:          allocation ratio to group 1 (i.e., an allocation ratio of group 1:group 2 = r:1, r > 0)
#  pi1:        allocation ratio of total sample size to stage 1
#              (i.e., N21 = ceiling(pi1 * N), and N11 = ceiling(r * N21), where N is the total sample size)
#  pi2:        allocation ratio of total sample size to stage 2
#              (i.e., N22 = ceiling(pi2 * N), and N12 = ceiling(r * N22), where N is the total sample size)
#  (Note) to control type I error rate, pi2 should be less than or equal to 1 - pi1
#  cutpoint:   a cutpoint of the risk difference (the value should be within -1 < cutpoint <= 1)
#  alpha:      one-sided level of significance
#  tar.power:  target power
#
## sample.size.2in1.binary returns the following result in R tibble format
#  p1:         true probability of responders for group 1
#  p2:         true probability of responders for group 2
#  r:          allocation ratio to group 1 (i.e., an allocation ratio of group 1:group 2 = r:1, r > 0)
#  pi1:        allocation ratio of total sample size to stage 1
#  pi2:        allocation ratio of total sample size to stage 2
#  alpha:      one-sided level of significance
#  tar.power:  target type II error rate
#  Test:       Statistical testing approach
#  cutpoint:   a cutpoint of the risk difference
#  N1.trad:    required sample size of group 1 for a traditional design
#  N2.trad:    required sample size of group 2 for a traditional design
#  N.trad:     required total sample size for a traditional design
#  N1.2in1:    required sample size of group 1 for a 2-in-1 adaptive phase 2/3 design
#  N2.2in1:    required sample size of group 2 for a 2-in-1 adaptive phase 2/3 design
#  N.2in1:     required total sample size for a 2-in-1 adaptive phase 2/3 design
#  ESS:        expected sample size 
#  Go.prob:    Go probability of expanding to a Phase 3 trial
#  Power.2:    power of tests for Phase 2 part
#  Power.3:    power of tests for Phase 3 part
#  Power.Total:total power of tests (i.e., sum of powers for Phase 2 and 3 parts)
######################################################################### How to use #########################################################################
# sample.size.2in1.binary(
#   p1 = 0.4, p2 = 0.2, r = 2, pi1 = 0.2, pi2 = 0.5, cutpoint = 0.2, alpha = 0.025, tar.power = 0.8
# )
# # A tibble: 5 × 20
#      p1    p2     r   pi1   pi2 alpha tar.power Test        cutpoint N1.trad N2.trad N.trad N1.2in1 N2.2in1 N.2in1   ESS Go.prob Power.2 Power.3 Power.Total
#   <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>     <dbl> <fct>          <dbl>   <dbl>   <dbl>  <dbl>   <dbl>   <dbl>  <dbl> <dbl>   <dbl>   <dbl>   <dbl>       <dbl>
# 1   0.4   0.2     2   0.2   0.5 0.025       0.8 Chisq            0.2     124      62    186     162      81    243  211.   0.537   0.292   0.514       0.806
# 2   0.4   0.2     2   0.2   0.5 0.025       0.8 Fisher           0.2     136      68    204     182      91    273  235.   0.517   0.314   0.500       0.813
# 3   0.4   0.2     2   0.2   0.5 0.025       0.8 Fisher-midP      0.2     124      62    186     162      81    243  211.   0.537   0.292   0.514       0.806
# 4   0.4   0.2     2   0.2   0.5 0.025       0.8 Z-pool           0.2     144      72    216     190      95    285  244.   0.517   0.306   0.499       0.805
# 5   0.4   0.2     2   0.2   0.5 0.025       0.8 Boschloo         0.2     124      62    186     162      81    243  211.   0.537   0.292   0.512       0.804
#
# sample.size.2in1.binary(
#   p1 = c(0.8, 0.7, 0.6), p2 = c(0.4, 0.3, 0.2), r = 1, pi1 = 0.3, pi2 = 0.3, cutpoint = 0.5, alpha = 0.025, tar.power = 0.9
# )
# # A tibble: 15 × 20
#      p1    p2     r   pi1   pi2 alpha tar.power Test        cutpoint N1.trad N2.trad N.trad N1.2in1 N2.2in1 N.2in1   ESS Go.prob Power.2 Power.3 Power.Total
#   <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>     <dbl> <fct>          <dbl>   <dbl>   <dbl>  <dbl>   <dbl>   <dbl>  <dbl> <dbl>   <dbl>   <dbl>   <dbl>       <dbl>
# 1   0.6   0.2     1   0.3   0.3 0.025       0.9 Chisq            0.5      30      30     60      47      47     94  69.4   0.276   0.630   0.275       0.905
# 2   0.6   0.2     1   0.3   0.3 0.025       0.9 Fisher           0.5      33      33     66      54      54    108  78.5   0.262   0.646   0.262       0.908
# 3   0.6   0.2     1   0.3   0.3 0.025       0.9 Fisher-midP      0.5      30      30     60      47      47     94  69.4   0.276   0.628   0.275       0.903
# 4   0.6   0.2     1   0.3   0.3 0.025       0.9 Z-pool           0.5      30      30     60      47      47     94  69.4   0.276   0.625   0.275       0.901
# 5   0.6   0.2     1   0.3   0.3 0.025       0.9 Boschloo         0.5      31      31     62      47      47     94  69.4   0.276   0.625   0.275       0.900
# 6   0.7   0.3     1   0.3   0.3 0.025       0.9 Chisq            0.5      31      31     62      51      51    102  77.1   0.344   0.580   0.344       0.924
# 7   0.7   0.3     1   0.3   0.3 0.025       0.9 Fisher           0.5      37      37     74      61      61    122  87.7   0.255   0.666   0.255       0.922
# 8   0.7   0.3     1   0.3   0.3 0.025       0.9 Fisher-midP      0.5      31      31     62      54      54    108  78.7   0.268   0.643   0.268       0.911
# 9   0.7   0.3     1   0.3   0.3 0.025       0.9 Z-pool           0.5      33      33     66      54      54    108  78.7   0.268   0.643   0.268       0.911
# 10  0.7   0.3     1   0.3   0.3 0.025       0.9 Boschloo         0.5      33      33     66      54      54    108  78.7   0.268   0.643   0.268       0.911
# 11  0.8   0.4     1   0.3   0.3 0.025       0.9 Chisq            0.5      30      30     60      47      47     94  69.4   0.276   0.630   0.275       0.905
# 12  0.8   0.4     1   0.3   0.3 0.025       0.9 Fisher           0.5      33      33     66      54      54    108  78.5   0.262   0.646   0.262       0.908
# 13  0.8   0.4     1   0.3   0.3 0.025       0.9 Fisher-midP      0.5      30      30     60      47      47     94  69.4   0.276   0.628   0.275       0.903
# 14  0.8   0.4     1   0.3   0.3 0.025       0.9 Z-pool           0.5      30      30     60      47      47     94  69.4   0.276   0.625   0.275       0.901
# 15  0.8   0.4     1   0.3   0.3 0.025       0.9 Boschloo         0.5      31      31     62      47      47     94  69.4   0.276   0.625   0.275       0.900
##############################################################################################################################################################
source('sample.size.binary.R')
source('power.2in1.binary.R')
sample.size.2in1.binary = function(p1, p2, r, pi1, pi2, cutpoint, alpha, tar.power) {
  # Step 0: sample size calculation for a traditional design of a phase 3 trial ensuring power = tar.power
  info.trad = tibble(
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
  N.min = info.trad %>% 
    filter(N == min(N)) %>% 
    slice(1) %>% 
    select(N1, N2, N)
  N2 = N.min[['N2']]
  N1 = N.min[['N1']]
  N21 = ceiling(pi1 * N2)
  N11 = ceiling(r * N21)
  N22 = ceiling(pi2 * N2)
  N12 = ceiling(r * N22)
  N23 = N2 - N21
  N13 = N1 - N11
  # Step 1: power calculation for a 2-in-1 adaptive phase 2/3 design given initial sample sizes
  info.2in1 = power.2in1.binary(p1, p2, N11, N21, N12, N22, N13, N23, cutpoint, alpha) %>% 
    mutate(
      N1 = N1,
      N2 = N2,
      N = N1 + N2
    )
  Power = info.2in1 %>% 
    pull(Power.Total) %>% 
    min()
  # Step 2: sample size will be increased until calculated power satisfying the target power
  while(Power %<<% tar.power) {
    N2 = N2 + 1
    N1 = ceiling(r * N2)
    N = N1 + N2
    N21 = ceiling(pi1 * N2)
    N11 = ceiling(r * N21)
    N22 = ceiling(pi2 * N2)
    N12 = ceiling(r * N22)
    N23 = N2 - N21
    N13 = N1 - N11
    info.2in1 = info.2in1 %>% 
      bind_rows(
        power.2in1.binary(p1, p2, N11, N21, N12, N22, N13, N23, cutpoint, alpha) %>% 
          mutate(
            N1 = N1,
            N2 = N2,
            N = N1 + N2
          )
      )
    Power = info.2in1 %>% 
      filter(N == .env$N) %>% 
      pull(Power.Total) %>% 
      min()
  }
  # Step 3: print summary result
  result = info.2in1 %>% 
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