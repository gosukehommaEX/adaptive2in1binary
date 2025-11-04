############################################################# README #############################################################
## The code was written by Gosuke Homma (gosuke.homma@boehringer-ingelheim.com), and was tested using R version 4.4.2.
#
## The "test.size.2in1.binary" function aims to calculate test size (the maximum type I error rate over the entire range of theta) 
## using closed-form formulas for a 2-in-1 adaptive phase 2/3 design with binary endpoints.
#  The following five tests results can be shown simultaneously:
#  (1) The one-sided Pearson chi-squared test (Chisq)
#  (2) The Fisher exact test (Fisher)
#  (3) The Fisher mid-p test (Fisher-midP)
#  (4) The Z-pooled exact unconditional test (Z-pool)
#  (5) The Boschloo exact unconditional test (Boschloo)
#  Note: The function only covers one-sided tests.
#
## test.size.2in1.binary has the following arguments.
#  N11:        sample size for group 1 at stage 1
#  N21:        sample size for group 2 at stage 1
#  N12:        sample size for group 1 at stage 2 (i.e., sample size of a Phase 2 part is N11 + N12)
#  N22:        sample size for group 2 at stage 2 (i.e., sample size of a Phase 2 part is N21 + N22)
#  N13:        sample size for group 1 at stage 3 (i.e., sample size of a Phase 2 part is N11 + N13)
#  N23:        sample size for group 2 at stage 3 (i.e., sample size of a Phase 2 part is N21 + N23)
#  cutpoint:   a cutpoint of the risk difference (user can choose NA or a single value)
#              (Note 1) if user select NA, results with all possible cutpoint values will be returned.
#              (Note 2) if user select a specified value, the value should be within -1 < cutpoint <= 1.
#  alpha:      one-sided level of significance
# 
## test.size.2in1.binary returns the following result in R tibble format
#  cutpoint:   cutpoint(s) of the risk difference over x11(=0,...,N11) and x21(=0,...,N21)
#  Test:       names of the statistical tests
#  Test.size:  test size (the maximum type I error rate over the entire range of theta)
########################################################### How to use ###########################################################
# test.size.2in1.binary(
#   N11 = 40, N21 = 20, N12 = 50, N22 = 25, N13 = 100, N23 = 50, cutpoint = NA, alpha = 0.025
# )
# # A tibble: 405 × 3
#   cutpoint Test        Test.size
#      <dbl> <fct>           <dbl>
# 1   -1     Chisq          0.0402
# 2   -1     Fisher         0.0180
# 3   -1     Fisher-midP    0.0273
# 4   -1     Z-pool         0.0245
# 5   -1     Boschloo       0.0247
# 6   -0.975 Chisq          0.0402
# 7   -0.975 Fisher         0.0180
# 8   -0.975 Fisher-midP    0.0273
# 9   -0.975 Z-pool         0.0245
# 10  -0.975 Boschloo       0.0247
# # ℹ 395 more rows
# # ℹ Use `print(n = ...)` to see more rows
#
# test.size.2in1.binary(
#   N11 = 40, N21 = 20, N12 = 50, N22 = 25, N13 = 100, N23 = 50, cutpoint = 0.2, alpha = 0.025
# )
# # A tibble: 5 × 3
#   cutpoint Test        Test.size
#      <dbl> <fct>           <dbl>
# 1      0.2 Chisq          0.0411
# 2      0.2 Fisher         0.0143
# 3      0.2 Fisher-midP    0.0256
# 4      0.2 Z-pool         0.0203
# 5      0.2 Boschloo       0.0212
##################################################################################################################################
source('power.2in1.binary.R')
test.size.2in1.binary = function(N11, N21, N12, N22, N13, N23, cutpoint, alpha) {
  
  # Calculate test size
  test.size = power.2in1.binary(
    seq(0, 1, l = 100), seq(0, 1, l = 100), N11, N21, N12, N22, N13, N23, cutpoint, alpha
  ) %>% 
    select(cutpoint, Test, Power.Total) %>% 
    group_by(cutpoint, Test) %>%
    reframe(
      Test.size = max(Power.Total)
    )
  
  # Return output
  return(test.size)
}