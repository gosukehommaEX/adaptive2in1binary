# adaptive2in1binary

<!-- badges: start -->
<!-- badges: end -->

## Overview

This repository contains R functions and reproducible code for the paper:

**"A 2-in-1 adaptive design for binary endpoints"**  
Gosuke Homma and Takuma Yoshida  
*Statistics in Medicine* (under review)

The 2-in-1 adaptive design allows a clinical trial to maintain a small trial or expand to a large trial adaptively based on decisions made at an interim analysis. This package provides functions to analytically evaluate the type I error rate and power for the 2-in-1 adaptive design with binary endpoints, without requiring Monte Carlo simulations.

## Key Features

- **Exact analytical calculations**: Evaluate type I error rate and power without simulations
- **Multiple statistical tests**: Support for chi-square test, Fisher's exact test, Fisher's exact test with mid-P correction, Z-pooled test, and Boschloo's exact unconditional test
- **Sample size calculations**: Functions for determining required sample sizes
- **Reproducible research**: All code used to generate figures and tables in the manuscript

## Installation

You can install the development version from GitHub:

``` r
# install.packages("pak")
pak::pak("gosukehommaEX/adaptive2in1binary")
```

## Repository Structure

```
adaptive2in1binary/
├── R/                          # Core R functions
│   ├── power.2in1.binary.R     # Calculate power for 2-in-1 design
│   ├── test.size.2in1.binary.R # Calculate type I error rate
│   ├── sample.size.2in1.binary.R # Sample size calculations
│   └── ...
├── examples/                   # Scripts to reproduce manuscript results
│   ├── fig2.R                  # Generate Figure 2
│   ├── fig3.R                  # Generate Figure 3
│   ├── fig4.R                  # Generate Figure 4
│   ├── fig5.R                  # Generate Figure 5
│   ├── fig6.R                  # Generate Figure 6
│   ├── table1.R                # Generate Table 1
│   └── table2.R                # Generate Table 2
└── README.md
```

## Basic Usage

### Calculate power for a 2-in-1 design

```r
library(adaptive2in1binary)

# Design parameters
p1 <- 0.6  # Response rate in treatment group
p2 <- 0.4  # Response rate in control group
N11 <- 30  # Sample size for treatment at stage 1
N21 <- 30  # Sample size for control at stage 1
N12 <- 30  # Sample size for treatment at stage 2
N22 <- 30  # Sample size for control at stage 2
N13 <- 40  # Sample size for treatment at stage 3
N23 <- 40  # Sample size for control at stage 3
cutpoint <- 0.1  # Decision cutpoint
alpha2 <- 0.025  # Significance level at stage 2
alpha3 <- 0.025  # Significance level at stage 3

# Calculate power
result <- power.2in1.binary(p1, p2, N11, N21, N12, N22, N13, N23, 
                            cutpoint, alpha2, alpha3)
print(result)
```

### Calculate type I error rate

```r
# Calculate type I error rate (test size) under H0
test_size <- test.size.2in1.binary(N11, N21, N12, N22, N13, N23, 
                                   cutpoint, alpha2, alpha3)
print(test_size)
```

## Reproducing Manuscript Results

To reproduce the figures and tables from the manuscript:

```r
# Navigate to the examples folder and run individual scripts
source("examples/fig2.R")  # Generates Figure 2
source("examples/table1.R")  # Generates Table 1
# etc.
```

Note: The example scripts require `ggplot2` and `patchwork` packages for visualization.

## Citation

If you use this code in your research, please cite:

```
Homma, G. and Yoshida, T. (2025). A 2-in-1 adaptive design for binary endpoints. 
Statistics in Medicine (under review).
```

## Contact

- Gosuke Homma: gosuke.homma@astellas.com
- Takuma Yoshida: yoshida@sci.kagoshima-u.ac.jp

## License

MIT License - see [LICENSE.md](LICENSE.md) file for details
