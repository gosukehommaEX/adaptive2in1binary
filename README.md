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
# Install pak if you don't have it
# install.packages("pak")

# Install the package
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

### Method 1: Using the installed package (Recommended)

After installing the package, you can reproduce the manuscript results by cloning the repository and running the example scripts:

```r
# 1. Install the package
pak::pak("gosukehommaEX/adaptive2in1binary")

# 2. Install additional packages required for examples
install.packages(c("dplyr", "tidyr", "ggplot2", "patchwork", "ggh4x", "scales", "kableExtra"))

# 3. Clone or download the repository
# Option A: Using Git (in terminal/command prompt)
#   git clone https://github.com/gosukehommaEX/adaptive2in1binary.git
# Option B: Download ZIP file
#   Visit: https://github.com/gosukehommaEX/adaptive2in1binary
#   Click "Code" → "Download ZIP" → Extract the ZIP file

# 4. Set working directory to the examples folder (IMPORTANT!)
# Replace with your actual path
setwd("path/to/adaptive2in1binary/examples")
# Example on Windows: setwd("C:/Users/YourName/Downloads/adaptive2in1binary/examples")
# Example on Mac/Linux: setwd("~/Downloads/adaptive2in1binary/examples")

# 5. Load required packages (important!)
library(dplyr)
library(tidyr)
library(adaptive2in1binary)

# 6. Run individual scripts
source("fig2.R")    # Generates Figure 2 (fig2.eps)
source("fig3.R")    # Generates Figure 3 (fig3.eps)
source("fig4.R")    # Generates Figure 4 (fig4.eps)
source("fig5.R")    # Generates Figure 5 (fig5.eps)
source("fig6.R")    # Generates Figure 6 (fig6.eps)
source("table1.R")  # Generates Table 1 (table1.tex)
source("table2.R")  # Generates Table 2 (table2.tex)
```

**Note**: The scripts must be run from the `examples/` folder because they reference each other and generate output files in the current directory.

### Method 2: Direct execution from GitHub

Alternatively, you can download the repository as a ZIP file from GitHub, extract it, and run the scripts as described above.

### Required Packages for Examples

The example scripts require the following additional packages:
- **For data manipulation**: `dplyr`, `tidyr` (required for all examples)
- **For figures**: `ggplot2`, `patchwork`, `ggh4x`, `scales`
- **For tables**: `kableExtra`

Install them all at once:
```r
install.packages(c("dplyr", "tidyr", "ggplot2", "patchwork", "ggh4x", "scales", "kableExtra"))
```

**Important Note**: Although `dplyr` and `tidyr` are listed as dependencies of `adaptive2in1binary`, you must load them explicitly when running the example scripts because the scripts use these packages' functions directly (e.g., `%>%`, `tibble()`, `group_by()`).

### Output Files

Running the example scripts will generate:
- **Figures**: `fig2.eps`, `fig3.eps`, `fig4.eps`, `fig5.eps`, `fig6.eps` (800 dpi EPS format)
- **Tables**: `table1.tex`, `table2.tex` (LaTeX format)

## Citation

If you use this code in your research, please cite:

```
Homma, G. and Yoshida, T. (2025). A 2-in-1 adaptive design for binary endpoints. 
Statistics in Medicine (under review).
```

## License

MIT License - see [LICENSE.md](LICENSE.md) file for details
