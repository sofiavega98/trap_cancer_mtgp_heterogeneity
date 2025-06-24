# Lymphoma Subsets Analysis: Treatment Effects by Demographic Groups

This repository contains the complete analysis code for the paper "MTGP_Lymphoma_Strata_paper.pdf", which examines treatment effects on lymphoma incidence across different demographic subgroups using Multi-Task Gaussian Process (MTGP) models.

## Overview

This study analyzes the impact of a 1995 treatment intervention on lymphoma incidence rates in California and Connecticut counties, comparing effects across demographic subgroups including age, sex, race, Hispanic ethnicity, and lymphoma type.

### Key Features
- **Study Period**: 1988-2003 (16 years)
- **Treatment Year**: 1995
- **Age Range**: 0-29 years
- **Geographic Scope**: SEER registry counties
- **Demographic Subgroups**: Age, Sex, Race, Hispanic Ethnicity, Lymphoma Type

##️ Project Structure

```
├── README.md                           # This file
├── config.R                            # Configuration file with paths and parameters
├── run_all_analyses.R                  # Master script to run entire analysis
├── run_sensitivity_analysis.R          # Sensitivity analysis script
├── data_checks.R                       # Data quality validation and consistency checks
├── MTGP_Lymphoma_Strata_paper.pdf      # Main paper
├── mt_gp_pois_multigroup.stan          # Stan model specification
├── functions.R                         # Helper functions for ATT calculation
├── 0_cleaning_data_*.R                 # Data cleaning scripts by demographic subset
├── 0_create_descriptive_statistics.R   # Descriptive statistics and map (Table 1, Figure S2)
├── 1_model_run_GSC_*.R                 # GSC model fitting scripts
├── 2_model_run_MTGP_*.R                # MTGP model fitting scripts
├── 3_export_results.R                  # Results generation and table/figure creation
├── SA/                                 # Sensitivity analysis scripts
│   ├── 0_cleaning_data_*.R
│   ├── 1_model_run_MC_*.R
│   ├── 2_model_run_MTGP_*.R
│   └── 3_export_results.R
├── data/                               # Data directory (create and add your data here)
└── results/                            # Output directory (created when running)
    ├── tables/                         # LaTeX tables
    └── figures/                        # PNG figures
```

## Quick Start

### Prerequisites

1. **R** (version 4.0 or higher)
2. **Required R packages**:
   ```r
   install.packages(c("tidyverse", "ggplot2", "readr", "gsynth", "rstan", 
                      "kableExtra", "cowplot", "posterior", "usmap"))
   ```

3. **Stan** (for MTGP models):
   ```r
   install.packages("rstan")
   ```

### Data Requirements

Before running the analysis, ensure you have the following data files in the `data/` directory:

- `data/lymphoma_full.csv` - Raw SEER lymphoma data
- `data/seer_population.RData` - SEER population data

### Configuration

The project uses a centralized configuration system for easy setup and portability:

1. **Edit `config.R`** to set your base directory:
   ```r
   # In config.R, modify this line if needed:
   BASE_DIR <- getwd()  # Uses current working directory (recommended)
   ```

2. **All paths are automatically generated** relative to `BASE_DIR`:
   - `DATA_DIR` = `BASE_DIR/data`
   - `RESULTS_DIR` = `BASE_DIR/results`
   - `CODE_DIR` = `BASE_DIR/code`
   - `SA_DIR` = `BASE_DIR/SA`

3. **The configuration system ensures**:
   - All scripts use consistent paths
   - Easy portability across different systems
   - No hardcoded absolute paths in the codebase
   - Centralized parameter management

### Running the Complete Analysis

1. **Clone or download this repository**
2. **Set your working directory** to the repository root:
   ```r
   setwd("path/to/your/Lymphoma_Subsets/")
   ```
3. **Add your data files** to the `data/` directory
4. **Run the master script**:
   ```r
   source("run_all_analyses.R")
   ```

This will execute the complete analysis pipeline and generate all tables and figures for the paper.

## Analysis Pipeline

### Step 0: Data Quality Checks (`data_checks.R`)

**Purpose**: Validate data consistency and quality across all demographic subsets.

**Key Validations**:
- Cross-check of case counts across demographic dimensions
- Verification that demographic breakdowns sum to totals (e.g., young + old = total)
- Treatment/control group balance checks
- Population denominator validation
- Small population identification for sensitivity analysis
- Missingness assessment across demographic groups

**Outputs**:
- Case counts validation table: `results/tables/case_counts.tex`
- Grouped case counts table: `results/tables/case_counts_grouped.tex`
- Console output with data quality summaries

### Step 1: Data Cleaning (`0_cleaning_data_*.R`)

**Purpose**: Prepare SEER lymphoma data for analysis by demographic subsets.

**Scripts**:
- `0_cleaning_data_total.R` - Overall population
- `0_cleaning_data_age.R` - Age groups (0-19, 20-29)
- `0_cleaning_data_sex.R` - Sex (Male, Female)
- `0_cleaning_data_race.R` - Race (White, Black, Other)
- `0_cleaning_data_hisp.R` - Hispanic status (Hispanic, Non-Hispanic)
- `0_cleaning_data_type.R` - Lymphoma type (Hodgkin, Non-Hodgkin)

**Key Operations**:
- Filter to study period (1988-2003) and age range (0-29)
- Handle Hawaii county FIPS code changes
- Create treatment indicators (CA/CT counties after 1995)
- Merge with population data for rate calculations
- Aggregate cases by county and year
- Use configuration-based paths for all file operations

**Outputs**: Cleaned datasets saved as `.RData` files in `DATA_DIR`

### Step 2: Descriptive Statistics (`0_create_descriptive_statistics.R`)

**Purpose**: Create descriptive statistics and visualizations for the paper.

**Outputs**:
- **Table 1**: Aggregate lymphoma rates by demographic group, treatment status, and time period
- **Figure S2**: Map showing treated vs control counties
- Additional descriptive tables for reference

**Key Features**:
- Calculates aggregate rates pre- and post-treatment
- Separates results by treatment status (treated vs control)
- Creates publication-ready LaTeX tables and PNG figures
- Uses configuration-based paths for all outputs

### Step 3: GSC Model Fitting (`1_model_run_GSC_*.R`)

**Purpose**: Fit Generalized Synthetic Control models as baseline comparison.

**Model Specification**:
- Estimator: Interactive Fixed Effects (IFE)
- Outcome: Lymphoma incidence rate per 100,000 population
- Bootstrap samples: 1000
- Cross-validation: Yes (for rank selection)

**Outputs**: GSC model fits saved as `.RData` files in `RESULTS_DIR`

### Step 4: MTGP Model Fitting (`2_model_run_MTGP_*.R`)

**Purpose**: Fit Multi-Task Gaussian Process models for main treatment effect estimation.

**Model Specification**:
- Likelihood: Poisson (for count data)
- Prior: Multi-task Gaussian Process
- Latent functions: 15 for time trends, 10 for unit effects
- MCMC: 1000 iterations, 500 warmup, 1 chain
- Treatment: California and Connecticut counties after 1995

**Outputs**: Stan model fits saved as `.RData` files in `RESULTS_DIR`

### Step 5: Results Generation (`3_export_results.R`)

**Purpose**: Generate all tables and figures for the paper.

**Tables Generated**:
- **Table 2**: Main MTGP treatment effect estimates (`mtgp_table.tex`)
- **Table 3**: Treatment effect ratios (`mtgp_table_ratio.tex`)
- **Table 4**: Sum of treatment effects across treated units/time periods (`mtgp_table_sum_TE.tex`)
- **Table S4**: MCMC convergence diagnostics - R-hat values (`rhat_table.tex`)
- **Table S5**: Sensitivity analysis GSC results (`gsc_table.tex`)

**Figures Generated**:
- **Figure 1**: Treatment effects over time for all demographic subsets (`ATT_plot_overall_trt1995.png`)
- **Figures S4-S7**: Treatment effects by demographic subset (`ATT_plot_[subset]_trt1995.png`)

## Key Results

The analysis produces treatment effect estimates for:

1. **Overall Population**: Average treatment effect across all demographic groups
2. **Lymphoma Type**: Hodgkin vs Non-Hodgkin lymphoma
3. **Sex**: Male vs Female
4. **Race**: White, Black, Other
5. **Hispanic Ethnicity**: Hispanic vs Non-Hispanic
6. **Age Groups**: Ages 0-19 vs Ages 20-29

Each estimate includes:
- Average Treatment Effect (ATT)
- 95% confidence intervals
- Separate estimates for California, Connecticut, and overall

## Customization

### Modifying Analysis Parameters

Key parameters can be modified in `config.R`:

```r
# Study period
YEARS <- c(1988:2003)

# Treatment year
TREATED_YEAR <- 1995

# Analysis date (auto-generated)
DATE <- format(Sys.Date(), "%m%d%Y")

# Base directory (set automatically)
BASE_DIR <- getwd()
```

### Running Individual Components

You can run individual steps by sourcing specific scripts:

```r
# Load configuration first
source("config.R")

# Data quality checks
source("data_checks.R")

# Data cleaning only
source("0_cleaning_data_total.R")

# Descriptive statistics only
source("0_create_descriptive_statistics.R")

# GSC model only
source("1_model_run_GSC_total.R")

# MTGP model only
source("2_model_run_MTGP_total.R")

# Results generation only
source("3_export_results.R")
```

### Sensitivity Analysis

Run the sensitivity analysis separately:

```r
source("run_sensitivity_analysis.R")
```

The sensitivity analysis includes:
- Extended study period (1983-2003 vs 1988-2003)
- Earlier treatment year (1990 vs 1995)
- Matrix Completion (MC) models as alternative to MTGP
- Additional robustness checks

## Output Structure

After running the analysis, you'll find:

```
results/
├── tables/
│   ├── case_counts.tex                # Data quality validation
│   ├── case_counts_grouped.tex        # Grouped case counts
│   ├── table1.tex                     # Descriptive statistics
│   ├── mtgp_table.tex                # Main treatment effects
│   ├── mtgp_table_ratio.tex          # Treatment effect ratios
│   └── mtgp_table_sum_TE.tex         # Sum of treatment effects
├── figures/
│   ├── map_counties.png              # County map
│   ├── ATT_plot_overall_trt1995.png  # Overall treatment effects
│   └── ATT_plot_[subset]_trt1995.png # Subgroup treatment effects
└── [date]/                           # Model objects by date
    ├── total_Lymphoma_*.RData        # GSC model fits
    ├── mtgp_fit_*.RData              # MTGP model fits
    └── ...
```

## Contact

For questions about this analysis, please contact the author, Sofia Vega, at sofialeighvega@gmail.com.
