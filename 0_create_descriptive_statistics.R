################################################################################
# Descriptive Statistics and Visualizations
# 
# This script creates descriptive statistics and visualizations for the paper,
# including Table 1 (descriptive statistics) and Figure S2 (county map).
#
# INPUTS:
# - Cleaned datasets for all demographic subsets from data cleaning scripts
#
# OUTPUTS:
# - Table 1: "results/tables/table1.tex"
# - Map Figure: "results/figures/map_counties.png"
#
# This script generates the descriptive statistics presented in the paper,
# showing aggregate lymphoma rates by demographic group, treatment status,
# and time period.
################################################################################

# R Script to create descriptive statistics and map
# Date Created: 08/26/24
# Updated: 2024 - Added comprehensive documentation

# Load configuration
source("config.R")

# Load required libraries
library(tidyverse)
library(ggplot2)
library(kableExtra)
library(usmap)

# Define study parameters
years <- YEARS
treated_year <- TREATED_YEAR

cat("Creating descriptive statistics and visualizations...\n")
cat("Study period:", min(years), "-", max(years), "\n")
cat("Treatment year:", treated_year, "\n")

# Load all cleaned datasets
load(file.path(DATA_DIR, "NH_Lymphoma_1988_2003_0-29.RData"))
data_full_pop_NH <- data_full_pop

load(file.path(DATA_DIR, "H_Lymphoma_1988_2003_0-29.RData"))
data_full_pop_H <- data_full_pop

load(file.path(DATA_DIR, "male_Lymphoma_1988_2003_0-29.RData"))
data_full_pop_male <- data_full_pop

load(file.path(DATA_DIR, "female_Lymphoma_1988_2003_0-29.RData"))
data_full_pop_female <- data_full_pop

load(file.path(DATA_DIR, "white_Lymphoma_1988_2003_0-29.RData"))
data_full_pop_white <- data_full_pop

load(file.path(DATA_DIR, "black_Lymphoma_1988_2003_0-29.RData"))
data_full_pop_black <- data_full_pop

load(file.path(DATA_DIR, "other_Lymphoma_1988_2003_0-29.RData"))
data_full_pop_other <- data_full_pop

load(file.path(DATA_DIR, "nonHisp_Lymphoma_1988_2003_0-29.RData"))
data_full_pop_nonHisp <- data_full_pop

load(file.path(DATA_DIR, "Hisp_Lymphoma_1988_2003_0-29.RData"))
data_full_pop_Hisp <- data_full_pop

load(file.path(DATA_DIR, "young_Lymphoma_1988_2003_0-29.RData"))
data_full_pop_young <- data_full_pop

load(file.path(DATA_DIR, "old_Lymphoma_1988_2003_0-29.RData"))
data_full_pop_old <- data_full_pop

cat("✓ All demographic datasets loaded\n")

#######################################################
## Table 1: Aggregate Rates by Demographic Group ##
#######################################################
# This table displays aggregate lymphoma rates pre- and post-treatment
# for each demographic subgroup, separated by treatment status

# Function to calculate aggregate rates based on the sum of rates divided by the number of unique FIPS
aggregate_rate <- function(df){
  # Identify treated counties
  treated_vec <- df %>% filter(C==1) %>% pull(FIPS) %>% unique()
  
  # Add a period column to distinguish before and after 1995 and treated column
  df <- df %>%
    mutate(period = ifelse(YEAR_DX < 1995, "before", "after"),
           treated = ifelse(FIPS %in% treated_vec, "treated", "control"))
  
  # Calculate aggregate rates by grouping
  aggregate_rates <- df %>%
    group_by(treated, period) %>%
    summarise(
      aggregate_cases = sum(CL_CASES),
      aggregate_pop = sum(POP),
      .groups = 'drop'
    ) %>%
    mutate(aggregate_rate = aggregate_cases / aggregate_pop * 100000)
  
  return(aggregate_rates[,c(1,2,5)])
}

# Apply function to each demographic subset
agg_rate_NH <- aggregate_rate(data_full_pop_NH)
agg_rate_H <- aggregate_rate(data_full_pop_H)
agg_rate_male <- aggregate_rate(data_full_pop_male)
agg_rate_female <- aggregate_rate(data_full_pop_female)
agg_rate_white <- aggregate_rate(data_full_pop_white)
agg_rate_black <- aggregate_rate(data_full_pop_black)
agg_rate_other <- aggregate_rate(data_full_pop_other)
agg_rate_nonHisp <- aggregate_rate(data_full_pop_nonHisp)
agg_rate_Hisp <- aggregate_rate(data_full_pop_Hisp)
agg_rate_young <- aggregate_rate(data_full_pop_young)
agg_rate_old <- aggregate_rate(data_full_pop_old)

# Combine all results into a single data frame
agg_rates_list <- list(
  "Non-Hodgkins" = agg_rate_NH,
  "Hodgkins" = agg_rate_H,
  "Male" = agg_rate_male,
  "Female" = agg_rate_female,
  "White" = agg_rate_white,
  "Black" = agg_rate_black,
  "Other" = agg_rate_other,
  "Non-Hispanic" = agg_rate_nonHisp,
  "Hispanic" = agg_rate_Hisp,
  "Ages 0-19" = agg_rate_young,
  "Ages 20-29" = agg_rate_old
)

# Create a summary table with rates rounded to 2 decimal places
summary_table <- bind_rows(agg_rates_list, .id = "Group") %>%
  pivot_wider(names_from = c(treated, period), values_from = aggregate_rate, names_prefix = "") %>%
  rename("Control Before" = control_before,
         "Control After" = control_after,
         "Treated Before" = treated_before,
         "Treated After" = treated_after) %>%
  select(Group, `Control Before`, `Control After`, `Treated Before`, `Treated After`)%>%
  mutate(across(`Control Before`:`Treated After`, ~ round(.x, 2)))

print(summary_table)

# Create and save LaTeX table
table1 <- summary_table %>%
  kbl(format = "latex", booktabs = TRUE) %>%
  kable_classic(full_width = FALSE, font_size = 10, latex_options = "HOLD_position")

save_kable(table1, file = file.path(RESULTS_DIR, "tables", "table1.tex"))

cat("✓ Table 1 (Aggregate Rates) created and saved\n")

#####################################
## Additional Descriptive Tables ##
#####################################
# These tables provide additional context about case distributions

# Combine datasets for additional analysis
data_full_pop_NH$strata <- "Type"
data_full_pop_H$strata <- "Type"
data_full_pop_male$strata <- "Sex"
data_full_pop_female$strata <- "Sex"
data_full_pop_white$strata <- "Race"
data_full_pop_black$strata <- "Race"
data_full_pop_other$strata <- "Race"
data_full_pop_nonHisp$strata <- "Hispanic"
data_full_pop_Hisp$strata <- "Hispanic"
data_full_pop_young$strata <- "Age"
data_full_pop_old$strata <- "Age"
data_full_pop_NH$group <- "Non-Hodgkins"
data_full_pop_H$group <- "Hodgkins"
data_full_pop_male$group <- "Male"
data_full_pop_female$group <- "Female"
data_full_pop_white$group <- "White"
data_full_pop_black$group <- "Black"
data_full_pop_other$group <- "Other"
data_full_pop_nonHisp$group <- "Non-Hispanic"
data_full_pop_Hisp$group <- "Hispanic"
data_full_pop_young$group <- "Ages 0-19"
data_full_pop_old$group <- "Ages 20-29"

data_full <- rbind(data_full_pop_NH, data_full_pop_H, data_full_pop_male, data_full_pop_female,
                   data_full_pop_white, data_full_pop_black, data_full_pop_other,
                   data_full_pop_nonHisp, data_full_pop_Hisp, data_full_pop_young, data_full_pop_old)

# Calculate total cases for each strata
total_cases_per_strata <- data_full %>%
  group_by(strata) %>%
  summarise(total_cases = sum(CL_CASES), .groups = 'drop')

# Calculate the proportion of cases for each group within the strata
proportions <- data_full %>%
  group_by(strata, group) %>%
  summarise(group_cases = sum(CL_CASES), .groups = 'drop') %>%
  left_join(total_cases_per_strata, by = "strata") %>%
  mutate(proportion = (group_cases / total_cases) * 100) %>%
  mutate(proportion = sprintf("%.2f", proportion)) %>%
  select(strata, group, proportion)

# Define the order for the groups
group_order <- c("Non-Hodgkins", "Hodgkins", "Male", "Female",
                 "White", "Black", "Other", "Non-Hispanic", "Hispanic",
                 "Ages 0-19", "Ages 20-29")

# Reorder the table according to the specified group order
table2 <- proportions %>%
  mutate(group = factor(group, levels = group_order)) %>%
  arrange(group)

print(table2)

# Create comprehensive table with proportions by treatment period
# Identify treated counties
treated_vec <- data_full %>% filter(C==1) %>% pull(FIPS) %>% unique()

# Add a period column to distinguish before and after 1995
data_full <- data_full %>%
  mutate(period = ifelse(YEAR_DX < 1995, "before", "after"),
         treated = ifelse(FIPS %in% treated_vec, "treated", "control"))

# Calculate total cases for each strata and treatment period
total_cases_per_strata_period <- data_full %>%
  group_by(strata, treated, period) %>%
  summarise(total_cases = sum(CL_CASES), .groups = 'drop')

# Calculate the total cases for each group within each strata and treatment period
total_cases_per_group <- data_full %>%
  group_by(strata, group, treated, period) %>%
  summarise(group_cases = sum(CL_CASES), .groups = 'drop')

# Calculate the proportion of cases for each group within each strata and treatment period
proportions <- total_cases_per_group %>%
  left_join(total_cases_per_strata_period, by = c("strata", "treated", "period")) %>%
  mutate(proportion = (group_cases / total_cases) * 100) %>%
  mutate(proportion = sprintf("%.2f", proportion))

# Pivot the table to get the desired format
final_table <- proportions %>% select(!c("group_cases","total_cases")) %>%
  pivot_wider(names_from = c(treated, period), values_from = proportion,names_prefix = "") %>%
  arrange(factor(group, levels = group_order)) %>%
  rename("Strata" = strata,
         "Group" = group,
         "Control Before" = control_before,
         "Control After" = control_after,
         "Treated Before" = treated_before,
         "Treated After" = treated_after)%>%
  select(Strata, Group, `Control Before`, `Control After`, `Treated Before`, `Treated After`)

# Add overall column from table 2
table3 <- cbind(final_table, Overall = table2$proportion)

print(table3)

cat("✓ Additional descriptive tables created\n")

#########################################
## Figure S2: Map of Counties in Dataset ##
#########################################
# This map shows which counties are treated vs control

# Install and load required packages for mapping
if (!require(usmap)) install.packages("usmap")
if (!require(ggplot2)) install.packages("ggplot2")
library(usmap)
library(ggplot2)

# All feature datasets should have the same counties, so use NH as the data for the map
data_full <- data_full_pop_NH

# Identify treated counties
treated_FIPS <- data_full %>% 
  filter(C == 1) %>% 
  pull(FIPS) %>% 
  unique()

untreated_FIPS <- data_full %>% 
  filter(!(FIPS %in% treated_FIPS)) %>% 
  pull(FIPS) %>% 
  unique()

# Create map_data dataframe
map_data <- data.frame(
  fips = c(treated_FIPS, untreated_FIPS),
  Treatment = c(rep(1, length(treated_FIPS)), rep(0, length(untreated_FIPS)))
)

# Convert FIPS codes to character and ensure they are 5 digits long
map_data$fips <- sprintf("%05d", as.numeric(map_data$fips))

# Convert Treatment to a factor
map_data$Treatment <- as.factor(map_data$Treatment)

# Create a map showing treated vs control counties
map <- plot_usmap(
  data = map_data, 
  values = "Treatment", 
  regions = "counties",
  color = "gray",  # Set county border color
  linewidth = .05
  
) +
  scale_fill_manual(
    values = c("0" = "lightblue", "1" = "darkred"),
    name = "Treatment Status",
    labels = c("Control", "Treated"),
    na.value = "white"  # Set color for counties not in the dataset
  ) +
  theme(
    legend.position = "right",
    panel.background = element_rect(color = "white", fill = "white"),  # White background
    plot.background = element_rect(color = "white", fill = "white")  # White plot area
  ) +
  geom_polygon(data = usmapdata::us_map(regions = "states"),
               aes(x, y, group = group), fill = NA, linewidth = .08, color = "gray")

# Display the map
print(map)

# Save map
map %>%
  ggsave(filename = file.path(RESULTS_DIR, "figures", "map_counties.png"),
         height = 5,
         width = 6,
         units = "in")

cat("✓ Figure S2 (County Map) created and saved\n")

################################################################################
# COMPLETION
################################################################################

cat("\n=== DESCRIPTIVE STATISTICS COMPLETE ===\n")
cat("Generated outputs:\n")
cat("- Table 1: Aggregate rates by demographic group and treatment status\n")
cat("- Figure S2: Map of treated vs control counties\n")
cat("- Additional descriptive tables for reference\n\n")

cat("Files saved:\n")
cat("- results/tables/table1.tex\n")
cat("- results/figures/map_counties.png\n\n")

