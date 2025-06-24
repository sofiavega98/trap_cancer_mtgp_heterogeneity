# Script to get the population data from SEER

# Define the file path (replace with your actual file path)
file_paths <- c("/Users/sofiavega/Documents/Harvard/Lymphoma_Subsets/data/Population/ca.1990_2022.19ages.txt",
                "/Users/sofiavega/Documents/Harvard/Lymphoma_Subsets/data/Population/ct.1990_2022.19ages.txt",
                "/Users/sofiavega/Documents/Harvard/Lymphoma_Subsets/data/Population/ga.1990_2022.19ages.txt",
                "/Users/sofiavega/Documents/Harvard/Lymphoma_Subsets/data/Population/hi.1990_2022.19ages.txt",
                "/Users/sofiavega/Documents/Harvard/Lymphoma_Subsets/data/Population/ia.1990_2022.19ages.txt",
                "/Users/sofiavega/Documents/Harvard/Lymphoma_Subsets/data/Population/mi.1990_2022.19ages.txt",
                "/Users/sofiavega/Documents/Harvard/Lymphoma_Subsets/data/Population/nm.1990_2022.19ages.txt",
                "/Users/sofiavega/Documents/Harvard/Lymphoma_Subsets/data/Population/ut.1990_2022.19ages.txt",
                "/Users/sofiavega/Documents/Harvard/Lymphoma_Subsets/data/Population/wa.1990_2022.19ages.txt")
seer_data <- c()
for(file_path in file_paths){
  # Define the column widths based on the data dictionary
  widths <- c(4, 2, 2, 3, 2, 1, 1, 1, 2, 8)
  
  # Define column names
  col_names <- c("Year", "State", "StateFIPS", "CountyFIPS", "Registry", "Race", "Origin", "Sex", "Age", "Population")
  
  # Read the fixed-width file
  seer_data_temp <- read.fwf(file_path, 
                        widths = widths, 
                        col.names = col_names, 
                        colClasses = c("integer", "character", "character", "character", 
                                       "integer", "integer", "integer", "integer", 
                                       "integer", "integer"),
                        strip.white = TRUE)
  
  # Convert StateFIPS and CountyFIPS to character to preserve leading zeros
  seer_data_temp$StateFIPS <- sprintf("%02d", as.integer(seer_data_temp$StateFIPS))
  seer_data_temp$CountyFIPS <- sprintf("%03d", as.integer(seer_data_temp$CountyFIPS))
  
  # Combine states
  seer_data <- rbind(seer_data, seer_data_temp)

}

# Print the first few rows to verify
print(head(seer_data))

# Verify 9 unique states
length(unique(seer_data$State))

# Combine StateFIPS and CountyFIPS into a single FIPS column
seer_data$FIPS <- paste0(seer_data$StateFIPS, seer_data$CountyFIPS)

# Remove the separate StateFIPS and CountyFIPS columns 
seer_data$StateFIPS <- NULL
seer_data$CountyFIPS <- NULL

# Combine Race groups 3 and 4 into group 3
seer_data$Race[seer_data$Race == 4] <- 3

# Filter to age of interests 0-29 (groups 0-6)
seer_data <- seer_data %>% filter(Age %in% c(0:6))
#seer_data$Age <- NULL
source("config.R")


# Aggregate the data, summing the Population for each unique combination
seer_data_aggregated <- seer_data %>%
  group_by(Year, State, FIPS, Registry, Race, Origin, Sex, Age) %>%
  summarise(Population = sum(Population)) %>%
  ungroup()

# Optional: Check for any remaining duplicates
duplicates <- seer_data_aggregated %>%
  group_by(Year, State, FIPS, Registry, Race, Origin, Sex, Age) %>%
  filter(n() > 1)

if (nrow(duplicates) > 0) {
  print("Warning: There are still some duplicate rows:")
  print(duplicates)
} else {
  print("No duplicates found after aggregation.")
}

# Reorder the columns to put FIPS after State
seer_population <- seer_data_aggregated[, c("Year", "State", "FIPS", "Registry", "Race", "Origin", "Sex", "Age","Population")]

# Save data
save(seer_population , file = file.path(DATA_DIR, "seer_population.RData"))

