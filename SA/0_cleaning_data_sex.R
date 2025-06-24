## This script cleans SEER lymphoma incidence data 
source("config.R")



# Define years of interest
years <- c(1983:2003)

# Define treated year
treated_year <- 1990

# Load packages
library(readr)
library(tidyverse)

# Load cleaned data
SEER <- read.csv(file.path(DATA_DIR, "sex_1983_2003.csv"))

# Load population data
load(file.path(DATA_DIR, "seer_population.RData"))

# Recode data
SEER_1 <- SEER %>%
  mutate(cancer_type = recode(`Lymphoma`, "0" = "Lymphoma"),
         
         state_county_recode = recode(`State.county`,
                                      "0"="San Francisco-Oakland SMSA Registry",
                                      "1"="  CA: Alameda County (06001)",
                                      "2"="  CA: Contra Costa County (06013)",
                                      "3"="  CA: Marin County (06041)",
                                      "4"="  CA: San Francisco County (06075)",
                                      "5"="  CA: San Mateo County (06081)",
                                      "6"="Connecticut",
                                      "7"="  CT: Fairfield County (09001)",
                                      "8"="  CT: Hartford County (09003)",
                                      "9"="  CT: Litchfield County (09005)",
                                      "10"="  CT: Middlesex County (09007)",
                                      "11"="  CT: New Haven County (09009)",
                                      "12"="  CT: New London County (09011)",
                                      "13"="  CT: Tolland County (09013)",
                                      "14"="  CT: Windham County (09015)",
                                      "15"="  CT: Unknown (09999)",
                                      "16"="Atlanta (Metropolitan) Registry",
                                      "17"="  GA: Clayton County (13063)",
                                      "18"="  GA: Cobb County (13067)",
                                      "19"="  GA: DeKalb County (13089)",
                                      "20"="  GA: Fulton County (13121)",
                                      "21"="  GA: Gwinnett County (13135)",
                                      "22"="Hawaii",
                                      "23"="  HI: Hawaii County (15001) - 2000+",
                                      "24"="  HI: Honolulu County (15003) - 2000+",
                                      "25"="  HI: Kalawao County (15005) - 2000+",
                                      "26"="  HI: Kauai County (15007) - 2000+",
                                      "27"="  HI: Maui County (15009) - 2000+",
                                      "28"="  HI: Hawaii County (15911) - Cases Only - pre-2000",
                                      "29"="  HI: Honolulu County (15912) - Cases Only - pre-2000",
                                      "30"="  HI: Kalawao County (15913) - Cases Only - pre-2000",
                                      "31"="  HI: Kauai County (15914) - Cases Only - pre-2000",
                                      "32"="  HI: Maui County (15915) - Cases Only - pre-2000",
                                      "33"="  HI: Unknown (15999)",
                                      "34"="Iowa",
                                      "35"="  IA: Adair County (19001)",
                                      "36"="  IA: Adams County (19003)",
                                      "37"="  IA: Allamakee County (19005)",
                                      "38"="  IA: Appanoose County (19007)",
                                      "39"="  IA: Audubon County (19009)",
                                      "40"="  IA: Benton County (19011)",
                                      "41"="  IA: Black Hawk County (19013)",
                                      "42"="  IA: Boone County (19015)",
                                      "43"="  IA: Bremer County (19017)",
                                      "44"="  IA: Buchanan County (19019)",
                                      "45"="  IA: Buena Vista County (19021)",
                                      "46"="  IA: Butler County (19023)",
                                      "47"="  IA: Calhoun County (19025)",
                                      "48"="  IA: Carroll County (19027)",
                                      "49"="  IA: Cass County (19029)",
                                      "50"="  IA: Cedar County (19031)",
                                      "51"="  IA: Cerro Gordo County (19033)",
                                      "52"="  IA: Cherokee County (19035)",
                                      "53"="  IA: Chickasaw County (19037)",
                                      "54"="  IA: Clarke County (19039)",
                                      "55"="  IA: Clay County (19041)",
                                      "56"="  IA: Clayton County (19043)",
                                      "57"="  IA: Clinton County (19045)",
                                      "58"="  IA: Crawford County (19047)",
                                      "59"="  IA: Dallas County (19049)",
                                      "60"="  IA: Davis County (19051)",
                                      "61"="  IA: Decatur County (19053)",
                                      "62"="  IA: Delaware County (19055)",
                                      "63"="  IA: Des Moines County (19057)",
                                      "64"="  IA: Dickinson County (19059)",
                                      "65"="  IA: Dubuque County (19061)",
                                      "66"="  IA: Emmet County (19063)",
                                      "67"="  IA: Fayette County (19065)",
                                      "68"="  IA: Floyd County (19067)",
                                      "69"="  IA: Franklin County (19069)",
                                      "70"="  IA: Fremont County (19071)",
                                      "71"="  IA: Greene County (19073)",
                                      "72"="  IA: Grundy County (19075)",
                                      "73"="  IA: Guthrie County (19077)",
                                      "74"="  IA: Hamilton County (19079)",
                                      "75"="  IA: Hancock County (19081)",
                                      "76"="  IA: Hardin County (19083)",
                                      "77"="  IA: Harrison County (19085)",
                                      "78"="  IA: Henry County (19087)",
                                      "79"="  IA: Howard County (19089)",
                                      "80"="  IA: Humboldt County (19091)",
                                      "81"="  IA: Ida County (19093)",
                                      "82"="  IA: Iowa County (19095)",
                                      "83"="  IA: Jackson County (19097)",
                                      "84"="  IA: Jasper County (19099)",
                                      "85"="  IA: Jefferson County (19101)",
                                      "86"="  IA: Johnson County (19103)",
                                      "87"="  IA: Jones County (19105)",
                                      "88"="  IA: Keokuk County (19107)",
                                      "89"="  IA: Kossuth County (19109)",
                                      "90"="  IA: Lee County (19111)",
                                      "91"="  IA: Linn County (19113)",
                                      "92"="  IA: Louisa County (19115)",
                                      "93"="  IA: Lucas County (19117)",
                                      "94"="  IA: Lyon County (19119)",
                                      "95"="  IA: Madison County (19121)",
                                      "96"="  IA: Mahaska County (19123)",
                                      "97"="  IA: Marion County (19125)",
                                      "98"="  IA: Marshall County (19127)",
                                      "99"="  IA: Mills County (19129)",
                                      "100"="  IA: Mitchell County (19131)",
                                      "101"="  IA: Monona County (19133)",
                                      "102"="  IA: Monroe County (19135)",
                                      "103"="  IA: Montgomery County (19137)",
                                      "104"="  IA: Muscatine County (19139)",
                                      "105"="  IA: OBrien County (19141)",
                                      "106"="  IA: Osceola County (19143)",
                                      "107"="  IA: Page County (19145)",
                                      "108"="  IA: Palo Alto County (19147)",
                                      "109"="  IA: Plymouth County (19149)",
                                      "110"="  IA: Pocahontas County (19151)",
                                      "111"="  IA: Polk County (19153)",
                                      "112"="  IA: Pottawattamie County (19155)",
                                      "113"="  IA: Poweshiek County (19157)",
                                      "114"="  IA: Ringgold County (19159)",
                                      "115"="  IA: Sac County (19161)",
                                      "116"="  IA: Scott County (19163)",
                                      "117"="  IA: Shelby County (19165)",
                                      "118"="  IA: Sioux County (19167)",
                                      "119"="  IA: Story County (19169)",
                                      "120"="  IA: Tama County (19171)",
                                      "121"="  IA: Taylor County (19173)",
                                      "122"="  IA: Union County (19175)",
                                      "123"="  IA: Van Buren County (19177)",
                                      "124"="  IA: Wapello County (19179)",
                                      "125"="  IA: Warren County (19181)",
                                      "126"="  IA: Washington County (19183)",
                                      "127"="  IA: Wayne County (19185)",
                                      "128"="  IA: Webster County (19187)",
                                      "129"="  IA: Winnebago County (19189)",
                                      "130"="  IA: Winneshiek County (19191)",
                                      "131"="  IA: Woodbury County (19193)",
                                      "132"="  IA: Worth County (19195)",
                                      "133"="  IA: Wright County (19197)",
                                      "134"="New Mexico",
                                      "135"="  NM: Bernalillo County (35001)",
                                      "136"="  NM: Catron County (35003)",
                                      "137"="  NM: Chaves County (35005)",
                                      "138"="  NM: Colfax County (35007)",
                                      "139"="  NM: Curry County (35009)",
                                      "140"="  NM: De Baca County (35011)",
                                      "141"="  NM: Dona Ana County (35013)",
                                      "142"="  NM: Eddy County (35015)",
                                      "143"="  NM: Grant County (35017)",
                                      "144"="  NM: Guadalupe County (35019)",
                                      "145"="  NM: Harding County (35021)",
                                      "146"="  NM: Hidalgo County (35023)",
                                      "147"="  NM: Lea County (35025)",
                                      "148"="  NM: Lincoln County (35027)",
                                      "149"="  NM: Los Alamos County (35028)",
                                      "150"="  NM: Luna County (35029)",
                                      "151"="  NM: McKinley County (35031)",
                                      "152"="  NM: Mora County (35033)",
                                      "153"="  NM: Otero County (35035)",
                                      "154"="  NM: Quay County (35037)",
                                      "155"="  NM: Rio Arriba County (35039)",
                                      "156"="  NM: Roosevelt County (35041)",
                                      "157"="  NM: Sandoval County (35043)",
                                      "158"="  NM: San Juan County (35045)",
                                      "159"="  NM: San Miguel County (35047)",
                                      "160"="  NM: Santa Fe County (35049)",
                                      "161"="  NM: Sierra County (35051)",
                                      "162"="  NM: Socorro County (35053)",
                                      "163"="  NM: Taos County (35055)",
                                      "164"="  NM: Torrance County (35057)",
                                      "165"="  NM: Union County (35059)",
                                      "166"="  NM: Cibola/Valencia",
                                      "167"="    NM: Cibola County (35006) - 1982+",
                                      "168"="    NM: Valencia County (35061) - 1982+",
                                      "169"="    NM: Cibola/Valencia (35910) - pre-1982",
                                      "170"="  NM: Unknown (35999)",
                                      "171"="Utah",
                                      "172"="  UT: Beaver County (49001)",
                                      "173"="  UT: Box Elder County (49003)",
                                      "174"="  UT: Cache County (49005)",
                                      "175"="  UT: Carbon County (49007)",
                                      "176"="  UT: Daggett County (49009)",
                                      "177"="  UT: Davis County (49011)",
                                      "178"="  UT: Duchesne County (49013)",
                                      "179"="  UT: Emery County (49015)",
                                      "180"="  UT: Garfield County (49017)",
                                      "181"="  UT: Grand County (49019)",
                                      "182"="  UT: Iron County (49021)",
                                      "183"="  UT: Juab County (49023)",
                                      "184"="  UT: Kane County (49025)",
                                      "185"="  UT: Millard County (49027)",
                                      "186"="  UT: Morgan County (49029)",
                                      "187"="  UT: Piute County (49031)",
                                      "188"="  UT: Rich County (49033)",
                                      "189"="  UT: Salt Lake County (49035)",
                                      "190"="  UT: San Juan County (49037)",
                                      "191"="  UT: Sanpete County (49039)",
                                      "192"="  UT: Sevier County (49041)",
                                      "193"="  UT: Summit County (49043)",
                                      "194"="  UT: Tooele County (49045)",
                                      "195"="  UT: Uintah County (49047)",
                                      "196"="  UT: Utah County (49049)",
                                      "197"="  UT: Wasatch County (49051)",
                                      "198"="  UT: Washington County (49053)",
                                      "199"="  UT: Wayne County (49055)",
                                      "200"="  UT: Weber County (49057)",
                                      "201"="  UT: Unknown (49999)",
                                      "202"="Seattle (Puget Sound) Registry",
                                      "203"="  WA: Clallam County (53009)",
                                      "204"="  WA: Grays Harbor County (53027)",
                                      "205"="  WA: Island County (53029)",
                                      "206"="  WA: Jefferson County (53031)",
                                      "207"="  WA: King County (53033)",
                                      "208"="  WA: Kitsap County (53035)",
                                      "209"="  WA: Mason County (53045)",
                                      "210"="  WA: Pierce County (53053)",
                                      "211"="  WA: San Juan County (53055)",
                                      "212"="  WA: Skagit County (53057)",
                                      "213"="  WA: Snohomish County (53061)",
                                      "214"="  WA: Thurston County (53067)",
                                      "215"="  WA: Whatcom County (53073)"),
         fips = str_extract(state_county_recode, "([:digit:]{5})"),
         year_recode = `Year.of.diagnosis.subset` + 1983
         
  )


# Subset to be ages less than 29
SEER_final <- SEER_1 %>% select(fips, year_recode, `Age.subset`, Count, `Sex`) %>%
  filter(`Age.subset` <= 6,
         !is.na(fips),
         !is.na(year_recode)) %>%
  mutate(age_recode = recode(`Age.subset`,
                             "0"="00 years",
                             "1"="01-04 years",
                             "2"="05-09 years",
                             "3"="10-14 years",
                             "4"="15-19 years",
                             "5"="20-24 years",
                             "6"="25-29 years")
  ) %>%
  select(fips, year_recode, age_recode, Count, `Sex`) %>%
  group_by(fips, year_recode, `Sex`) %>%
  summarise(count = sum(Count)) %>%
  rename(year = year_recode) %>%
  ungroup()

# Rename columns
colnames(SEER_final) <- c("FIPS","YEAR_DX","SEX","CL_CASES") 

# Subset SEER to start at start year and end at end year
data_full <- SEER_final
data_full <- data_full %>% filter(!YEAR_DX < years[1] & !YEAR_DX > years[length(years)])

# For Hawaii counties, combine appropriate FIPS
HI_pre2000 <- c(15911, 15912, 15913, 15914, 15915)
HI_post2000 <- c(15001, 15003, 15005, 15007, 15009)

for (i in 1:length(HI_pre2000)){
  # Find the indices prior to 2000 in post HI FIPS to be replaced
  ind_rep_post2000 <- which((data_full$FIPS == HI_post2000[i]) & (data_full$YEAR_DX < 2000))
  # Find the indices prior 2000 in the pre HI FIPS to replace 
  ind_rep_pre2000 <- which((data_full$FIPS == HI_pre2000[i]) & (data_full$YEAR_DX < 2000))
  
  print(data_full[ind_rep_pre2000,2:3])
  print(data_full[ind_rep_post2000,2:3])
  
  # Replace every column but FIPS column (this keeps current FIPS)
  data_full[ind_rep_post2000,2:3] <- data_full[ind_rep_pre2000,2:3]
  print(data_full[ind_rep_post2000,2:3])
}

# Remove old HI counties 
data_full <- data_full %>% filter(!(FIPS %in% HI_pre2000))

# Remove 35910 (SEER combined 35006 and 35061 for 1969-1981 so does not apply for our study period)
# https://seer.cancer.gov/seerstat/variables/countyattribs/ruralurban.html
data_full <- data_full %>% filter((FIPS !=35910 ))

## Set Up Data

# Specify treated CT and CA Counties
data_full <- data_full %>% mutate(
  # Specify CT as treated after 2003
  C = case_when(
    str_starts(FIPS, "09") & YEAR_DX >= treated_year ~ 1,
    TRUE ~ 0),
  # Specify CA as treated after 2003
  C = case_when(
    str_starts(FIPS, "06") & YEAR_DX >= treated_year ~ 1, 
    TRUE ~ C)
)

# Combine with SEER subgroup data

# Subset to years of interest
seer_population_0_29 <- seer_population %>% filter(Year %in% years)

# Remove variables not considered
seer_population_0_29$Registry <- NULL
seer_population_0_29$Race <- NULL
seer_population_0_29$Origin <- NULL
seer_population_0_29$State <- NULL
seer_population_0_29$Age <- NULL

# Aggregate the data, summing the Population for each unique combination
seer_data_aggregated <- seer_population_0_29 %>%
  group_by(Year, FIPS, Sex) %>%
  summarise(Population = sum(Population)) %>%
  ungroup()

## dont have census data from 1983 to 1989 so copy 1990
data_full_hisp <- left_join(data_full,seer_data_aggregated, by = c("FIPS" = "FIPS", "YEAR_DX" = "Year", "SEX" = "Sex"))

# Remove Sex ==0
data_full_hisp <- data_full_hisp %>% filter(SEX != 0)

# Fill years prior to 1990
data_full_hisp <- data_full_hisp %>%
  group_by(FIPS, SEX) %>%
  arrange(FIPS, SEX, YEAR_DX) %>%
  fill(Population, .direction = "up") %>%
  ungroup()

# Assign column names
colnames(data_full_hisp) <- c("FIPS","YEAR_DX","SEX","CL_CASES","C",
                              "POP")
# FIPS 15005 was incorporated into 15009 (HI) I'm going to remove
data_full_hisp <- data_full_hisp %>% filter(!(FIPS==15005))

# Remove unknown counties from the dataset
data_full_hisp<- data_full_hisp %>% filter(! (FIPS %in% c("09999", "15999", "35999", "49999")) )

# Identify if a row (county) has all 0s
#data_full_hisp %>%
#  group_by(FIPS) %>%
#  summarise(total_cases = sum(CL_CASES)) %>%
#  filter(total_cases == 0) %>%
#  pull(FIPS) -> zero_case_FIPS

# Remove counties with all 0s
#data_full_hisp_noZero <- data_full_hisp %>%
#  filter(! FIPS %in% zero_case_FIPS)

data_full_hisp_noZero <- data_full_hisp
# Save final data
save(data_full_hisp_noZero, file = file.path(DATA_DIR, paste0('sex_Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData')))

# Print summary statistics
cat("Sensitivity analysis sex data cleaning completed!\n")
cat("Final dataset dimensions:", dim(data_full_hisp_noZero), "\n")
cat("Sex groups:", unique(data_full_hisp_noZero$SEX), "\n")
cat("Years:", range(data_full_hisp_noZero$YEAR_DX), "\n")
cat("Total counties:", length(unique(data_full_hisp_noZero$FIPS)), "\n")
cat("Total cases:", sum(data_full_hisp_noZero$CL_CASES), "\n")
cat("Treated observations:", sum(data_full_hisp_noZero$C), "\n")



