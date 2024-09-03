## NOTE: Data operations on specific countries' data have been removed from this script for confidentiality

library(tidyverse)
library(magrittr)
library(lubridate)
library(here)
library(readxl)

# Load grey lit data file
grey_new_raw <- haven::read_dta(here::here("data/update-2024/round1_round2_appended_12022024.dta")) %>% 
  mutate(prefix = str_extract(idno, "^[A-Z]+")) |>
  filter(prefix == "GRY") |>
  select(-prefix)
  
stopifnot(all(grey_new_raw$dataset == "new"))

grey_old_raw <- haven::read_dta(here("data/current/raw/grey/grey-combined-2020-04-18.dta"))

# Check for differences in column names
colnames(grey_old_raw) |>
  (\(.x).x[!.x %in% colnames(grey_new_raw)])()

colnames(grey_new_raw) |>
  (\(.x).x[!.x %in% colnames(grey_old_raw)])()

# Make columns align better where possible
grey_old <- grey_old_raw |>
  rename(
    icd_guess1 = icdguess1,
    icd_guess2 = icdguess2,
    icd_guess3 = icdguess3,
    icd_guess4 = icdguess4,
    note = notes
  ) |>
  mutate(dataset = "old") |>
  select(
    # neither of these columns are actually used in the grey data
    -denom_other,
    -biblio
  ) 

# The new data and old data are not coded in the same way for some variables
# Make tables to map the new values to old codings

# agestrat (manual -- labels are different text)
# check with attributes(grey_new_raw$agestrat); attributes(grey_old_raw$agestrat)
agestrat_recode <- tibble(
  new_code = c(3, 4),
  old_code = c(1, 0)
)

# area (automatic -- text labels are the same, but the numerical codings are different)
new_area_coding <- grey_new_raw$area |> attr("labels")

area_recode <- tibble(
  new_code = new_area_coding,
  old_code = (grey_old_raw$area |> attr("labels"))[names(new_area_coding)]
)
  
# timeframe (automatic-ish -- small fix to make the text labels align )
new_timeframe_coding <- grey_new_raw$timeframe |> 
  attr("labels")

names(new_timeframe_coding) <- names(new_timeframe_coding) |>
  str_replace("Late only", "Late only (43-1y)")

timeframe_recode <- tibble(
  new_code = new_timeframe_coding,
  old_code = (grey_old_raw$timeframe |> attr("labels"))[names(new_timeframe_coding)]
)  

# va (automatic -- text labels are the same, but numerical codings are different)
new_va_coding <- grey_new_raw$va |> attr("labels")

va_recode <- tibble(
  new_code = new_va_coding,
  old_code = (grey_old_raw$va |> attr("labels"))[names(new_va_coding)]
)

# denom (automatic -- text labels are the same, but numerical codings are different)
# some of the new codes don't have equivalents in the old coding system, but doesn't matter for grey lit
new_denom_coding <- grey_new_raw$denom |> attr("labels")

denom_recode <- tibble(
  new_code = new_denom_coding,
  old_code = (grey_old_raw$denom |> attr("labels"))[names(new_denom_coding)]
)

## cod_source (this variable is not used)
new_cod_source_coding <- grey_new_raw$cod_source |> attr("labels")

cod_source_recode <- tibble(
  new_code = new_cod_source_coding,
  old_code = (grey_old_raw$cod_source |> attr("labels"))[names(new_cod_source_coding)]
)

grey_new <- grey_new_raw |>
  select(-whoregion_name, -whoregion) |> # can add the regions back later
  select(-total_late_md) |> # this is not used in the grey data
  mutate(across(c("mmr", "num_t1", "num_t2", "num_t3"), as.numeric)) |>
  mutate(
    agestrat = plyr::mapvalues(agestrat, agestrat_recode$new_code, agestrat_recode$old_code),
    area = plyr::mapvalues(area, area_recode$new_code, area_recode$old_code), 
    timeframe = plyr::mapvalues(timeframe, timeframe_recode$new_code, timeframe_recode$old_code),
    va = plyr::mapvalues(va, va_recode$new_code, va_recode$old_code),
    denom = plyr::mapvalues(denom, denom_recode$new_code, denom_recode$old_code),
    cod_source = plyr::mapvalues(cod_source, cod_source_recode$new_code, cod_source_recode$old_code)
  ) |>
  mutate(across(c("agestrat", "area", "timeframe", "va", "denom", "cod_source"), as.integer))

# Check that there are no NAs in the variables that are used
stopifnot(all(!is.na(grey_new$agestrat)))
stopifnot(all(!is.na(grey_new$area)))
stopifnot(all(!is.na(grey_new$timeframe)))
stopifnot(all(!is.na(grey_new$denom)))

combined_grey <- bind_rows(grey_old, grey_new) 
  
stopifnot(nrow(combined_grey) == sum(nrow(grey_old), nrow(grey_new)))

# Take care of swapped start and end days in Botswana
# (GRY-9907203-BWA)
grey_data <- combined_grey %>% 
  select(-start_temp) 

# Make dates
grey_data <- grey_data |> 
  mutate(
    start_date = parse_date_time(datestart, "m y"),
     end_date = parse_date_time(dateend, "m y")
  ) |>
  # Take end of month for end dates
  mutate(end_date = ceiling_date(end_date, unit = "month", change_on_boundary = TRUE)) |>
  # Calculate study duration: 
  mutate(duration_days = end_date - start_date) |>
  # Take date midpoint:
  mutate(date = start_date + (end_date - start_date)/2) |> 
  mutate(year = lubridate::year(date)) |> 
  select(-date)

# Code national/subnational geographic level
grey_data <- grey_data |> 
  mutate(source = case_when( 
    area == 0 ~ "grey_national",
    area == 1 ~ "grey_ADM1+",
    area == 2 ~ "grey_below_ADM1")) |>
  select(-area)

# Check that all the age-stratified observations have non-age-stratified versions
grey_data |>
  group_by(idno_iso, start_date, end_date) |>
  summarise(
    has_agestrat = any(agestrat == 1), 
    has_not_agestrat = any(agestrat == 0),
    problem = has_agestrat & !has_not_agestrat
  ) |>
  filter(problem) |>
  (\(.x) nrow(.x) == 0)() |>
  stopifnot()
  

# Remove late deaths from total deaths in studies of timeframe==1
grey_data <- grey_data |> 
  mutate(num = if_else(timeframe==1 & !is.na(num_t3), num-num_t3, num)) 

# Remove timeframe==2 studies altogether
grey_data <- grey_data |> 
  filter(timeframe!=2)

# For a vector of date intervals,
# return TRUE iff any 2 UNIQUE intervals overlap
any_overlaps <- function(intervals) {
  interval_list <- as.list(intervals)
  
  for (i in 1:length(intervals)) {
    flag <- plyr::laply(
      interval_list,
      function(x) {int_overlaps(x, intervals[i]) & (x != intervals[i])}
    ) %>%
      any()
    
    if (flag) {
      return(TRUE)
    }
  }
  return(FALSE)
}


overlaps <- grey_data %>%
  mutate(date_interval = interval(start_date, (end_date - as.difftime(1, units = "secs")))) %>%
  distinct(idno, idno_iso, cod, start_date, end_date, date_interval) |> 
  # Flag studies that have overlapping date intervals
  group_by(idno_iso, cod) %>%
  summarise(any_overlaps = any_overlaps(date_interval)) 

if (sum(overlaps$any_overlaps) != 0) {
  stop("Unchecked overlapping date ranges!")
}

# Compute death counts where possible ----------
## Using pct and total_md
grey_data <- grey_data |> 
  mutate(prop_cause = sapply(
    prop_cause,
    function(x) eval(parse(text = x))
  )) %>%
  mutate(num = ifelse(
    is.na(num),
    round(0.01 * pct * total_md * prop_cause),
    num
  ))

## If still NA 
grey_data <- grey_data |> 
  mutate(num = ifelse(
    is.na(num),
    round(prop_cause*total_md*mmr/total_denom_mmr),
    num
  ))

# Deal with cause of death -----------

# read in clean icd-to-main cause file
icd <- read_csv(here("data/current/clean/icd10_to_main_cause_correspond_clean.csv"))

# clean up the cod and ICD guess columns
grey_data <- grey_data %>% 
  mutate(across(
    c("cod", "icd_guess1", "icd_guess2", "icd_guess3", "icd_guess4"),
    function(x) {str_remove(x, "\\.")}
  )) 

grey_data %<>%
  separate(
    col = cod,
    sep = ",",
    into = c("cod_1", "cod_2"),
    remove = FALSE,
    extra = "merge",
    fill = "right"
  ) %>%
  mutate(cause = cod_1)

## Split the data based on ICD_GUESS
grey_no_guesses <- grey_data %>%
  filter(!cod == "ICD_GUESS") %>% 
  rename(deaths_age_15_to_49 = num)

grey_guesses_only <- grey_data %>%
  filter(cod == "ICD_GUESS")


grey_guesses_only <- grey_guesses_only |> 
  pivot_longer(icd_guess1:icd_guess4, names_to = "guess", values_to = "cod_guess") |>
  separate(
    col = cod_guess,
    sep = " \\(",
    into = c("cod_code", "confidence"),
    remove = FALSE,
    extra = "merge",
    fill = "right"
  ) %>%
  # Remove things after the comma 
  mutate(cod_code = str_remove(cod_code, ",.*")) |>
  # retrieve the % "guess"
  mutate(percent = 0.01*as.numeric(str_remove(confidence, "%\\)"))) |> 
  # apportion caounts accordingly
  mutate(num = round(num * percent, 1)) |>
  # for the fractions
  mutate(halves = num %% 1 == 0.5) |>
  # for causes split 3 ways equally
  mutate(thirds = percent == 0.33) |>
  # deal with the halves (round up for guess 1 down for rest)
  mutate(num_new = case_when(
    !halves ~ round(num), # 
    guess == "icd_guess1" ~ ceiling(num),
    TRUE ~ floor(num)
  )) |>
  mutate(num_new = case_when(
    !thirds ~ num_new,
    guess == "icd_guess1" ~ ceiling(num),
    guess == "icd_guess2" & num - floor(num) > 0.5 ~ ceiling(num),
    TRUE ~ floor(num)
  )) |> 
  filter(!is.na(num_new))


## Clean a few things up
grey_guesses_only %<>%
  mutate(
    cause = cod_code,
    deaths_age_15_to_49 = num_new) %>%
  select(
    -guess,
    -cod_guess,
    -cod_code,
    -confidence,
    -percent,
    -num,
    -num_new,
    -halves,
    -thirds
  )

## Put split data sets back together

grey_all <- grey_no_guesses %>% 
  bind_rows(grey_guesses_only)
## Account for recent changes to some Ocodes
## (see data issues and changes Rmd)
grey_all %<>% 
  mutate(cause = case_when(cause == "O71" ~ "O719",
                           cause == "O99" ~ "Group 7",
                           TRUE ~ cause))



#  Assign main causes
grey_all %<>% 
  left_join(icd %>% select(cause,analysis_group),
            "cause") 



# For HIV assign NA not IND
grey_all %<>% 
  mutate(analysis_group = ifelse(cause == "O987", NA, analysis_group))


# Remove unwanted columns, clean names 
grey_all %<>% 
  select(
    source,
    idno_iso,
    iso,
    year,
    start_date,
    end_date,
    duration_days,
    cause,
    analysis_group,
    deaths_age_15_to_49,
    dataset
  ) %>% 
  rename(
    study_id = idno_iso
  )

# Save cleaned data
write_csv(grey_all, "data/current/clean/grey_data_clean.csv")

