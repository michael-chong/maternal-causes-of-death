library(readxl)
library(tidyverse)
library(magrittr)
library(here)

vr_new <- read_csv(here("data/update-2024/vr_mat-causes_21Feb2024.csv")) %>% 
  janitor::clean_names()

vr_old <- read_xls(here("data/current/raw/vr/Copy of maternal_causes_19mar2020.xls"),
                   guess_max = 30000) %>% 
    janitor::clean_names()

vr_old_clean <- read_csv(here("data/current/clean/vr_data_clean.csv"))


# The new VR file (21 Feb 2024) does not have country name or ISO, only 4 digit codes
# Pulling in an old country info file to try to get the right ISOs
country_info <- read_csv(here("data/current/raw/info/country_info_20180216.csv")) %>% 
  janitor::clean_names() %>%
  select(iso_code,
         country_name,
         country_name_long,
         whoname,
         whocode,
         whocod2)

countries <- read_csv(here("data/current/clean/country_list.csv"))

vr_new <- vr_new %>% 
  left_join(country_info, by = c("country"="whocod2"))

# isolate problematic numeric codes which cannot be matched
problem_country_codes <- vr_new %>% filter(is.na(iso_code)) %>%  distinct(country) 

vr_new <- vr_new %>% 
  mutate(source = "VR")

vr_data <- vr_new %>% 
  rename(who_code2 = country,
         iso = iso_code) %>% 
  select(- c(x1, admin1, sub_div, list,frmat, im_frmat, sex)) %>% 
  filter(year > 2008) %>% 
  select(source, iso, year, cause, tot, starts_with("a")) %>% 
  filter(!is.na(iso)) 

# Fix the age of death labelling
clean_names <- 
  readr::read_csv(here("data/current/raw/other/age_labelling_correspondence.csv")) %>% 
  rename(age_group = aug_dataset) %>% 
  mutate(age_group = tolower(age_group)) %>%
  select(-oct_dataset)

vr_temp <- vr_data %>% 
  pivot_longer(cols = c(tot, starts_with("a")),
               names_to = "age_group",
               values_to = "death_count") %>% 
  ## old age groups "translation" file has a95ov group which new VR data does not (it has a95)
  ## NOTE: making assumption that a95ov in previous iterations refers to same age group as a95 in new VR data
  mutate(age_group = case_when(
    age_group == "a95" ~ "a95ov",
    TRUE ~ age_group)
    ) %>% 
  left_join(clean_names, by = "age_group") %>% 
  select(-age_group) %>% 
  rename(age_group = new) %>% 
  pivot_wider(id_cols = c(source, iso, year, cause),
              names_from = age_group,
              values_from = death_count)
  

# Make main age_group of interest
vr_tidied <- vr_temp %>% 
  mutate(deaths_age_15_to_49 = rowSums(across(deaths_age_15_to_19:deaths_age_45_to_49), na.rm=T)) %>% 
  mutate(deaths_age_50_over = rowSums(across(deaths_age_50_to_54:deaths_age_95_over), na.rm=T)) %>% 
  select(source,
         iso,
         year,
         cause,
         deaths_all_ages,
         deaths_age_10_to_14,
         deaths_age_15_to_49, 
         deaths_age_50_over)

vr_tidied %>% 
  rowwise() %>% 
  mutate(flag = if_else(sum(deaths_age_10_to_14,
                            deaths_age_15_to_49, 
                            deaths_age_50_over,
                            na.rm = TRUE) > deaths_all_ages,
                        1,
                        0
  )) %>% 
  filter(flag == 1)

# Small fixes
vr_tidied <- 
  vr_tidied %>% 
  arrange(iso, year, cause) %>% 
  mutate(study_id = "VR") %>% 
  select(source, study_id, iso, year, cause, everything())

icd <- read_csv(here("data/current/clean/icd10_to_main_cause_correspond_clean.csv"))

vr_tidied <- vr_tidied %>% 
mutate(cause = case_when(cause == "O71" ~ "O719",
                           cause == "O99" ~ "Group 7",
                           TRUE ~ cause)) %>% 
  # Assign main causes (don't filter NA and "other" categories)
  left_join(icd %>% select(cause,analysis_group),
            by = "cause") %>% 
  # For HIV assign NA not IND
  mutate(analysis_group = ifelse(cause == "O987", NA, analysis_group))


vr_tidied <- 
  vr_tidied %>%
  select(
    source, 
    study_id,
    iso, 
    year, 
    cause, 
    analysis_group,
    deaths_age_15_to_49
  )


# Save the data ----------
write_csv(vr_tidied, "data/current/clean/vr_data_clean.csv")

