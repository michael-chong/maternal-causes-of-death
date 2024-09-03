library(tidyverse)
library(here)
library(caret)
library(readxl)


# Import data ------------------------------------------------------------------
# vr
df_vr <- read_csv(here("data/current/clean/vr_data_clean.csv"))
df_studies <- read_csv(here("data/current/clean/studies_data_clean.csv"))
df_grey <- read_csv(here("data/current/clean/grey_data_clean.csv"))

df_all<- bind_rows(
  df_vr,
  df_studies,
  df_grey %>% select(-dataset))

d_cod <-  df_all %>%
  transmute(source = source,
            study_id = study_id,
            iso = iso,
            #country = country_name,
            year = year,
            main_group = analysis_group,
            deaths = deaths_age_15_to_49,
            group7 = cause == "Group 7",
            o987 = cause == "O987",
            o98 = cause == "O98",
            o99 = (cause == "Group 7" & source == "VR"),
            o95 = cause == "O95",
            deaths_group7 = deaths*group7,
            deaths_o98 = deaths*o98,
            deaths_o95 = deaths*o95,
            deaths_o987 = deaths*o987,
            start_date, 
            end_date) %>%
  filter(year >= 2009, year <= 2020) %>% 
  group_by(source, study_id, iso, year, main_group, start_date, end_date) %>%
  summarise(deaths = sum(deaths),
            group7 = any(group7),
            o987 = any(o987),
            deaths_group7 = sum(deaths_group7),
            deaths_o98 = sum(deaths_o98),
            deaths_o987 = sum(deaths_o987),
            deaths_o95 = sum(deaths_o95),
            o98 = any(o98),
            o99 = any(o99)) %>%
  group_by(iso, study_id, year) %>%
  mutate(group7 = any(group7),
         o987 = any(o987),
         o98 = any(o98),
         o99 = any(o99),
         deaths_o98 = sum(deaths_o98),
         deaths_group7 = sum(deaths_group7),
         deaths_o987 = sum(deaths_o987),
         deaths_o95 = sum(deaths_o95)) %>%
  ungroup() %>%
  filter(!is.na(main_group)) %>% 
  filter(main_group != "other") %>%
  spread(key = main_group, value = deaths)%>%
  rowwise() %>%
  mutate(TOT = sum(ABO,DIR,EMB,HEM,SEP,IND,HYP, na.rm = TRUE))  

# Import country list
d_countries <- read_csv(here("data/current/clean/country_list.csv")) %>%
  # clean up some column names
  transmute(iso = iso,
            name_official = country,
            name = name_clean, 
            include = include,
            region_mdg = as_factor(region_mdg),
            region_sdg1 = as_factor(region_sdg1),
            region_sdg2 = as_factor(region_sdg2),
            region_sdg3 = as_factor(region_sdg3),
            region_who = as_factor(region_who),
            region_wb = as_factor(region_wb),
            region_alkema = as_factor(region_alkema),
            region_gbd = as_factor(region_gbd),
            region_epi = as_factor(region_epi))


# Import MMR, aids deaths, and maternal deaths data
# MMR
d_mmr <- read_csv(here("data/update-2024/other_data_and_info/estimates.csv")) %>% 
  rename(mmr = `0.5`,
         iso = iso_alpha_3_code,
         year = year_mid) %>%
  select(iso, year, parameter, mmr) %>% 
  filter(parameter == "mmr") %>% 
  select(-parameter) %>% 
  filter(between(year, 2009, 2020)) %>% 
  mutate(log_mmr_standard = as.numeric(scale(log(mmr)))) %>% 
  group_by(iso) %>% 
  mutate(country_mean_log_mmr_stnd = mean(log_mmr_standard),
         country_mean_mmr = mean(mmr)) %>% 
  ungroup()
  
# Maternal deaths data
d_matdeaths <- read_csv(here("data/update-2024/other_data_and_info/estimates.csv")) %>% 
  rename(matdeaths = `0.5`,
         iso = iso_alpha_3_code,
         year = year_mid) %>%
  select(iso, year, parameter, matdeaths) %>% 
  filter(parameter == "maternal_deaths") %>% 
  select(-parameter) 
  

# Import HIV/AIDS estimates
d_aids <- read_csv(here("data/update-2024/other_data_and_info/estimates.csv")) %>% 
  rename(aidsdeaths = `0.5`,
         iso = iso_alpha_3_code,
         year = year_mid) %>%
  select(iso, year, parameter, aidsdeaths) %>% 
  filter(parameter == "hiv_related_indirect_maternal_deaths") %>% 
  select(-parameter) 

# Subtract AIDS deaths from total maternal deaths
d_matdeaths <- d_matdeaths %>% 
  left_join(d_aids) %>%
  mutate(matdeaths = matdeaths - aidsdeaths) %>%
  select(-aidsdeaths)

# Import WPP estimates
wpp <- readxl::read_excel(here::here("data/update-2024/other_data_and_info/WPP2022_GEN_F01_DEMOGRAPHIC_INDICATORS_COMPACT_REV1.xlsx"), 
                          skip = 16, 
                          col_types = c("ISO3 Alpha-code"= "text"))  %>%
  mutate(country = `Region, subregion, country or area *`,
         iso = `ISO3 Alpha-code`,
         year = as.numeric(Year),
         deaths = 1000* as.numeric(`Female Deaths (thousands)`)) %>%
  mutate(country = stringi::stri_trans_general(country, id = "Latin-ASCII")) %>%
  rename(deaths_wpp = deaths) %>% 
  ## select only country estimates
  filter(!is.na(iso)) %>% 
  select(country, iso, year, deaths_wpp)

# Deal with contributing causes, usability, and type  ---------------------------------
icd_issues <- read_csv(here("data/current/clean/icd10_codes_vr_quality_issues.csv"))
contributing <- icd_issues %>% 
  filter(reason == "contributory") %>%
  distinct()

df_cont <- df_all %>% 
  left_join( (contributing %>% select(cause, reason)), by = "cause") %>% 
  left_join(d_countries %>% select(iso, name), by = "iso") %>% 
  filter(!is.na(analysis_group)) %>% 
  filter( analysis_group != "other") %>%
  mutate(reason = replace_na(reason, "underlying")) %>% 
  group_by(iso,  year, source, study_id) %>% 
  summarise(prop_contributory =
              sum(deaths_age_15_to_49*(reason == "contributory")) / sum(deaths_age_15_to_49), 
            total_counts = sum(deaths_age_15_to_49) ) %>% 
  ungroup()


# Import all-cause and ill-defined counts, calculate coverage
# load (older) countries info file to match "mystery" numeric country codes to appropriate country iso
country_info <- read_csv(here("data/current/raw/info/country_info_20180216.csv")) %>% 
  janitor::clean_names() %>%
  select(iso_code,
         country_name,
         country_name_long,
         whoname,
         whocode,
         whocod2)

d_u <- read_csv(here("data/update-2024/vr_21_Feb_2024.csv")) %>%
  janitor::clean_names() %>% 
  select(-c(x1,admin1, sub_div, list)) %>% 
  left_join(country_info %>% select(iso_code, whocod2), by = c("country"="whocod2")) %>% 
  rename(iso = iso_code) %>% 
  left_join(d_countries %>% select(iso, name), by = "iso") %>% 
  filter(!is.na(iso)) %>% 
  distinct(iso, year, cause, .keep_all = TRUE) %>% 
  # only keep all deaths and ill defined
  filter(cause %in% c("all","ill")) %>%
  rename(deaths = tot) %>% 
  select(country, iso, name, year, cause, deaths) %>% 
  pivot_wider(names_from = cause, values_from=deaths) %>% 
  left_join(wpp %>% select(-country), by = c("iso", "year")) %>%
  mutate(ill_prop = ill/all) %>% # calculate proportion of ill-defined
  mutate(coverage = pmin(all/deaths_wpp, 1)) %>% # calculate coverage
  mutate(u = coverage*(1 - ill_prop)) %>% # calculate usability
  filter(year >= 2007) %>% # this is the most we need to "reach back" to determine the type for 2009
  left_join((df_cont %>%
               filter( source == "VR") %>% 
               select(iso, year, prop_contributory, total_counts)),
            by = c("iso" = "iso", "year" = "year")) %>% 
  # define new usability index based on prop_contributory:
  mutate(u_new = ifelse(
    (total_counts <= 5 | is.na(total_counts)), 
    u,
    (1 - prop_contributory)*coverage*(1 - ill_prop))) %>% 
  #for consistency with code downstream
  # remove country numeric code
  select(-country) %>% 
  # rename name into country
  rename(country = name)


# make new variables where years are lagged
to_join <- d_u %>%
  transmute(
    iso = iso, 
    year_m1 = year - 1, 
    year_m2 = year - 2,
    year_p1 = year + 1,
    year_p2 = year + 2,
    u = u,
    u_new = u_new
  )

# join using lagged years
d_t <- left_join(d_u, to_join, by = c("year" = "year_m1", "iso" = "iso"), suffix = c("", "_m1")) %>%
  left_join(to_join, by = c("year" = "year_m2", "iso" = "iso"), suffix = c("", "_m2")) %>%
  left_join(to_join, by = c("year" = "year_p1", "iso" = "iso"), suffix = c("", "_p1")) %>%
  left_join(to_join, by = c("year" = "year_p2", "iso" = "iso"), suffix = c("", "_p2")) %>%
  select(country, iso, year, starts_with("u")) %>%
  # Indicate whether surrounding years have u > 0.6:
  mutate_at(vars(contains("u_p")|contains("u_m")), ~ (. > .6) ) %>%
  mutate_at(vars(contains("new_p")|contains("new_m")), ~ (. > .65) ) %>% 
  rowwise() %>%
  # Check whether ANY consecutive 3 years have u > 0.6:
  mutate(t1 = prod(u, u_m1, u_m2),
         t2 = prod(u, u_m1, u_p1),
         t3 = prod(u, u_p1, u_p2),
         t1_new = prod(u_new, u_new_m1, u_new_m2),
         t2_new = prod(u_new, u_new_m1, u_new_p1),
         t3_new = prod(u_new, u_new_p1, u_new_p2)
  ) %>%
  mutate(t = max(0, t1, t2, t3, na.rm = TRUE), 
         t_new = max(0, t1_new, t2_new, t3_new, na.rm = TRUE)) %>%
  # Assign types as appropriate:
  mutate(type = 
           case_when(
             t > 0.8 ~ 1,
             t > 0.6 ~ 2,
             TRUE ~ 3
           ),
         type_new = 
           case_when(
             t_new > 0.85 ~ 2,
             t_new > 0.65 ~ 3,
             TRUE ~ 4
           )) %>%
  filter(year >= 2009) %>%
  ungroup() %>%
  mutate(source = "VR") %>%
  select(source, country, iso, year, u, u_new, type, type_new) 


# Calculate maternal coverage -------------------------------------
library(lubridate)

# (objects to use in calculation later)
all_days <- seq(ymd("2009-01-01"), ymd("2020-12-31"), by = "day")
matdeaths_daily <- data.frame(
  day = all_days
) %>%
  mutate(year = year(day))%>%
  group_by(year) %>%
  mutate(n_days = n()) %>%
  ungroup() %>%
  left_join(d_matdeaths, by = "year") %>%
  rename(md_iso= iso) %>%
  mutate(matdeath_daily = matdeaths/n_days)

aidsdeaths_daily <- data.frame(
  day = all_days
) %>%
  mutate(year = year(day))%>%
  group_by(year) %>%
  mutate(n_days = n()) %>%
  ungroup() %>%
  left_join(d_aids, by = "year") %>%
  rename(aids_iso= iso) %>%
  mutate(aidsdeath_daily = aidsdeaths/n_days)

# First join matdeaths for VR:
d_cod <- left_join(d_cod, d_matdeaths %>% mutate(source = "VR"), by = c("iso", "year", "source")) %>%
  rowwise() %>%
  # Then use start and end dates to get matdeaths for studies:
  mutate(matdeaths = ifelse(is.na(matdeaths),
                            sum(
                              matdeaths_daily %>%
                                filter(md_iso == iso, day >= start_date, day < end_date) %>%
                                pull(matdeath_daily)), matdeaths)) %>%
  # calculate mat coverage
  # (add back O95 and O98.7 deaths to match their inclusion in envelope estimates)
  mutate(mat_coverage = min((TOT + deaths_o987 + deaths_o95)/matdeaths, 1)) %>% 
  group_by(iso) %>%
  mutate(max_mat_coverage = min(max(mat_coverage), 1)) %>%
  ungroup()


# Subtract HIV/AIDS estimate ---------------------------------------------------
# Subtract when Group 7 appears in that study-year or O98 code was reported
# remove all O987 deaths (if any) from the adjusted AIDS deaths before subtracting from IND
# First join aidsdeaths for VR:
d_cod <- left_join(d_cod, d_aids %>% mutate(source = "VR"), by = c("iso", "year", "source")) %>%
  rowwise() %>%
  mutate(aidsdeaths = if_else(is.na(aidsdeaths),
                              sum(
                                aidsdeaths_daily %>%
                                  filter(aids_iso == iso, day >= start_date, day < end_date) %>%
                                  pull(aidsdeath_daily)),
                              aidsdeaths)) %>%
  mutate(
    possibly_aids = ifelse(
      source == "VR", 
      pmin(sum(deaths_group7, deaths_o98), mat_coverage*aidsdeaths),
      pmin(deaths_o98, mat_coverage*aidsdeaths)
    ),
    aids_unaccounted_for = pmax(0, possibly_aids - deaths_o987)
  ) %>%
  mutate(IND = IND - aids_unaccounted_for) %>%
   mutate(TOT = sum(ABO,DIR,EMB,HEM,SEP,IND,HYP, na.rm = TRUE))


# Put data frames together -----------------------------------------------------
d <- d_cod %>%
  left_join(d_countries, by = "iso") %>%
  left_join(d_mmr, by = c("iso", "year")) %>% 
  left_join(select(d_t, -country), by = c("iso", "year", "source")) %>%
  ungroup() %>%
  # set PSE to be Type 2 because we don't have all-cause deaths here
  mutate(type = ifelse(iso == "PSE", 2, type),
         type_new = ifelse(iso == "PSE", 3, type_new)) %>% 
  mutate(type = ifelse(source == "VR", type, 3) %>% factor(),
         type_new = case_when(
           source == "VR" ~ type_new,
           source == "grey_national" ~ 1,
           TRUE ~ 4) %>% factor()) %>%
  # adjust down type_new from 1 to 2 for following grey
  mutate(type_new = ifelse(study_id %in% c("GRY-9900401-AFG", "GRY-9905001-BGD", "GRY-9928801-GHA"), 2, type_new) %>% factor()) %>% 
  filter(include) %>%
  filter(TOT > 0) %>%
  #filter(!(iso == "SVK" & study_id == "VR" & year == 2010)) %>%
  filter(year < 2021) %>% 
  #change existing country name to be name_official
  #select(-country) %>% 
  rename(country=name_official) %>% 
  select(source:iso, country, everything()) %>%
  select(-start_date, -end_date)

## Make the subcauses data frames ----------------------------------------------

# make df_all compatible with d
df_all <- df_all %>% 
  left_join((d_countries %>% select(iso, name_official, include)), by = "iso") %>% 
  filter(include) %>% 
  filter(year < 2021) %>% 
  filter(!is.na(analysis_group)) %>% 
  filter(analysis_group != "other") %>% 
  rename(deaths = deaths_age_15_to_49, 
         country = name_official) %>% 
  select(source, study_id, iso, country, everything())

# DIR subcauses data
df_dir <- read_xlsx(here("data/current/raw/info/other_DIR.xlsx"))
df_dir <- df_dir %>% 
  rename(dir_sub = `other-DIR`) %>% 
  mutate(dir_sub = case_when(
    dir_sub == "other" ~ "oth",
    dir_sub=="obtrauma_other" ~ "obt",
    dir_sub=="obstructedlabour" ~ "obs",
    dir_sub=="anesthesia"| dir_sub=="anaesthesia" ~ "ane",
    TRUE~"NA"
  )) %>% 
  mutate(subanalysis_group = paste(analysis_group, dir_sub, sep = "_"))
## Manually add missing subanalysis group for ICD codes that were previously not there
# (refer to email from Sahar Ahmed on 26-April-2024)
add_dir <- tribble(
  ~analysis_group, ~cause, ~dir_sub, ~subanalysis_group,
  "DIR", "O611", NA, NA,
  "DIR", "O650", "obs", "DIR_obs",
  "DIR", "O294", "ane", "DIR_ane",
  "DIR", "O332", "obt", "DIR_obt",
  "DIR", "O262", "oth", "DIR_oth",
  "DIR", "O893", "ane", "DIR_ane",
  "DIR", "O331", "obt", "DIR_obt",
  "DIR", "O755", "oth", "DIR_oth",
  "DIR", "O668", "obs", "DIR_obs"
)

df_dir <- bind_rows(df_dir, add_dir)

# Data on timing of death subcauses for HEM and SEP
df_timing <- read_xlsx(here("data/current/raw/info/ocodes-timing.xlsx"))
df_timing <- df_timing %>% 
  select(analysis_group, cause, timing) %>% 
  mutate(timing = recode(timing, 
                         antepartum = "ante",
                         intrapartum = "intra",
                         postpartum = "post",
                         unknown = "unkwn")) %>% 
  mutate(subanalysis_group = paste(analysis_group, timing, sep="_")) %>% 
  filter(analysis_group %in% c("HEM", "SEP"))

# Manually adding SEP code whith missing timing assignment (O910)
add_sep <- tribble(
  ~analysis_group, ~cause, ~timing, ~subanalysis_group,
  "SEP", "O910", "unkwn", "SEP_unkwn"
)
  
df_timing <- bind_rows(df_timing, add_sep)

# DIR subcauses
d_DIR <- df_all  %>%
  filter(analysis_group == "DIR") %>% 
  left_join(df_dir %>% select(cause, subanalysis_group), by = "cause") %>% 
  filter(!is.na(subanalysis_group)) %>% 
  group_by(source, study_id,country, iso, year, subanalysis_group, start_date, end_date) %>%  
  summarise(sub_deaths = sum(deaths)) %>% 
  pivot_wider(values_from = sub_deaths, names_from = subanalysis_group) %>% 
  rowwise() %>% 
  mutate(TOT_DIR = sum(DIR_obs, DIR_ane, DIR_obt, DIR_oth, na.rm = T)) %>% 
  ungroup() %>%
  left_join(select(d_countries, iso, starts_with("region"))) %>%
  filter(TOT_DIR > 0)

# HEM subcauses 
d_HEM <- df_all %>% 
  filter(analysis_group == "HEM") %>% 
  left_join((df_timing %>% filter(analysis_group == "HEM") %>% select(cause, subanalysis_group)), by = "cause") %>% 
  group_by(source, study_id,country, iso, year, subanalysis_group, start_date, end_date) %>%  
  summarise(sub_deaths = sum(deaths)) %>% 
  pivot_wider(values_from = sub_deaths, names_from = subanalysis_group) %>% 
  rowwise() %>% 
  mutate(TOT_HEM = sum(HEM_ante, HEM_intra, HEM_post, na.rm = T)) %>% 
  ungroup() %>%
  left_join(select(d_countries, iso, starts_with("region"))) %>%
  filter(TOT_HEM > 0) %>%
  # Calculate proportion of unknown HEM
  mutate(prop_known_HEM = TOT_HEM/(TOT_HEM + pmax(HEM_unkwn, 0, na.rm = TRUE))) %>%
  group_by(iso) %>%
  mutate(max_prop_known_HEM = max(prop_known_HEM)) %>%
  ungroup()

# SEP subcauses 
d_SEP <- df_all %>% 
  filter(analysis_group == "SEP") %>% 
  left_join((df_timing %>% filter(analysis_group == "SEP") %>% select(cause, subanalysis_group)), by = "cause") %>% 
  filter(!is.na(subanalysis_group)) %>% 
  group_by(source, study_id,country, iso, year, subanalysis_group, start_date, end_date) %>%  
  summarise(sub_deaths = sum(deaths), .groups = "keep") %>% 
  pivot_wider(values_from = sub_deaths, names_from = subanalysis_group) %>% 
  rowwise() %>% 
  mutate(TOT_SEP = sum(SEP_ante, SEP_intra, SEP_post, na.rm = T)) %>% 
  ungroup()  %>%
  select(source:end_date, SEP_unkwn, SEP_ante, SEP_intra, SEP_post, TOT_SEP) %>%
  left_join(select(d_countries, iso, starts_with("region"))) %>%
  filter(TOT_SEP > 0) %>%
  # Calculate proportion of unknown SEP
  mutate(prop_known_SEP = TOT_SEP / (TOT_SEP + pmax(SEP_unkwn,0, na.rm = TRUE))) %>%
  group_by(iso) %>%
  mutate(max_prop_known_SEP = max(prop_known_SEP)) %>%
  ungroup()

# Write useable dataframe
#(don't save separately here, use d in objects for modeling instead)
#write_rds(d, here::here("data/current/clean/main_causes_analysis_ready.rds"))

# Other objects ---------------------------------------------------------------
# * Hard code relevant years -----
years <- 2009:2020
nyears <- length(years)

# * Full data frame with covariates to generate predictions ---------
d_complete <- d_mmr %>%
  filter(year >= 2009 & year <= 2020) %>%
  left_join(d_countries, by = c("iso")) %>%
  left_join(d_matdeaths, by = c("iso", "year")) %>%
  left_join(d_aids, by = c("iso", "year")) %>% 
  left_join((d %>% select(iso, max_mat_coverage) %>% distinct), by = "iso") %>% 
  mutate(max_mat_coverage = replace_na(max_mat_coverage, 0)) %>% 
  mutate(type = factor(1, levels = c(1, 2, 3))) 


country_mat <- dummyVars(~iso, levelsOnly = TRUE, data = d_complete) %>% 
  predict(d_complete) %>% 
  t()
# * Region model matrix for region estimates ----
region_mdg_mat <- dummyVars(~region_mdg, levelsOnly = TRUE, data = d_complete) %>% 
  predict(d_complete) %>% 
  t()

# * Region-year model matrix for region-year estimates ----
region_mdg_year_mat <- dummyVars(~1 + region_mdg:year, levelsOnly = TRUE, data = d_complete %>% mutate(year = as_factor(year))) %>%
  predict(d_complete %>% mutate(year = as_factor(year))) %>%
  t()

# * SDG region matrix for region estimates
region_sdg_mat <- dummyVars(~region_sdg1, levelsOnly = TRUE, data = d_complete) %>% 
  predict(d_complete) %>% 
  t()

# * SDG Region-year model matrix for region-year estimates ----
region_sdg_year_mat <- dummyVars(~1 + region_sdg1:year, levelsOnly = TRUE, data = d_complete %>% mutate(year = as_factor(year))) %>%
  predict(d_complete %>% mutate(year = as_factor(year))) %>%
  t()

# * Labels and factors -----------
# For use in mapping to and from factor levels and numbers
region_mdg_labels <- data.frame(name = levels(d$region_mdg), number = 1:nlevels(d$region_mdg), stringsAsFactors = FALSE)
region_sdg_labels <- data.frame(name = levels(d$region_sdg1), number = 1:nlevels(d$region_sdg1), stringsAsFactors = FALSE)

# Give an error if these are not the same as in region_mat
if (!all(rownames(region_mdg_mat) == region_mdg_labels$name)) {
  stop("MDG region labels inconsistent")
}
if (!all(rownames(region_sdg_mat) == region_sdg_labels$name)) {
  stop("SDG region labels inconsistent")
}

region_alkema_labels <- data.frame(name = levels(d$region_alkema), number = 1:nlevels(d$region_alkema), stringsAsFactors = FALSE)
region_wb_labels <- data.frame(name = levels(d$region_wb), number = 1:nlevels(d$region_wb), stringsAsFactors = FALSE)
region_gbd_labels <- data.frame(name = levels(d$region_gbd), number = 1:nlevels(d$region_gbd), stringsAsFactors = FALSE)
region_epi_labels <- data.frame(name = levels(d$region_epi), number = 1:nlevels(d$region_epi), stringsAsFactors = FALSE)

# For use in mapping to and from countries
country_labels <- data.frame(
  name = unique(d_complete$iso), 
  number = 1:length(unique(d_complete$iso)), stringsAsFactors = FALSE) %>%
  left_join(select(d_countries, iso, starts_with("region")), by = c("name" = "iso"))
ncountry <- nrow(country_labels)

# Hard code cause names
causes <- c("ABO", "DIR", "EMB", "HEM", "SEP", "IND", "HYP")
ncause <- length(causes)
cause_labels <- data.frame(name = causes, number = 1:ncause, stringsAsFactors = FALSE)


subcauses_DIR <- c("DIR_obs", "DIR_ane", "DIR_obt", "DIR_oth")
nsubcause_DIR <- length(subcauses_DIR)
subcauses_DIR_labels <- data.frame(name = subcauses_DIR, number = 1:nsubcause_DIR, stringsAsFactors = FALSE)

subcauses_HEM <- c("HEM_ante", "HEM_intra", "HEM_post")
nsubcause_HEM <- length(subcauses_HEM)
subcauses_HEM_labels <- data.frame(name = subcauses_HEM, number = 1:nsubcause_HEM, stringsAsFactors = FALSE)


subcauses_SEP <- c("SEP_ante", "SEP_intra", "SEP_post")
nsubcause_SEP <- length(subcauses_SEP)
subcauses_SEP_labels <- data.frame(name = subcauses_SEP, number = 1:nsubcause_SEP, stringsAsFactors = FALSE)

# * Observation dataframe as matrix to feed into Stan ----
d_mat <- d %>%
  select(all_of(causes)) %>%
  as.matrix()

d_mat_DIR <- d_DIR %>% 
  select(all_of(subcauses_DIR)) %>% 
  as.matrix()

d_mat_HEM <- d_HEM %>% 
  select(all_of(subcauses_HEM)) %>% 
  as.matrix()

d_mat_SEP <- d_SEP %>% 
  select(all_of(subcauses_SEP)) %>% 
  as.matrix()

# save locations of true zeroes
true_zeroes <- d_mat == 0
true_zeroes[is.na(true_zeroes)] <- FALSE


true_zeroes_DIR <- d_mat_DIR == 0
true_zeroes_DIR[is.na(true_zeroes_DIR)] <- FALSE

true_zeroes_HEM <- d_mat_HEM == 0
true_zeroes_HEM[is.na(true_zeroes_HEM)] <- FALSE

true_zeroes_SEP <- d_mat_SEP == 0
true_zeroes_SEP[is.na(true_zeroes_SEP)] <- FALSE

# Put zeroes in place of NAs in Stan's matrix
d_mat[is.na(d_mat)] <- 0

d_mat_DIR[is.na(d_mat_DIR)] <- 0
d_mat_HEM[is.na(d_mat_HEM)] <- 0
d_mat_SEP[is.na(d_mat_SEP)] <- 0

# zero_mat == 1 where observation is > 0, or is a true zero
# dictates which values contribute to likelihood
zero_mat <- 1*(d_mat > 0 | true_zeroes)

zero_mat_DIR <- 1*(d_mat_DIR > 0 | true_zeroes_DIR)
zero_mat_HEM <- 1*(d_mat_HEM > 0 | true_zeroes_HEM)
zero_mat_SEP <- 1*(d_mat_SEP > 0 | true_zeroes_SEP)

d_complete <- d_complete %>%
  left_join(d_HEM %>% select(iso, max_prop_known_HEM) %>% distinct()) %>%
  left_join(d_SEP %>% select(iso, max_prop_known_SEP) %>% distinct()) %>%
  mutate(max_prop_known_HEM = replace_na(max_prop_known_HEM, 0)) %>%
  mutate(max_prop_known_SEP = replace_na(max_prop_known_SEP, 0))

# Save objects to use in model scripts ------------------------
save(
  d,
  d_DIR,
  d_HEM,
  d_SEP,
  d_mat,
  d_mat_DIR,
  d_mat_HEM,
  d_mat_SEP,
  d_complete, 
  zero_mat,
  zero_mat_DIR,
  zero_mat_HEM,
  zero_mat_SEP,
  country_mat,
  region_mdg_mat, 
  region_mdg_year_mat, 
  region_mdg_labels,
  region_sdg_mat,
  region_sdg_year_mat,
  region_sdg_labels,
  region_alkema_labels,
  region_wb_labels,
  region_gbd_labels,
  country_labels,
  ncountry,
  causes,
  ncause,
  cause_labels,
  subcauses_DIR,
  nsubcause_DIR,
  subcauses_DIR_labels,
  subcauses_HEM,
  nsubcause_HEM,
  subcauses_HEM_labels,
  subcauses_SEP,
  nsubcause_SEP,
  subcauses_SEP_labels,
  years,
  nyears,
  file = here::here("data/current/clean/objects_for_modelling.RData")
)

