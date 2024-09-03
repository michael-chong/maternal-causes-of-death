library(tidyverse)
library(lubridate)
library(magrittr)
library(here)
library(haven)

# Read in the studies data
studies_new <- read_dta(here("data/update-2024/round1_round2_appended_12022024.dta")) %>% 
  mutate(prefix = str_extract(idno, "^[A-Z]+")) %>% 
  filter(!prefix == "GRY")

studies_old <- read_csv(here("data/current/clean/studies_data_clean.csv"))
studies_old_raw <- studies_data <- read_dta(here("data/current/raw/studies/journal-2020-05-16.dta"))

# Deal with dates ---------------------------
studies_new <- studies_new %>%
  # if month not specified, assign Jan if start date, Dec if end date
  mutate(datestart = str_replace(datestart, "^99", "01"),
         dateend = str_replace(dateend, "^99", "12")) %>%
  mutate(start_date = parse_date_time(datestart, "m y"),
         end_date = parse_date_time(dateend, "m y")) %>%
  # Take end of month for end dates
  mutate(end_date = ceiling_date(end_date, unit = "month", change_on_boundary = TRUE)) %>%
  # Calculate study duration: 
  mutate(duration_days = end_date - start_date) %>%
  # Take date midpoint:
  mutate(date = start_date + (end_date - start_date)/2) %>%
  mutate(year = lubridate::year(date)) %>%
  select(-date)


any_overlaps <- function(intervals) {
  interval_list <- as.list(intervals)
  
  for (i in 1:length(intervals)) {
    flag <- plyr::laply(interval_list,
                        function(x) {int_overlaps(x, intervals[i]) & (x != intervals[i])}
    ) %>%
      any()
    
    if (flag) {
      return(TRUE)
    }
  }
  return(FALSE)
}

# check for overlapping date ranges
overlaps <- studies_new %>%
  mutate(date_interval = interval(start_date, (end_date - as.difftime(1, units = "secs")))) %>%
  # Flag studies that have overlapping date intervals
  group_by(idno_iso, cod) %>%
  summarise(any_overlaps = any_overlaps(date_interval)) 

if (sum(overlaps$any_overlaps) != 0) {
 stop("Unchecked overlapping date ranges!")
}
## ^ !!Overlapping intervals; Deal with them further below

attr(studies_new$area, "labels")

# for some variables value labels have changed since last version of the data 
# but since this round has both old and new studies (unlike grey)
# okay to just change the codes directly; no real need for mapping 

# Recode national/subnational
studies_new <- studies_new %>%
  mutate(source = case_when(
    area %in% c(1,5) ~ "studies_national",
    area == 2 ~ "studies_ADM1+",
    area == 3 ~ "studies_below_ADM1"
  )) %>%
  select(-area)


attr(studies_new$agestrat, "labels")
studies_new %>% select(agestrat) %>%  table(.)

# problematic studies: COC-12633-CMR; COC-15148-NGA;
 
studies_new <- studies_new %>% 
  mutate(flag_agestrat = ifelse(agestrat == 2, TRUE, FALSE))
## NOTE: two studies have age-"stratified" counts
## COC-12633-CMR only reports on adolescents deaths (see issues Rmd) and N is missing
## COC-15148-NGA has aggregated and 15-19 death counts; for those 15-19 N is missing and only % of cause is given
## FIX: keep all Cameroon observations; keep only 15-49 for Nigeria 
studies_new <- studies_new %>% 
  filter(!(idno_iso=="COC-15148-NGA" & ageend==19))

# Tidy up MMR variable
studies_new <- studies_new %>% 
  # first remove all parenthesis and everything inside them
  # this takes care of the "(SE)" part in the variable
  mutate(mmr = gsub("\\s*\\([^\\)]+\\)",
                    "",
                    as.character(mmr))) %>% 
  # next remove the "MMR = " part
  mutate(mmr = as.numeric(str_replace(mmr, "MMR = ", ""))) 

# Impute the missing num
#
studies_new <- studies_new %>%
  mutate(prop_cause = ifelse(
    (prop_cause != "Unknown"),
    prop_cause,
    NA
  )) %>%
  mutate(prop_cause = sapply(
    prop_cause,
    function(x) eval(parse(text = x))
  )) 

stopifnot(all(studies_new$prop_cause <= 1))

studies_new <- studies_new %>% 
  mutate(prop_cause = case_when(
    #prop_cause is wrong for all observations in both studies
    idno_iso %in% c("COC-18027-IND", "COC-2182-CIV") ~ 1,
    TRUE ~ prop_cause
  ))

{if(!all(studies_new$prop_cause<=1, na.rm = T))stop("Not all prop_cause are <= 1")} 

studies_new <- studies_new %>%
  # if N not available and pct (% deaths in COD category) available then
  mutate(num = ifelse(
    is.na(num),
    round(0.01 * pct * total_md * prop_cause),
    num
  )) %>%
  # if N (num) not available and pct not available, then use mmr
  mutate(num = ifelse(
    is.na(num),
    round(prop_cause * total_md * mmr / total_denom_mmr),
    num
  ))

# Fill in the remaining num either using MMR and live births
# or MMR and total births

studies_new <- studies_new %>%
  mutate(num = ifelse(
    is.na(num),
    (ifelse(!is.na(prop_cause),
            round(prop_cause * mmr * total_denom / 100000),
            round(mmr * total_denom / 100000))),
    num
  )) %>%
  # Deal with the remaining num  still left as NA
  mutate(num = ifelse(
    (is.na(num) & is.na(prop_cause) & is.na(total_denom)),
    round(total_md * mmr / total_denom_mmr),
    num
  ))


{if(any(is.na(studies_new$num))) stop("there are still NA values in num")}
# Deal with cause of death ---------------------------

# clean up the ICD guess columns
# read in clean icd-to-main cause file
icd <- read_csv(here("data/current/clean/icd10_to_main_cause_correspond_clean.csv"))

studies_new <- studies_new %>% 
  mutate_at(vars(cod, icd_guess1, icd_guess2, icd_guess3, icd_guess4), list(~str_remove(., "\\."))) %>%
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
studies_no_guesses <- studies_new %>%
  filter(!str_detect(cod, "ICD_GUESS"))

studies_guesses_only <- studies_new %>%
  filter(str_detect(cod, "ICD_GUESS"))


studies_guesses_only <- studies_guesses_only %>% 
  # not much we can do about ICD_GUESS_R when only organ is given (id_rorgan)
  # so filter out if icd_guess1 is left empty
  filter(!(cod == "ICD_GUESS_R" & icd_guess1=="")) %>% 
  pivot_longer(icd_guess1:icd_guess4, names_to = "guess", values_to = "cod_guess") %>%
  separate(
    col = cod_guess,
    sep = " \\(",
    into = c("cod_g_1", "cod_g_2"),
    remove = FALSE,
    extra = "merge",
    fill = "right"
  ) %>%
  ## clean up cod_g_1
  mutate(cod_g_1 = gsub(",.*", "", cod_g_1)) %>%
  ## retrieve the % "guess"
  mutate(cod_g_2 = str_remove(cod_g_2, "%\\)")) %>%
  mutate(cod_g_2 = 0.01 * as.numeric(cod_g_2)) %>%
  ## Change the icd guess % for COC-55348-TZA manually (see data issues Rmd)
  mutate(cod_g_2 = ifelse(idno_iso == "COC-55348-TZA" & cod_g_1 == "O85", 0.17, cod_g_2)) %>% 
  mutate(cod_g_2 = ifelse(idno_iso == "COC-55348-TZA" & cod_g_1 == "EXCLUDE", 0.24, cod_g_2))%>% 
  ## apportion counts accordingly
  mutate(num = round(num * cod_g_2, 1)) %>%
  ## deal with fractions
  mutate(halves = num %% 1 == 0.5) %>%
  ## When halves true and  num is odd, 
  ## split the counts
  ## by rounding up the for icd_guess1
  ## and rounding down for icd_guess2
  ## Making anote here that Tanzania (COC-55348-TZA) may be problematic
  mutate(num = case_when(
    !halves ~ round(num),
    halves & guess == "icd_guess1" ~ ceiling(num),
    halves & !guess == "icd_guess1" ~ floor(num)
  )) %>% 
  filter(!is.na(num))

## Clean a few things up
studies_guesses_only %<>%
  mutate(cause = cod_g_1) %>%
  select(
    -guess,
    -cod_guess,
    -cod_g_1,
    -cod_g_2,
    -halves,
  )


## Put split data sets back together
setdiff(names(studies_guesses_only), names(studies_no_guesses))
setdiff(names(studies_no_guesses), names(studies_guesses_only))

studies_all <- studies_no_guesses %>% 
  select(-c(icd_guess1, icd_guess2, icd_guess3, icd_guess4)) %>% 
  bind_rows(studies_guesses_only)

## Account for recent changes to some Ocodes
## (see data issues and changes Rmd 
studies_all <- studies_all %>% 
  mutate(cause = case_when(cause == "O71" ~ "O719",
                           cause == "O99" ~ "Group 7",
                           TRUE ~ cause)) %>% 
  ## Assign main causes
  left_join(icd %>% select(cause,analysis_group),
            by = "cause") 

# there are 4 observations with non-standard cod: "R", "OXX", "Other" and "Other causes"
# for now "Other" and "Other causes" can be manually places in the "other" analysis group.
# "R" and "OXX" get assigned an NA (regardless, none of these make it to the analysis stage)
studies_all <- studies_all %>% 
  mutate(analysis_group = case_when(
    cause == "Other" ~ "other",
    cause == "Other causes" ~ "other",
    TRUE ~ analysis_group
  ))

##
## ---------------------- Clean up nested intervals ------------------------
## (i.e. cases where counts for each cause are reported for 
##  the full study range AND some or all individual study years)

any_nested <- function(intervals) {
  interval_list <- as.list(intervals)
  
  for (i in 1:length(intervals)) {
    flag_nest <- plyr::laply(interval_list,
                             function(x) { (x %within% intervals[i]) & (x != intervals[i])}
    ) %>%
      any()
    
    if (flag_nest) {
      return(TRUE)
    }
  }
  return(FALSE)
}

studies_with_overlaps <- overlaps %>% 
  filter(any_overlaps) %>% distinct(idno_iso, .keep_all = T)

problematic_studies <- c("COC-1563-GMB",
                         "COC-2547-PAK",  
                         "COC-26734-RWA",
                         "COC-28991-IND",
                         "COC-36570-IND",
                         "COC-52934-IND",
                         "COC-56010-MWI" )


studies_all <-  studies_all %>% 
  arrange(idno_iso, 
          cause, 
          start_date,
          end_date
  ) %>% 
  mutate(make_int = interval(start_date, end_date)) %>% 
  group_by(idno_iso, cause) %>% 
  mutate(nest = any_nested(make_int)) %>%
  mutate(find_longest = ifelse(
    (idno_iso %in% studies_with_overlaps$idno_iso & n()>1 & !(idno_iso %in% problematic_studies) ),
    max(duration_days), 
    NA
  )) %>% 
  mutate(find_shortest = ifelse(
    (idno_iso %in% studies_with_overlaps$idno_iso & n()>1 & !(idno_iso %in% problematic_studies)), 
    min(duration_days), 
    NA
  )) %>% 
  mutate(flag_interval = 
           (find_longest == duration_days) & (find_longest != find_shortest)
  ) %>% 
  ungroup()

## For the non-problematic studies remove the longest interval  
## which ends up being just a duplicate of available individual ones
## (non-problematic: COC-14008-GHA and COC-14166-IND)

studies_all <- studies_all %>% 
  mutate(remove = ifelse(
    (nest & flag_interval), T, NA
  )) %>% 
  filter(is.na(remove))


## NOTE (17Jul20-MP): For COC-14166-IND the aggregated intervals are duplicates of the
## disaggregated ones, except for one observation (4 deaths of cause O994 (IND) for 01-01-2011 to 01-01-2013)
## This causes problems when fed to the model (duplicate country-year observation because of different end/start date)
## Change start date for just this observation to 2012 assuming the disaggregated year count would have been the same 
## and is just missing 
## (use if_else to preserve date type)
studies_all <- studies_all %>% 
  mutate(
    start_date = as_datetime(
      if_else(
        idno_iso == "COC-14166-IND" & end_date==ymd("2013-01-01") & start_date==ymd("2011-01-01"),
        ymd(start_date) + years(1), 
        ymd(start_date))),
    duration_days = if_else(
      idno_iso == "COC-14166-IND",
      end_date - start_date,
      duration_days
    )) 

## Handle COC-1563-GMB by hand: remove shorter nested period;
## 2011-2015 and 2007-2015 are the only two periods appearing 
## for every cause available in this study
studies_all <- studies_all %>% 
  filter(!(idno_iso == "COC-1563-GMB" & start_date == as.Date("2011-01-01")))

## Handle Malawi study COC-56010-MWI by hand; 
## Remove the three month period (march to june)
studies_all <- studies_all %>% 
  filter(!(idno_iso == "COC-56010-MWI" & end_date == as.Date("2016-06-01")))


## For the following studies the longest/aggregated time period observations 
## can be safely taken out: COC-26734-RWA, COC-36570-IND, COC-52934-IND
## Counts in available disaggregated periods add up (give or take)
## to the count in the aggregated interval

studies_all <- studies_all %>% 
  group_by(idno_iso, cause) %>% 
  mutate(find_longest = ifelse(
    (idno_iso %in% c("COC-26734-RWA",
                     "COC-36570-IND",
                     "COC-52934-IND")& n()>1),
    max(duration_days), 
    NA
  )) %>% 
  mutate(find_shortest = ifelse(
    (idno_iso %in% c("COC-26734-RWA",
                     "COC-36570-IND",
                     "COC-52934-IND")& n()>1), 
    min(duration_days), 
    NA
  )) %>% 
  mutate(remove_interval = ifelse(
    (find_longest == duration_days) & (find_longest != find_shortest), T, NA)
  ) %>% 
  ungroup() %>%
  filter(is.na(remove_interval))


## Studies COC-28991-IND and COC-2547-PAK are a bit more problematic
## Decision for the moment is to remove the observations for a cause 
## where the time interval is aggregated (i.e. full study period)
studies_all <- studies_all %>% 
  group_by(idno_iso, cause) %>%
  mutate(period_type = ifelse(
    idno_iso %in% c("COC-28991-IND","COC-2547-PAK"), ifelse(
      duration_days == max(duration_days) &
        duration_days != min(duration_days) &
        # add a fix for leap years (366 days)
        #!leap_year(start_date),
        max(duration_days) - min(duration_days) > 1,
      "aggregated", "disaggregated"
    ), NA
  )) %>% 
  ungroup()

studies_all <- studies_all %>% 
  filter((is.na(period_type) | period_type == "disaggregated")) 

## Fix the iso in idno_iso == COC-58954-NGA
## This should be a Niger study (sub NGA with NER)
studies_all <- studies_all %>% 
  mutate(idno_iso = ifelse(
    idno_iso == "COC-58954-NGA",
    str_replace(idno_iso,"NGA", "NER"),
    idno_iso))  

## For HIV assign NA not IND
studies_all <- studies_all %>% 
  mutate(analysis_group = ifelse(cause == "O987", NA, analysis_group))


# A final check for overlapping date ranges
overlaps_final <- studies_all %>%
  mutate(date_interval = interval(start_date, (end_date - as.difftime(1, units = "secs")))) %>%
  # Flag studies that have overlapping date intervals
  group_by(idno_iso, cod) %>%
  summarise(any_overlaps = any_overlaps(date_interval)) 

if (sum(overlaps_final$any_overlaps) != 0) {
  stop("Unchecked overlapping date ranges!")
}

## Adjust the counts for COC-55348-TZA manually (see data issues rmd)
## Add 215 ABO deaths and remove 10 deaths from HEM and 205 from SEP.

# Remove HEM and SEP death
studies_all <- studies_all %>% 
  mutate(num = ifelse(idno_iso == "COC-55348-TZA" & cause == "Hemorrhage NOS",
                                      num - 10, num),
         num = ifelse(idno_iso == "COC-55348-TZA" & cause == "O85",
                                      num-205, num))
## Make a row with the ABO deaths (very hacky)
tza_abo_obs <- studies_all %>% 
  filter(idno_iso == "COC-55348-TZA") %>%
  slice(1) %>% 
  mutate(cause = "Group 1",
         analysis_group = "ABO", 
         num = 215
  )

studies_all <- studies_all %>% 
  bind_rows(tza_abo_obs) %>% 
  arrange(iso, idno_iso)

studies_final <- studies_all %>% 
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
    num
  ) %>% 
  rename(
    #iso = country_iso,
    study_id = idno_iso,
    deaths_age_15_to_49 = num
  )

write_csv(studies_final, "data/current/clean/studies_data_clean.csv")    

