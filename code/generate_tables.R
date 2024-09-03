library(tidyverse)
library(here)

load(here("data/current/clean/objects_for_modelling.RData"))
result <- read_rds(here("analysis/full_model.rds"))

# Global main cause -----
# Write global results
result_global <- result$summary |>
  filter(str_detect(variable, "p_global")) |>
  mutate(cause = colnames(result$data$Y)) |>
  select(cause, q2.5, q10, mean, median = q50, q90, q97.5)

result_count_global <- result$summary |>
  filter(str_detect(variable, "count_global")) |>
  mutate(cause = colnames(result$data$Y)) |>
  select(cause, q2.5, q10, mean, median = q50, q90, q97.5)

bind_rows(
  result_global |> mutate(type = "prop"),
  result_count_global |> mutate(type = "count"),
  ) |> 
  arrange(cause) |> 
  select(cause, type, everything()) |>
  write_csv(here("deliverables/global_update-2024.csv"))

# SDG Region main cause ------
# Write SDG region results
result_sdg <- result$summary |>
  filter(str_detect(variable, "p_region_sdg")) |>
  separate(variable, into = c("row", "cause_number"), sep = ",") |>
  group_by(cause_number) |>
  mutate(
    region_sdg = rownames(result$data$region_sdg_mat)
  ) |>
  group_by(region_sdg) |>
  mutate(cause = colnames(result$data$Y)) |>
  select(region_sdg, cause, q2.5, q10, mean, median = q50, q90, q97.5)

result_count_sdg <- result$summary |>
  filter(str_detect(variable, "count_region_sdg")) |>
  separate(variable, into = c("row", "cause_number"), sep = ",") |>
  group_by(cause_number) |>
  mutate(
    region_sdg = rownames(result$data$region_sdg_mat)
  ) |>
  group_by(region_sdg) |>
  mutate(cause = colnames(result$data$Y)) |>
  select(region_sdg, cause, q2.5, q10, mean, median = q50, q90, q97.5)

bind_rows(
  result_sdg |> mutate(type = "prop"),
  result_count_sdg |> mutate(type = "count")
  ) |>
  arrange(region_sdg, cause) |> 
  select(region_sdg, cause, type, everything()) |>
  write_csv(here("deliverables/sdg_update-2024.csv"))
  
# Country-level main cause ------
# Write country-level results
result_estimates <- result$summary %>%
  filter(str_detect(variable, "p_country"))|>
  separate(variable, into = c("country", "cause_number"), sep = ",") |>
  group_by(cause_number) |>
  mutate(country = rownames(result$data$country_mat) |> str_remove("iso")) |>
  group_by(country) |>
  mutate(cause = colnames(result$data$Y)) |>
  ungroup() |>
  select(country, cause, q2.5, q10, mean, median = q50, q90, q97.5)

result_count_estimates <- result$summary %>%
  filter(str_detect(variable, "count_country"))|>
  separate(variable, into = c("country", "cause_number"), sep = ",") |>
  group_by(cause_number) |>
  mutate(country = rownames(result$data$country_mat) |> str_remove("iso")) |>
  group_by(country) |>
  mutate(cause = colnames(result$data$Y)) |>
  ungroup() |>
  select(country, cause, q2.5, q10, mean, median = q50, q90, q97.5)

bind_rows(
  result_estimates |> mutate(type = "prop"),
  result_count_estimates |> mutate(type = "count")
  ) |>
  arrange(country, cause) |> 
  select(country, cause, type, everything()) |>
  write_csv(here("deliverables/country_update-2024.csv"))


# DIR subcauses --------
# Global subDIR proportions
result_DIR_prop_global <- result$summary |>
  filter(str_detect(variable, "p_subDIR_global")) |>
  mutate(subcause = str_extract(variable, "[0-9]+") %>% as.numeric()) |>
  mutate(subcause = as_factor(plyr::mapvalues(
    subcause,
    subcauses_DIR_labels$number,
    subcauses_DIR_labels$name
  ))) |>
  mutate(region = "Global") |>
  select(region, subcause, q2.5, q10, mean, median = q50, q90, q97.5) 

# SDG subDIR proportions
result_DIR_prop_sdg <- result$summary |>
  filter(str_detect(variable, "p_subDIR_region_sdg")) |>
  separate(variable, sep = ",", into = c("region", "subcause")) %>%
  mutate(region = str_extract(region, "[0-9]+") |> as.numeric()) |>
  mutate(region = plyr::mapvalues(
    str_extract(region, "[0-9]+"), 
    region_sdg_labels$number, 
    region_sdg_labels$name)) %>%
  mutate(subcause = str_extract(subcause, "[0-9]+") |>  as.numeric()) |>
  mutate(subcause = as_factor(plyr::mapvalues(
    subcause,
    subcauses_DIR_labels$number,
    subcauses_DIR_labels$name
  ))) |>
  select(region, subcause, q2.5, q10, mean, median = q50, q90, q97.5) 

# Global subDIR counts
result_DIR_count_global <- result$summary |>
  filter(str_detect(variable, "count_subDIR_global")) |>
  mutate(subcause = str_extract(variable, "[0-9]+") %>% as.numeric()) |>
  mutate(subcause = as_factor(plyr::mapvalues(
    subcause,
    subcauses_DIR_labels$number,
    subcauses_DIR_labels$name
  ))) |>
  mutate(region = "Global") |>
  select(region, subcause, q2.5, q10, mean, median = q50, q90, q97.5) 

# SDG subDIR counts 
result_DIR_count_sdg <- result$summary |>
  filter(str_detect(variable, "count_subDIR_region_sdg")) |>
  separate(variable, sep = ",", into = c("region", "subcause")) %>%
  mutate(region = str_extract(region, "[0-9]+") |> as.numeric()) |>
  mutate(region = plyr::mapvalues(
    str_extract(region, "[0-9]+"), 
    region_sdg_labels$number, 
    region_sdg_labels$name)) %>%
  mutate(subcause = str_extract(subcause, "[0-9]+") |>  as.numeric()) |>
  mutate(subcause = as_factor(plyr::mapvalues(
    subcause,
    subcauses_DIR_labels$number,
    subcauses_DIR_labels$name
  ))) |>
  select(region, subcause, q2.5, q10, mean, median = q50, q90, q97.5) 


bind_rows(
  result_DIR_prop_global |> mutate(type = "prop"), 
  result_DIR_prop_sdg |> mutate(type = "prop"),
  result_DIR_count_global |> mutate(type = "count"),
  result_DIR_count_sdg |> mutate(type = "count")
) |>
  select(region, subcause, type, everything())|>
  write_csv(here("deliverables/subDIR_update-2024.csv"))


# HEM subcauses --------
# Global subHEM proportions
result_HEM_prop_global <- result$summary |>
  filter(str_detect(variable, "p_subHEM_global")) |>
  mutate(subcause = str_extract(variable, "[0-9]+") %>% as.numeric()) |>
  mutate(subcause = as_factor(plyr::mapvalues(
    subcause,
    subcauses_HEM_labels$number,
    subcauses_HEM_labels$name
  ))) |>
  mutate(region = "Global") |>
  select(region, subcause, q2.5, q10, mean, median = q50, q90, q97.5) 

# SDG subHEM proportions
result_HEM_prop_sdg <- result$summary |>
  filter(str_detect(variable, "p_subHEM_region_sdg")) |>
  separate(variable, sep = ",", into = c("region", "subcause")) %>%
  mutate(region = str_extract(region, "[0-9]+") |> as.numeric()) |>
  mutate(region = plyr::mapvalues(
    str_extract(region, "[0-9]+"), 
    region_sdg_labels$number, 
    region_sdg_labels$name)) %>%
  mutate(subcause = str_extract(subcause, "[0-9]+") |>  as.numeric()) |>
  mutate(subcause = as_factor(plyr::mapvalues(
    subcause,
    subcauses_HEM_labels$number,
    subcauses_HEM_labels$name
  ))) |>
  select(region, subcause, q2.5, q10, mean, median = q50, q90, q97.5) 

# Global subHEM counts
result_HEM_count_global <- result$summary |>
  filter(str_detect(variable, "count_subHEM_global")) |>
  mutate(subcause = str_extract(variable, "[0-9]+") %>% as.numeric()) |>
  mutate(subcause = as_factor(plyr::mapvalues(
    subcause,
    subcauses_HEM_labels$number,
    subcauses_HEM_labels$name
  ))) |>
  mutate(region = "Global") |>
  select(region, subcause, q2.5, q10, mean, median = q50, q90, q97.5) 

# SDG subHEM counts 
result_HEM_count_sdg <- result$summary |>
  filter(str_detect(variable, "count_subHEM_region_sdg")) |>
  separate(variable, sep = ",", into = c("region", "subcause")) %>%
  mutate(region = str_extract(region, "[0-9]+") |> as.numeric()) |>
  mutate(region = plyr::mapvalues(
    str_extract(region, "[0-9]+"), 
    region_sdg_labels$number, 
    region_sdg_labels$name)) %>%
  mutate(subcause = str_extract(subcause, "[0-9]+") |>  as.numeric()) |>
  mutate(subcause = as_factor(plyr::mapvalues(
    subcause,
    subcauses_HEM_labels$number,
    subcauses_HEM_labels$name
  ))) |>
  select(region, subcause, q2.5, q10, mean, median = q50, q90, q97.5) 

bind_rows(
  result_HEM_prop_global |> mutate(type = "prop"), 
  result_HEM_prop_sdg |> mutate(type = "prop"),
  result_HEM_count_global |> mutate(type = "count"),
  result_HEM_count_sdg |> mutate(type = "count")
) |>
  select(region, subcause, type, everything()) |>
  write_csv(here("deliverables/subHEM_update-2024.csv"))

# SEP subcauses --------
# Global subSEP proportions
result_SEP_prop_global <- result$summary |>
  filter(str_detect(variable, "p_subSEP_global")) |>
  mutate(subcause = str_extract(variable, "[0-9]+") %>% as.numeric()) |>
  mutate(subcause = as_factor(plyr::mapvalues(
    subcause,
    subcauses_SEP_labels$number,
    subcauses_SEP_labels$name
  ))) |>
  mutate(region = "Global") |>
  select(region, subcause, q2.5, q10, mean, median = q50, q90, q97.5) 

# SDG subSEP proportions
result_SEP_prop_sdg <- result$summary |>
  filter(str_detect(variable, "p_subSEP_region_sdg")) |>
  separate(variable, sep = ",", into = c("region", "subcause")) %>%
  mutate(region = str_extract(region, "[0-9]+") |> as.numeric()) |>
  mutate(region = plyr::mapvalues(
    str_extract(region, "[0-9]+"), 
    region_sdg_labels$number, 
    region_sdg_labels$name)) %>%
  mutate(subcause = str_extract(subcause, "[0-9]+") |>  as.numeric()) |>
  mutate(subcause = as_factor(plyr::mapvalues(
    subcause,
    subcauses_SEP_labels$number,
    subcauses_SEP_labels$name
  ))) |>
  select(region, subcause, q2.5, q10, mean, median = q50, q90, q97.5) 

# Global subSEP counts
result_SEP_count_global <- result$summary |>
  filter(str_detect(variable, "count_subSEP_global")) |>
  mutate(subcause = str_extract(variable, "[0-9]+") %>% as.numeric()) |>
  mutate(subcause = as_factor(plyr::mapvalues(
    subcause,
    subcauses_SEP_labels$number,
    subcauses_SEP_labels$name
  ))) |>
  mutate(region = "Global") |>
  select(region, subcause, q2.5, q10, mean, median = q50, q90, q97.5) 

# SDG subSEP counts 
result_SEP_count_sdg <- result$summary |>
  filter(str_detect(variable, "count_subSEP_region_sdg")) |>
  separate(variable, sep = ",", into = c("region", "subcause")) %>%
  mutate(region = str_extract(region, "[0-9]+") |> as.numeric()) |>
  mutate(region = plyr::mapvalues(
    str_extract(region, "[0-9]+"), 
    region_sdg_labels$number, 
    region_sdg_labels$name)) %>%
  mutate(subcause = str_extract(subcause, "[0-9]+") |>  as.numeric()) |>
  mutate(subcause = as_factor(plyr::mapvalues(
    subcause,
    subcauses_SEP_labels$number,
    subcauses_SEP_labels$name
  ))) |>
  select(region, subcause, q2.5, q10, mean, median = q50, q90, q97.5) 

bind_rows(
  result_SEP_prop_global |> mutate(type = "prop"), 
  result_SEP_prop_sdg |> mutate(type = "prop"),
  result_SEP_count_global |> mutate(type = "count"),
  result_SEP_count_sdg |> mutate(type = "count")
) |>
  select(region, subcause, type, everything()) |>
  write_csv(here("deliverables/subSEP_update-2024.csv"))


