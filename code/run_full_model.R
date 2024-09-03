
library(tidyverse)
library(cmdstanr)
library(here)
library(caret)
library(posterior)

options(mc.cores = 4)

load(here("data/current/clean/objects_for_modelling.RData"))

studies_durations <- read_csv(here::here("data/current/clean/studies_data_clean.csv")) %>%
  distinct(study_id, year, duration_days) %>%
  mutate(duration_years = pmax(duration_days/365, 1))

grey_durations <- read_csv(here::here("data/current/clean/grey_data_clean.csv")) %>%
  distinct(study_id, year, duration_days) %>%
  mutate(duration_years = pmax(duration_days/365, 1))

durations <- bind_rows(studies_durations, grey_durations) %>%
  select(-duration_days)

# Turn year into a factor with specified levels
d <- d %>%
  left_join(durations) %>%
  mutate(year = factor(year, levels = years)) %>%
  replace_na(list(duration_years = 1)) %>%
  filter(!source == "studies_below_ADM1") %>% 
  rowwise() %>% 
  mutate_at(causes, ~ifelse(
    TOT/duration_years > 5, 
    round(./duration_years),
    round(.))) %>%
  mutate(TOT = sum(ABO, DIR, EMB, HEM, HYP, IND, SEP, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(quality = case_when(
    source == "grey_national" ~ 0,
    source == "VR" ~ 1 - u_new,
    TRUE ~ 1
  ))

d_complete <- d_complete %>%
  mutate(year = factor(year, levels = years))

# Recreate some other objects: -----------------------
## * Observation dataframe as matrix to feed into Stan ----
d_mat <- d %>%
  select(all_of(causes)) %>%
  as.matrix()

# save locations of true zeroes
true_zeroes <- d_mat == 0
true_zeroes[is.na(true_zeroes)] <- FALSE

# Put zeroes in place of NAs in Stan's matrix
d_mat[is.na(d_mat)] <- 0

# zero_mat == 1 where observation is > 0, or is a true zero
# dictates which values contribute to likelihood
zero_mat <- 1*(d_mat > 0 | true_zeroes)

# Change zero_mat
# If an observation is type 1 and that source has ever recorded a death from that cause, then treat it as a true zero
cause_ever_present <- d %>% 
  group_by(study_id, iso) %>%
  mutate_at(tidyselect::all_of(causes), function(x) {any(x > 0, na.rm = TRUE)}) %>%
  ungroup() %>%
  select(all_of(causes))

# Type 1 observation, and analysis group has previously been observed
condition1 <- diag(d$type_new %in% 1:2) %*% as.matrix(cause_ever_present)

# Type 1 observation, predicted deaths < 7
condition2 <- diag(d$type_new %in% 1:2) %*% (matrix(d$matdeaths < 7, nrow = nrow(d), ncol = 7))

new_zero_mat <- 1*(condition1 | condition2 | zero_mat)

# Make model matrix with tools from caret package
model_setup <- dummyVars(~1 + region_modelling, data = d)

X <- model_setup %>% predict(d)
X_complete <- model_setup %>% predict(d_complete)

X_subDIR <- model_setup %>% predict(d_DIR) 
X_subHEM <- model_setup %>% predict(d_HEM) 
X_subSEP <- model_setup %>% predict(d_SEP)

# Check that columns match
if (!identical(colnames(X), colnames(X_complete))) {
  stop("Column names of X and X_complete do not match")
}

# Make country labels for each observation
c_labels <- plyr::mapvalues(
  d$iso,
  from = country_labels$name,
  to = country_labels$number
) %>% as.numeric()

# Make country labels for predictions
c_labels_complete <- plyr::mapvalues(
  d_complete$iso,
  from = country_labels$name,
  to = country_labels$number
) %>% as.numeric()

# Make region labels for each country
r_labels <- plyr::mapvalues(
  country_labels$region_modelling,
  from = region_modelling_labels$name,
  to = region_modelling_labels$number
) %>% as.numeric()

# Make year labels for each observation
year_labels <- plyr::mapvalues(
  d$year,
  from = years,
  to = 1:nyears
) %>% as.numeric()

# Make year labels for prediction
year_labels_complete <- plyr::mapvalues(
  d_complete$year,
  from = years,
  to = 1:nyears
) %>% as.numeric()

# Boolean for when an observation is not Type 1
not_type_1 <- d %>% 
  mutate(x = !(as.numeric(type_new) == 1)) %>%
  pull(x)

# Make country labels for DIR subcause observations
c_labels_subDIR <- plyr::mapvalues(
  d_DIR$iso,
  from = country_labels$name,
  to = country_labels$number
) %>% as.numeric()

# Make country labels for HEM subcause observations
c_labels_subHEM <- plyr::mapvalues(
  d_HEM$iso,
  from = country_labels$name,
  to = country_labels$number
)  %>% as.numeric()

# Make country labels for SEP subcause observations
c_labels_subSEP <- plyr::mapvalues(
  d_SEP$iso,
  from = country_labels$name,
  to = country_labels$number
) %>% as.numeric()


# Put data in a list with appropriate names
data_list <- list(
  # Data objects for main model
  N = nrow(d),
  ncat = ncause,
  Y = d_mat,
  tot = d$TOT,
  K = ncol(X),
  X = X,
  zero_mat = new_zero_mat,
  nyear = nyears,
  year = year_labels,
  nregion = nrow(region_modelling_labels),
  region = r_labels,
  ncountry = ncountry,
  c = c_labels,
  not_type_1 = not_type_1,
  quality = d$quality,
  delta_scale = 1,
  # specific to DIR subcause model
  Y_subDIR = d_mat_DIR, 
  N_subDIR = nrow(d_mat_DIR),
  ncat_subDIR = nsubcause_DIR,
  X_subDIR = X_subDIR,
  c_subDIR = c_labels_subDIR,
  zero_mat_subDIR = zero_mat_DIR,
  index_DIR = which(colnames(d_mat) == "DIR"),
  # specific to HEM subcause model
  Y_subHEM = d_mat_HEM, 
  N_subHEM = nrow(d_mat_HEM),
  ncat_subHEM = nsubcause_HEM,
  X_subHEM = X_subHEM,
  c_subHEM = c_labels_subHEM,
  zero_mat_subHEM = zero_mat_HEM,
  index_HEM = which(colnames(d_mat) == "HEM"),
  # specific to SEP subcause model
  Y_subSEP = d_mat_SEP, 
  N_subSEP = nrow(d_mat_SEP),
  ncat_subSEP = nsubcause_SEP,
  X_subSEP = X_subSEP,
  c_subSEP = c_labels_subSEP,
  zero_mat_subSEP = zero_mat_SEP,
  index_SEP = which(colnames(d_mat) == "SEP"),
  # Prediction
  N_complete = nrow(X_complete),
  X_complete = X_complete,
  c_complete = c_labels_complete,
  year_complete = year_labels_complete,
  matdeaths = d_complete$matdeaths,
  region_sdg_mat = region_sdg_mat,
  region_sdg_year_mat = region_sdg_year_mat,
  nregion_sdg = nrow(region_sdg_labels),
  country_mat = country_mat,
  # Coverage calcluations
  w = d_complete$max_mat_coverage, 
  # HIV/AIDS
  aidsdeaths = d_complete$aidsdeaths
)

attributes(data_list)$d_modified <- d

file <- here("code/models/model81_full_continuous-var_test.stan")
mod <- cmdstan_model(file)

time_to_run <- system.time(
  fit <- mod$sample(
    data = data_list, 
    seed = 2023,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 3000,
    iter_sampling = 2000,
    refresh = 100,
    max_treedepth = 13,
    adapt_delta = 0.99,
    output_dir = here("analysis")
  )
)

variables_to_save <- c(
  # Parameters
  "alpha",
  "beta", 
  "sd_beta",
  "sd_delta",
  # Diagnostics
  "log_lik",
  "Y_rep",
  # Main model
  "Sigma",
  "p_country",
  "p_complete",
  "p_region_sdg",
  "p_global",
  "count_country",
  "count_region_sdg",
  "count_global",
  # DIR 
  "Sigma_subDIR",
  "p_subDIR_complete",
  "p_subDIR_region_sdg",
  "p_subDIR_global",
  "count_subDIR_region_sdg",
  "count_subDIR_global",
  # HEM
  "Sigma_subHEM",
  "p_subHEM_complete",
  "p_subHEM_region_sdg",
  "p_subHEM_global",
  "count_subHEM_region_sdg",
  "count_subHEM_global",
  # SEP
  "Sigma_subSEP",
  "p_subSEP_complete",
  "p_subSEP_region_sdg",
  "p_subSEP_global",
  "count_subSEP_region_sdg",
  "count_subSEP_global"
)

q_probs <- c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975)

summary_table <- fit$summary(
  variables = variables_to_save,
  mean, sd, rhat, ~quantile2(.x, probs = q_probs), ess_bulk, ess_tail                   
)

to_save <- list(
  model = "full_model",
  data = data_list,
  time = time_to_run,
  fit = fit,
  variables = variables_to_save,
  summary = summary_table
)

write_rds(to_save, here("analysis/full_model.rds"))

