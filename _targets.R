library(targets)
library(tarchetypes)
library(tidyverse)
library(rlang)

source("R/spectral_helper.R")
source("R/targets_helper.R")
source("R/brms_helper.R")

# brms will save models here to avoid recompiling
options(cmdstanr_write_stan_file_dir="models/cmdstanr_models")
# run stan with 4 parallel chains
# this option is not invoked when running targets with crew in parallel
options(mc.cores=4)

library(crew)
tar_option_set(
  #controller = crew_controller_local(workers = 4)
)

tar_option_set(packages = c(
  "tidyverse", # data processing
  "cmdstanr", # bayesian modeling
  "scales", # variable transformations
  "mdatools", # PLSR
  #"glmnet", # ridge regression
  "tidybayes",
  "resemble",
  #"rstanarm",
  "tidymodels",
  #"rwavelet",
  "brms",
  "rlang"#,
  #"sf",
  #"soilDB"
  ),
  workspace_on_error=FALSE)

sites = read_csv("data/hudson_sites.csv")
HUDSON_SITES = sites$site_id 
HUDSON_FARM_IDS = sites$farm_id %>% unique

hudson_targets = list(
  # spectroscopy scans
  tar_target(hudson_moist_file,
             "data/HudsonEDF_459set_FieldMoist_NSscans_avg_meta.csv",
             format="file"),
  tar_target(
    hudson_NIR,
    read_csv(hudson_moist_file) %>%
      select(., soil_id=sample_id, sort_numeric_colnames(get_numeric_colnames(.)))
  ),
  
  # soil carbon concentration data
  tar_target(hudson_lab_file,
             "data/HudsonEDF_459set_Cper_BD.csv",
             format="file"),
  tar_target(
    hudson_lab,
    read_csv(hudson_lab_file) %>%
      select(soil_id=sample_id, eoc_tot_c=Cper)
  ),
  
  # soil sample level data
  tar_target(
    hudson_soils,
    read_csv(hudson_moist_file) %>%
      select(sample_id) %>%
      mutate(layer=str_sub(sample_id, -1),
             site_id = str_split_i(sample_id, "\\.", 2),
             location_id = str_sub(sample_id, 1, -4)) %>%
      mutate(top=recode(layer, `1`=0, `2`=15, `3`=30),
             bottom=recode(layer, `1`=15, `2`=30,`3`=60)) %>%
      select(-layer) %>%
      rename(soil_id=sample_id)
  ),
  
  # site level data
  tar_target(hudson_sites_file,
             "data/hudson_sites.csv",
             format="file"),
  tar_target(
    hudson_sites,
    read_csv(hudson_sites_file)
  )
)

spectra_targets = list(
  tar_target(hudson_NIR_snv,
             get_NIR_snv(hudson_NIR)),
  tar_target(hudson_NIR_snv_64,
             interpolate_spectra(hudson_NIR_snv, length.out=64)),
  # very coarse spectrum used for quickly testing models
  tar_target(hudson_NIR_snv_16,
             interpolate_spectra(hudson_NIR_snv, length.out=16))
)
  
# define training and test splits
split_targets = list( 
  tar_map(
    values=expand_grid(
      test_farm=HUDSON_FARM_IDS,
      n_train = c(0, 3),
      k=1:3) %>%
      filter(k == 1 | n_train > 0),
    names=c(test_farm, p, n_train, k),
    tar_target(hudson_lfarmo_split,
               hudson_soils %>%
                 select(location_id, site_id) %>%
                 distinct %>%
                 inner_join(hudson_sites) %>%
                 group_by(site_id) %>%
                 # a hacky way to select n_train samples from all the sites in the test_farm
                 # as well as all the other farms
                 mutate(t = ifelse(farm_id==test_farm, n_train, n())) %>%
                 mutate(train = row_number() %in% sample(n(), first(t))) %>%
                 ungroup %>%
                 expand_grid(layer=c(1,2,3)) %>%
                 mutate(soil_id = str_glue("{location_id}.S{layer}")) %>%
                 select(-t)
    ),
    tar_target(hudson_lfarmo_train,
               hudson_lfarmo_split %>%
                 filter(train) %>%
                 select(-train, -site_id)),
    tar_target(hudson_lfarmo_test,
               hudson_lfarmo_split %>%
                 filter(!train) %>%
                 select(-train, -site_id))
  )
  
)

# define training and test sets
cv_values = bind_rows(
  expand_grid(site = HUDSON_FARM_IDS %>% head(1),
              n_train=c(0, 3) %>% head(1),
              k=1:1) %>%
    filter(k == 1 | n_train > 0) %>%
    rowwise() %>%
    transmute(train_sites = str_glue("hudson_lfarmo_train_{site}_{n_train}_{k}"),
              test_sites = str_glue("hudson_lfarmo_test_{site}_{n_train}_{k}"),
              args=list(list(group=str_glue("hudson_lfarmo"), cluster=site, n=n_train, k=k))),

) %>%
  expand_grid(NIR_transform=c(
    "snv",
    "snv_16"
  )) %>%
  mutate(train_sites=syms(train_sites),
         test_sites=syms(test_sites),
         NIR_data=syms(str_glue("hudson_NIR_{NIR_transform}")),
         NIR_data_test=syms(str_glue("hudson_NIR_{NIR_transform}"))
         )

cv_targets = tar_map(
  cv_values,
  names=-c(NIR_data, NIR_data_test, args),
  tar_target(
    std_rec,
    recipe(train_sites %>% 
             # get rid of any extra cols in train_sites, otherwise the recipe will require them...
             select(all_of(intersect(colnames(train_sites), colnames(NIR_data)))) %>% 
             inner_join(NIR_data)) %>%
      update_role(get_numeric_colnames(NIR_data), new_role="predictor") %>%
      step_center(all_predictors()) %>% 
      step_scale(all_predictors()) %>%
      prep(training=train_sites %>% inner_join(NIR_data))
  ),
  tar_target(
    train, 
    train_sites %>% 
      inner_join(hudson_soils) %>% 
      inner_join(hudson_lab) %>%
      inner_join(bake(std_rec, NIR_data))),
  tar_target(
    test, 
    test_sites %>% 
      inner_join(hudson_soils) %>% 
      inner_join(bake(std_rec, NIR_data_test))),
  
  # we actually mostly don't want to standardize...
  tar_target(
    train_nostd, 
    train_sites %>% 
      inner_join(hudson_soils) %>% 
      inner_join(hudson_lab) %>%
      inner_join(NIR_data)),
  tar_target(
    test_nostd, 
    test_sites %>% 
      inner_join(hudson_soils) %>% 
      inner_join(NIR_data_test)),
  
  tar_target(
    truth, 
    test %>%
      select(soil_id) %>%
      inner_join(hudson_lab)),
  tar_target(
    truth_sitelayer,
    test %>%
      select(soil_id, site_id, top) %>%
      inner_join(hudson_lab) %>%
      select(-soil_id) %>%
      group_by(site_id, top) %>%
      summarize(across(everything(), mean))
  )
)

model_values = bind_rows(
  tibble_row(model="plsr", y=list(list(eoc_tot_c="log1p"))),
  tibble_row(model="plsr_stratified", y=list(list(eoc_tot_c="log1p"))),
  tibble_row(model="qrf", y=list(list(eoc_tot_c="log1p"))),
  tibble_row(model="mbl", y=list(list(eoc_tot_c="log1p"))),
  #tibble_row(model="cubist", y=list(list(eoc_tot_c="log1p"))),
  tibble_row(model="cubist_conformal", y=list(list(eoc_tot_c="log1p"))),

  tibble_row(model="lmer_SSURGO",
             y=list(list(eoc_tot_c="log1p")),
             fit_args=list(list(group="farm_id + site_id*top",
                                SSURGO=c("clay", "sand", "silt", "pH")))),
  # 
  # Supplementary models
  # Cubist with SSURGO
  tibble_row(model="cubist_SSURGO", y=list(list(eoc_tot_c="log1p")),
           fit_args = list(list(extraFactors=c("farm_id", "top", "site_id"),
                                SSURGO=c("silt", "sand", "clay", "pH")))),
  # 
  # # lmer no SSURGO
  # tibble_row(model="lmer_SSURGO_mean2.brms", 
  #            y=list(list(eoc_tot_c="log1p")),
  #            fit_args=list(list(group="farm_id + site_id*top",
  #                               SSURGO=c(
  #                               ) ))),
  # # lmer no varying slopes
  # tibble_row(model="lmer_SSURGO_mean2.brms",
  #            y=list(list(eoc_tot_c="log1p")),
  #            fit_args=list(list(group="farm_id + site_id*top",
  #                               SSURGO=c("clay", "sand", "silt", "pH"),
  #                               varying_slopes=FALSE))),
  # 
  # 
  # 
) %>%
  # make sure there is a fit_args list in the table
  bind_rows(tibble_row(fit_args=list(list()))) %>%
  head(-1) %>% # then remove the empty row
  mutate(fit_function = syms(paste0("fit_", model))) %>%
  mutate(predict_function=syms(paste0("predict_", model))) %>%
  rowwise() %>%
  mutate(
    y_transforms = list(to_transforms(y)), # replace "log1p" with transform_log1p
    y_str = named_list_to_str2(y),
    fit_args_str = named_list_to_str2(fit_args) %>%
           str_replace("^_$", "none") %>%
           str_replace("\\)$", "")) %>% # if there are no args call it none
  ungroup

fit_values = tar_add_steps_to_values(cv_values, cv_targets) %>%
  expand_grid(model_values) %>%
  mutate(args_str = named_list_to_str2(args)) %>%
  # don't fit correlation matrix models with more than 64 wavelengths
  filter(xor(NIR_transform == "snv_16", !grepl("lmer", model) ))


fit_targets = tar_map(
  fit_values,
  names=c(train, model, fit_args_str, y_str),
  tar_target(
    fit,
    do.call(fit_function, c(list(
      y=names(y),
      train=train %>%
        sample_frac(0.1) %>%
        transform_cols(y_transforms),
      test=test),
      fit_args) ) ),
  tar_target(
    predict,
    predict_function(fit,
            newdata=test,
            y=names(y)) %>%
      inv_transform_preds(y_transforms) %>%
      ensure_row_in_preds
    ),
  tar_target(
    predict_soil,
    predict %>%
      group_by(.row) %>%
      summarize_preds %>%
      bind_cols(truth) %>%
      mutate(provenance)
  ),
  tar_target(
    provenance,
    tibble(train=as_name(quo(train)),
           model=model,
           fit_args=fit_args_str,
           NIR_transform=NIR_transform,
           y=list(y),
           as_tibble(args))
  ),
  tar_target(
    predict_sitelayer,
    predict %>%
      inner_join(test %>% select(site_id, top) %>% mutate(.row=row_number())) %>%
      group_by(.draw, site_id, top) %>%
      aggregate_preds %>%
      group_by(site_id, top) %>%
      summarize_preds %>%
      inner_join(truth_sitelayer, by=c("site_id", "top")) %>%
      mutate(provenance)
  )
)

c(hudson_targets,
  spectra_targets,
  split_targets,
  cv_targets,
  fit_targets,
  tar_combine(
    predict_soil_combined,
    fit_targets$predict_soil,
    command=bind_rows(!!!.x, .id="name")
  ),
  tar_combine(
    predict_sitlayer_combined,
    fit_targets$predict_sitelayer,
    command=bind_rows(!!!.x, .id="name")
  )
)
