get_numeric_colnames = function(neospectra) {
  c = colnames(neospectra)
  suppressWarnings(c[!is.na(as.numeric(c))])
}

sort_numeric_colnames = function(numeric_colnames) {
  numeric_colnames[order(as.numeric(numeric_colnames))]
}

get_nonnumeric_colnames = function(neospectra) {
  c = colnames(neospectra)
  suppressWarnings(c[is.na(as.numeric(c))])
}

# brms doesn't like numeric colnames, e.g. 1350 so rename to x_1350
prefix_numeric_colnames = function(df, prefix="x_") {
  numeric_colnames = get_numeric_colnames(df)
  
  if(length(numeric_colnames > 0)) {
    df %>%
      rename_with(~paste0(prefix, .x), all_of(numeric_colnames))
  } else {
    df
  }
}

# select columns with numeric names and convert to matrix
get_NIR_matrix = function(df) {
  df %>%
    select(any_of(get_numeric_colnames(df))) %>%
    as.matrix
}
# 
# # select columns that dont have numeric names
# get_NIR_metadata = function(df) {
#   df %>%
#     select(!any_of(get_numeric_colnames(df)))
# }
# 
# get_neospectra_NIR = function(neospectra) {
#   neospectra %>%
#     select(kssl_id, scanner_SerialNo, Lab, SCAN_ID, any_of(sort_numeric_colnames(get_numeric_colnames(neospectra))))
# }

get_NIR_snv = function(neospectra_NIR, center=TRUE, scale=TRUE) {
  NIR_cols = get_numeric_colnames(neospectra_NIR)

  m = neospectra_NIR %>%
    select(all_of(NIR_cols)) %>%
    as.matrix

  m = t(scale(t(m), center=center, scale=scale))

  neospectra_NIR %>%
    select(-all_of(NIR_cols)) %>%
    bind_cols(m)
}

fit_plsr = function(train, y, ...) {
  if(length(y) != 1) {
    stop("plsr needs exactly one y variable")
  }

  NIR_cols = get_numeric_colnames(train)
  x = train[,NIR_cols] %>% as.matrix
  y = train %>% pull(y) %>% as.matrix

  set.seed(123)
  mdatools::pls(x = x, y = y,
      ncomp = 20,
      center = TRUE, scale = TRUE,
      cv = min(nrow(x), 10), lim.type = "ddmoments", cv.scope = 'local')
}

# train plsr only on the texture classes that appear in the test set
fit_plsr_stratified = function(train, y, test, ...) {
  target_texture_classes = test$texture_class %>% unique
  train_stratum = train %>% filter(texture_class %in% target_texture_classes)
  fit_plsr(train_stratum, y)
}

predict_plsr = function(model, newdata, y) {
  ncomp <- model$ncomp.selected
  Xnew = newdata[,get_numeric_colnames(newdata)] %>% as.matrix
  pred <- predict(model, Xnew)
  y_hat <- pred$y.pred[, ncomp, 1]
  
  # Faber and Kowalski (1996) uncertainty quantification
  n <- model$res$cal$y.pred %>% nrow
  conf_level = 0.95
  
  T_cal <- model$res$cal$xdecomp$scores[, 1:ncomp, drop = FALSE]

  T_new <- pred$xdecomp$scores[, 1:ncomp, drop = FALSE]

  T_cal_ssq <- colSums(T_cal^2)
  ho <- rowSums(sweep(T_new^2, 2, T_cal_ssq, "/"))

  df_naive <- ncomp + 1
  rmsec <- model$res$cal$rmse[ncomp]
  sigma <- rmsec * sqrt(n / (n - df_naive))

  s <- sigma * sqrt(1 + 1/n + ho)

  alpha <- 1 - conf_level
  t_val <- qt(1 - alpha/2, n - df_naive)

  lower_ci <- y_hat - (t_val * s)
  upper_ci <- y_hat + (t_val * s)

  # Return data frame
  tibble(
    ".pred_{y}_point" := y_hat,
    ".pred_{y}_q025" := lower_ci,
    ".pred_{y}_q975" := upper_ci
  )
}
predict_plsr_stratified = predict_plsr

fit_mbl = function(train, y, test, ...) {
  if(length(y) != 1) {
    stop("plsr needs exactly one y variable")
  }

  NIR_cols = get_numeric_colnames(train)
  x = train[,NIR_cols] %>% as.matrix
  y = train %>% pull(y) %>% as.matrix
  x_pred = test[,NIR_cols] %>% as.matrix

  my_diss <- "cor"
  my_ks <- seq(nrow(x),21, by = -100)
  ignore_diss <- "none"
  my_waplsr <- fit_wapls(min_ncomp = 3, max_ncomp = 20)
  nnv_val_control <- mbl_control(validation_type = "NNv")

  my_ks <- seq(nrow(x),21, by = -100)
  ignore_diss <- "none"
  my_wapls <- fit_wapls(min_ncomp = 3, max_ncomp = 20, scale = TRUE)
  nnv_val_control <- mbl_control()

  local_ciso <- mbl(
    Xr = x,
    Yr = y,
    Xu = x_pred,
    neighbors=neighbors_k(my_ks),
    diss_usage="predictors",
    spike=which(train$site_id %in% unique(c(test$site_id)))
  )
}

predict_mbl = function(fit, newdata, y) {
  k_opt = fit$validation_results$nearest_neighbor_validation %>%
    arrange(rmse) %>%
    pull(k) %>%
    first

  colname = str_glue("k_{k_opt}")

  p=get_predictions(fit) %>%
    as_tibble %>%
    select(all_of(colname))
  colnames(p) = sprintf(".pred_%s_point", y)

  p
}

fit_cubist = function(train, y, ...) {
  if(length(y) != 1) {
    stop("plsr needs exactly one y variable")
  }

  NIR_cols = get_numeric_colnames(train)
  x = train[,NIR_cols] %>% as.matrix
  y = train %>% pull(y) %>% as.vector

  train_control <- caret::trainControl(
    method = "cv", number = 5, savePredictions = "final")
  tune_grid <- expand.grid(committees = c(1,5,10,15,20),
                           neighbors=0)

  #library(doParallel)
  #cl <- makePSOCKcluster(5) # use 5 cores for 5 fold CV
  #registerDoParallel(cl)

  cubist_model <- caret::train(
    x=x, y=y,
    method = "cubist",
    trControl = train_control,
    #metric="RMSE", maximize=FALSE, # minimize RMSE
    tuneGrid = tune_grid
  )

  #stopCluster(cl)

  cubist_model
}

fit_cubist_conformal = function(train, y, ...) {
  cubist_model = fit_cubist(train, y)

  # get CV residuals from main model
  abs_resid = cubist_model$pred %>%
    arrange(rowIndex) %>%
    mutate(abs_resid = abs(pred - obs)) %>%
    pull(abs_resid)

  # train model to predict abs(resid)
  train = cubist_model$trainingData %>% select(-.outcome)
  conformal_model = Cubist::cubist(
    x=train,
    y=abs_resid,
    committees=cubist_model$bestTune$committees,
    neighbors=cubist_model$bestTune$neighbors)

  # get alpha, i.e. sample predictions from conformal model
  pred_abs_resid = predict(conformal_model, newdata = train)
  alpha_scores = abs_resid/pred_abs_resid

  # save conformal model and alpha_scores as attributes of main model
  cubist_model$conformal_model = conformal_model
  cubist_model$alpha_scores = alpha_scores

  cubist_model
}

predict_cubist_conformal = function(fit, newdata, y) {
  # first get point predictions
  preds = predict_cubist(fit, newdata, y)

  # next get predicted abs resid
  preds$abs_resid_pred = predict(fit$conformal_model,
                                 newdata = newdata[,get_numeric_colnames(newdata)] %>% as.matrix)

  # multiply by appropriate quantile from training set...
  target_coverage = 0.95
  n_calibration <- nrow(fit$trainingData)
  corrected_quantile <- ((n_calibration + 1)*(target_coverage))/n_calibration
  alpha_corrected_95 <- quantile(fit$alpha_scores, corrected_quantile)

  # Final output
  preds %>%
    mutate(".pred_{y}_q025" := get(str_glue(".pred_{y}_point")) - abs_resid_pred*alpha_corrected_95,
           ".pred_{y}_q975" := get(str_glue(".pred_{y}_point")) + abs_resid_pred*alpha_corrected_95) %>%
    select(-abs_resid_pred)
}

fit_cubist_SSURGO = function(train, y, extraFactors, SSURGO, ...) {
  if(length(y) != 1) {
    stop("plsr needs exactly one y variable")
  }
  SSURGO = paste0(c("SSURGO_site_"), SSURGO, "_r")

  train = train %>%
    mutate(across(all_of(extraFactors), factor))
  NIR_cols = get_numeric_colnames(train)
  x = train[,c(NIR_cols, extraFactors, SSURGO)] %>% as.data.frame
  y = train %>% pull(y) %>% as.vector

  train_control <- caret::trainControl(method = "cv", number = 5)
  tune_grid <- expand.grid(committees = c(1,5,10,15,20),
                           neighbors=0)

  cubist_model <- caret::train(
    x=x, y=y,
    method = "cubist",
    trControl = train_control,
    tuneGrid = tune_grid
  )

  cubist_model
}

predict_cubist_SSURGO = function(fit, newdata, y) {
  # get factor variables from training data
  extraFactors = names(select_if(fit$trainingData, is.factor))
  # convert test data to factors with levels from training data
  newdata = newdata %>%
    mutate(across(all_of(extraFactors),
                  ~ factor(., levels=levels(fit$trainingData[, cur_column()]))))
  library(Cubist)
  pred = predict(fit$finalModel,
                 newdata = newdata[,fit$trainingData %>% select(-.outcome) %>% colnames] %>% as.data.frame)

  p = tibble(pred=pred)
  colnames(p) = sprintf(".pred_%s_point", y)

  p
}

predict_cubist = function(fit, newdata, y) {
  library(Cubist)
  pred = predict(fit$finalModel,
                 newdata = newdata[,get_numeric_colnames(newdata)] %>% as.matrix)

  p = tibble(pred=pred)
  colnames(p) = sprintf(".pred_%s_point", y)

  p
}

fit_qrf = function(train, y,...) {
  if(length(y) != 1) {
    stop("qrf needs exactly one y variable")
  }

  NIR_cols = get_numeric_colnames(train)
  train = train %>%
    select(all_of(c(y, NIR_cols))) %>%
    prefix_numeric_colnames()

  f = as.formula(sprintf("%s ~ .", y))
  quantregRanger::quantregRanger(
    f,
    data=train,
    params.ranger=list(min.bucket=1)
  )
}

predict_qrf = function(fit, newdata, y) {
  library(quantregRanger)
  library(ranger)

  pred = predict(fit,
                 data = newdata %>% prefix_numeric_colnames,
                 quantiles = c(.025, 0.25, .5, 0.75, .975))

  p = as_tibble(pred)
  colnames(p) = sprintf(".pred_%s_%s", y, c("q025", "q25", "point", "q75", "q975"))
  p
}

to_spectra_long = function(df) {
  df %>%
    pivot_longer(get_numeric_colnames(df),
               names_to="wavelength",
               values_to="reflectance",
               names_transform = as.numeric)
}

named_list_to_str2 = function(l) {
  paste0(names(l), "_", l, collapse = "_")
}

interpolate_spectra = function(df, ...) {
  index_colnames = get_nonnumeric_colnames(df)
  wavelengths = get_numeric_colnames(df) %>%
    as.numeric %>%
    sort %>%
    round

  x = do.call(seq, modifyList(list(from=wavelengths %>% first, to=wavelengths %>% last), list(...)))

  df %>%
    to_spectra_long() %>%
    nest_by(!!!syms(index_colnames)) %>%
    mutate(f = list(splinefun(data$wavelength, data$reflectance))) %>%
    mutate(data2 = list(tibble(wavelength=x, reflectance=f(x)))) %>%
    select(data2) %>%
    unnest(data2) %>%
    to_spectra_wide
}


to_spectra_wide = function(df) {
  df %>%
    pivot_wider(names_from="wavelength", values_from="reflectance") %>%
    ungroup
}

# transforms is a list of colname=transform
# where transforms follow the structure of transforms in the scales package
# e.g. list(eoc_tot_c=transform_log1p)
# drop infinite and NA values
transform_cols = function(df, transforms) {
  for(col in names(transforms)) {
    df[,col] = transforms[[col]]()$transform(df %>% pull(col))
    df = df[is.finite(df %>% pull(col)),]
  }
  df
}

# invert the above transform
# the prediction columns are expected to start with .pred_colname
# e.g. the point prediction might be .pred_eoc_tot_c
# while the prediction sd might be .pred_eoc_tot_c_sd
inv_transform_preds = function(df, transforms) {
  for(col in names(transforms)) {
    cols = colnames(df)[grepl(sprintf(".pred_%s", col), colnames(df))]
    if(length(cols) > 0) {
      for (i in 1:length(cols)) {
        c = cols[[i]]
        df[,c] = transforms[[col]]()$inverse(df %>% pull(c))
      }
    }
  }
  df
}

# for methods that only give point predictions (no .draw column)
# give them a .draw = 1 and .row column
ensure_row_in_preds = function(preds) {
  if(!".row" %in% colnames(preds)) {
    preds %>%
      mutate(.row = row_number(), .draw=NA)
  } else {
    preds
  }
}

# expects a grouped data frame
# .draw should be one of the grouping variables
aggregate_preds = function(preds) {
  if(n_distinct(preds$.draw) == 1) {
    # with uncertainty on the individual predictions we can naively get uncertainty on the aggregate
    if(".pred_eoc_tot_c_q975" %in% colnames(preds)) {
      preds %>%
        mutate(se=(.pred_eoc_tot_c_q975 - .pred_eoc_tot_c_q025)/1.96/2) %>%
        summarize(.pred_eoc_tot_c_se = sqrt(sum(se^2)/n()^2),
                  .pred_eoc_tot_c_point = mean(.pred_eoc_tot_c_point),
                  .pred_eoc_tot_c_q25 = .pred_eoc_tot_c_point - 0.67*.pred_eoc_tot_c_se,
                  .pred_eoc_tot_c_q75 = .pred_eoc_tot_c_point + 0.67*.pred_eoc_tot_c_se,
                  .pred_eoc_tot_c_q025 = .pred_eoc_tot_c_point - 1.96*.pred_eoc_tot_c_se,
                  .pred_eoc_tot_c_q975 = .pred_eoc_tot_c_point + 1.96*.pred_eoc_tot_c_se)
    # without individual uncertainty just return a point estimate
    } else {
      preds %>%
        summarize(across(.pred_eoc_tot_c_point, mean))
    }
  # with many draws just aggregate them
  } else {
    preds %>%
      summarize(across(.pred_eoc_tot_c, mean))
  }
}

# preds is a grouped tibble with columns
# .draw
# .pred_{property}
# for example it might be grouped by .row for sample summaries
# or by site_id for site_id summaries
summarize_preds = function(preds) {
  if(n_distinct(preds$.draw) == 1) {
    preds %>%
      summarize(across(starts_with(".pred_"), first))
  } else {
    preds %>%
      summarize(across(starts_with(".pred_"), list(
        point=median,
        sd=sd,
        `q025`=~quantile(.x, .025),
        `q16`=~quantile(.x, .16),
        `q25`=~quantile(.x, .25),
        `q75`=~quantile(.x, .75),
        `q84`=~quantile(.x, .84),
        `q975`=~quantile(.x, .975))))
  }
}

# turn named list of transform strings into transforms
# e.g. to_transforms(list(eoc_tot_c="log1p")) = list(eoc_tot_c=transform_log1p)
to_transforms = function(str_list) {
  sapply(str_list, function(s) syms(paste0("transform_", s)))
}
