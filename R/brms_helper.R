# 
# BRMS_SAVE_DIR="models/brms_output"
# 

# this runs before fitting a brms model
# prefix NIR colnames, e.g. 1350 -> x_1350
#   because brms doesn't like nuermic colnames
# standardize outcome variable(s) and save their 
#   mean and sd so that transform can be inverted when predicting
# generate a formula with all the NIR variables on the RHS
fit.brms_helper_pre = function(train, y) {
  NIR_cols = get_numeric_colnames(train)

  train = train %>%
    prefix_numeric_colnames

  y_stats = train %>%
    select(all_of(y)) %>%
    pivot_longer(everything()) %>%
    group_by(name) %>%
    summarize(mean=mean(value), sd=sd(value))

  # standardize outcome
  train = train %>%
    mutate(across(all_of(y), function(x) (x - mean(x))/sd(x)))

  list(y_stats=y_stats,
       train=train,
       NIR_cols=NIR_cols,
       RHS = paste(paste0("x_", NIR_cols), collapse="+"))
}

# fit a brms model with specific arguments
fit.brms_helper_fit = function(formula, pre,
                               ...) {
  default_args = list(formula=formula, data=pre$train,
                      backend="cmdstanr",
                      iter=200, chains=4,
                      stan_model_args=list(stanc_options = list("O1")))
  args = modifyList(default_args, list(...))
  fit = do.call(brm, args)
  fit$y_stats = pre$y_stats

  fit
}

# for each column {SSURGO} in the character vector SSURGO, e.g. SSURGO_site_clay_r
# and each spectral column x_{wavelength}
# generate the interaction column {SSURGO}_x_{wavelength}
add_SSURGO_interactions = function(df, SSURGO) {
  interaction_fns = map(SSURGO, function(x) function(y) y*df[[x]])
  names(interaction_fns) = SSURGO

  df %>%
    mutate(across(starts_with("x_"), interaction_fns, .names="{.fn}_{.col}"))
}
# 
# add_SSURGO_interactions2 = function(df, SSURGO) {
#   if(length(SSURGO) > 1) {
#     interaction_fns = map(SSURGO, function(x) function(y) y*df[[x]])
#     names(interaction_fns) = SSURGO
#     
#     df %>%
#       mutate(across(starts_with("x_"), interaction_fns, .names="{.fn}_{.col}"))
#   } else {
#     df
#   }
# }
# 

fit_lmer_SSURGO.brms = function(train, y, group, SSURGO,
                                varying_slopes=TRUE,
                                sigmabeta=TRUE,
                                cor=FALSE, ...) {

  pre = fit.brms_helper_pre(train, y)
  SSURGO = str_glue("SSURGO_site_{SSURGO}_r")
  pre$train = add_SSURGO_interactions(pre$train, SSURGO)

  # scale SSURGO and interactions
  recipe_SSURGO = pre$train %>%
    recipe() %>%
    update_role(contains(SSURGO), new_role="predictor") %>%
    step_center(all_predictors()) %>%
    step_scale(all_predictors()) %>%
    prep(training=pre$train) %>%
    # do not require other columns at bake
    update_role_requirements(NA, bake=FALSE)

  pre$train = bake(recipe_SSURGO, pre$train)

  NIR_cols = str_split(pre$RHS, "\\+")[[1]]
  SSURGO_int = expand_grid(SSURGO, NIR_cols) %>%
    mutate(int=str_glue("{SSURGO}_{NIR_cols}")) %>%
    pull(int)

  SSURGO_formula = paste(SSURGO, collapse="+")
  SSURGO_int_formula = paste(SSURGO_int, collapse="+")
  # add priors, one set for SSURGO, one for SSURGO interactions

  f_template = "{y} ~ {pre$RHS}"
  if(varying_slopes) {
    f_template = paste0(f_template, " + ({pre$RHS}|{group})")
  } else {
    f_template = paste0(f_template, " + (1|{group})")
  }

  if(length(SSURGO) > 0) {
    f_template = paste0(f_template, " + {SSURGO_formula}")
    if(varying_slopes) {
      f_template = paste0(f_template, " + {SSURGO_int_formula}")
    }
  }

  f = as.formula(str_glue(f_template))

  # get groups by calling get_prior
  groups = get_prior(formula=f, data=pre$train) %>%
    filter(group != "") %>% pull(group) %>% unique
  # clean up names
  group_vars = 1:length(groups)

  prior = c(
    set_prior("normal(0,2.5)", class="Intercept"),
    set_prior("exponential(1)", class="sigma"))

  if(length(SSURGO) > 0) {
    # std normal prior for the SSURGO slopes (if any)
    prior = prior +
      set_prior("normal(0,1)", class="b")
  }

  stanvars=c(
    # main NIR columns shrinkage
    stanvar(scode="real<lower=0> sigmabeta;", name="sigmabeta", block="parameters"),
    stanvar(scode="sigmabeta ~ std_normal();", block="model"))

  if(sigmabeta) {
    prior = prior +
      Reduce("+",
        # main NIR cols
        map(NIR_cols, function(c)
          set_prior(
            str_glue("normal(0,sigmabeta)"),
            class = "b", coef = c
          )))
  }

  # Varying intercepts prior
  prior = prior + Reduce(
    "+", map2(groups, group_vars, function(group, group_var)
      set_prior(
        str_glue("constant(sd_intercept_{group_var})"),
        class = "sd", coef = "Intercept", group = group
      )))

  stanvars = stanvars + Reduce("+", c(
    map2(groups, group_vars, function(group, group_var)
      stanvar(scode=str_glue("real<lower=0> sd_intercept_{group_var};"),
              name=str_glue("sd_intercept_{group_var}"), block="parameters")),

    map2(groups, group_vars, function(group, group_var)
      stanvar(scode=str_glue("sd_intercept_{group_var} ~ normal(0, 1);"), block="model"))
  ))

  # whether to make random slopes (and intercepts) correlated
  if(varying_slopes) {
    if(cor) {
      prior = prior + set_prior("lkj(2)", class="cor")
    } else {
      prior = prior + prior(constant(Id, broadcast=TRUE), class="cor")
      stanvars = stanvars + stanvar(diag(1, length(pre$NIR_cols)+1), name="Id")
    }

    prior = prior + Reduce(
      "+",
      map2(groups, group_vars, function(group, group_var)
        set_prior(
          str_glue("constant(sd_slope_{group_var})"),
          class = "sd", group = group
        )))

    stanvars = stanvars + Reduce("+", c(
      map2(groups, group_vars, function(group, group_var)
        stanvar(scode=str_glue("real<lower=0> sd_slope_{group_var};"),
                name=str_glue("sd_slope_{group}"), block="parameters")),
      map2(groups, group_vars, function(group, group_var)
        stanvar(scode=str_glue("sd_slope_{group_var} ~ normal(0, 1);"), block="model"))
    ))

    # interaction NIR-SSURGO cols
    if(length(SSURGO) > 0) {
      prior = prior + Reduce(
        "+", map(SSURGO_int, function(c) set_prior(
          str_glue("normal(0,sigmabetaint)"), class = "b", coef = c
          )))

      stanvars = stanvars +
        stanvar(scode="real<lower=0> sigmabetaint;", name="sigmabetaint", block="parameters") +
        stanvar(scode="sigmabetaint ~ std_normal();", block="model")
    }
  }

  f = fit.brms_helper_fit(
    formula=f, pre=pre,
    prior=prior,
    stanvars=stanvars
  )

  f$recipe_SSURGO = recipe_SSURGO
  f$SSURGO_cols = SSURGO
  f
}
# 
# predict_lmer_vivs2_cor_SSURGO_sigmabeta.brms = function(fit, newdata, y) {
#   newdata = newdata %>%
#     prefix_numeric_colnames
#   newdata = add_SSURGO_interactions(newdata, fit$SSURGO_cols)
#   
#   newdata = bake(fit$recipe_SSURGO, newdata)
#   predict_lm.brms(fit, newdata, y)
# }
# 
# 
# # add site level means2
# predict_lmer_vivs2_cor_SSURGO2.brms2 = function(fit, newdata, y) {
#   
#   newdata = newdata %>% inner_join(get_SSURGO_site_means2(newdata), by="site_id")
#   newdata = newdata %>%
#     prefix_numeric_colnames
#   newdata = add_SSURGO_interactions(newdata, fit$SSURGO_cols)
#   
#   newdata = bake(fit$recipe_SSURGO, newdata)
#   
#   predict_lm.brms(fit, newdata, y)
# }
# 
# 
# predict_lmer_SSURGO_mean2.brms = predict_lmer_vivs2_cor_SSURGO2.brms2
# 
# predict_lm.brms = function(fit, newdata, y) {
#   p1 = add_predicted_draws(fit, 
#                           newdata=newdata %>% prefix_numeric_colnames,
#                           allow_new_levels=TRUE,
#                           sample_new_levels="uncertainty") %>%
#     ungroup
#   
#   p2 = add_epred_draws(fit, 
#                        newdata=newdata %>% prefix_numeric_colnames,
#                        allow_new_levels=TRUE,
#                        sample_new_levels="uncertainty") %>%
#     rename(.prediction=.epred)
#     ungroup
#   
#   trans = function(p, suffix="") {
#     p %>%
#       mutate(.prediction = .prediction * fit$y_stats$sd[[1]] + fit$y_stats$mean[[1]]) %>%
#       select(.row, .draw, .prediction) %>%
#       rename(".pred_{y}{suffix}" := .prediction)
#   }
#   
#   inner_join(trans(p1),
#              trans(p2, "_epred"))
# }