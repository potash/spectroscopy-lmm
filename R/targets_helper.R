# # note that the names need to be the last names used in the old targets!
# tar_combine_sub = function(new_values, old_targets, names) {
#   names_quosure <- rlang::enquo(names)
#   names <- eval_tidyselect(names_quosure, base::names(new_values))
#   suffix = tar_map_produce_suffix(new_values, names)
#   
#   # create an empty dataframe with colnames equal to old target names
#   vec = names(old_targets)
#   old_vars <- bind_rows(setNames(rep("", length(vec)), vec))[0, ]
#   
#   new_values %>%
#     mutate(.suffix=suffix) %>%
#     bind_rows(old_vars) %>%
#     rowwise() %>%
#     mutate(across(all_of(vec), 
#                   ~ list(tar_filter(old_targets, cur_column(), .suffix)))) %>%
#     select(-.suffix)
# }
# 
# given the values passed to tar_map, and the result (unlist=FALSE)
# create an expanded values that adds the mapped steps as columns to do another map that references them
tar_add_steps_to_values = function(values, targets) {
  new_values = as_tibble(values)

  for(name in names(targets)) {
    new_values[[name]] = tar_syms(targets[[name]])
  }
  new_values
}
# 
# tar_cross = function(new_values, old_targets, old_values, names) {
#   # create an empty dataframe with colnames equal to old target names
#   vec = names(old_targets)
#   old_vars <- bind_rows(setNames(rep("", length(vec)), vec))[0, ]
#   
#   new_values %>%
#     bind_rows(old_vars) %>%
#     mutate(across(all_of(vec), ~ rlang::syms(paste0(cur_column(), "_", get(names)))))
# }
# 
# tar_select_syms = function(targets, match) {
#   list(rlang::syms(
#     tar_select_names(
#       targets,
#       matches(match)
#     )))
# }
# 
tar_syms = function(targets) {
  targets::tar_assert_target_list(targets)
  names_chr <- map_chr(targets, ~.x$settings$name)
  names_sym <- lapply(names_chr, as.symbol)
  names(names_sym) <- names_chr
  names_sym
}
# 
# tar_filter = function(old_targets, name, suffix) {
#   targets = old_targets[[name]]
#   syms = tar_syms(targets)
#   syms[str_ends(names(syms), suffix)]
# }
# 
# eval_tidyselect <- function(names_quosure, choices) {
#   if (is.null(rlang::quo_squash(names_quosure)) || !length(choices)) {
#     return(NULL)
#   }
#   names(choices) <- choices
#   out <- tidyselect::eval_select(names_quosure, data = choices, strict = FALSE)
#   out <- names(out)
#   out
# }
# 
# tar_map_produce_suffix <- function(values, names) {
#   data <- values[names] %||% tar_map_default_suffixes(values)
#   data <- map(data, ~as.character(unlist(.x)))
#   out <- apply(as.data.frame(data), 1, paste, collapse = "_")
#   out <- gsub("'", "", out)
#   out <- gsub("\"", "", out)
#   make.unique(out, sep = "_")
# }

# bind 
bind_rows_with_name = function(..., .names_to="name") {
  args = list(...)
  dots = match.call(expand.dots = FALSE)$...
  names = sapply(dots, deparse)
  for(i in 1:length(args)) {
    args[[i]][.names_to] = names[[i]]
  }
  df = bind_rows(args)
  df
}
