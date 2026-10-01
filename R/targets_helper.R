# given the values passed to tar_map, and the result (unlist=FALSE)
# create an expanded values that adds the mapped steps as columns to do another map that references them
tar_add_steps_to_values = function(values, targets) {
  new_values = as_tibble(values)

  for(name in names(targets)) {
    new_values[[name]] = tar_syms(targets[[name]])
  }
  new_values
}

tar_syms = function(targets) {
  targets::tar_assert_target_list(targets)
  names_chr <- map_chr(targets, ~.x$settings$name)
  names_sym <- lapply(names_chr, as.symbol)
  names(names_sym) <- names_chr
  names_sym
}
