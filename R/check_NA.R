## After using BayesX, it checks if some of the models return NA values
## for the estimates

check_NA <- function(list_folders){
  as.logical(sapply(list_folders, function(a){
    fixed_effects_files <- grepl("_FixedEffects[0-9]+.res",
                                 list.files(paste0(a, '/')))

    files_results <- list.files(paste0(a, '/'))[fixed_effects_files]

    if(sum(fixed_effects_files) > 1){
      all_files <- lapply(files_results, function(aa){
        utils::read.table(paste0(a, '/', aa), head = TRUE)
      })
      fixedEffects <- do.call(rbind, all_files)
    } else {
      info <- utils::read.table(paste0(a, '/',
                                       files_results), head = TRUE)
      fixedEffects <- info
    }

    spline_files <- grepl("spline.res", list.files(paste0(a, '/')))
    splines_results <- list.files(paste0(a, '/'))[spline_files]
    info_splines <- utils::read.table(paste0(a, '/', splines_results),
                                      head = TRUE)

    if (
      any( is.na(fixedEffects$pmean) ) | any( is.na(info_splines$pmean) )
    ) TRUE
    else FALSE

  }))
}