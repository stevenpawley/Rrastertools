#' Hyperparameter tuning using the out-of-bag error of ranger models
#'
#' @param object partial function of ranger
#' @param formula model formula
#' @param param_grid data.frame of hyperparameters, with names of
#' hyperparameters as columns and each trial as a row
#' @param data data.frame of training data
#' @param ... additional named arguments to pass onto ranger function
#'
#' @return list
#' best_model : model refitted using best hyperparameters
#' param_grid :  data.frame of param_grid with additional score column
#' @export
#' @importFrom magrittr "%>%"
ranger_tuned <- function(object, formula, param_grid, data, ...) {
  
  kwargs <- list(...)
  
  # iterate through each hyperparameter set and score oobe
  scores <- 
    iterators::iter(param_grid, by = 'row') %>%
    sapply(function(pm) {
      pm <- c(pm, kwargs)
      m <- do.call(pryr::partial(
        object, data = data, formula = formula), pm)
      m$prediction.error
      })
  
  # select best hyperparameter combination
  best_params <- param_grid[which.min(scores), ]
  
  # refit model (slower but saves memory)
  best_model <- do.call(pryr::partial(object, data = data, formula = formula),
                        best_params)
  
  # insert param grid and scores into the ranger class object
  best_model[['best_params']] <- cbind(param_grid, score = scores)
  
  return(best_model)
}


#' K-Fold cross-validation scheme
#'
#' @param object partial function of ranger
#' @param formula model forumula
#' @param data data.frame of training data
#' @param param_grid pass a function that generates an expand.grid of
#' hyperparameters
#' @param n_splits integer, number of cross validation folds
#' @param keep_models logical, save each model to the tibble per cross
#' validation fold. Note that this can consume a lot of memory for large
#' datasets
#' @param ... additional named arguments to pass onto ranger function
#' @return tibble containing cross validation results
#' @export
#' @importFrom magrittr "%>%"
cross_validate <- function(object, formula, data, param_grid = NULL,
                           n_splits = 10, keep_models = FALSE, ...) {
  
  me <- function(truth, estimate)
    return(mean(truth - estimate))
  
  # use tuned forest if param_grid passed
  if (!missing(param_grid)) {
    estimator <- pryr::partial(
      ranger_tuned, object = object, param_grid = param_grid)
    
  } else {
    estimator <- object
  }
  
  # fit models on analysis partitions and predict assessment partitions
  if (keep_models == TRUE) {
    
    folds = rsample::vfold_cv(data, v = n_splits) %>%
      dplyr::mutate(models = purrr::map(
        .$splits, ~ estimator(formula = formula, data = rsample::analysis(x))))
    
    folds <- folds %>% dplyr::mutate(
      test_predictions = purrr::map2(
        .$splits, .$models, ~ predict(.y, rsample::assessment(.x))$predictions)
    )
    
  } else {
    
    folds = rsample::vfold_cv(data, v = n_splits) %>% 
      dplyr::mutate(test_predictions = purrr::map(.$splits, ~ {
        m <- estimator(formula = formula, data = rsample::analysis(.x))
        predict(m, rsample::assessment(.x))$predictions
      }))
  }
  
  # calculate metrics
  folds = folds %>%
    dplyr::mutate(
      me = purrr::map2_dbl(
        .x = .$splits,
        .y = .$test_predictions,
        .f = ~ me(truth = rsample::assessment(.x)[[formula.tools::lhs(formula)]],
                  estimate = .y)),
      mae = purrr::map2_dbl(
        .x = .$splits,
        .y = .$test_predictions,
        .f = ~ yardstick::mae_vec(truth = rsample::assessment(.x)[[formula.tools::lhs(formula)]],
                                  estimate = .y)),
      rmse = purrr::map2_dbl(
        .x = .$splits,
        .y = .$test_predictions,
        .f = ~ yardstick::rmse_vec(truth = rsample::assessment(.x)[[formula.tools::lhs(formula)]],
                                   estimate = .y))
    )
  
  return(folds)
}


#' Apply trained ranger model to a RasterStack or RasterBrick
#'
#' @param model fitted ranger object
#' @param data passed from raster::predict function
#' @param factor_ldf optional, named list of dataframes containing
#' reclassification for factors in the raster data.  Each data.frame must
#' contain 3 columns, with column 1 containing raster values (ID),
#' column 2 containing the original levels that are associated with each raster 
#' value, and column 3 containing the collapsed levels
#' @param ... additional named arguments to pass to the predict method
#'
#' @return RasterLayer object
#' @export
pred_fun_ranger <- function(model, data, factor_ldf = NULL, ...) {
  
  kwargs <- list(...)
  
  if (!is.null(factor_ldf)) {
    
    factor_var_names <- names(factor_ldf)
    
    for (name in factor_var_names) {
      data_for_level <- data[, which(colnames(data) == name)]
      data[, which(colnames(data) == name)] <- 
        factor_ldf[[name]][data_for_level, 'collapsed']
    }
  }
  
  do.call(pryr::partial(stats::predict, object = model, data = data), 
          kwargs)$predictions
}
