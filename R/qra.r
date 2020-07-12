##' Helper function to convert a data frame of forecasts into a list of matrices
##'
##' @param x input data frame
##' @return list of matrices
##' @author Sebastian Funk
##' @importFrom dplyr select
##' @importFrom tidyr unite spread
##' @keywords internal
to_matrix <- function(x) {
  x %>%
    dplyr::select(-model, -horizon) %>%
    tidyr::unite(prediction_date, creation_date, value_date) %>%
    tidyr::spread(quantile, value) %>%
    dplyr::select(starts_with("0")) %>%
    as.matrix()
}

##' Helper function to estimate weights for QRA.
##'
##' for a data frame
##'
##' @param x input data frame containing \code{model}, \code{quantile}, \code{boundary},
##' \code{value}, \code{interval} columns.
##' @return data frame with weights per quantile (which won't vary unless
##' \code{per_quantile_weights} is set to TRUE), per model
##' @importFrom quantgen combine_into_array quantile_ensemble
##' @importFrom dplyr mutate select distinct group_split 
##' @importFrom tidyr expand_grid unite
##' @keywords internal
qra_estimate_weights <-
  function(x, per_quantile_weights, enforce_normalisation, ...) {

    pred_matrices <- x %>%
      dplyr::select(-data) %>%
      dplyr::group_split(model) %>%
      lapply(to_matrix) %>%
      quantgen::combine_into_array()

    data <- x %>%
      tidyr::unite(prediction_date, creation_date, value_date) %>%
      dplyr::group_by_at(dplyr::vars(-model, -horizon, -quantile, -value, -data)) %>%
      dplyr::summarise(data = unique(data)) %>%
      .$data

    tau <- x %>%
      dplyr::select(quantile) %>%
      dplyr::distinct() %>%
      .$quantile

    if (per_quantile_weights) {
      tau_groups <- seq_along(tau)
    } else {
      tau_groups <- rep(1, length(tau))
    }

    qe <-
      quantgen::quantile_ensemble(pred_matrices, data, tau,
                                  tau_groups = tau_groups,
                                  nonneg = enforce_normalisation,
                                  unit_sum = enforce_normalisation, 
                                  time_limit = 60, ...)

    ## retrieve weights from optimisation
    if (per_quantile_weights) {
      weights <- c(t(qe$alpha))
    } else {
      weights <- rep(qe$alpha, each = length(unique(x$quantile)))
    }

    ## create return tibble
    ret <- tidyr::expand_grid(model = unique(sort(x$model)),
                              quantile = unique(x$quantile)) %>%
      dplyr::mutate(weight = weights)

    return(tibble(weights = list(ret), res = list(qe)))
}

##' @name qra
##' @title Quantile Regression Average
##' Calculates a quantile regression average for forecasts.
##' @param forecasts data frame with forecasts; this is expected to have columns
##' \code{model} (a character string), \code{creation_date} (the date at which the
##' forecast was made), \code{value_date} (the date for which a forecsat was
##' created), \code{quantile} (the quantile level, between 0 and
##' 1) and \code{value} (the forecast at the quantile level)
##' @param data data frame with a \code{value} column, and otherwise matching columns
##' to \code{forecasts} (especially \code{value_date})
##' @param target_date the date for which to create a QRA; by default, will use
##' the latest \\code{creation_date} in \\code{forecasts}
##' @param min_date the minimum creation date for past forecasts to be included
##' @param max_date the maximum creation date for past forecasts to be included
##' @param history number of historical forecasts to include
##' @param pool any columns to pool as a list of character vectors (e.g.,
##' "horizon", "geography_scale", etc.) indicating columns in the \code{forecasts}
##' and \code{data} data frames; by default, will not pool across anything
##' @param quantiles Numeric - which quantiles to consider; by default will
##' consider the maximum spanning set. Options are determined by data but will be between
##' 0 and 1.
##' @param max_future Numeric - the maximum number of days of forecast to consider
##' @importFrom dplyr filter arrange desc inner_join mutate rename select bind_rows group_by_at starts_with
##' @importFrom tidyr gather complete nest spread
##' @importFrom rlang !!! syms
##' @importFrom readr parse_number
##' @importFrom tidyselect all_of
##' @importFrom purrr map
##' @inheritParams qra_estimate_weights
##' @return a list of \code{ensemble}, a data frame similar to the input forecast,
##' but with \\code{model} set to "Quantile regression average" and the values
##' set to the weighted averages; and \code{weights}, a data frame giving the weights
##' @export
qra <- function(forecasts, data, target_date, min_date, max_date, history,
                pool, quantiles, max_future = Inf,
                per_quantile_weights = FALSE, enforce_normalisation = TRUE,
                ...) {

  ## set target date to last forecast date if missing
  if (missing(target_date)) {target_date <- max(forecasts$creation_date)}
  ## initialise pooling to empty vector if not given
  if (missing(pool)) {pool <- c()}

  if (!missing(history) && history > 0 && (!missing(min_date) || !missing(max_date))) {
    stop("If 'history' is given and > 0, 'min_date' and 'max_date' can't be." )
  }

  forecasts <- forecasts %>%
    dplyr::arrange(dplyr::desc(creation_date), model, quantile)

  ## data frame with the forecasts that are being combined
  latest_forecasts <- forecasts %>%
    dplyr::filter(creation_date == target_date)

  ## prepare data frame containing data and predictions
  obs_and_pred <- forecasts %>%
    dplyr::filter(creation_date < target_date) %>%
    dplyr::mutate(horizon = value_date - creation_date)

  creation_dates <- unique(obs_and_pred$creation_date)

  ## determine dates for training data set
  if (!missing(min_date)) {
    creation_dates <- creation_dates[creation_dates >= min_date]
  }
  if (!missing(max_date)) {
    creation_dates <- creation_dates[creation_dates <= max_date]
  }
  if (!missing(history) && history > 0) {
    if (history <= length(creation_dates)) {
      creation_dates <-
        creation_dates[seq_len(history)]
    } else {
      creation_dates <- c()
    }
  }

  ## create helper vars for creating a complete set of models to be used for
  ## training below
  grouping_vars <-
    setdiff(colnames(obs_and_pred),
            c("creation_date", "value_date", "value", "model", "data", "quantile", pool))
  pooling_vars <-
    c("model", "creation_date", pool)

  present_models <- latest_forecasts %>%
    dplyr::select_at(tidyselect::all_of(c("model", grouping_vars))) %>%
    dplyr::distinct()

  ## create training data set
  obs_and_pred <- obs_and_pred %>%
    dplyr::filter(creation_date %in% creation_dates) %>%
    ## select present models
    dplyr::inner_join(present_models, by = colnames(present_models))

  ## maximum horizon by grouping variables - the last date on which all models
  ## are available
  max_horizons <- obs_and_pred %>%
    filter(horizon <= max_future) %>%
    dplyr::group_by_at(
             tidyselect::all_of(c(setdiff(grouping_vars, "horizon"),
                                  "creation_date", "model"))) %>%
    dplyr::mutate(max_horizon = max(horizon)) %>%
    dplyr::group_by_at(tidyselect::all_of(c(grouping_vars, "creation_date"))) %>%
    dplyr::summarise(max_horizon = min(max_horizon)) %>%
    dplyr::ungroup()

  obs_and_pred <- obs_and_pred %>%
    ## filter <= max horizon
    dplyr::left_join(max_horizons, by = c(grouping_vars, "creation_date")) %>%
    dplyr::filter(horizon <= max_horizon) %>%
    dplyr::select(-max_horizon) %>%
    ## join data
    dplyr::left_join(data,
               by = setdiff(colnames(data), c("value"))) %>%
    dplyr::rename(value = value.x, data = value.y)

  ## check if only specific quantiles are to be used
  if (!missing(quantiles)) {
    obs_and_pred <- obs_and_pred %>%
      dplyr::filter(quantile %in% quantiles)
  }

  ## require a complete set of forecasts to be include in QRA
  complete_set <- obs_and_pred %>%
    dplyr::group_by_at(tidyselect::all_of(c(grouping_vars, "horizon"))) %>%
    ## create complete tibble of all combinations of creation date, model
    ## and pooling variables
    tidyr::complete(!!!syms(pooling_vars)) %>%
    ## check if anything is missing and filter out
    dplyr::group_by_at(
             tidyselect::all_of(c(grouping_vars, "model"))) %>%
    dplyr::mutate(any_na = any(is.na(value))) %>%
    dplyr::filter(!any_na) %>%
    dplyr::select(-any_na) %>%
    ungroup()

  ## perform QRA
  ensemble <- complete_set %>%
    filter(!is.na(data)) %>%
    tidyr::nest(test_data = c(-setdiff(grouping_vars, "creation_date"))) %>%
    dplyr::mutate(weights =
                    purrr::map(test_data, qra_estimate_weights,
                                      per_quantile_weights, enforce_normalisation,
                               ...)) %>%
    tidyr::unnest(weights) %>%
    select(-test_data)

  weights <- ensemble %>%
    dplyr::select(-res) %>%
    tidyr::unnest(weights)

  if (nrow(weights) > 0) {
    pred <- latest_forecasts %>%
      mutate(creation_date = target_date) %>%
      ## only keep value dates which have all models present
      dplyr::group_by_at(
               tidyselect::all_of(c(grouping_vars, "value_date", "quantile"))) %>%
      dplyr::mutate(n = n()) %>%
      dplyr::group_by_at(tidyselect::all_of(grouping_vars)) %>%
      dplyr::filter(n == max(n)) %>%
      dplyr::select(-n) %>%
      ## join weights
      dplyr::mutate(horizon = value_date - creation_date) %>%
      dplyr::inner_join(weights, by = c(setdiff(grouping_vars, "creation_date"),
                                        "model", "quantile")) %>%
      dplyr::select(-horizon) %>%
      ## weigh quantiles
      dplyr::group_by_at(dplyr::vars(-model, -weight, -value)) %>%
      dplyr::summarise(value = sum(weight * value, na.rm = TRUE)) %>%
      dplyr::ungroup() %>%
      ## give model a name
      dplyr::mutate(model = "Quantile regression average")
  } else {
    pred <- latest_forecasts %>%
      filter(FALSE)
  }

  return(list(ensemble = pred, weights = weights))
}
