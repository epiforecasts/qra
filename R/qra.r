##' Helper function to calculate a weighted average interval score for a set of
##' model predictions and data and using a given set of weights
##'
##' @param weights weights given to each model, as real vector
##' @param x input data frame, containing the columns \code{model}, \code{quantile} and
##' \code{value}
##' @param enforce_normalisation if TRUE, normalisation (weights >0 and sum to
##' 1) is enforced by adding a penalty to the score if the weights are not
##' normalised
##' @param per_quantile_weights if TRUE, separate weights are calculated for
##' each quantiles
##' @importFrom dplyr rowwise summarise ungroup mutate group_by_at vars
##' @importFrom tidyr expand_grid
##' @importFrom scoringutils interval_score
##' @return an average interval score
##' @keywords internal
qra_weighted_average_interval_score <-
  function(weights, x, per_quantile_weights) {
    ## build tibble holding the per-model, (possibly) per-quantile weights
    models <- unique(x$model)
    quantiles <-unique(x$quantile)
    mw <- tidyr::expand_grid(model = models, quantile = quantiles)

    ## set weights
    if (!per_quantile_weights) {
      weights <- rep(weights, each = length(quantiles))
    }

    mw <- mw %>%
      dplyr::mutate(weight = weights)

    ## join weights with input data frame
    y <- x %>%
      dplyr::left_join(mw, by = c("model", "quantile")) %>%
      dplyr::group_by(quantile) %>%
      ungroup()

    ## calculate mean score
    mean_score <- y %>%
      dplyr::group_by_at(dplyr::vars(-model, -value, -quantile, -weight)) %>%
      dplyr::summarise(value = sum(value * weight) / sum(weight)) %>%
      dplyr::ungroup() %>%
      tidyr::spread(boundary, value) %>%
      dplyr::rowwise() %>%
      dplyr::summarise(score = scoringutils::interval_score(data, lower, upper,
                                       100 * interval, weigh = TRUE)) %>%
      dplyr::ungroup() %>%
      .$score %>%
      mean

  return(mean_score)
}

##' Helper function to estimate weights for QRA.
##'
##' for a data frame
##'
##' @param x input data frame containing \code{model}, \code{quantile}, \code{boundary},
##' \code{value}, \code{interval} columns.
##' @return data frame with weights per quantile (which won't vary unless
##' \code{per_quantile_weights} is set to TRUE), per model
##' @inheritParams qra_weighted_average_interval_score
##' @importFrom nloptr sbplx
##' @importFrom dplyr mutate
##' @importFrom tidyr expand_grid
##' @keywords internal
qra_estimate_weights <-
  function(x, per_quantile_weights, enforce_normalisation) {

  ## change x to only have present models
  x <- x %>%
    mutate(model = factor(x$model))
  ## number of models
  nmodels <- length(unique(x$model))
  init_weight <- 1 / nmodels
  ## number of quantiles
  nquant <- length(unique(x$quantile))

  if (per_quantile_weights) {
    nweights <- nmodels * nquant
  } else {
    nweights <- nmodels
  }

  ## initial weights
  init_weights <- rep(init_weight, nweights)
  ## set up optimisation
  slsqp_opts <-
    list(x0 = init_weights, fn = qra_weighted_average_interval_score,
         x = x, per_quantile_weights = per_quantile_weights)
  if (enforce_normalisation) {
    slsqp_opts[["lower"]] <- rep(0, length(init_weights))
    slsqp_opts[["upper"]] <- rep(1, length(init_weights))
    heq <- function(x) sum(x) - 1
    slsqp_opts[["heq"]] <- heq
    slsqp_opts[["heqjac"]] <- function(x) nl.jacobian(x, heq)
  }
  res <- do.call(nloptr::slsqp, slsqp_opts)

  ## retrieve weights from optimisation
  if (per_quantile_weights) {
    weights <- res$par
  } else {
    weights <- rep(res$par, each = nquant)
  }

  ## create return tibble
  ret <- tidyr::expand_grid(model = unique(x$model),
                            quantile = unique(x$quantile)) %>%
    dplyr::mutate(weight = weights)

  return(ret)
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
#' the latest \\code{creation_date} in \\code{forecasts}
##' @param min_date the minimum creation date for past forecasts to be included
##' @param max_date the maximum creation date for past forecasts to be included
##' @param history number of historical forecasts to include
##' @param pool any columns to pool as a list of character vectors (e.g.,
##' "horizon", "geography_scale", etc.) indicating columns in the \code{forecasts}
##' and \code{data} data frames; by default, will not pool across anything
##' @param intervals Numeric - which central intervals to consider; by default will
##' consider the maximum spanning set. Options are determined by data but will be between
##' 0 and 1.
##' @param max_future Numeric - the maximum number of days of forecast to consider
##' @importFrom dplyr filter arrange desc inner_join mutate rename select bind_rows group_by_at starts_with
##' @importFrom tidyr gather complete nest spread
##' @importFrom rlang !!! syms
##' @importFrom readr parse_number
##' @importFrom tidyselect all_of
##' @importFrom future plan
##' @importFrom furrr future_map
##' @inheritParams qra_weighted_average_interval_score
##' @return a list of \code{ensemble}, a data frame similar to the input forecast,
##' but with \\code{model} set to "Quantile regression average" and the values
##' set to the weighted averages; and \code{weights}, a data frame giving the weights
##' @export
qra <- function(forecasts, data, target_date, min_date, max_date, history,
                pool, intervals, max_future = Inf,
                per_quantile_weights = FALSE, enforce_normalisation = TRUE) {

  ## set target date to last forecast date if missing
  if (missing(target_date)) {target_date <- max(forecasts$creation_date)}
  ## initialise pooling to empty vector if not given
  if (missing(pool)) {pool <- c()}

  if (!missing(history) && (!missing(min_date) || !missing(max_date))) {
    stop("If 'history' is given, 'min_date' and 'max_date' can't be." )
  }

  ## data frame with the forecasts that are being combined
  latest_forecasts <- forecasts %>%
    dplyr::filter(creation_date == target_date)

  ## prepare data frame containing data and predictions
  obs_and_pred <- forecasts %>%
    dplyr::filter(creation_date < target_date) %>%
    dplyr::arrange(dplyr::desc(creation_date)) %>%
    dplyr::mutate(horizon = value_date - creation_date)

  creation_dates <- unique(obs_and_pred$creation_date)

  ## determine dates for training data set
  if (!missing(min_date)) {
    creation_dates <- creation_dates[creation_dates >= min_date]
  }
  if (!missing(max_date)) {
    creation_dates <- creation_dates[creation_dates <= max_date]
  }
  if (!missing(history)) {
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
            c("value_date", "value", "model", "data",
              "quantile", "boundary", "interval", pool))
  pooling_vars <-
    c("creation_date", "model", pool)

  present_models <- latest_forecasts %>%
    dplyr::select_at(tidyselect::all_of(c("model", grouping_vars))) %>%
    select(-creation_date) %>%
    distinct()

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
             tidyselect::all_of(c(setdiff(grouping_vars, "horizon"), "model"))) %>%
    dplyr::mutate(max_horizon = max(horizon)) %>%
    dplyr::group_by_at(tidyselect::all_of(c(grouping_vars))) %>%
    dplyr::summarise(max_horizon = min(max_horizon)) %>%
    dplyr::ungroup()

  obs_and_pred <- obs_and_pred %>%
    ## filter <= max horizon
    dplyr::left_join(max_horizons, by = grouping_vars) %>%
    dplyr::filter(horizon <= max_horizon) %>%
    dplyr::select(-max_horizon) %>%
    ## join data
    dplyr::inner_join(data,
               by = setdiff(colnames(data), c("value"))) %>%
    dplyr::rename(value = value.x, data = value.y) %>%
    ## work out intervals and boundaries
    dplyr::mutate(interval = round(2 * abs(quantile - 0.5), 2),
                  boundary =
                    dplyr::if_else(quantile < 0.5, "lower", "upper"))

  ## check if only specific intervals are to be used
  if (!missing(intervals)) {
    obs_and_pred <- obs_and_pred %>%
      dplyr::filter(interval %in% intervals)
  }

  ## duplicated median (quantila == 0.5) as it's both the lower and upper bound
  ## of the corresponding interval
  alpha_lower <- obs_and_pred %>%
    dplyr::filter(quantile == 0.5) %>%
    dplyr::mutate(boundary = "lower")

  obs_and_pred_double_alpha <- obs_and_pred %>%
    dplyr::bind_rows(alpha_lower) %>%
    select(-value_date)

  ## require a complete set of forecasts to be include in QRA
  complete_set <- obs_and_pred_double_alpha %>%
    dplyr::group_by_at(
             tidyselect::all_of(
                           c(grouping_vars, "interval"))) %>%
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

  future::plan(multiprocess)

  ## perform QRA
  weights <- complete_set %>%
    tidyr::nest(test_data = c(-setdiff(grouping_vars, "creation_date"))) %>%
    dplyr::mutate(weights =
                    furrr::future_map(test_data, qra_estimate_weights,
                                      per_quantile_weights, enforce_normalisation)) %>%
    tidyr::unnest(weights) %>%
    dplyr::select(-test_data)

  if (nrow(weights) > 0) {
    ensemble <- latest_forecasts %>%
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

    weights_ret <- weights %>%
      dplyr::arrange(quantile) %>%
      dplyr::mutate(quantile =
                      paste0("percentile_", sprintf("%.2f", quantile))) %>%
      tidyr::spread(quantile, weight)
  } else {
    ensemble <- latest_forecasts %>%
      filter(FALSE)
    weights_ret <- NULL
  }

  return(list(ensemble = ensemble, weights = weights_ret))
}
