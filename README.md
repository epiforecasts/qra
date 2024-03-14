
## Quantile regression average

NB: This is a transient package that will probably be merged into the
[stackr package](https://github.com/nikosbosse/stackr).

## Prerequisites

This depends on the
[scoringutils](https://github.com/epiforecasts/scoringutils) package,
which can be installed with

``` r
remotes::install_github("epiforecasts/scoringutils")
```

The code itself can be tested by installing the corresponding package:

``` r
remotes::install_github("epiforecasts/qra")
```

### Create the forecast/data structure

We will create 3 weeks of daily (toy) forecasts, produced at 3 different
dates in May. There will be (toy) forecasts at a regional and national
level. national level

``` r
library("qra")
library("dplyr")
library("tidyr")
library("readr")
df <- tidyr::expand_grid(
                      value_type = c("cases", "deaths"),
                      geography = c(paste("region", 1:3), "country"),
                      creation_date = as.Date(c("2020-05-11",
                                                "2020-05-18",
                                                "2020-05-25")),
                      horizon = 1:21) %>%
  dplyr::mutate(value_date = creation_date + horizon) %>%
  dplyr::select(-horizon) %>%
  dplyr::mutate(geography_scale =
           dplyr::if_else(grepl("region", geography), "region", "nation"))
```

### create toy “forecasts” (draws from negative binomial distributions)

``` r
mean <- c(10L, 20L, 30L, 40L, 50L)
k <- c(0.5, 1, 1.5, 2, 3)
quantile_levels <- seq(0.05, 0.95, by = 0.05)

flist <- lapply(seq_along(mean), function(x) {
  df %>%
    rowwise() %>%
    dplyr::mutate(model = paste("model", x),
                  quantiles = list(as_tibble(t(setNames(
                    qnbinom(quantile_levels, size = 1/k[x], mu = mean[x]),
                    paste0("quantile_", quantile_levels)))))) %>%
    tidyr::unnest(quantiles) %>%
    tidyr::pivot_longer(names_to = "quantile", starts_with("quantile_")) %>%
    dplyr::mutate(quantile = readr::parse_number(quantile))
})

forecasts <- flist %>%
  dplyr::bind_rows()
```

### create toy “data”

``` r
true_mean <- 25L
true_k <- 2
data <- df %>%
  select(value_type, geography, value_date) %>%
  distinct() %>%
  mutate(value = rnbinom(n(), true_mean, 1/true_k))
```

### calculate QRA

Forecasts are pooled by forecast horizon and geography, use last \<14
days of forecasts for optimising the weights.

``` r
res <- qra::qra(forecasts,
                data, pool = c("horizon", "geography"),
                min_date = max(forecasts$creation_date) - 13)
res
```

    ## $ensemble
    ## # A tibble: 3,192 × 10
    ##    value_type geography_scale intercepts creation_date geography value_date
    ##    <chr>      <chr>           <list>     <date>        <chr>     <date>    
    ##  1 cases      nation          <tibble>   2020-05-25    country   2020-05-26
    ##  2 cases      nation          <tibble>   2020-05-25    country   2020-05-27
    ##  3 cases      nation          <tibble>   2020-05-25    country   2020-05-28
    ##  4 cases      nation          <tibble>   2020-05-25    country   2020-05-29
    ##  5 cases      nation          <tibble>   2020-05-25    country   2020-05-30
    ##  6 cases      nation          <tibble>   2020-05-25    country   2020-05-31
    ##  7 cases      nation          <tibble>   2020-05-25    country   2020-06-01
    ##  8 cases      nation          <tibble>   2020-05-25    country   2020-06-02
    ##  9 cases      nation          <tibble>   2020-05-25    country   2020-06-03
    ## 10 cases      nation          <tibble>   2020-05-25    country   2020-06-04
    ## # ℹ 3,182 more rows
    ## # ℹ 4 more variables: horizon <int>, quantile <dbl>, value <dbl>, model <chr>
    ## 
    ## $weights
    ## # A tibble: 380 × 5
    ##    value_type geography_scale model   quantile weight
    ##    <chr>      <chr>           <chr>      <dbl>  <dbl>
    ##  1 cases      nation          model 1     0.05      0
    ##  2 cases      nation          model 1     0.1       0
    ##  3 cases      nation          model 1     0.15      0
    ##  4 cases      nation          model 1     0.2       0
    ##  5 cases      nation          model 1     0.25      0
    ##  6 cases      nation          model 1     0.3       0
    ##  7 cases      nation          model 1     0.35      0
    ##  8 cases      nation          model 1     0.4       0
    ##  9 cases      nation          model 1     0.45      0
    ## 10 cases      nation          model 1     0.5       0
    ## # ℹ 370 more rows
    ## 
    ## $intercepts
    ## # A tibble: 76 × 4
    ##    value_type geography_scale quantile intercept
    ##    <chr>      <chr>              <dbl>     <dbl>
    ##  1 cases      nation              0.05         0
    ##  2 cases      nation              0.1          0
    ##  3 cases      nation              0.15         0
    ##  4 cases      nation              0.2          0
    ##  5 cases      nation              0.25         0
    ##  6 cases      nation              0.3          0
    ##  7 cases      nation              0.35         0
    ##  8 cases      nation              0.4          0
    ##  9 cases      nation              0.45         0
    ## 10 cases      nation              0.5          0
    ## # ℹ 66 more rows
