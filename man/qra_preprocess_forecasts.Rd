% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qra.r
\name{qra_preprocess_forecasts}
\alias{qra_preprocess_forecasts}
\title{Preprocess forecasts for QRA}
\usage{
qra_preprocess_forecasts(forecast)
}
\arguments{
\item{forecast}{a data.table representing forecast; this is expected to
have been created using \code{\link[scoringutils:as_forecast]{scoringutils::as_forecast()}}}
}
\value{
A list with three elements: \code{predictions} (a data.table only
containing the predictions), \code{data} (a data.table only containing the data)
and \code{models} (a data.table listing the models and whether they are
included)
}
\description{
This splits forecasts into data and forecast, and excludes any models that
have missing forecasts
}
\keyword{internal}
