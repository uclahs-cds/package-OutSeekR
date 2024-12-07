% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/identify.bic.optimal.data.distribution.R
\name{identify.bic.optimal.data.distribution}
\alias{identify.bic.optimal.data.distribution}
\title{Identify optimal distribution of data}
\usage{
\method{identify}{bic.optimal.data.distribution}(x)
}
\arguments{
\item{x}{A numeric vector.}
}
\value{
A numeric code representing which distribution optimally fits \code{x}.  Possible values are
\itemize{
\item 1 = normal,
\item 2 = log-normal,
\item 3 = exponential, and
\item 4 = gamma.
}
}
\description{
Identify which of four distributions---normal, log-normal, exponential, or gamma---best fits the given data according to BIC.
}
\examples{
# Generate fake data.
set.seed(1234);
x <- rgamma(
    n = 20,
    shape = 2,
    scale = 2
    );
identify.bic.optimal.data.distribution(
    x = x
    );
}
