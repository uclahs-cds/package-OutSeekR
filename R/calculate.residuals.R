#' Calculate residuals
#'
#' Calculate residuals between quantiles of the input and quantiles of one of four distributions: normal, log-normal, exponential, or gamma.
#'
#' @param x A numeric vector.
#' @param distribution A number corresponding to the optimal distribution of `x` as returned by, e.g., `identify.bic.optimal.data.distribution()`.  One of
#' * 1 = normal,
#' * 2 = log-normal,
#' * 3 = exponential, and
#' * 4 = gamma.
#'
#' @return A numeric vector of the same length as `x`.  Names are not retained.
#' @export
#' @examples
#' # Generate fake data.
#' set.seed(1234);
#' x <- rgamma(
#'     n = 20,
#'     shape = 2,
#'     scale = 2
#'     );
#' names(x) <- paste(
#'     'Sample',
#'     seq_along(x),
#'     sep = '.'
#'     );
#' calculate.residuals(
#'     x = x,
#'     distribution = 4
#'     );
#'
calculate.residuals <- function(x, distribution) {
    add.minimum.value <- least.significant.digit(x);
    # Add a minimum value to ensure the values in `x` are strictly
    # positive.  Do the same for a 5% trimmed sample of the data.  The
    # trimmed data will be used to calculate the parameters of the
    # distribution, while the untrimmed data will be used to generate
    # residuals.
    x.nozero <- x + add.minimum.value;
    x.trim <- trim.sample(
        x = x,
        trim = 0.05
        );
    x.nozero.trim <- x.trim + add.minimum.value;

    # Generate the percentiles at which theoretical quantiles will be
    # computed for the specified distribution.
    p <- stats::ppoints(
        n = x.nozero
        );
    # Get the quantiles of `x`.
    x.quantiles <- stats::quantile(
        x = x.nozero,
        probs = p
        );
    # Get the quantiles of the specified distribution.
    if (1 == distribution) {
        norm.mean <- mean(x.nozero.trim);
        norm.sd <- stats::sd(x.nozero.trim);
        theoretical.quantiles <- stats::qnorm(
            p = p,
            mean = norm.mean,
            sd = norm.sd
            );
        }
    else if (2 == distribution) {
        mean.log <- mean(x.nozero.trim);
        sd.log <- stats::sd(x.nozero.trim);
        m2 <-  log(mean.log^2 / sqrt(sd.log^2 + mean.log^2));
        sd2 <- sqrt(log(1 + (sd.log^2 / mean.log^2)));
        theoretical.quantiles <- stats::qlnorm(
            p = p,
            meanlog = m2,
            sdlog = sd2
            );
        }
    else if (3 == distribution) {
        exp.rate <- 1 / mean(x.nozero.trim);
        theoretical.quantiles <- stats::qexp(
            p = p,
            rate = exp.rate
            );
        }
    else if (4 == distribution) {
        mean.gamma <- mean(x.nozero.trim);
        sd.gamma <- stats::sd(x.nozero.trim);
        gamma.shape <- (mean.gamma / sd.gamma)^2;
        gamma.rate <- mean.gamma / (sd.gamma^2);
        theoretical.quantiles <- stats::qgamma(
            p = p,
            shape = gamma.shape,
            rate = gamma.rate
            );
        }
    # Calculate residuals between the quantiles of `x` and the
    # quantiles of the specified distribution.
    residuals <- x.quantiles - theoretical.quantiles;
    residuals;
    }
