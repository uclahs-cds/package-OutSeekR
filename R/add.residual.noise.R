#' Add noise to residuals
#'
#' Add additional noise to residuals based on a specified null distribution.
#'
#' @param x A numeric vector.
#' @param distribution  A numeric code corresponding to the optimal distribution of `x` as returned by `identify.bic.optimal.residuals.distribution()`.  Possible values are
#' * 1 = normal,
#' * 2 = log-normal,
#' * 3 = exponential, and
#' * 4 = gamma.
#'
#' @return A numeric vector of the same length as `x`.
#'
#' @noRd
add.residual.noise <- function(x, distribution) {
    # Add a minimum value to ensure the values in `x` are strictly
    # positive.
    add.minimum.value <- least.significant.digit(x)
    if (min(x) < 0) {
        x.nozero <- x - min(x) + add.minimum.value;
        }
    else {
        x.nozero <- x + add.minimum.value;
        }
    # Generate null values from the distribution coded by
    # `distribution`.
    if (1 == distribution) {
        norm.mean <- mean(x.nozero);
        norm.sd <- stats::sd(x.nozero);
        simulated.noise <- truncnorm::rtruncnorm(
            n = length(x),
            mean = norm.mean,
            sd = norm.sd,
            a = 0
            );
        }
    else if (2 == distribution) {
        mean.log <- mean(x.nozero);
        sd.log <- stats::sd(x.nozero);
        m2 <- log(mean.log^2 / sqrt(sd.log^2 + mean.log^2));
        sd2 <- sqrt(log(1 + (sd.log^2 / mean.log^2)));
        simulated.noise <- stats::rlnorm(
            n = length(x),
            meanlog = m2,
            sdlog = sd2
            );
        }
    else if (3 == distribution) {
        exp.rate <- 1 / mean(x.nozero);
        simulated.noise <- stats::rexp(
            n = length(x),
            rate = exp.rate
            );
        }
    else if (4 == distribution) {
        mean.gamma <- mean(x.nozero);
        sd.gamma <- stats::sd(x.nozero);
        gamma.shape <- (mean.gamma / sd.gamma)^2;
        gamma.rate <- mean.gamma / (sd.gamma^2);
        simulated.noise <- stats::rgamma(
            n = length(x),
            shape = gamma.shape,
            rate = gamma.rate
            );
        }
    # Add the simulated noise to `x`, but make sure the result is
    # positive.
    print(simulated.noise);
    abs(x + simulated.noise);
    }
