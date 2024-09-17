#' Simulate from a null distribution
#'
#' Simulate transcripts from a specified null distribution.
#'
#' @param x A numeric vector of transcripts.
#' @param x.distribution A numeric code corresponding to the optimal distribution of `x` as returned by `identify.bic.optimal.data.distribution()`.  Possible values are
#' * 1 = normal,
#' * 2 = log-normal,
#' * 3 = exponential, and
#' * 4 = gamma.
#' @param r A numeric vector of residuals calculated for this transcript.
#' @param r.distribution A numeric code corresponding to the optimal distribution of `x` as returned by `identify.bic.optimal.residuals.distribution()`.  Possible values are the same as those for `x.distribution`.
#'
#' @return A numeric vector of the same length as `x`.  Names are not retained.
#' @export simulate.null
#' @examples
#' # Prepare fake data.
#' set.seed(1234);
#' x <- rgamma(
#'     n = 20,
#'     shape = 2,
#'     scale = 2
#'     );
#' names(x) <- paste('Sample', seq_along(x), sep = '.');
#' x.dist <- identify.bic.optimal.data.distribution(
#'     x = x
#'     );
#' r <- calculate.residuals(
#'     x = x,
#'     distribution = x.dist
#'     );
#' r.trimmed <- trim.sample(
#'     x = r
#'     );
#' r.dist <- identify.bic.optimal.residuals.distribution(
#'     x = r.trimmed
#'     );
#' null <- simulate.null(
#'     x = x,
#'     x.distribution = x.dist,
#'     r = r.trimmed,
#'     r.distribution = r.dist
#'     );
simulate.null <- function(
    x,
    x.distribution,
    r,
    r.distribution
    ) {
    #
    # Simulate transcripts
    #
    # Ensure the values in `x` are strictly positive.
    add.minimum.value <- least.significant.digit(x);
    x.nozero <- x + add.minimum.value;
    # Apply 5% trimming.
    # x.trim <- trim.sample(x);
    x.nozero.trim <- trim.sample(x.nozero);
    # Generate null values according to the optimal
    # distribution for this transcript.
    if (1 == x.distribution) {
        norm.mean <- mean(x.nozero.trim);
        norm.sd <- stats::sd(x.nozero.trim);
        simulated.null <- truncnorm::rtruncnorm(
            n = length(x),
            mean = norm.mean,
            sd = norm.sd,
            a = 0
            );
        }
    else if (2 == x.distribution) {
        mean.log <- mean(x.nozero.trim);
        sd.log <- stats::sd(x.nozero.trim);
        m2 <-  log(mean.log^2 / sqrt(sd.log^2 + mean.log^2));
        sd2 <- sqrt(log(1 + (sd.log^2 / mean.log^2)));
        simulated.null <- stats::rlnorm(
            n = length(x),
            meanlog = m2,
            sdlog = sd2
            );
        }
    else if (3 == x.distribution) {
        exp.rate <- 1 / mean(x.nozero.trim);
        simulated.null <- stats::rexp(
            n = length(x),
            rate = exp.rate
            );
        }
    else if (4 == x.distribution) {
        mean.gamma <- mean(x.nozero.trim);
        sd.gamma <- stats::sd(x.nozero.trim);
        gamma.shape <- (mean.gamma / sd.gamma)^2;
        gamma.rate <- mean.gamma / (sd.gamma^2);
        simulated.null <- stats::rgamma(
            n = length(x),
            shape = gamma.shape,
            rate = gamma.rate
            );
        }
    #
    # Simulate noise
    #
    # Ensure the values in `r` are strictly positive.
    add.minimum.value <- least.significant.digit(r)
    if (min(r) < 0) {
        r.nozero <- r - min(r) + add.minimum.value;
        }
    else {
        r.nozero <- r + add.minimum.value;
        }
    # Generate null values from the distribution coded by
    # `distribution`.
    if (1 == r.distribution) {
        norm.mean <- mean(r.nozero);
        norm.sd <- stats::sd(r.nozero);
        simulated.noise <- truncnorm::rtruncnorm(
            n = length(x),
            mean = norm.mean,
            sd = norm.sd,
            a = 0
            );
        }
    else if (2 == r.distribution) {
        mean.log <- mean(r.nozero);
        sd.log <- stats::sd(r.nozero);
        m2 <- log(mean.log^2 / sqrt(sd.log^2 + mean.log^2));
        sd2 <- sqrt(log(1 + (sd.log^2 / mean.log^2)));
        simulated.noise <- stats::rlnorm(
            n = length(x),
            meanlog = m2,
            sdlog = sd2
            );
        }
    else if (3 == r.distribution) {
        exp.rate <- 1 / mean(r.nozero);
        simulated.noise <- stats::rexp(
            n = length(x),
            rate = exp.rate
            );
        }
    else if (4 == r.distribution) {
        mean.gamma <- mean(r.nozero);
        sd.gamma <- stats::sd(r.nozero);
        gamma.shape <- (mean.gamma / sd.gamma)^2;
        gamma.rate <- mean.gamma / (sd.gamma^2);
        simulated.noise <- stats::rgamma(
            n = length(x),
            shape = gamma.shape,
            rate = gamma.rate
            );
        }
    if (min(r) < 0) {
        simulated.noise <- simulated.noise + min(r) - add.minimum.value;
        }
    else {
        simulated.noise <- simulated.noise - add.minimum.value;
        }
    # Add the simulated noise to the simulated transcript.
    abs(simulated.null + simulated.noise)
    }
