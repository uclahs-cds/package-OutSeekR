#' Simulate from a null distribution
#'
#' Simulate transcripts from a specified null distribution.
#'
#' @param data A matrix or data frame from which to sample transcripts.
#' @param distributions A numeric vector consisting of codes corresponding to the optimal distribution of `x` as returned by, e.g., `identify.bic.optimal.data.distribution()`.  The value of `distributions[i]` corresponds to the code for the optimal distribution of the transcript in `data[i, ]`.  Possible values are
#' * 1 = normal,
#' * 2 = log-normal,
#' * 3 = exponential, and
#' * 4 = gamma.
#' @param num.null The number of simulated transcripts to generate.
#'
#' @return A matrix of simulated transcripts with `num.null` rows and `ncol(data)` columns.  Row names are retained from `data` and correspond to the row names of the transcripts used to simulate each row.  Column names are not retained.
#'
#' @examples
#' # Prepare fake data.
#' set.seed(1234);
#' n <- 25;
#' x1 <- rnorm(
#'     n = n,
#'     mean = 5,
#'     sd = 2
#'     );
#' x2 <- rlnorm(
#'     n = n,
#'     meanlog = 0,
#'     sdlog = 0.5
#'     );
#' x3 <- rexp(
#'     n = n,
#'     rate = 0.2
#'     );
#' x4 <- rexp(
#'     n = n,
#'     rate = 0.1
#'     );
#' x5 <- rgamma(
#'     n = n,
#'     shape = 2,
#'     scale = 2
#'     );
#' x6 <- rgamma(
#'     n = n,
#'     shape = 4,
#'     scale = 0.5
#'     );
#' x <- rbind(x1, x2, x3, x4, x5, x6);
#' rownames(x) <- letters[seq_len(nrow(x))];
#' colnames(x) <- paste(
#'      'Sample',
#'      seq_len(ncol(x)),
#'      sep = '.'
#'      );
#' 
#' # Identify optimal distributions for each transcript.
#' optimal.distribution.x <- apply(
#'     X = x,
#'     MARGIN = 1,
#'     FUN = identify.bic.optimal.data.distribution
#'     );
#' 
#' # Simulate from null distributions.
#' simulated.x <- simulate.null(
#'     data = x,
#'     distributions = optimal.distribution.x,
#'     num.null = 100
#'     );
#'
#' @noRd
simulate.null <- function(
    data,
    distributions,
    num.null = 10000
    ) {
    # Sample `num.null` indices from `data`.  The corresponding
    # transcripts will constitute our simulated null data set.
    sampled.indices <- sample(
        x = nrow(data),
        size = num.null,
        replace = TRUE
        );
    # For each sampled index, generate a null transcript.
    simulated.null <- future.apply::future_lapply(
        X = sampled.indices,
        FUN = function(i) {
            # Select the corresponding transcript and the code for its
            # optimal distribution.
            x <- as.numeric(data[i, ]);
            distribution <- distributions[i];
            # Ensure the values in `x` are strictly positive.
            add.minimum.value <- least.significant.digit(x);
            x.nozero <- x + add.minimum.value;
            # Apply 5% trimming.
            x.trim <- trim.sample(x);
            x.nozero.trim <- trim.sample(x.trim);
            # Generate null values according to the optimal
            # distribution for this transcript.
            if (1 == distribution) {
                norm.mean <- mean(x.nozero.trim);
                norm.sd <- stats::sd(x.trim);
                simulated.null <- truncnorm::rtruncnorm(
                    n = length(x),
                    mean = norm.mean,
                    sd = norm.sd,
                    a = 0
                    );
                }
            else if (2 == distribution) {
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
            else if (3 == distribution) {
                exp.rate <- 1 / mean(x.nozero.trim);
                simulated.null <- stats::rexp(
                    n = length(x),
                    rate = exp.rate
                    );
                }
            else if (4 == distribution) {
                mean.gamma <- mean(x.nozero.trim);
                sd.gamma <- stats::sd(x.nozero.trim);
                gamma.shape <- (mean.gamma/sd.gamma)^2;
                gamma.rate <- mean.gamma/(sd.gamma^2);
                simulated.null <- stats::rgamma(
                    n = length(x),
                    shape = gamma.shape,
                    rate = gamma.rate
                    );
                }
            },
        future.seed = TRUE
        );
    simulated.null <- do.call(
        what = rbind,
        args = simulated.null
        );
    rownames(simulated.null) <- rownames(data)[sampled.indices];
    simulated.null;
    }
