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
#' @return A matrix of simulated transcripts with `num.null` rows and `ncol(data)` columns.  Row names are retained, but column names are not.
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
            x <- data[i, ];
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
                norm.sd <- sd(x.trim);
                simulated.null <- rtnorm(
                    n = length(x),
                    mean = norm.mean,
                    sd = norm.sd,
                    a = 0
                    );
                }
            else if (2 == distribution) {
                mean.log <- mean(x.nozero.trim);
                sd.log <- sd(x.nozero.trim);
                m2 <-  log(mean.log^2 / sqrt(sd.log^2 + mean.log^2));
                sd2 <- sqrt(log(1 + (sd.log^2 / mean.log^2)));
                simulated.null <- rlnorm(
                    n = length(x),
                    meanlog = m2,
                    sdlog = sd2
                    );
                }
            else if (3 == distribution) {
                exp.rate <- 1 / mean(x.nozero.trim);
                simulated.null <- rexp(
                    n = length(x),
                    rate = exp.rate
                    );
                }
            else if (4 == distribution) {
                mean.gamma <- mean(x.nozero.trim);
                sd.gamma <- sd(x.nozero.trim);
                gamma.shape <- (mean.gamma/sd.gamma)^2;
                gamma.rate <- mean.gamma/(sd.gamma^2);
                simulated.null <- rgamma(
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

