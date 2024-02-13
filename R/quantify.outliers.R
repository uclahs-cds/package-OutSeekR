#' Compute quantities for outlier detection
#'
#' Compute quantities for use in the detection of outliers.  Specifically, compute z-scores based on the mean / standard deviation, the trimmed mean / trimmed standard deviation, or the median / median absolute deviation, or the cluster assignment from k-means with two clusters.
#'
#' @param x A numeric vector.
#' @param method A string indicating the quantities to be computed.  Possible values are
#' * `'mean'` : z-scores based on mean and standard deviation or trimmed mean and trimmed standard deviation if `trim > 0`,
#' * `'median'` : z-scores based on median and median absolute deviation, or
#' * `'kmeans'` : cluster assignment from k-means with two clusters.
#' The default is z-scores based on the mean and standard deviation.
#' @param trim A number, the fraction of observations to be trimmed from each end of `x`.  Default is no trimming.
#' @param nstart A number, for k-means clustering, the number of random initial centers for the clusters.  Default is `1`.  See [stats::kmeans()] for further information.
#' @param exclude.zero A logical, whether zeros should be excluded (`TRUE`) or not excluded (`FALSE`, the default) from computations.  For `method = 'mean'` and `method = 'median'`, this means zeros will not be included in computing the summary statistics; for `method = 'kmeans'`, this means zeros will be placed in their own cluster, coded `0`.
#'
#' @return A numeric vector the same size as `x` whose values are the requested quantities computed on the corresponding elements of `x`.
#'
#' @examples
#' # Generate fake data.
#' set.seed(1234);
#' x <- rgamma(
#'     n = 20,
#'     shape = 2,
#'     scale = 2
#'     );
#' # Add missing values and zeros for demonstration.  Missing values are
#' # ignored, and zeros can be ignored with `exclude.zeros = TRUE`.
#' x[1:5] <- NA;
#' x[6:10] <- 0;
#'
#' # Compute z-scores based on mean and standard deviation.
#' quantify.outliers(
#'     x = x,
#'     method = 'mean',
#'     trim = 0
#'     );
#' # Exclude zeros from the calculation of the mean and standard
#' # deviation.
#' quantify.outliers(
#'     x = x,
#'     method = 'mean',
#'     trim = 0,
#'     exclude.zero = TRUE
#'     );
#'
#' # Compute z-scores based on the 5% trimmed mean and 5% trimmed
#' # standard deviation.
#' quantify.outliers(
#'     x = x,
#'     method = 'mean',
#'     trim = 0.05
#'     );
#'
#' # Compute z-scores based on the median and median absolute deviation.
#' quantify.outliers(
#'     x = x,
#'     method = 'median'
#'     );
#'
#' # Compute cluster assignments using k-means with k = 2.
#' quantify.outliers(
#'     x = x,
#'     method = 'kmeans'
#'     );
#' # Try different initial cluster assignments.
#' quantify.outliers(
#'     x = x,
#'     method = 'kmeans',
#'     nstart = 10
#'     );
#' # Assign zeros to their own cluster.
#' quantify.outliers(
#'     x = x,
#'     method = 'kmeans',
#'     exclude.zero = TRUE
#'     );
#'
#' @noRd
quantify.outliers <- function(x, method = 'mean', trim = 0, nstart = 1, exclude.zero = FALSE) {
    x.na <- na.omit(as.numeric(x));
    if ('median' == method) {
        if (exclude.zero) {
            x.nonzero <- x.na[0 != x.na];
            data.median <- median(x.nonzero);
            data.mad <- mad(x.nonzero);
            }
        else {
            data.median <- median(x.na);
            data.mad <- mad(x.na);
            }
        result.na <- (x.na - data.median) / data.mad;
        x[which(!is.na(x))] <- result.na;
        x;
        }
    else if ('kmeans' == method) {
        if (exclude.zero) {
            if (1 == length(unique(x.na))) {
                kmeans.matrix <- rep(NA, length(x.na));
                names(kmeans.matrix) <- names(x.na);
                }
            else {
                data.order <- sort(x.na, decreasing = TRUE);
                non.zero <- data.order[data.order > 0];
                if (length(unique(non.zero)) <= 2) {
                    na.matrix <- rep(NA, length(non.zero));
                    cluster.zero <- c(na.matrix, rep(0, length(x.na[0 == x.na])));
                    kmeans.matrix <- cluster.zero[match(x.na, data.order)];
                    names(kmeans.matrix) <- names(x.na);
                    }
                else {
                    kmeans <- kmeans(non.zero, 2, nstart = nstart);
                    cluster <- kmeans$cluster;
                    cluster.zero <- c(cluster, rep(0, length(x[0 == x])));
                    kmeans.matrix <- cluster.zero[match(x.na, data.order)];
                    names(kmeans.matrix) <- names(x.na);
                    }
                }
            }

        else {
            if (1 == length(unique(x.na))) {
                kmeans.matrix <- rep(NA, length(x.na));
                names(kmeans.matrix) <- names(x.na);
                }
            else {
                kmeans <- kmeans(x.na, 2, nstart = nstart);
                cluster <- kmeans$cluster;
                kmeans.matrix <- cluster;
                names(kmeans.matrix) <- names(x.na);
                }
            }
        result.na <- kmeans.matrix;
        x[which(!is.na(x))] <- result.na;
        x;
        }
    else if ('mean' == method) {
        gene.order <- x.na[order(x.na, decreasing = TRUE)];
        if (exclude.zero) {
            gene.order.nonzero <- gene.order[0 != gene.order];
            top.patient <- round(
                x = length(gene.order.nonzero) * trim,
                digits = 0
                );
            low.patient <- round(
                x = length(gene.order.nonzero) * (1 - trim),
                digits = 0
                );
            data.mean <- mean(
                x = gene.order.nonzero,
                trim = trim
                );
            data.sd <- sd(gene.order.nonzero[(top.patient + 1):(low.patient)]);
            }
        else {
            top.patient <- round(
                x = length(x.na) * trim,
                digits = 0
                );
            low.patient <- round(
                x = length(x.na) * (1 - trim),
                digits = 0
                );
            data.mean <- mean(
                x = gene.order,
                trim = trim
                );
            data.sd <- sd(gene.order[(top.patient + 1):(low.patient)]);
            }
        result.na <- (x.na - data.mean) / data.sd;
        x[which(!is.na(x))] <- result.na;
        x;
        }
    }

#' Cosine similarity
#'
#' Compute the cosine similarity between pairs of values.
#'
#' @param x A number.
#' @param y A number.
#' @param large.value.percent The percentage of values to use in the calculation.  Values are selected starting with the maximums.  `large.value.percent = 0` compares only the maximums.
#'
#' @return A number.
#'
#' #' # Generate input data.
#' set.seed(1234);
#' x <- rgamma(
#'     n = 20,
#'     shape = 2,
#'     scale = 2
#'     );
#' y <- seq(
#'     from = 0.5,
#'     to = 20,
#'     by = 0.5
#'     );
#' cosine.similarity.large.value.percent(
#'     x = x,
#'     y = y,
#'     large.value.percent = 10
#'     );
#'
#' @noRd
cosine.similarity.large.value.percent <- function(x, y, large.value.percent) {
    if (0 == large.value.percent) {
        large.value.number.integer <- 1;
        }
    else {
        large.value.number <- length(x) * (large.value.percent / 100);
        large.value.number.integer <- round(large.value.number);
        }
    # Subset to the largest values.
    indices.for.comparison <- (length(x) - large.value.number.integer + 1):length(x);
    observed.value <- sort(y);
    theoretical.value <- sort(x);
    mid.value <- c(1, 1);
    value.x.y <- data.frame(theoretical.value, observed.value);

    # Calculate cosine similarity.
    cosine.large.value <- future.apply::future_sapply(
        X = indices.for.comparison,
        FUN = function(i) {
            lsa::cosine(
                x = as.numeric(value.x.y[i, ]), mid.value
                );
            }
        );
    cosine.large.value;
    }

#' Compute cosine similarity for detection of outliers
#'
#' Generate theoretical quantiles based on the optimal distribution of the data, and compute cosine similarity between the data and quantiles.
#' .
#' @param x A numeric vector.
#' @param distribution A numeric code corresponding to the optimal distribution of `x` as returned by `identify.bic.optimal.data.distribution()`.
#' @param value.portion The percentage of values to use in the calculation of cosine similarity.  Values are selected starting with the maximum.  `value.portion = 0` uses only the maximum.  The default is to use the largest 1% of values.
#'
#' @return A numeric vector containing the cosine similarity value(s) along with the distribution code.
#'
#' @examples
#' # Generate fake data.
#' set.seed(1234);
#' x <- rgamma(
#'     n = 20,
#'     shape = 2,
#'     scale = 2
#'     );
#' outlier.detection.cosine(
#'     x = x,
#'     distribution = 4,
#'     value.portion = 0
#'     );
#'
#' @noRd
outlier.detection.cosine <- function (x, distribution, value.portion = 1) {
    # Define a minimum value to ensure the values in `x` are strictly
    # positive.
    add.minimum.value <- least.significant.digit(
        x = x
        );
    x.nozero <- x + add.minimum.value;
    x.trim <- trim.sample(
        x = x,
        trim = 0.05
        );
    x.nozero.trim <- x.trim + add.minimum.value;

    # Generate the percentiles at which theoretical quantiles will be
    # computed for the optimal distribution of the data.
    p <- stats::ppoints(
        n = x.nozero
        );
    # Generate quantiles.
    if (1 == distribution) {
        norm.mean <- mean(x.nozero.trim);
        norm.sd <- stats::sd(x.nozero.trim);
        # Use truncated norm
        norm.quantiles <- truncnorm::qtruncnorm(
            p,
            a = 0,
            b = Inf,
            mean = norm.mean,
            sd = norm.sd
            );
        obs.quantile.norm <- stats::quantile(
            x = x.nozero,
            probs = p
            );
        last.cos <- cosine.similarity.large.value.percent(
            x = norm.quantiles,
            y = obs.quantile.norm,
            large.value.percent = value.portion
            );
        }
    else if (2 == distribution) {
        mean.log <- mean(x.nozero.trim);
        sd.log <- stats::sd(x.nozero.trim);
        m2 <-  log(mean.log^2 / sqrt(sd.log^2 + mean.log^2));
        sd2 <- sqrt(log(1 + (sd.log^2 / mean.log^2)));
        lnorm.quantile <- stats::qlnorm(
            p = p,
            meanlog = m2,
            sdlog = sd2
            );
        obs.quantile.lnorm <- stats::quantile(
            x = x.nozero,
            probs = p
            );
        last.cos <- cosine.similarity.large.value.percent(
            x = lnorm.quantile,
            y = obs.quantile.lnorm,
            large.value.percent = value.portion
            );
        }
    else if (3 == distribution) {
        exp.rate <- 1 / mean(x.nozero.trim);
        exp.quantile <- stats::qexp(
            p = p,
            rate = exp.rate
            );
        obs.quantile.exp <- stats::quantile(
            x = x.nozero,
            probs = p
            );
        last.cos <- cosine.similarity.large.value.percent(
            x = exp.quantile,
            y = obs.quantile.exp,
            large.value.percent = value.portion
            );
        }
    else if (4 == distribution) {
        mean.gamma <- mean(x.nozero.trim);
        sd.gamma <- stats::sd(x.nozero.trim);
        gamma.shape <- (mean.gamma / sd.gamma)^2;
        gamma.rate <- mean.gamma / (sd.gamma^2);
        gamma.quantile <- stats::qgamma(
            p = p,
            shape = gamma.shape,
            rate = gamma.rate
            );
        obs.quantile.gamma <- stats::quantile(
            x = x.nozero,
            probs = p
            );
        last.cos <- cosine.similarity.large.value.percent(
            x = gamma.quantile,
            y = obs.quantile.gamma,
            large.value.percent = value.portion
            );
        }

    cosine.sum.distribution.fit <- c(last.cos, distribution);
    cosine.sum.distribution.fit;
    }
