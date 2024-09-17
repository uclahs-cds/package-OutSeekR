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
#' @export
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
quantify.outliers <- function(x, method = 'mean', trim = 0, nstart = 1, exclude.zero = FALSE) {
    x.na <- stats::na.omit(x);
    if ('median' == method) {
        if (exclude.zero) {
            x.nonzero <- x.na[0 != x.na];
            data.median <- stats::median(x.nonzero);
            data.mad <- stats::mad(x.nonzero);
            }
        else {
            data.median <- stats::median(x.na);
            data.mad <- stats::mad(x.na);
            }
        if (0 == data.mad || is.na(data.mad) || is.na(data.median)) {
            return(rep(NA, length(x)));
            }
        result.na <- (x.na - data.median) / data.mad;
        x[which(!is.na(x))] <- result.na;
        return(x);
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
                    kmeans <- stats::kmeans(
                        x = non.zero,
                        centers = 2,
                        nstart = nstart
                        );
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
                kmeans <- stats::kmeans(
                    x = x.na,
                    centers = 2,
                    nstart = nstart
                    );
                cluster <- kmeans$cluster;
                kmeans.matrix <- cluster;
                names(kmeans.matrix) <- names(x.na);
                }
            }
        result.na <- kmeans.matrix;
        x[which(!is.na(x))] <- result.na;
        return(x);
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
            data.sd <- stats::sd(gene.order.nonzero[(top.patient + 1):(low.patient)]);
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
            data.sd <- stats::sd(gene.order[(top.patient + 1):(low.patient)]);
            }
        if (0 == data.sd || is.na(data.sd) || is.na(data.mean)) {
            return(rep(NA, length(x)));
            }
        result.na <- (x.na - data.mean) / data.sd;
        x[which(!is.na(x))] <- result.na;
        return(x);
        }
    }

#' Cosine similarity
#'
#' Compute cosine similarity for detection of outliers.  Generate theoretical quantiles based on the optimal distribution of the data, and compute cosine similarity between a point made up of the largest observed quantile and the largest theoretical quantile and a point on the line y = x.
#' .
#' @param x A numeric vector.
#' @param distribution A numeric code corresponding to the optimal distribution of `x` as returned by `identify.bic.optimal.data.distribution()`.
#'
#' @return A number.
#' @export
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
#'     distribution = 4
#'     );
outlier.detection.cosine <- function(x, distribution) {
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
        # For a normal distribution, generate quantiles from a
        # truncated normal distribution.
        theoretical.quantiles <- truncnorm::qtruncnorm(
            p = p,
            a = 0,
            b = Inf,
            mean = norm.mean,
            sd = norm.sd
            );
        observed.quantiles <- stats::quantile(
            x = x.nozero,
            probs = p
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
        observed.quantiles <- stats::quantile(
            x = x.nozero,
            probs = p
            );
        }
    else if (3 == distribution) {
        exp.rate <- 1 / mean(x.nozero.trim);
        theoretical.quantiles <- stats::qexp(
            p = p,
            rate = exp.rate
            );
        observed.quantiles <- stats::quantile(
            x = x.nozero,
            probs = p
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
        observed.quantiles <- stats::quantile(
            x = x.nozero,
            probs = p
            );
        }
    cosine.similarity <- lsa::cosine(
        x = c(theoretical.quantiles[length(p)], observed.quantiles[length(p)]),
        y = c(1, 1)
        );
    cosine.similarity[1, 1];
    }
