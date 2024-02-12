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
            top.patient <- round(length(gene.order.nonzero) * trim, digit = 0);
            low.patient <- round(length(gene.order.nonzero) * (1 - trim), digit = 0);
            data.mean <- mean(gene.order.nonzero, trim = trim);
            data.sd <- sd(gene.order.nonzero[(top.patient + 1):(low.patient)]);
            }
        else {
            top.patient <- round(length(x.na) * trim, digit = 0);
            low.patient <- round(length(x.na) * (1 - trim), digit = 0);
            data.mean <- mean(gene.order, trim = trim);
            data.sd <- sd(gene.order[(top.patient + 1):(low.patient)]);
            }
        result.na <- (x.na - data.mean) / data.sd;
        x[which(!is.na(x))] <- result.na;
        x;
        }
    }
