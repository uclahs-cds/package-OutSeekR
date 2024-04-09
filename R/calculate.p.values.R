#' Calculate p-values
#'
#' Calculate p-values for a single transcript.
#'
#' @return A list consisting of the following entries:
#' * `p.values`: a vector of p-values for the outlier test run on each sample (up until the p-value exceeds `p.value.threshold`); and
#' * `outlier.statistics.list`, a list of vectors containing the values of the outlier statistics calculated from the remaining samples.  The list will be of length equal to one plus the total number of outliers (i.e., the number of samples with an outlier test p-value less than `p.value.threshold`) and will contain entries `outlier.statistics.N`, where `N` is between zero and the total number of outliers.  `outlier.statistics.N` is the vector of outlier statistics after excluding the `N`th outlier sample, with `outlier.statistics.0` being for the complete transcript.
#'
#' @noRd
calculate.p.values <- function(
    x,
    x.distribution,
    null.data,
    null.distributions,
    p.value.threshold = 0.05,
    kmeans.nstart = 1
    ) {
    p.values <- rep(NA, length(x));
    names(p.values) <- names(x);
    # Get the name and index of the most abundant sample in `x`.
    index.most.abundant.sample <- which.max(x);
    most.abundant.sample <- names(x)[index.most.abundant.sample];
    # Initialize a list to store the outlier statistics for `x` as
    # well as subsequent copies of `x` as each outlier is excluded.
    outlier.statistics.list <- list();

    data <- rbind(x, null.data);
    # Compute quantities for outlier detection: (1) z-scores based on
    # the mean / standard deviation, (2) z-scores based on the trimmed
    # mean / trimmed standard deviation, (3) z-scores based on the
    # median / median absolute deviation, and (4) the cluster
    # assignment from k-means with two clusters.
    zrange.mean <- future.apply::future_apply(
        X = data,
        MARGIN = 1,
        FUN = quantify.outliers,
        method = 'mean'
        );
    zrange.median <- future.apply::future_apply(
        X = data,
        MARGIN = 1,
        FUN = quantify.outliers,
        method = 'median'
        );
    zrange.trimmean <- future.apply::future_apply(
        X = data,
        MARGIN = 1,
        FUN = quantify.outliers,
        method = 'mean',
        trim = 0.05
        );
    fraction.kmeans <- future.apply::future_apply(
        X = data,
        MARGIN = 1,
        FUN = quantify.outliers,
        method = 'kmeans',
        nstart = kmeans.nstart,
        future.seed = TRUE
        );
    # Compute the ranges of the z-score statistics.
    zrange.mean <- future.apply::future_apply(
        X = zrange.mean,
        MARGIN = 2,
        FUN = zrange
        );
    zrange.median <- future.apply::future_apply(
        X = zrange.median,
        MARGIN = 2,
        FUN = zrange
        );
    zrange.trimmean <- future.apply::future_apply(
        X = zrange.trimmean,
        MARGIN = 2,
        FUN = zrange
        );
    # Compute the k-means fraction.
    fraction.kmeans <- future.apply::future_apply(
        X = fraction.kmeans,
        MARGIN = 2,
        FUN = kmeans.fraction
        );
    # Compute the cosine similarity.
    cosine.similarity <- future.apply::future_sapply(
        X = seq_len(nrow(data)),
        FUN = function(i) {
            if (1 == i) {
                outlier.detection.cosine(
                    x = x,
                    distribution = x.distribution
                    );
                }
            else {
                outlier.detection.cosine(
                    x = data[i, ],
                    distribution = null.distributions[i - 1]
                    );
                }
            }
        );
    # Record the outlier statistics for `x`.
    outlier.statistics.list[['outlier.statistics.0']] <- c(
        zrange.mean[1],
        zrange.median[1],
        zrange.trimmean[1],
        fraction.kmeans[1],
        cosine.similarity[1]
        );
    names(outlier.statistics.list[['outlier.statistics.0']]) <- c(
        'zrange.mean',
        'zrange.median',
        'zrange.trimmean',
        'fraction.kmeans',
        'cosine.similarity'
        );
    # Assign ranks within each method.
    rank.zrange.mean <- outlier.rank(
        outlier.statistic = zrange.mean,
        method = 'zrange.mean'
        );
    rank.zrange.median <- outlier.rank(
        outlier.statistic = zrange.median,
        method = 'zrange.median'
        );
    rank.zrange.trimmean <- outlier.rank(
        outlier.statistic = zrange.trimmean,
        method = 'zrange.trimmean'
        );
    rank.fraction.kmeans <- outlier.rank(
        outlier.statistic = fraction.kmeans,
        method = 'fraction.kmeans'
        );
    rank.cosine.similarity <- outlier.rank(
        outlier.statistic = cosine.similarity,
        method = 'cosine.similarity'
        );
    # Compute the rank product for each transcript.
    rank.product <- outlier.rank.product(
        zrange.mean = rank.zrange.mean,
        zrange.median = rank.zrange.median,
        zrange.trimmean = rank.zrange.trimmean,
        fraction.kmeans = rank.fraction.kmeans,
        cosine.similarity = rank.cosine.similarity
        );
    # Calculate the p-value associated with the most abundant sample
    # in `x`.
    p.value <- (sum(rank.product[1] >= rank.product[-1]) + 1) / nrow(data);
    p.values[most.abundant.sample] <- p.value;

    while (p.value < p.value.threshold) {
        # Remove the most abundant sample from `x`.
        x <- x[-index.most.abundant.sample];
        # Get the name and index of the next most abundant sample.
        index.most.abundant.sample <- which.max(x);
        most.abundant.sample <- names(x)[index.most.abundant.sample];
        # Remove a random column from `null.data`.
        column.to.remove <- sample(
            x = ncol(null.data),
            size = 1
            );
        null.data <- null.data[, -column.to.remove];

        data <- rbind(x, null.data);
        zrange.mean <- future.apply::future_apply(
            X = data,
            MARGIN = 1,
            FUN = quantify.outliers,
            method = 'mean'
            );
        zrange.median <- future.apply::future_apply(
            X = data,
            MARGIN = 1,
            FUN = quantify.outliers,
            method = 'median'
            );
        zrange.trimmean <- future.apply::future_apply(
            X = data,
            MARGIN = 1,
            FUN = quantify.outliers,
            method = 'mean',
            trim = 0.05
            );
        fraction.kmeans <- future.apply::future_apply(
            X = data,
            MARGIN = 1,
            FUN = quantify.outliers,
            method = 'kmeans',
            nstart = kmeans.nstart,
            future.seed = TRUE
            );
        # Compute the ranges of the z-score statistics.
        zrange.mean <- future.apply::future_apply(
            X = zrange.mean,
            MARGIN = 2,
            FUN = zrange
            );
        zrange.median <- future.apply::future_apply(
            X = zrange.median,
            MARGIN = 2,
            FUN = zrange
            );
        zrange.trimmean <- future.apply::future_apply(
            X = zrange.trimmean,
            MARGIN = 2,
            FUN = zrange
            );
        # Compute the k-means fraction.
        fraction.kmeans <- future.apply::future_apply(
            X = fraction.kmeans,
            MARGIN = 2,
            FUN = kmeans.fraction
            );
        # Compute the cosine similarity.
        cosine.similarity <- future.apply::future_sapply(
            X = seq_len(nrow(data)),
            FUN = function(i) {
                if (1 == i) {
                    outlier.detection.cosine(
                        x = x,
                        distribution = x.distribution
                        );
                    }
                else {
                    outlier.detection.cosine(
                        x = data[i, ],
                        distribution = null.distributions[i - 1]
                        );
                    }
                }
            );
        # Record the outlier statistics for `x`.
        next.entry.name <- paste0(
            'outlier.statistics.',
            length(outlier.statistics.list)
            );
        outlier.statistics.list[[next.entry.name]] <- c(
            zrange.mean[1],
            zrange.median[1],
            zrange.trimmean[1],
            fraction.kmeans[1],
            cosine.similarity[1]
            );
        names(outlier.statistics.list[[next.entry.name]]) <- c(
            'zrange.mean',
            'zrange.median',
            'zrange.trimmean',
            'fraction.kmeans',
            'cosine.similarity'
            );
        # Assign ranks within each method.
        rank.zrange.mean <- outlier.rank(
            outlier.statistic = zrange.mean,
            method = 'zrange.mean'
            );
        rank.zrange.median <- outlier.rank(
            outlier.statistic = zrange.median,
            method = 'zrange.median'
            );
        rank.zrange.trimmean <- outlier.rank(
            outlier.statistic = zrange.trimmean,
            method = 'zrange.trimmean'
            );
        rank.fraction.kmeans <- outlier.rank(
            outlier.statistic = fraction.kmeans,
            method = 'fraction.kmeans'
            );
        rank.cosine.similarity <- outlier.rank(
            outlier.statistic = cosine.similarity,
            method = 'cosine.similarity'
            );
        # Compute the rank product for each transcript.
        rank.product <- outlier.rank.product(
            zrange.mean = rank.zrange.mean,
            zrange.median = rank.zrange.median,
            zrange.trimmean = rank.zrange.trimmean,
            fraction.kmeans = rank.fraction.kmeans,
            cosine.similarity = rank.cosine.similarity
            );
        # Calculate the p-value associated with the most abundant sample
        # in `x`.
        p.value <- (sum(rank.product[1] >= rank.product[-1]) + 1) / nrow(data);
        p.values[most.abundant.sample] <- p.value;
        }

    # After the while loop terminates, all remaining p-values are
    # greater than `p.value.threshold`.  Replace any missing values or
    # values greater than `p.value.threshold` in `p.values` with
    # `p.value.threshold`.
    p.values[is.na(p.values) | p.values > p.value.threshold] <- p.value.threshold;

    list(
        p.values = p.values,
        outlier.statistics.list = outlier.statistics.list
        );
    }
