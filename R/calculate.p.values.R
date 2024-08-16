#' Calculate p-values
#'
#' Calculate p-values for each sample of a single transcript.
#'
#' @param x A numeric vector of values for an observed transcript.
#' @param x.distribution A numeric code corresponding to the optimal distribution of `x` as returned by `identify.bic.optimal.data.distribution()`.
#' @param x.zrange.mean A number, the range of the z-scores calculated using the mean and standard deviation of `x`.
#' @param x.zrange.median A number, the range of the z-scores calculated using the median and median absolute deviation of `x`.
#' @param x.zrange.trimmean A number, the range of the z-scores calculated using the trimmed mean and trimmed standard deviation of `x`.
#' @param x.fraction.kmeans A number, the k-means fraction of `x`.
#' @param x.cosine.similarity A number, the cosine similarity of `x`.
#' @param null.zrange.mean A numeric vector, the ranges of the z-scores calculated using the mean and standard deviation of each transcript in the null data.
#' @param null.zrange.median A numeric vector, the ranges of the z-scores calculated using the median and median absolute deviation of each transcript in the null data.
#' @param null.zrange.trimmean A numeric vector, the ranges of the z-scores calculated using the trimmed mean and trimmed standard deviation of each transcript in the null data.
#' @param null.fraction.kmeans A numeric vector, the k-means fraction of each transcript in the null data.
#' @param null.cosine.similarity A numeric vector, the cosine similarity of each transcript in the null data.
#' @param p.value.threshold The p-value threshold for the outlier test; default is 0.05.  Once the p-value for a sample exceeds `p.value.threshold`, testing ceases, and all samples (including the one that initially exceeded `p.value.threshold`) are assigned p-values equal to `p.value.threshold`.
#' @param kmeans.nstart The number of random starts when computing k-means fraction; default is 1.  See `?stats::kmeans` for further details.
#'
#' @return A list consisting of the following entries:
#' * `p.values`: a vector of p-values for the outlier test run on each sample (up until the p-value exceeds `p.value.threshold`); and
#' * `outlier.statistics.list`, a list of vectors containing the values of the outlier statistics calculated from the remaining samples.  The list will be of length equal to one plus the total number of outliers (i.e., the number of samples with an outlier test p-value less than `p.value.threshold`) and will contain entries `outlier.statistics.N`, where `N` is between zero and the total number of outliers.  `outlier.statistics.N` is the vector of outlier statistics after excluding the `N`th outlier sample, with `outlier.statistics.0` being for the complete transcript.
#'
#' @noRd
#' @export
calculate.p.values <- function(
    x,
    x.distribution,
    x.zrange.mean,
    x.zrange.median,
    x.zrange.trimmean,
    x.fraction.kmeans,
    x.cosine.similarity,
    null.zrange.mean,
    null.zrange.median,
    null.zrange.trimmean,
    null.fraction.kmeans,
    null.cosine.similarity,
    p.value.threshold = 0.05,
    kmeans.nstart = 1
    ) {
    p.values <- rep(NA, length(x));
    names(p.values) <- names(x);
    # Get the name and index of the most abundant sample in `x`.
    index.most.abundant.sample <- which.max(x);
    most.abundant.sample <- names(x)[index.most.abundant.sample];
    # Initialize a list to store the results of the outlier testing at
    # each round.
    results.list <- list();

    # Combine the outlier statistics from this observed transcript
    # with those of the null data.
    zrange.mean <- c(x.zrange.mean, null.zrange.mean);
    zrange.median <- c(x.zrange.median, null.zrange.median);
    zrange.trimmean <- c(x.zrange.trimmean, null.zrange.trimmean);
    fraction.kmeans <- c(x.fraction.kmeans, null.fraction.kmeans);
    cosine.similarity <- c(x.cosine.similarity, null.cosine.similarity);
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
    p.value <- (sum(rank.product[1] >= rank.product[-1]) + 1) / length(rank.product);
    p.values[most.abundant.sample] <- p.value;
    # Record results for the original `x`.
    results.list[['round1']] <- data.frame(
        sample = most.abundant.sample,
        zrange.mean = x.zrange.mean,
        zrange.median = x.zrange.median,
        zrange.trimmean = x.zrange.trimmean,
        fraction.kmeans = x.fraction.kmeans,
        cosine.similarity = x.cosine.similarity,
        rank.product = rank.product[1],
        p.value = p.value
        );

    while (p.value < p.value.threshold) {
        # Remove the most abundant sample from `x`.
        x <- x[-index.most.abundant.sample];
        # Get the name and index of the next most abundant sample.
        index.most.abundant.sample <- which.max(x);
        most.abundant.sample <- names(x)[index.most.abundant.sample];

        # Compute quantities for outlier detection in the remaining
        # samples of the observed transcript.
        x.zrange.mean <- quantify.outliers(
            x = x,
            method = 'mean'
            );
        x.zrange.median <- quantify.outliers(
            x = x,
            method = 'median'
            );
        x.zrange.trimmean <- quantify.outliers(
            x = x,
            method = 'mean',
            trim = 0.05
            );
        x.fraction.kmeans <- quantify.outliers(
            x = x,
            method = 'kmeans',
            nstart = kmeans.nstart
            );
        # Compute the ranges of the z-score statistics.
        x.zrange.mean <- zrange(
            x = x.zrange.mean
            );
        x.zrange.median <- zrange(
            x = x.zrange.median
            );
        x.zrange.trimmean <- zrange(
            x = x.zrange.trimmean
            );
        # Compute the k-means fraction.
        x.fraction.kmeans <- kmeans.fraction(
            x = x.fraction.kmeans
            );
        # Compute the cosine similarity.
        x.cosine.similarity <- outlier.detection.cosine(
            x = x,
            distribution = x.distribution
            );

        # Combine the recomputed outlier statistics from the remaining
        # samples of this observed transcript with those of the null
        # data.
        zrange.mean <- c(x.zrange.mean, null.zrange.mean);
        zrange.median <- c(x.zrange.median, null.zrange.median);
        zrange.trimmean <- c(x.zrange.trimmean, null.zrange.trimmean);
        fraction.kmeans <- c(x.fraction.kmeans, null.fraction.kmeans);
        cosine.similarity <- c(x.cosine.similarity, null.cosine.similarity);
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
        p.value <- (sum(rank.product[1] >= rank.product[-1]) + 1) / length(rank.product);
        p.values[most.abundant.sample] <- p.value;
        # Record results for the remaining samples in `x`.
        next.entry.name <- paste0(
            'round',
            length(results.list) + 1
            );
        results.list[[next.entry.name]] <- data.frame(
            sample = most.abundant.sample,
            zrange.mean = x.zrange.mean,
            zrange.median = x.zrange.median,
            zrange.trimmean = x.zrange.trimmean,
            fraction.kmeans = x.fraction.kmeans,
            cosine.similarity = x.cosine.similarity,
            rank.product = rank.product[1],
            p.value = p.value
            );
        }

    list(
        p.values = p.values,
        results.list = results.list
        );
    }
