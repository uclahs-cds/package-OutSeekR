#' Detect outliers
#'
#' Detect outliers in normalized RNA-seq data.
#'
#' @param data A matrix or data frame of normalized RNA-seq data, organized with transcripts on rows and samples on columns.  Transcript identifiers should be stored as `rownames(data)`.
#' @param num.null The number of transcripts to generate when simulating from null distributions; default is one million.
#' @param initial.screen.method The statistical criterion for initial gene selection; valid options are 'p-value' and 'FDR'.
#' @param p.value.threshold The p-value threshold for the outlier test; default is 0.05.  Once the p-value for a sample exceeds `p.value.threshold`, testing for that transcript ceases, and all remaining samples will have p-values equal to `NA`.
#' @param fdr.threshold The false discovery rate (FDR)-adjusted p-value threshold for determining the final count of outliers; default is 0.01.
#' @param kmeans.nstart The number of random starts when computing k-means fraction; default is 1.  See `?stats::kmeans` for further details.
#'
#' @return A list consisting of the following entries:
#' * `p.values`: a matrix of unadjusted p-values for the outlier test run on each transcript in `data`.
#' * `fdr`: a matrix of FDR-adjusted p-values for the outlier test run on each transcript in `data`.
#' * `num.outliers`: a vector giving the number of outliers detected for each transcript based on the threshold.
#' * `outlier.test.results.list`: a list of length `max(num.outliers) + 1` containing entries `roundN`, where `N` is between one and `max(num.outliers) + 1`.  `roundN` is the data frame of results for the outlier test after excluding the (N-1)th outlier sample, with `round1` being for the original data set (i.e., before excluding any outlier samples).
#' * `distributions`: a numeric vector indicating the optimal distribution for each transcript.  Possible values are 1 (normal), 2 (log-normal), 3 (exponential), and 4 (gamma).
#' * `initial.screen.method`: Specifies the statistical criterion for initial feature selection. Valid options are 'p-value' and 'FDR'.
#' @export
detect.outliers <- function(
        data,
        num.null = 1e6,
        initial.screen.method = c('p.value', 'fdr'),
        p.value.threshold = 0.05,
        fdr.threshold = 0.01,
        kmeans.nstart = 1
    ) {
    # Match the initial screening method
    initial.screen.method <- match.arg(initial.screen.method);

    # Determine which of the normal, log-normal, exponential, or gamma
    # distributions provides the best fit to each row of values in
    # `data`.
    optimal.distribution.data <- future.apply::future_apply(
        X = data,
        MARGIN = 1,
        FUN = identify.bic.optimal.data.distribution,
        future.seed = TRUE
        );
    # Calculate residuals of the observed data with respect to the
    # optimal distribution.  (We use `as.numeric()` on the input to
    # `calculate.residuals()` in order to handle data frame input.
    observed.residuals <- future.apply::future_lapply(
        X = seq_len(nrow(data)),
        FUN = function(i) {
            calculate.residuals(
                x = as.numeric(data[i, ]),
                distribution = optimal.distribution.data[i]
                );
            }
        );
    observed.residuals <- do.call(
        what = rbind,
        args = observed.residuals
        );
    rownames(observed.residuals) <- rownames(data);
    # Apply 5% trimming to each row of residuals.
    observed.residuals <- future.apply::future_apply(
        X = observed.residuals,
        MARGIN = 1,
        FUN = trim.sample
        );
    observed.residuals <- t(observed.residuals);
    # Determine which of the normal, log-normal, exponential, or gamma
    # distributions provides the best fit to each row of residuals.
    optimal.distribution.residuals <- future.apply::future_apply(
        X = observed.residuals,
        MARGIN = 1,
        FUN = identify.bic.optimal.residuals.distribution
        );
    # Compute quantities for outlier detection: (1) z-scores based on
    # the mean / standard deviation, (2) z-scores based on the trimmed
    # mean / trimmed standard deviation, (3) z-scores based on the
    # median / median absolute deviation, and (4) the cluster
    # assignment from k-means with two clusters.
    data.zrange.mean <- future.apply::future_apply(
        X = data,
        MARGIN = 1,
        FUN = quantify.outliers,
        method = 'mean'
        );
    data.zrange.median <- future.apply::future_apply(
        X = data,
        MARGIN = 1,
        FUN = quantify.outliers,
        method = 'median'
        );
    data.zrange.trimmean <- future.apply::future_apply(
        X = data,
        MARGIN = 1,
        FUN = quantify.outliers,
        method = 'mean',
        trim = 0.05
        );
    data.fraction.kmeans <- future.apply::future_apply(
        X = data,
        MARGIN = 1,
        FUN = quantify.outliers,
        method = 'kmeans',
        nstart = kmeans.nstart,
        future.seed = TRUE
        );
    # Compute the ranges of the z-score statistics.
    data.zrange.mean <- future.apply::future_apply(
        X = data.zrange.mean,
        MARGIN = 2,
        FUN = zrange
        );
    data.zrange.median <- future.apply::future_apply(
        X = data.zrange.median,
        MARGIN = 2,
        FUN = zrange
        );
    data.zrange.trimmean <- future.apply::future_apply(
        X = data.zrange.trimmean,
        MARGIN = 2,
        FUN = zrange
        );
    # Compute the k-means fraction.
    data.fraction.kmeans <- future.apply::future_apply(
        X = data.fraction.kmeans,
        MARGIN = 2,
        FUN = kmeans.fraction
        );
    # Compute the cosine similarity.
    data.cosine.similarity <- future.apply::future_sapply(
        X = seq_len(nrow(data)),
        FUN = function(i) {
            outlier.detection.cosine(
                x = as.numeric(data[i, ]),
                distribution = optimal.distribution.data[i]
                );
            }
        );
    # Generate a matrix of null transcripts by simulating from their
    # respective optimal distributions.
    sampled.indices <- sample(
        x = nrow(data),
        size = num.null,
        replace = TRUE
        );
    null.data <- future.apply::future_lapply(
        X = sampled.indices,
        FUN = function(i) {
            simulate.null(
                x = as.numeric(data[i, ]),
                x.distribution = optimal.distribution.data[i],
                r = as.numeric(observed.residuals[i, ]),
                r.distribution = optimal.distribution.residuals[i]
                );
            },
        future.seed = TRUE
        );
    null.data <- do.call(
        what = rbind,
        args = null.data
        );
    rownames(null.data) <- rownames(data)[sampled.indices];
    # Determine which of the normal, log-normal, exponential, or gamma
    # distributions provides the best fit to each row of values in
    # `null.data`.
    optimal.distribution.null.data <- future.apply::future_apply(
        X = null.data,
        MARGIN = 1,
        FUN = identify.bic.optimal.data.distribution,
        future.seed = TRUE
        );
    null.zrange.mean <- future.apply::future_apply(
        X = null.data,
        MARGIN = 1,
        FUN = quantify.outliers,
        method = 'mean'
        );
    null.zrange.median <- future.apply::future_apply(
        X = null.data,
        MARGIN = 1,
        FUN = quantify.outliers,
        method = 'median'
        );
    null.zrange.trimmean <- future.apply::future_apply(
        X = null.data,
        MARGIN = 1,
        FUN = quantify.outliers,
        method = 'mean',
        trim = 0.05
        );
    null.fraction.kmeans <- future.apply::future_apply(
        X = null.data,
        MARGIN = 1,
        FUN = quantify.outliers,
        method = 'kmeans',
        nstart = kmeans.nstart,
        future.seed = TRUE
        );
    # Compute the ranges of the z-score statistics.
    null.zrange.mean <- future.apply::future_apply(
        X = null.zrange.mean,
        MARGIN = 2,
        FUN = zrange
        );
    null.zrange.median <- future.apply::future_apply(
        X = null.zrange.median,
        MARGIN = 2,
        FUN = zrange
        );
    null.zrange.trimmean <- future.apply::future_apply(
        X = null.zrange.trimmean,
        MARGIN = 2,
        FUN = zrange
        );
    # Compute the k-means fraction.
    null.fraction.kmeans <- future.apply::future_apply(
        X = null.fraction.kmeans,
        MARGIN = 2,
        FUN = kmeans.fraction
        );
    # Compute the cosine similarity.
    null.cosine.similarity <- future.apply::future_sapply(
        X = seq_len(nrow(null.data)),
        FUN = function(i) {
            outlier.detection.cosine(
                x = as.numeric(null.data[i, ]),
                distribution = optimal.distribution.null.data[i]
                );
            }
        );
    # Calculate p-values.  The result is a list of length equal to
    # `nrow(data)`, with each sublist containing the results of
    # `calculate.p.values()` for a transcript in the observed data.
    # See the documentation for `calculate.p.values()` for a
    # description of its return value.
    outlier.test.results <- future.apply::future_lapply(
        X = seq_len(nrow(data)),
        FUN = function(i) {
            x <- as.numeric(data[i, ]);
            names(x) <- colnames(data);
            calculate.p.values(
                x = x,
                x.distribution = optimal.distribution.data[i],
                x.zrange.mean = data.zrange.mean[i],
                x.zrange.median = data.zrange.median[i],
                x.zrange.trimmean = data.zrange.trimmean[i],
                x.fraction.kmeans = data.fraction.kmeans[i],
                x.cosine.similarity = data.cosine.similarity[i],
                null.zrange.mean = null.zrange.mean,
                null.zrange.median = null.zrange.median,
                null.zrange.trimmean = null.zrange.trimmean,
                null.fraction.kmeans = null.fraction.kmeans,
                null.cosine.similarity = null.cosine.similarity,
                kmeans.nstart = kmeans.nstart
                );
            },
        future.seed = TRUE
        );
    outlier.test.results.list <- list();
    k <- 1;
    # Record results for the remaining samples in `x`.
    next.name <- paste0(
        'round',
        k
        );
    # Combine results from all genes
    outlier.test.results.list[[next.name]] <- do.call(rbind, outlier.test.results)
    # Calculate FDR
    outlier.test.results.list[[next.name]]$fdr <- stats::p.adjust(
        p = outlier.test.results.list[[next.name]]$p.value,
        method = 'fdr'
        );

    if (initial.screen.method == 'p.value') {
        while (sum(stats::na.omit(outlier.test.results.list[[next.name]]$p.value) < p.value.threshold) > 0) {
            k <- k + 1;
            next.name <- paste0('round', k);
            outlier.test.results.iter <- future.apply::future_lapply(
                X = seq_len(nrow(data)),
                FUN = function(i) {
                    x <- as.numeric(data[i, ]);
                    names(x) <- colnames(data);
                    # Sort x in decreasing order
                    sorted.indices <- order(x, decreasing = TRUE);
                    # Remove the top k-1 values (k-1 because we're keeping the k-th value and below)
                    x <- x[sorted.indices[k:length(x)]];
                    # The most abundant sample in this round is the k-th highest overall
                    most.abundant.sample <- names(x)[1];
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
                        distribution = optimal.distribution.data[i]
                        );
                    calculate.p.values(
                        x = x,
                        x.distribution = optimal.distribution.data[i],
                        x.zrange.mean = x.zrange.mean,
                        x.zrange.median = x.zrange.median,
                        x.zrange.trimmean = x.zrange.trimmean,
                        x.fraction.kmeans = x.fraction.kmeans,
                        x.cosine.similarity = x.cosine.similarity,
                        null.zrange.mean = null.zrange.mean,
                        null.zrange.median = null.zrange.median,
                        null.zrange.trimmean = null.zrange.trimmean,
                        null.fraction.kmeans = null.fraction.kmeans,
                        null.cosine.similarity = null.cosine.similarity,
                        kmeans.nstart = kmeans.nstart
                        );
                    },
                future.seed = TRUE
                )
            # Combine results from all genes
            outlier.test.results.list[[next.name]] <- do.call(rbind, outlier.test.results.iter);
            # Calculate FDR
            outlier.test.results.list[[next.name]]$fdr <- stats::p.adjust(
                p = outlier.test.results.list[[next.name]]$p.value,
                method = 'fdr'
                );
            }
        }
    else {
        while (sum(stats::na.omit(outlier.test.results.list[[next.name]]$fdr) < fdr.threshold) > 0) {
            k <- k + 1;
            next.name <- paste0('round', k);
            outlier.test.results.iter <- future.apply::future_lapply(
                X = seq_len(nrow(data)),
                FUN = function(i) {
                    x <- as.numeric(data[i, ]);
                    names(x) <- colnames(data);
                    # Sort x in decreasing order
                    sorted.indices <- order(x, decreasing = TRUE);
                    # Remove the top k-1 values (k-1 because we're keeping the k-th value and below)
                    x <- x[sorted.indices[k:length(x)]];
                    # The most abundant sample in this round is the k-th highest overall
                    most.abundant.sample <- names(x)[1];
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
                        distribution = optimal.distribution.data[i]
                        );
                    calculate.p.values(
                        x = x,
                        x.distribution = optimal.distribution.data[i],
                        x.zrange.mean = x.zrange.mean,
                        x.zrange.median = x.zrange.median,
                        x.zrange.trimmean = x.zrange.trimmean,
                        x.fraction.kmeans = x.fraction.kmeans,
                        x.cosine.similarity = x.cosine.similarity,
                        null.zrange.mean = null.zrange.mean,
                        null.zrange.median = null.zrange.median,
                        null.zrange.trimmean = null.zrange.trimmean,
                        null.fraction.kmeans = null.fraction.kmeans,
                        null.cosine.similarity = null.cosine.similarity,
                        kmeans.nstart = kmeans.nstart
                        );
                    },
                future.seed = TRUE
                );
            # Combine results from all genes
            outlier.test.results.list[[next.name]] <- do.call(rbind, outlier.test.results.iter);
            # Calculate FDR
            outlier.test.results.list[[next.name]]$fdr <- stats::p.adjust(
                p = outlier.test.results.list[[next.name]]$p.value,
                method = 'fdr'
                );
            }
        }
    # Give rownames
    for (i in 1:length(outlier.test.results.list)) {
        rownames(outlier.test.results.list[[i]]) <- rownames(data);
        }
    #
    # Prepare output
    #
    # Assemble a matrix of p-values for each transcript and sample.
    # Currently, the p-values are stored in the 'p.values' entry of
    # each sublist of `outlier.test.results`, with each sublist
    # corresponding to a different transcript.  We extract these
    # entries and then `rbind()` them to form a matrix.
    p.values <- do.call(
        what = cbind,
        args = lapply(
            X = outlier.test.results.list,
            FUN = function(x) {
                getElement(
                    object = x,
                    name = 'p.value'
                    );
                }
            )
        );
    rownames(p.values) <- rownames(data);
    # Assemble a matrix of FDR-adjusted p-values for each transcript
    # and sample.
    fdr <- do.call(
        what = cbind,
        args = lapply(
            X = outlier.test.results.list,
            FUN = function(x) {
                getElement(
                    object = x,
                    name = 'fdr'
                    );
                }
            )
        );
    rownames(fdr) <- rownames(data);
    # # Get counts of the number of outliers per transcript based on the
    # # FDR-adjusted p-values.
    # num.outliers.adjusted <- apply(
    #     X = fdr,
    #     MARGIN = 1,
    #     FUN = function(x) sum(x < fdr.threshold, na.rm = TRUE)
    # );
    # Calculate outliers based on the chosen method
    if (initial.screen.method == 'p.value') {
        num.outliers <- apply(
            X = p.values,
            MARGIN = 1,
            FUN = function(x) sum(x < p.value.threshold, na.rm = TRUE)
            )
        }
    else {  # initial.screen.method == "fdr"
        num.outliers <- apply(
            X = fdr,
            MARGIN = 1,
            FUN = function(x) sum(x < fdr.threshold, na.rm = TRUE)
            );
        }
    list(
        p.values = p.values,
        fdr = fdr,
        num.outliers = num.outliers,
        outlier.test.results.list = outlier.test.results.list,
        distributions = optimal.distribution.data,
        initial.screen.method = initial.screen.method
        );
    }
