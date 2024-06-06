#' Detect outliers
#'
#' Detect outliers in normalized RNA-seq data.
#'
#' @param data A matrix or data frame of normalized RNA-seq data, organized with transcripts on rows and samples on columns.  Transcript identifiers should be stored as `rownames(data)`.
#' @param num.null The number of transcripts to generate when simulating from null distributions.
#' @param p.value.threshold The p-value threshold for the outlier test; default is 0.05.  Once the p-value for a sample exceeds `p.value.threshold`, testing for that transcript ceases, and all remaining samples will have p-values equal to `NA`.
#' @param fdr.threshold The false discovery rate (FDR)-adjusted p-value threshold for determining the final count of outliers; default is 0.01.
#' @param kmeans.nstart The number of random starts when computing k-means fraction; default is 1.  See `?stats::kmeans` for further details.
#'
#' @return A list consisting of the following entries:
#' * `p.values`: a matrix of unadjusted p-values for the outlier test run on each sample for each transcript in `data`.
#' * `fdr`: a matrix of FDR-adjusted p-values for the outlier test run on each sample for each transcript in `data`.
#' * `num.outliers.unadjusted`: a vector giving the number of samples for which the outlier test yielded a p-value less than `p.value.threshold` for each transcript.
#' * `num.outliers.adjusted`: a vector giving the number of samples for which the outlier test yielded an FDR-adjusted p-value less than `fdr.threshold` for each transcript.
#' * `outlier.test.results.list`: a list of length `max(num.outliers.unadjusted) + 1` containing entries `roundN`, where `N` is between one and `max(num.outliers.unadjusted) + 1`.  `roundN` is the data frame of results for the outlier test after excluding the (N-1)th outlier sample, with `round1` being for the original data set (i.e., before excluding any outlier samples).  A transcript will only appear in the data frames up to the total number of outliers in that transcript plus one; for example, a transcript with one outlier will appear in `round1` and `round2`.  All transcripts appear in `round1`.
#' * `distributions`: a numeric vector indicating the optimal distribution for each transcript.  Possible values are 1 (normal), 2 (log-normal), 3 (exponential), and 4 (gamma).
#'
#' @export
detect.outliers <- function(
    data,
    num.null,
    p.value.threshold = 0.05,
    fdr.threshold = 0.01,
    kmeans.nstart = 1
    ) {
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
                p.value.threshold = p.value.threshold,
                kmeans.nstart = kmeans.nstart
                );
            },
        future.seed = TRUE
        );
    #
    # Prepare output
    #
    # Assemble a matrix of p-values for each transcript and sample.
    # Currently, the p-values are stored in the 'p.values' entry of
    # each sublist of `outlier.test.results`, with each sublist
    # corresponding to a different transcript.  We extract these
    # entries and then `rbind()` them to form a matrix.
    p.values <- do.call(
        what = rbind,
        args = lapply(
            X = outlier.test.results,
            FUN = function(x) {
                getElement(
                    object = x,
                    name = 'p.values'
                    );
                }
            )
        );
    rownames(p.values) <- rownames(data);
    # Get counts of the number of outliers per transcript based on the
    # unadjusted p-values from the outlier test.
    num.outliers.unadjusted <- apply(
        X = p.values,
        MARGIN = 1,
        FUN = function(x) sum(x < p.value.threshold, na.rm = TRUE)
        );
    # Create a list of outlier test results data frames.  (See the
    # documentation for the return value of this function for details
    # of its structure.)
    outlier.test.results.list <- list();
    for (i in 1:(max(num.outliers.unadjusted) + 1)) {
        next.name <- paste0('round', i);
        outlier.test.results.list[[next.name]] <- list();
        for (j in seq_along(outlier.test.results)) {
            outlier.test.results.list[[next.name]][[j]] <- getElement(
                object = outlier.test.results[[j]],
                name = c(
                    'results.list',
                    paste0('round', i)
                    )
                );
            }
        # Convert the list of one-row data frames to a data frame.
        outlier.test.results.list[[next.name]] <- do.call(
            what = rbind,
            args = outlier.test.results.list[[next.name]]
            );
        # Adjust the p-values from this round of results using the
        # false discovery rate (FDR).
        outlier.test.results.list[[next.name]]$fdr <- stats::p.adjust(
            p = outlier.test.results.list[[next.name]]$p.value,
            method = 'fdr'
            );
        # Add transcript name to the data frame as the row names and
        # as a dedicated column.
        outlier.test.results.list[[next.name]]$transcript <- names(num.outliers.unadjusted[num.outliers.unadjusted >= (i - 1)]);
        outlier.test.results.list[[next.name]] <- outlier.test.results.list[[next.name]][, c('transcript', colnames(outlier.test.results.list[[next.name]])['transcript' != colnames(outlier.test.results.list[[next.name]])])];
        rownames(outlier.test.results.list[[next.name]]) <- names(num.outliers.unadjusted[num.outliers.unadjusted >= (i - 1)]);
        }
    # Assemble a matrix of FDR-adjusted p-values for each transcript
    # and sample.
    fdr <- matrix(
        data = NA,
        nrow = nrow(data),
        ncol = ncol(data),
        dimnames = list(
            rownames(data),
            colnames(data)
            )
        );
    # Combine the data frames of outlier test results across rounds.
    # Keep only the transcript names (the rows of the matrix `fdr`),
    # the sample names (the columns of the matrix `fdr`), and the
    # FDR-adjusted p-values.
    all.results <- do.call(
        what = rbind,
        args = outlier.test.results.list
        );
    all.results <- all.results[, c('transcript', 'sample', 'fdr')];
    # Reshape from long to wide format.
    temp.fdr <- stats::reshape(
        data = all.results,
        direction = 'wide',
        idvar = 'transcript',
        timevar = 'sample'
        );
    rownames(temp.fdr) <- temp.fdr$transcript;
    temp.fdr$transcript <- NULL;
    # Drop the prefix 'fdr.' from the column names.
    colnames(temp.fdr) <- sub(
        pattern = 'fdr.',
        replacement = '',
        x = colnames(temp.fdr),
        fixed = TRUE
        );
    temp.fdr <- as.matrix(temp.fdr);
    # Fill in `fdr` with the FDR-adjusted p-values.
    fdr[rownames(temp.fdr), colnames(temp.fdr)] <- temp.fdr;
    # Get counts of the number of outliers per transcript based on the
    # FDR-adjusted p-values.
    num.outliers.adjusted <- apply(
        X = fdr,
        MARGIN = 1,
        FUN = function(x) sum(x < fdr.threshold, na.rm = TRUE)
        );

    list(
        p.values = p.values,
        fdr = fdr,
        num.outliers.unadjusted = num.outliers.unadjusted,
        num.outliers.adjusted = num.outliers.adjusted,
        outlier.test.results.list = outlier.test.results.list,
        distributions = optimal.distribution.data
        );
    }
