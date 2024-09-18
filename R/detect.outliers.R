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
#' * `initial.screen.method`: Specifies the statistical criterion for initial feature selection. Valid options are 'p-value' and 'FDR' (p-value used by default).
#' @export
#' @examples
#' data(outliers);
#' outliers.subset <- outliers[1:10,];
#' results <- detect.outliers(
#'    data = outliers.subset,
#'    num.null = 10
#'    );
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
    obs.data <- future.apply::future_lapply(
        X = seq_len(nrow(data)),
        FUN = function(i) {
            x <- as.numeric(data[i, ]);
            zrange.mean <- zrange(
                quantify.outliers(
                    x,
                    method = 'mean'
                    )
                );
            zrange.median <- zrange(
                quantify.outliers(
                    x,
                    method = 'median'
                    )
                );
            zrange.trimmean <- zrange(
                quantify.outliers(
                    x,
                    method = 'mean',
                    trim = 0.05
                    )
                );
            fraction.kmeans <- kmeans.fraction(
                quantify.outliers(
                    x,
                    method = 'kmeans',
                    nstart = kmeans.nstart
                    )
                );
            cosine.similarity <- outlier.detection.cosine(
                x = x,
                distribution = optimal.distribution.data[i]
                );
            list(
                zrange.mean = zrange.mean,
                zrange.median = zrange.median,
                zrange.trimmean = zrange.trimmean,
                fraction.kmeans = fraction.kmeans,
                cosine.similarity = cosine.similarity
                );
            },
        future.seed = TRUE
        );
    data.zrange.mean <- future.apply::future_sapply(
        X = obs.data,
        FUN = function(x) {
            x$zrange.mean
            }
        );
    data.zrange.median <- future.apply::future_sapply(
        X = obs.data,
        FUN = function(x) {
            x$zrange.median
            }
        );
    data.zrange.trimmean <- future.apply::future_sapply(
        X = obs.data,
        FUN = function(x) {
            x$zrange.trimmean
            }
        );
    data.fraction.kmeans <- future.apply::future_sapply(
        X = obs.data,
        FUN = function(x) {
            x$fraction.kmeans
            }
        );
    data.cosine.similarity <- future.apply::future_sapply(
        X = obs.data,
        FUN = function(x) {
            x$cosine.similarity
            }
        );
    names(data.zrange.mean) <- rownames(data);
    names(data.zrange.median) <- rownames(data);
    names(data.zrange.trimmean) <- rownames(data);
    names(data.fraction.kmeans) <- rownames(data);
    names(data.cosine.similarity) <- rownames(data);
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
            null.row <- simulate.null(
                x = as.numeric(data[i, ]),
                x.distribution = optimal.distribution.data[i],
                r = as.numeric(observed.residuals[i, ]),
                r.distribution = optimal.distribution.residuals[i]
                );
            # Determine which of the normal, log-normal, exponential, or gamma
            # distributions provides the best fit to each row of values in
            # `null.data`.
            optimal.distribution.null.data <- identify.bic.optimal.data.distribution(null.row);
            zrange.mean <- zrange(
                quantify.outliers(
                    null.row,
                    method = 'mean'
                    )
                );
            zrange.median <- zrange(
                quantify.outliers(
                    null.row,
                    method = 'median'
                    )
                );
            zrange.trimmean <- zrange(
                quantify.outliers(
                    null.row,
                    method = 'mean',
                    trim = 0.05
                    )
                );
            fraction.kmeans <- kmeans.fraction(
                quantify.outliers(
                    null.row,
                    method = 'kmeans',
                    nstart = kmeans.nstart
                    )
                );
            cosine.similarity <- outlier.detection.cosine(
                x = as.numeric(null.row),
                distribution = optimal.distribution.null.data
                );
            list(
                zrange.mean = zrange.mean,
                zrange.median = zrange.median,
                zrange.trimmean = zrange.trimmean,
                fraction.kmeans = fraction.kmeans,
                cosine.similarity = cosine.similarity
                );
            },
        future.seed = TRUE
        );
    null.zrange.mean <- future.apply::future_sapply(
        X = null.data,
        FUN = function(x) {
            x$zrange.mean
            }
        );
    null.zrange.median <- future.apply::future_sapply(
        X = null.data,
        FUN = function(x) {
            x$zrange.median
            }
        );
    null.zrange.trimmean <- future.apply::future_sapply(
        X = null.data,
        FUN = function(x) {
            x$zrange.trimmean
            }
        );
    null.fraction.kmeans <- future.apply::future_sapply(
        X = null.data,
        FUN = function(x) {
            x$fraction.kmeans
            }
        );
    null.cosine.similarity <- future.apply::future_sapply(
        X = null.data,
        FUN = function(x) {
            x$cosine.similarity
            }
        );
    names(null.zrange.mean) <- rownames(data)[sampled.indices];
    names(null.zrange.median) <- rownames(data)[sampled.indices];
    names(null.zrange.trimmean) <- rownames(data)[sampled.indices];
    names(null.fraction.kmeans) <- rownames(data)[sampled.indices];
    names(null.cosine.similarity) <- rownames(data)[sampled.indices];
    # Calculate p-values.  The result is a list of length equal to
    # `nrow(data)`, with each sublist containing the results of
    # `calculate.p.values()` for a transcript in the observed data.
    # See the documentation for `calculate.p.values()` for a
    # description of its return value.
    # ### save example data for testing calcualte.p.values
    # browser();
    # example.data.for.calculate.p.values <- list(
    #     data = data,
    #     x.distribution = optimal.distribution.data,
    #     x.zrange.mean = data.zrange.mean,
    #     x.zrange.median = data.zrange.median,
    #     x.zrange.trimmean = data.zrange.trimmean,
    #     x.fraction.kmeans = data.fraction.kmeans,
    #     x.cosine.similarity = data.cosine.similarity,
    #     null.zrange.mean = null.zrange.mean,
    #     null.zrange.median = null.zrange.median,
    #     null.zrange.trimmean = null.zrange.trimmean,
    #     null.fraction.kmeans = null.fraction.kmeans,
    #     null.cosine.similarity = null.cosine.similarity,
    #     kmeans.nstart = kmeans.nstart
    #     );
    # usethis::use_data(example.data.for.calculate.p.values);
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
                    x.zrange.mean <- zrange(
                        quantify.outliers(
                            x = x,
                            method = 'mean'
                            )
                        );
                    x.zrange.median <- zrange(
                        quantify.outliers(
                            x = x,
                            method = 'median'
                            )
                        );
                    x.zrange.trimmean <- zrange(
                        quantify.outliers(
                            x = x,
                            method = 'mean',
                            trim = 0.05
                            )
                        );
                    x.fraction.kmeans <- kmeans.fraction(
                        quantify.outliers(
                            x = x,
                            method = 'kmeans',
                            nstart = kmeans.nstart
                            )
                        );
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
                    x.zrange.mean <- zrange(
                        quantify.outliers(
                            x = x,
                            method = 'mean'
                            )
                        );
                    x.zrange.median <- zrange(
                        quantify.outliers(
                            x = x,
                            method = 'median'
                            )
                        );
                    x.zrange.trimmean <- zrange(
                        quantify.outliers(
                            x = x,
                            method = 'mean',
                            trim = 0.05
                            )
                        );
                    x.fraction.kmeans <- kmeans.fraction(
                        quantify.outliers(
                            x = x,
                            method = 'kmeans',
                            nstart = kmeans.nstart
                            )
                        );
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
