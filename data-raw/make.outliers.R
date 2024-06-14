simulate.data <- function(n, distr, parameters) {
    if (1 == distr) {
        x <- truncnorm::rtruncnorm(
            n = n,
            mean = parameters['norm.mean'],
            sd = parameters['norm.sd'],
            a = 0
            );
        }
    else if (2 == distr) {
        x <- stats::rlnorm(
            n = n,
            meanlog = parameters['m2'],
            sdlog = parameters['sd2']
            );
        }
    else if (3 == distr) {
        x <- stats::rexp(
            n = n,
            rate = parameters['exp.rate']
            );
        }
    else if (4 == distr) {
        x <- stats::rgamma(
            n = n,
            shape = parameters['gamma.shape'],
            rate = parameters['gamma.rate']
            );
        }
    x;
    }

generate.outliers <- function(data, freq, mult) {
    row.order <- apply(
        X = data,
        MARGIN = 1,
        FUN = order,
        decreasing = TRUE
        );
    row.order <- t(row.order);
    row.num.outliers <- sample(
        x = seq_along(freq) - 1,
        size = nrow(data),
        replace = TRUE,
        prob = freq
        );
    for (row in seq_len(nrow(data))) {
        for (i in seq_len(row.num.outliers[row])) {
            col <- row.order[row, ][i];
            data[row, col] <- data[row, col] * mult[i];
            }
        }
    data;
    }

set.seed(122085);

prep.outliers <- lapply(
    X = seq_len(length(outliers.distributions)),
    FUN = function(i) {
        simulate.data(
            n = 50,
            distr = outliers.distributions[i],
            parameters = outliers.params[[i]]
            );
        }
    );
prep.outliers <- do.call(
    what = 'rbind',
    args = prep.outliers
    );
rownames(prep.outliers) <- sprintf('T%03d', 1:nrow(prep.outliers));
colnames(prep.outliers) <- sprintf('S%02d', 1:ncol(prep.outliers));

outliers <- generate.outliers(
    data = prep.outliers,
    freq = c(0.50, 0.20, 0.20, 0.10),
    mult = c(2.5, 5, 10)
    );
outliers <- as.data.frame(outliers);

usethis::use_data(outliers, overwrite = TRUE);
