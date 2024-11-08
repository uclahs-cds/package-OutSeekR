# OutSeekR

1. [Introduction](#instroduction)
2. [Methodology](#methodology)
3. [Installation](#installation)
4. [Parallelization](#parallelization)
5. [Output](#output)
6. [Example](#example)
7. [Discussions](#discussions)
8. [License](#license)


## Introduction

![OutSeekR](https://github.com/user-attachments/assets/3569bdb4-f03a-4616-9e79-e577aa51c1d0)

The **OutSeekR** package implements a statistical approach to detecting outliers in RNA-seq and related data.

---

## Methodology

An analysis using **OutSeekR** should start with a matrix or data frame of normalized RNA-seq data, e.g., fragments per kilobase of transcript per million fragments mapped (FPKM), or similar data; in the case of RNA-seq data, we require normalized values rather than the counts themselves because the statistics we calculate assume continuous data. For simplicity, we'll refer to the rows of the input data set as 'transcripts' and the columns of the input data set as 'samples'. Transcript identifiers should be stored as the row names of the matrix or data frame, and sample identifiers should be stored as the column names.

The statistical approach of **OutSeekR** centers around the use of five statistics calculated for each transcript in the observed data:

- 1) the range of standard scores (z-scores) computed using the mean and standard deviation;
- 2) the range of z-scores computed using the median and median absolute deviation (MAD);
- 3) the range of z-scores computed using the 5%-trimmed mean and 5%-trimmed standard deviation;
- 4) the fraction of observations assigned to the smaller cluster based on K-means clustering with K = 2 clusters; and
- 5) the cosine similarity between the most extreme observed value and the largest quantile of a representative theoretical distribution (see [Simulating null data]).

Specifically, it uses the five statistics calculated on the observed data and compares the distributions of these statistics with counterparts calculated using simulated null data. Observed data yielding statistics more extreme than those of the null data suggest the presence of outliers.

Note: Input data with low variation (e.g., where 50% or more of the values are identical) may cause issues in calculating the MAD or performing K-means clustering. It is recommended to exclude such data from the analysis to ensure reliable results.
**OutSeekR** first identifies which transcripts/features contain outlier values.  For a given transcript, the null hypothesis is there are no outlier values (*i.e.* all values come from some unknown population distribution) while the alternative hypothesis is that the transcript contains 1 or more outlier values from another distribution.  A small p-value or FDR-value gives evidence for the presence of 1 or more outliers.

After determining which transcripts contain outliers, an iterative procedure is used to identify the exact number of outliers each transcript has.

---

## Installation

To install the latest public release of OutSeekR from CRAN:

```
install.packages("OutSeekR")
```

Or to install the latest development version from GitHub:
```
# install.packages("devtools")

devtools::install_github("uclahs-cds/package-OutSeekR")
```

---

## Parallelization in **OutSeekR**

The statistical approach to outlier detection implemented by **OutSeekR** is time-consuming. To reduce runtime, the code in **OutSeekR** has been written to allow parallelization. Parallelization is supported through the use of the package [**future.apply**](https://future.apply.futureverse.org/). **future.apply** is built on top of the future framework implemented by the [**future**](https://future.futureverse.org/) package, which provides a uniform way to parallelize code independent of operating system or computing environment. In particular, **OutSeekR** uses the `future_*apply()` functions defined in **future.apply** rather than their base R equivalents in order to take advantage of a parallelization strategy set with `future::plan()` by the user.

**Note:** Depending on the size of the input data set and the number of null transcripts requested, the objects created by **OutSeekR** may exceed a limit defined in **future** to prevent exporting very large objects, which will trigger an error with a message referring to `future.globals.maxSize`. This can be avoided by setting `future.globals.maxSize` to a large value, e.g., `options(future.globals.maxSize = 8000 * 1024^2);  # = 8 GB`. See [the documentation of `future.options`](https://search.r-project.org/CRAN/refmans/future/html/future.options.html) for further details.

---

## Output

The output of the main function in **OutSeekR**, `detect.outliers()`, is a named list with the following components:

- `p.values`: a matrix of unadjusted p-values for the outlier test run on each transcript in the observed data.
- `fdr`: a matrix of FDR-adjusted p-values for the outlier test run on each transcript in the observed data.
- `num.outliers`: a vector giving the number of outliers detected for each transcript based on the threshold.
- `outlier.test.results.list`: a list of length `max(num.outliers) + 1` containing entries `roundN`, where `N` is between one and `max(num.outliers) + 1`.  `roundN` is the data frame of results for the outlier test after excluding the $(N-1)$th outlier sample, with `round1` being for the original data set (i.e., before excluding any outlier samples).
- `distributions`: a numeric vector indicating the optimal distribution for each transcript. Possible values are 1 (normal), 2 (log-normal), 3 (exponential), and 4 (gamma).
- `initial.screen.method`: Specifies the statistical criterion for initial feature selection. Valid options are 'p-value' and 'FDR'.

---

## Example

**OutSeekR** includes a sample data set called `outliers` that we shall use as an example. It is a data frame consisting of 50 samples across 500 transcripts. Note that transcript identifiers and sample identifiers are stored as `rownames(outliers)` and `colnames(outliers)`, respectively.

```{r view-data}
str(outliers, list.len = 5);
outliers[1:6, 1:6];
outliers[495:500, 45:50];
```

An analysis in **OutSeekR** is run using the function `detect.outliers()`. With default settings, it can be run by simply passing the input data set via the `data` parameter. Other parameters include the following:

- `num.null`: the number of transcripts to generate when simulating from null distributions; default is 1000.
- `initial.screen.method`: the statistical criterion for initial gene selection.
- `p.value.threshold`: the p-value threshold for the outlier test; default is 0.05.
- `fdr.threshold`: the false discovery rate (FDR)-adjusted p-value threshold for determining the final count of outliers; default is 0.01.
- `kmeans.nstart`: the number of random starts when computing k-means fraction; default is 1. See `?stats::kmeans` for further details.

In our first example, we pass the `outliers` data frame and use the defaults for most other parameters. To keep runtime manageable, we limit the number of null transcripts generated to 1,000. Prior to running `detect.outliers()`, we set up parallel processing. See [Parallelization in **OutSeekR**] for further details.

```{r run-1}
# Set random seed for reproducibility.
set.seed(371892);

# Set up parallel processing.
future::plan(future::multisession);

outlier.test.run.1 <- detect.outliers(
    data = outliers,
    num.null = 1e3
    );

str(outlier.test.run.1, max.level = 2);

# Restore sequential processing.
future::plan(future::sequential);
```

Unadjusted p-values for the outlier test on each transcript are contained in the `p.values` entry of the object returned by `detect.outliers()`; FDR-adjusted p-values are contained in the `fdr` entry. 

```{r examine-p-values}
head(outlier.test.run.1$p.values);
head(outlier.test.run.1$fdr);
```

As explained in [Overview], an iterative procedure is used to determine the exact number of outliers per transcript.  The iterations stop after the p-value or FDR-value exceeds the specified threshold.

The entries `num.outliers` give transcript-specific counts of the number of outliers based on the threshold of specified `initial.screen.method`.

```{r transcript-level-outlier-counts}
head(outlier.test.run.1$num.outliers);
```

We can get a quick view of the distribution of the number of outliers across transcripts using `table()`:

```{r table}
table(outlier.test.run.1$num.outliers);
```

More granular results can be found in the `outlier.test.results.list` entry of the object returned by `detect.outliers()`. It is a named list of data frames, where each reflects the results of a different 'round' of outlier testing; by 'round', we mean the process of excluding the most outlying sample (for all rounds except the first), calculating the five outlier statistics, ranking them alongside the statistics of the null data, computing the rank product, and computing the unadjusted p-value. Each data frame in `outlier.test.results.list` contains the transcript and sample identifiers, the five statistics calculated on the observed transcript, the rank product, and the unadjusted and FDR-adjusted p-value for the outlier test.

```{r view-outlier-test-results-list}
str(outlier.test.run.1$outlier.test.results.list);
```

It may be of use to combine the results across all rounds of outlier testing into a single data frame. One way to do this (while retaining the round information) using `lapply()` and `rbind()` is shown below:

```{r collapse-rounds}
outlier.test.results.combined <- lapply(
    X = seq_along(outlier.test.run.1$outlier.test.results.list),
    FUN = function(i) {
        df <- outlier.test.run.1$outlier.test.results.list[[i]];
        df$round <- i;
        df <- df[, c(
            'round',
            colnames(outlier.test.run.1$outlier.test.results.list[[i]])
            )];
        }
    );
outlier.test.results.combined <- do.call(
    what = 'rbind',
    args = outlier.test.results.combined
    );
# Combining the data frames produces duplicates in the row names.  R
# will de-duplicate them, but as all the necessary information is
# included in the columns of the data frame (specifically, 'round' and
# 'transcript'), we'll simply discard the row names.
rownames(outlier.test.results.combined) <- NULL;
head(outlier.test.results.combined);
```

We can observe the optimal theoretical distribution for each transcript and the frequencies across all transcripts.

```{r distributions}
head(outlier.test.run.1$distributions);
table(outlier.test.run.1$distributions);
```

The results are numeric codes, with 1 corresponding to the normal distribution, 2 the log-normal distribution, 3 the exponential distribution, and 4 the gamma distribution.

Lastly, we show a second run of the algorithm where we adjust the p-value and FDR thresholds:

```{r run-2}
# Set up parallel processing.
future::plan(future::multisession);

outlier.test.run.2 <- detect.outliers(
    data = outliers,
    num.null = 1e3,
    initial.screen.method = 'fdr',
    p.value.threshold = 0.25,
    fdr.threshold = 0.05
    );

# Restore sequential processing.
future::plan(future::sequential);

str(outlier.test.run.2, max.level = 2);

# Examine p-value and FDR matrices.
head(outlier.test.run.2$p.values);
head(outlier.test.run.2$fdr);

# Check the distribution of number of outliers detected.
table(outlier.test.run.2$num.outliers);
```

---

## Discussions

- [Issue tracker](https://github.com/uclahs-cds/package-OutSeekR/issues) to report errors and enhancement ideas.
- Discussions can take place in [package-OutSeekR Discussions](https://github.com/uclahs-cds/package-OutSeekR/discussions)
- [package-OutSeekR pull requests](https://github.com/uclahs-cds/package-OutSeekR/pulls) are also open for discussion

---

## License

Author: Jee Yun Han(jyhan@mednet.ucla.edu), John M Sahrmann(jsahrmann@mednet.ucla.edu)

This project is licensed under the GNU General Public License version 2. See the file LICENSE.md for the terms of the GNU GPL license.

Copyright (C) 2024 University of California Los Angeles ("Boutros Lab") All rights reserved.

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
