# Unreleased

# OutSeekR 1.0.0 - 2024-11-15

## Added
* Implementation of core *Outlier Detection Algorithm*, a statistical approach for detecting transcript-level outliers in RNA-seq or related data types, leveraging normalized data (e.g., FPKM) and several statistical metrics.
* Five distinct statistics for robustly assessing outliers:
  * Z-scores using mean and standard deviation.
  * Z-scores using median and median absolute deviation.
  * Z-scores with 5%-trimmed mean and standard deviation.
  * Fraction of observations in the smaller cluster from K-means (K=2).
  * Cosine similarity between extreme observed values and theoretical distribution quantiles.
* Comprehensive null simulation functionality. Generates null datasets mimicking the observed data distribution (without outliers) through generalized additive modeling of four potential distributions.
* Outlier p-value calculation by comparing rank products from observed and null data across multiple rounds, refining the detection by iteratively removing the most extreme outliers.
* Support for *false discovery rate* (FDR) correction to control for multiple testing.
* Optimization for high-performance analysis using `future.apply` to enable parallelization, compatible with various computing environments.
* Sample `outliers` data and usage demonstration.
