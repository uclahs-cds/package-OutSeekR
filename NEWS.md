# OutSeekR 1.0.0 - 2024-11-07 

## Features

### Core Functionality
- **Outlier Detection Algorithm**: Implements a statistical approach for detecting transcript-level outliers in RNA-seq or related data types, leveraging normalized data (e.g., FPKM) and several statistical metrics.
- **Multiple Statistics for Robust Analysis**: Utilizes five distinct statistics for each transcript to robustly assess outliers:
  - Z-scores using mean and standard deviation.
  - Z-scores using median and median absolute deviation.
  - Z-scores with 5%-trimmed mean and standard deviation.
  - Fraction of observations in the smaller cluster from K-means (K=2).
  - Cosine similarity between extreme observed values and theoretical distribution quantiles.

### Null Data Simulation
- **Comprehensive Null Simulation**: Generates null datasets mimicking the observed data distribution, without outliers, through generalized additive modeling of four potential distributions.

### Iterative Outlier Testing
- **P-value and Rank Product Calculation**: Calculates outlier p-values by comparing rank products from observed and null data across multiple rounds, refining the detection by iteratively removing the most extreme outliers.
- **FDR Adjustment**: Applies False Discovery Rate (FDR) correction to control for multiple testing.

### Performance Optimization
- **Parallel Processing Support**: Optimized for high-performance analysis using `future.apply` to enable parallelization, compatible with various computing environments.

### Output
- **Results**: Provides  p-values, FDR values, the count of detected outliers, and the preferred statistical threshold for initial selection (`p-value` or `FDR`).

### Example Data and Parameters
- **Sample Data and Example Workflow**: Includes `outliers` sample data and demonstrates usage with the `detect.outliers()` function.
- **Adjustable Parameters**: Allows customization of p-value and FDR thresholds, the number of null data iterations, and the initial screen method for flexible analysis.

