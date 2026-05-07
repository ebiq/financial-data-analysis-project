# Financial Data Analysis Project
# Kalman Filter Based Dynamic Beta Modeling for Financial Time Series

## Project Overview

This project investigates whether dynamic beta models can better capture time-varying market risk sensitivity than traditional static beta models. Based on the CAPM framework, I model the relationship between industry excess returns and market excess returns, with a focus on estimating time-varying beta using state-space models and Kalman filtering.

The project compares static OLS, rolling OLS, and three Kalman filter-based models to evaluate their performance in fitting, prediction, and cross-sectional ranking of industry risk sensitivity.

## Research Motivation

Traditional CAPM assumes that beta is constant over time. However, financial markets often exhibit structural changes, volatility clustering, non-normal returns, and time-varying risk exposure. These characteristics suggest that a static beta may fail to accurately capture the dynamic relationship between industry returns and market returns.

Through exploratory analysis and diagnostic testing, the data shows evidence of:

- Fat-tailed return distributions
- Volatility clustering
- Time-varying variance
- Strong but potentially unstable industry-market correlations
- Parameter instability indicated by CUSUMSQ tests

These findings motivate the use of dynamic beta models based on state-space modeling and Kalman filtering.

## Data

The project uses weekly excess return data for 8 pan-European industries from 2010 to 2024.

The raw data is not included in this repository due to data access restrictions. The repository focuses on the modeling pipeline, methodology, and reproducible code structure.

## Methodology

The project compares the following models:

### 1. Static OLS

A baseline CAPM regression model that assumes beta is constant over the full sample period.

### 2. Rolling OLS

A rolling-window regression model using a fixed 52-week window to estimate time-varying beta.

### 3. Kalman Filter Models

Three state-space models are implemented:

- **Random Walk Model (RW)**  
  Assumes beta evolves as a random walk over time.

- **Mean-Reverting Model (MR)**  
  Assumes beta fluctuates around a long-term mean with mean-reversion behavior.

- **Moving Mean-Reverting Model (MMR)**  
  Extends the MR model by allowing the long-term mean itself to evolve over time, providing greater flexibility under structural changes.

## Technical Implementation

The project was implemented in Python using the following tools:

- `Pandas` for data cleaning, time alignment, and structured data management
- `NumPy` for matrix operations and vectorized numerical computation
- `Statsmodels` for OLS regression and statistical diagnostics
- `Scikit-learn` for baseline regression modeling
- `pykalman` for Kalman filtering implementation
- `SciPy` for likelihood optimization using BFGS
- `Matplotlib` for visualization of beta paths, diagnostics, and model comparison results

## Workflow

The project follows the workflow below:

1. Data loading and cleaning
2. Weekly return and excess return calculation
3. Exploratory data analysis
4. Statistical diagnostics
5. Static OLS estimation
6. Rolling OLS estimation with a 52-week window
7. Kalman filter model implementation
8. Maximum likelihood estimation of model parameters
9. In-sample and out-of-sample forecasting
10. Model evaluation and visualization

## Model Evaluation

The models are evaluated using multiple metrics:

- **MAE**: Measures average prediction error and is relatively robust to extreme values.
- **MSE**: Penalizes large errors more heavily and is useful for identifying poor performance during volatile periods.
- **Spearman Rank Correlation**: Evaluates cross-sectional ranking consistency, which is relevant for industry allocation and risk comparison.
- **Beta Path Visualization**: Compares the smoothness, stability, and responsiveness of estimated beta paths.
- **Boxplots**: Used to assess the distribution and dispersion of beta estimates across models.

## Key Findings

The main findings are:

- Dynamic Kalman filter models generally outperform static OLS and rolling OLS in residual behavior, model fit, and forecasting performance.
- The MR and MMR models show the strongest overall performance, especially during periods of elevated volatility.
- Static beta can still serve as a long-term benchmark, but it is less adaptive during short-term market fluctuations and structural changes.
- Modeling time-varying beta provides a more flexible framework for understanding dynamic market risk sensitivity.
- The state-space modeling framework can potentially be extended to other time series applications, such as electricity price forecasting, risk monitoring, and industrial data modeling.

