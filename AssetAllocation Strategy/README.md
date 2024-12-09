# Asset Allocation Strategy

This folder contains projects from the **Investment Strategies** course, focusing on momentum-based asset allocation.  
The strategies and examples presented here are for informational purposes only and serve as simple illustrations of the concepts taught in the course.

## Content

- **Momentum-Based Strategies**: 
  - Projects exploring how momentum indicators can be applied to asset allocation.
  - Examples include ranking assets based on their historical performance and constructing portfolios accordingly.
  - Combining momentum Strategies with traditional allocation schemes:
     - Equal weighting.
     - Risk parity.
     - Mean-variance optimization.
     - Volatility Parity.
     - Minimum Variance.

- **Analysis and Reports**:
  - Visualizations of results and performance metrics.
 
The files required to replicate the analysis are:

ES3_final → Code

Functions for strategies: (EquallyWeightedMomentum, MeanVariance_optLong, MeanVarianceMomLong, MinVariance, MinVarianceMomentum, RiskParity, RiskParityMomentum, VolParity, VolParityMomentum)

ETF Dataset (2) → Dataset

backtester_(FullSample, InSample, Oos) → Backtester object, for performance measures (lines 196-198), in case  want to avoid minutes of runtime in the code.

ES 3 Report strategie (1) → Report

In-Sample Test Check:
find_optimal_comb → Contains the code for testing and finding the optimal combination of strategies for the adaptive strategy.
backtester_opt_adaptive → Backtester object, in case you want to avoid minutes of runtime for backtesting all possible strategy combinations and directly check (line 83).
calib_Nmomentum → File used to determine the number of assets to select in the momentum strategies.

