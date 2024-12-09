clc
clear all

%% Test for optimal loockback
%con il dataset del professore
load("DatasetTable.mat");
%load("assetClass_category.mat");

%category= categorical(Category.assetClass);


%Periodo di formazione per il momentum
%formationPeriod=84
formationPeriod=126
%formationPeriod=252

%Conversione "table" in "timetable" e definizione dei periodi di Backtest
Dataset = table2timetable(DatasetETF2S3);


%prices = Dataset(timerange("1999-01-05","2017-06-30"), :); % Full sample
%prices = Dataset(timerange("1999-01-05","2009-12-31"), :); % In sample
prices = Dataset(timerange("2010-01-01","2018-01-01"), :); % Out of sample



disp("prova")

numAssets=size(prices.Variables,2);
numrows=size(prices.Variables,1);

disp("prova")

%Partiamo dall'ipotesi di base (EW)
initialWeights=ones(1,numAssets)/numAssets;


% Warm-up partition of data set timetable.
formationPeriodTT = prices(1:formationPeriod,:);

% Compute the initial portfolio weights for each strategy.
% Classical allocation strategies
start_risk_parity = RiskParity(initialWeights,formationPeriodTT);
start_equally_weighted = EW(initialWeights,formationPeriodTT);
start_markowitz = MeanVariance_optLong(initialWeights,formationPeriodTT);
start_markowitz=start_markowitz';

start_min_variance = MinVariance(initialWeights,formationPeriodTT);
start_min_variance=start_min_variance';

start_volatility_parity = VolParity(initialWeights,formationPeriodTT);

disp(sum(start_risk_parity))
disp(sum(start_equally_weighted))
disp(sum(start_markowitz))
disp(sum(start_min_variance))
disp(sum(start_volatility_parity))

% Strategie Momentum 
mom_equallyWeighted_start = EquallyWeightedMomentum(initialWeights,formationPeriodTT);
mom_riskParity_start = RiskParityMomentum(initialWeights,formationPeriodTT);
mom_minVariance_start = MinVarianceMomentum(initialWeights,formationPeriodTT);
mom_markowitz_start = MeanVarianceMomLong(initialWeights,formationPeriodTT);
mom_volParity_start = VolParityMomentum(initialWeights,formationPeriodTT);


% Strategia Adaptive 
start_adaptive_strategy = AaaStrategy(initialWeights, formationPeriodTT);% Compute the initial portfolio weights for each strategy.
disp("prova1")%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Come da paper, ribilanciamenti mensili, ipotizzando 252 giorni (annuali) di
%negoziazioni: approx (252 / 12 = 21)
rebalancing= 21;

% Set the rolling lookback window to be at most 252 days (coherenthly with the
% period of necessarly data needed for the definition of momentum.
% The minimum is set in a way to avoid an ill conditionet Var-Cov matrix.
% we have choose to set it as a proportional number to the number of assets
lookbackWindow  = [2*numAssets formationPeriod]; 
%lookback  = [126 252*4]; 


transactionsFee = 0.002; 

%Per la strategia EW abbiamo impostato la finestra di lookback a 0 in quanto non richiede informazioni pregresse per determinare lo schema di allocazione 

strategy_EW = backtestStrategy('Equally Weighted', @EW, ...
    'RebalanceFrequency', rebalancing, ...
    'LookbackWindow', 0, ...
    'TransactionCosts', transactionsFee, ...
    'InitialWeights', start_equally_weighted);

strategy_mw_long = backtestStrategy('Markowitz Long', @MeanVariance_optLong, ...
    'RebalanceFrequency', rebalancing, ...
    'LookbackWindow', [formationPeriod, 252*4], ...
    'TransactionCosts', transactionsFee, ...
    'InitialWeights', start_markowitz);

strategy_minvar = backtestStrategy('Minimum Variance', @MinVariance, ...
    'RebalanceFrequency', rebalancing, ...
    'LookbackWindow', lookbackWindow, ...
    'TransactionCosts', transactionsFee, ...
    'InitialWeights', start_min_variance);

strategy_riskp = backtestStrategy('Risk Parity', @RiskParity, ...
    'RebalanceFrequency', rebalancing, ...
    'LookbackWindow', lookbackWindow, ...
    'TransactionCosts', transactionsFee, ...
    'InitialWeights', start_risk_parity);

strategy_volp = backtestStrategy('Volatility Parity', @VolParity, ...
    'RebalanceFrequency', rebalancing, ...
    'LookbackWindow', lookbackWindow, ...
    'TransactionCosts', transactionsFee, ...
    'InitialWeights', start_volatility_parity);



strategy_mom_EW = backtestStrategy('Momentum Equally Weighted', @EquallyWeightedMomentum, ...
    'RebalanceFrequency', rebalancing, ...
    'LookbackWindow', lookbackWindow, ...
    'TransactionCosts', transactionsFee, ...
    'InitialWeights', mom_equallyWeighted_start);

strategy_mom_mw_long = backtestStrategy('Momentum Markowitz Long', @MeanVarianceMomLong, ...
    'RebalanceFrequency', rebalancing, ...
    'LookbackWindow', [formationPeriod, 252*4], ...
    'TransactionCosts', transactionsFee, ...
    'InitialWeights', mom_markowitz_start);

strategy_mom_minvar = backtestStrategy('Momentum Minimum Variance', @MinVarianceMomentum, ...
    'RebalanceFrequency', rebalancing, ...
    'LookbackWindow', lookbackWindow, ...
    'TransactionCosts', transactionsFee, ...
    'InitialWeights', mom_minVariance_start);


strategy_mom_riskp = backtestStrategy('Momentum Risk Parity', @RiskParityMomentum, ...
    'RebalanceFrequency', rebalancing, ...
    'LookbackWindow', lookbackWindow, ...
    'TransactionCosts', transactionsFee, ...
    'InitialWeights', mom_riskParity_start);


strategy_mom_volp = backtestStrategy('Momentum Volatility Parity', @VolParityMomentum, ...
    'RebalanceFrequency', rebalancing, ...
    'LookbackWindow', lookbackWindow, ...
    'TransactionCosts', transactionsFee, ...
    'InitialWeights', mom_volParity_start);

strategy_adaptive = backtestStrategy('Adaptive Strategy', @AaaStrategy, ...
    'RebalanceFrequency', rebalancing, ...
    'LookbackWindow', lookbackWindow, ...
    'TransactionCosts', transactionsFee, ...
    'InitialWeights', start_adaptive_strategy);

% Aggregate the strategy objects into an array.
%strategies = [strat_mom_vp, strat_mom_ew, strat_as];
%strategies = [strat_mom_vp, strat_mom_ew strat_as];
strategies=[strategy_EW, strategy_mw_long, strategy_minvar, strategy_riskp, strategy_volp, strategy_mom_EW, strategy_mom_mw_long, strategy_mom_minvar, strategy_mom_riskp, strategy_mom_volp, strategy_adaptive];

%backtester = backtestEngine(strategies); % This creat a backtesting engine object

disp("prova3")
%backtester = runBacktest(backtester, prices, 'Start', formationPeriod);

backtester_Oos = backtestEngine(strategies)
backtester_Oos = runBacktest(backtester_Oos, prices, 'Start', formationPeriod);
disp("prova4")

figure(1)
equityCurve(backtester_Oos)

%save("backtester_FullSample.mat", "backtester_FullSample" )

%save("backtester_inSample.mat", "backtester_inSample" )

%save("backtester_Oos.mat", "backtester_Oos" )