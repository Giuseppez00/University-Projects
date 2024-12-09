clc;
clear all;

% Carica il dataset
load("DatasetTable.mat");

% Definisci i parametri
formationPeriod = 252;  % Setta il periodo di formazione per il momentum
Dataset = table2timetable(DatasetETF2S3);

% Setta il periodo in-sample
prices = Dataset(timerange("1999-01-05", "2009-12-31"), :);

numAssets = size(prices.Variables, 2);
initialWeights = ones(1, numAssets) / numAssets;

% Warm-up partition of data set timetable
formationPeriodTT = prices(1:formationPeriod, :);

% Step 2: Compute Initial Portfolio Weights
start_risk_parity = RiskParity(initialWeights, formationPeriodTT);
start_equally_weighted = EW(initialWeights, formationPeriodTT);
start_markowitz = MeanVariance_optLong(initialWeights, formationPeriodTT)';
start_min_variance = MinVariance(initialWeights, formationPeriodTT)';
start_volatility_parity = VolParity(initialWeights, formationPeriodTT);

% Step 3: Generate All Possible Combinations
strategies = {@RiskParity, @EW, @MeanVariance_optLong, @MinVariance, @VolParity};
numStrategies = length(strategies);

% Genera tutte le combinazioni di strategie
combinations = nchoosek(1:numStrategies, 2);
comb = combinations;
comb(:, [1, 2]) = comb(:, [2, 1]);
combinations = [combinations; comb];

% Step 4: Backtest Each Combination
rebalancing = 21;  % Frequenza di ribilanciamento
transactionsFee = 0.002;  %

sharpeRatios = [];  

strategy1Index = zeros(size(combinations, 1), 1);
strategy2Index = zeros(size(combinations, 1), 1);

% Ottieni gli indici delle strategie per tutte le combinazioni
for i = 1:size(combinations, 1)
    strategy1Index(i) = combinations(i, 1);
    strategy2Index(i) = combinations(i, 2);
end

% Inizializza variabili per i risultati
adaptive_strat = cell(size(combinations, 1), 1);
backtester = cell(size(combinations, 1), 1);
Returns = cell(size(combinations, 1), 1);
monthly_ret = cell(size(combinations, 1), 1);
monthly_stdev = zeros(size(combinations, 1), 1);
monthly_ExpectedRet = zeros(size(combinations, 1), 1);
sharpe_ratio = zeros(size(combinations, 1), 1);

% Ciclo su tutte le combinazioni di strategie
for i = 1:size(combinations, 1)
    % Passa le funzioni strategiche alle funzioni anonime
    strategy_func1 = strategies{strategy1Index(i)};
    strategy_func2 = strategies{strategy2Index(i)};
    
    % Calcola la strategia adattiva utilizzando AaaStrategy2
    adaptive_strat{i} = @(current_weights, prices) AaaStrategy2(current_weights, prices, ...
        strategy_func1, strategy_func2);
    
    % Esegui il backtest
    strategy_adaptive = backtestStrategy('Adaptive', adaptive_strat{i}, ...
        'RebalanceFrequency', rebalancing, ...
        'LookbackWindow', [2 * numAssets, formationPeriod], ...
        'TransactionCosts', transactionsFee, ...
        'InitialWeights', initialWeights);

    % Crea il motore di backtesting
    backtester{i} = backtestEngine(strategy_adaptive);
    backtester{i} = runBacktest(backtester{i}, prices, 'Start', formationPeriod);
end

load("backtester_opt_adaptive.mat")
%backtester_opt_adaptive=backtester
%save("backtester_opt_adaptive", "backtester_opt_adaptive")

Returns = cell(20, 1);  % Inizializza Returns come cell array




for i = 1:20
    Returns{i} = backtester_opt_adaptive{i}.Returns;  % Assegna ogni timetable a una cella
    monthly_ret{i} = convert2monthly(Returns{i}, "Aggregation", @compound);
    monthly_stdev(i) = std(monthly_ret{i}.Adaptive);
    monthly_ExpectedRet(i) = mean(monthly_ret{i}.Adaptive, 1);
end

% statistiche annualizzate
SR = (monthly_ExpectedRet.*(12))./(monthly_stdev.*(12^0.5));

%backtester_opt_adaptive{1,1}.Returns;
% Trova il massimo rapporto di Sharpe e la combinazione corrispondente
[maxSharpe, bestIndex] = max(SR);
bestCombination = combinations(bestIndex, :);
disp('Best Combination of Strategies:');
disp(bestCombination);
disp('Maximum Sharpe Ratio:');
disp(maxSharpe);

[s,b]= sort(SR);
bestCombination2 = combinations(b==2, :);
bestCombination3 = combinations(b==3, :);

