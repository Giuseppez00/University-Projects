clc
clear all

% Carica il dataset
load("DatasetTable.mat");

% Periodo di formazione per il momentum
formationPeriod = 252;

% Conversione "table" in "timetable" e definizione dei periodi di Backtest
Dataset = table2timetable(DatasetETF2S3);
prices = Dataset(timerange("1999-01-05","2009-12-31"), :); % In sample

numAssets = size(prices.Variables, 2);
numrows = size(prices.Variables, 1);

initialWeights=ones(1,numAssets)/numAssets;


% Warm-up partition of data set timetable.
formationPeriodTT = prices(1:formationPeriod,:);


% Inizializza l'array per memorizzare gli Sharpe Ratios
N_values = 1:6; % Range per il numero di asset da considerare
sharpe_ratios = zeros(size(N_values)); % Array per memorizzare i valori di Sharpe Ratio
avg_return_y=zeros(size(N_values));
std_y=zeros(size(N_values));



% Simulazione per ogni valore di N
for i = 1:length(N_values)
    momassets = N_values(i); % Numero di asset da considerare

    % Calcola i pesi ottimizzati usando la strategia adattiva
    %new_weights = MeanVarianceMomLong_(initialWeights, prices, momassets); % Modifica il parametro per ogni strategia
    %new_weights1 = EquallyWeightedMomentum_(initialWeights, prices, momassets);
    new_weights2 = MinVarianceMomentum_(initialWeights, prices, momassets);
    %new_weights3 = RiskParityMomentum_(initialWeights, prices, momassets);  
    %new_weights4 = VolParityMomentum_(initialWeights, prices, momassets);

    rebalancing = 21;

    % Imposta la finestra di lookback
    lookbackWindow = [126 formationPeriod*4]; 

    % Costi di transazione
    transactionsFee = 0.002; % Commissione fissa

    strategy_mom_minv=backtestStrategy('Momentum Minimum Variance', @MinVarianceMomentum, ...
            'RebalanceFrequency', rebalancing, ...
            'LookbackWindow', lookbackWindow, ...
            'TransactionCosts', transactionsFee, ...
            'InitialWeights', new_weights2);
    strategies =strategy_mom_minv

    % Crea l'engine per il backtest
    backtester = backtestEngine(strategies);
    backtester = runBacktest(backtester, prices, 'Start', formationPeriod);

    figure(i)
    equityCurve(backtester)

    % Calcola i rendimenti mensili usando la tua funzione convert2monthly
    monthly_ret = table2array(convert2monthly(backtester.Returns, "Aggregation", @compound));
    
    % Annualizza i rendimenti e la deviazione standard
    avg_return_annual = mean(monthly_ret) .* 12; % Rendimento medio annualizzato
    avg_return_y(i)=avg_return_annual;
    std_dev_annual = std(monthly_ret) .* sqrt(12); % Deviazione standard annualizzata
    std_y(i)=std_dev_annual;

    % Calcola lo Sharpe Ratio (assumendo un tasso privo di rischio di 0)
    sharpe_ratio = avg_return_annual ./ std_dev_annual; % Calcola lo Sharpe Ratio
    sharpe_ratios(i) = sharpe_ratio; % Memorizza il valore di Sharpe Ratio
end

% Trova il valore ottimo di N che massimizza lo Sharpe Ratio
[~, optimal_index] = max(sharpe_ratios); % Indice del massimo
optimal_N = N_values(optimal_index); % Ottimo N
disp(['N ottimo: ', num2str(optimal_N)]);
disp(['Sharpe Ratio con N ottimo: ', num2str(sharpe_ratios(optimal_index))]);


% Backtest delle strategie
    strategies = [
        backtestStrategy('Momentum Equally Weighted', @EquallyWeightedMomentum, ...
            'RebalanceFrequency', rebalancing, ...
            'LookbackWindow', lookbackWindow, ...
            'TransactionCosts', transactionsFee, ...
            'InitialWeights', new_weights1);
        
        backtestStrategy('Momentum Markowitz Long', @MeanVarianceMomLong, ...
            'RebalanceFrequency', rebalancing, ...
            'LookbackWindow', [126, 252*4], ...
            'TransactionCosts', transactionsFee, ...
            'InitialWeights', new_weights);
        
        backtestStrategy('Momentum Minimum Variance', @MinVarianceMomentum, ...
            'RebalanceFrequency', rebalancing, ...
            'LookbackWindow', lookbackWindow, ...
            'TransactionCosts', transactionsFee, ...
            'InitialWeights', new_weights2);
        
        backtestStrategy('Momentum Risk Parity', @RiskParityMomentum, ...
            'RebalanceFrequency', rebalancing, ...
            'LookbackWindow', lookbackWindow, ...
            'TransactionCosts', transactionsFee, ...
            'InitialWeights', new_weights3);
        
        backtestStrategy('Momentum Volatility Parity', @VolParityMomentum, ...
            'RebalanceFrequency', rebalancing, ...
            'LookbackWindow', lookbackWindow, ...
            'TransactionCosts', transactionsFee, ...
            'InitialWeights', new_weights4);

             backtestStrategy('Momentum Risk Parity', @RiskParityMomentum, ...
            'RebalanceFrequency', rebalancing, ...
            'LookbackWindow', lookbackWindow, ...
            'TransactionCosts', transactionsFee, ...
            'InitialWeights', new_weights3);
    ];
