clc
clear all

%con il dataset del professore
load("DatasetTable.mat");
%load("assetClass_category.mat");

%category= categorical(Category.assetClass);


%Periodo di formazione per il momentum
%formationPeriod=84
%formationPeriod=126
formationPeriod=252

%Conversione "table" in "timetable" e definizione dei periodi di Backtest
Dataset = table2timetable(DatasetETF2S3);


prices = Dataset(timerange("1999-01-05","2017-06-30"), :); % Full sample
%prices = Dataset(timerange("1999-01-05","2009-12-31"), :); % In sample
%prices = Dataset(timerange("2010-01-01","2017-06-30"), :); % Out of sample

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

disp(sum(mom_equallyWeighted_start))
disp(sum(mom_riskParity_start))
disp(sum(mom_minVariance_start))
disp(sum(mom_markowitz_start))
disp(sum(mom_volParity_start))

% Strategia Adaptive 
start_adaptive_strategy = AaaStrategy(initialWeights, formationPeriodTT);% Compute the initial portfolio weights for each strategy.
disp("prova1")%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Come da paper, ribilanciamenti mensili, ipotizzando 252 giorni (annuali) di
%negoziazione: approx (252 / 12 = 21)
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


disp("prova2")%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Aggregate the strategy objects.

%strategies = [strategy_EW, strategy_mw_long, strategy_minvar,strategy_riskp, strategy_volp]; %strategie senza momentum
%strategies = [strategy_mom_EW, strategy_mom_mw_long, strategy_mom_minvar,strategy_mom_riskp, strategy_mom_volp]; % solo strategie momentum
%strategies=[strategy_EW, strategy_mw_long, strategy_minvar, strategy_riskp, strategy_volp, strategy_mom_EW, strategy_mom_mw_long, strategy_mom_minvar, strategy_mom_riskp, strategy_mom_volp, strategy_adaptive]; %confronto tra tutte le strategie
strategies = [strategy_mom_EW, strategy_mom_mw_long, strategy_mom_minvar,strategy_mom_riskp, strategy_mom_volp,strategy_adaptive] %momentum +adaptive

backtester = backtestEngine(strategies); % This creat a backtesting engine object

disp("prova3")
backtester = runBacktest(backtester, prices, 'Start', formationPeriod);

disp("prova4")

figure(1)
equityCurve(backtester)

%backtester_Oos=backtester
%backtester_inSample=backtester
%backtester_FullSample=backtester
%save("backtester_FullSample.mat", "backtester_FullSample" )

%save("backtester_inSample.mat", "backtester_inSample" )

%save("backtester_Oos.mat", "backtester_Oos" )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Misure di Performance %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load("backtester_FullSample.mat")
load("backtester_inSample.mat")
load("backtester_Oos.mat")

summaries = summary(backtester)

%backtester_FullSample=backtester;
%backtester_inSample=backtester;
%backtester_Oos=backtester;

num_strategy= width(backtester.Strategies);
%load("Strategies.mat")

% Store the strategy Name for each strategy in an 1*num_strategy string
% variable
for i=1:num_strategy
    Strategy_name(i) = strategies(1, i).Name;
end

%  for each strategy, retrieve the returns from the backtester object.
for i=1:num_strategy
    ret(:,i) = backtester.Returns.(Strategy_name(i));
end

% equity curves
for i=1:num_strategy
    equity_curve(:,i) = sum(table2array(backtester.Positions.(Strategy_name(i))), 2);
end

% Maximum Drawdown
Max_dd = maxdrawdown(equity_curve);

% Retrive returns from backtester object and convert it to monthly returns
% to calculate performance measures.
Returns= backtester.Returns;
Turnover = backtester.Turnover;

% Monthly performance measures
monthly_ret = convert2monthly(Returns,"Aggregation", @compound);
monthly_stdev = std(monthly_ret);
monthly_ExpectedRet = mean(monthly_ret, 1);

SR = monthly_ExpectedRet./monthly_stdev.*(12^0.5);

% CAGR
CAGR = (equity_curve(end, :)./equity_curve(1, :)).^(252/(size(equity_curve, 1)))-1

% Value at Risk (VaR) using portvrisk
VaR_5_pct = portvrisk(table2array(monthly_ExpectedRet), table2array(monthly_stdev))'*(12^0.5);

%best month
best_month = max(monthly_ret)

%worst month
worst_month = min(monthly_ret)

% Profitable Months
perc_sopra_benchmark= sum(monthly_ret>0, 1)./size(monthly_ret, 1) %[   %%%% conta quante volte la strategia > benchmark

%Monthly average turnover
monthly_avg_turnover = mean(convert2monthly(Turnover, "Aggregation", "sum"));

%Test diagnostici delle strategie
%Testare le performance (Annualizzati)
dev_std= monthly_stdev.* (12^0.5); %dev std annuale sqrt(12) .
mean_ret =monthly_ExpectedRet.*12; %media dei rendimenti mensili del vettore con tutti i rendimenti (però sharpe ratio, indicatori di rischio e rendimento vanno esplicitati annualizzati, sono piu comprensibili)

efficiency = mean_ret./dev_std; 

% stats table 
summary_stats = table();

summary_stats.Variables = CAGR; 
summary_stats.Properties.VariableNames=worst_month.Properties.VariableNames;
summary_stats(2,:) =dev_std; % monthly_STD.*sqrt(12);
summary_stats(3,:) = SR;
summary_stats(4,:) = array2table(VaR_5_pct);
summary_stats(5,:) = array2table(Max_dd);
summary_stats(6,:) = worst_month;
summary_stats(7,:) = best_month;
summary_stats(8,:) = perc_sopra_benchmark;
summary_stats(9, :) = monthly_avg_turnover;

summary_stats.Properties.RowNames=["CAGR"; "Ann. STD"; "Sharpe Ratio"; "VaR at 5%"; "Maximum Drawdown"; 
    "Worst Month Return"; "Best Month Return"; "Profitable Months"; "Monthly Average Turnover"];


%% Full Sample, In Sample & Out of Sample

monthly_ret_fullSample = convert2monthly(backtester.Returns, "Aggregation",@compound);
monthly_ExpectedRet_fullSample = mean(monthly_ret_fullSample, 1).*12;
monthly_STD_fullSample = std(monthly_ret_fullSample).*(12^0.5);
efficienctRatio_fullSample = monthly_ExpectedRet_fullSample./monthly_STD_fullSample.*(12^0.5);

monthly_ret_inSample = convert2monthly(backtester_inSample.Returns, "Aggregation",@compound);
monthly_ExpectedRet_inSample = mean(monthly_ret_inSample, 1).*12;
monthly_STD_inSample = std(monthly_ret_inSample).*(12^0.5);
efficienctRatio_inSample = monthly_ExpectedRet_inSample./monthly_STD_inSample.*(12^0.5);

monthly_ret_Oos = convert2monthly(backtester_Oos.Returns, "Aggregation",@compound);
monthly_ExpectedRet_Oos = mean(monthly_ret_Oos, 1).*12;
monthly_STD_Oos = std(monthly_ret_Oos).*(12^0.5);
efficienctRatio_outSample = monthly_ExpectedRet_Oos./monthly_STD_Oos.*(12^0.5);


%% PLOT
% Risk-Reward plot

shortName = ["EW","MWL","MINV","RP","VP","EW-M","MWL-M","MINV-M","RP-M","VP-M","AS"];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Scatter Plot only momentum strategies %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure()
col = lines(num_strategy-6);
scatter(table2array(monthly_STD_fullSample(:,6:10)), table2array(monthly_ExpectedRet_fullSample(:,6:10)), [], col, "h", "filled")
hold on
scatter(table2array(monthly_STD_inSample(:,6:10)), table2array(monthly_ExpectedRet_inSample(:,6:10)), [], col, "o", "filled")
scatter(table2array(monthly_STD_Oos(:,6:10)), table2array(monthly_ExpectedRet_Oos(:,6:10)), [], col, "^","filled")

% To avoid overlapping between data and text 
delta_y = 0.002;
delta_x = 0.001/2;  
text(table2array(monthly_STD_fullSample(:,6:10))+delta_x, table2array(monthly_ExpectedRet_fullSample(:,6:10))+delta_y, shortName(:,6:10));
text(table2array(monthly_STD_inSample(:,6:10))+delta_x, table2array(monthly_ExpectedRet_inSample(:,6:10))+delta_y, shortName(:,6:10));
text(table2array(monthly_STD_Oos(:,6:10))+delta_x, table2array(monthly_ExpectedRet_Oos(:,6:10))+delta_y, shortName(:,6:10));
legend("Full Sample", "In Sample", "Out of Sample", "Location","northeast")
xlabel("\sigma", "FontSize", 20)
ylabel("\mu", "FontSize",20, Rotation=1)
title("Risk-Reward Scatter Plot")
subtitle("Momentum Strategies")

% Set Y-axis limits
ylim([0 0.1])
hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Scatter Plot only classic strategies %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure()
col = lines(num_strategy-6);
scatter(table2array(monthly_STD_fullSample(:,1:5)), table2array(monthly_ExpectedRet_fullSample(:,1:5)), [], col, "h", "filled")
hold on
scatter(table2array(monthly_STD_inSample(:,1:5)), table2array(monthly_ExpectedRet_inSample(:,1:5)), [], col, "o", "filled")
scatter(table2array(monthly_STD_Oos(:,1:5)), table2array(monthly_ExpectedRet_Oos(:,1:5)), [], col, "^","filled")

% To avoid overlapping between data and text 
delta_y = 0.002;
delta_x = 0.001/2;  
text(table2array(monthly_STD_fullSample(:,1:5))+delta_x, table2array(monthly_ExpectedRet_fullSample(:,1:5))+delta_y, shortName(:,1:5));
text(table2array(monthly_STD_inSample(:,1:5))+delta_x, table2array(monthly_ExpectedRet_inSample(:,1:5))+delta_y, shortName(:,1:5));
text(table2array(monthly_STD_Oos(:,1:5))+delta_x, table2array(monthly_ExpectedRet_Oos(:,1:5))+delta_y, shortName(:,1:5));
legend("Full Sample", "In Sample", "Out of Sample", "Location","northeast")
xlabel("\sigma", "FontSize", 20)
ylabel("\mu", "FontSize",20, Rotation=1)
title("Risk-Reward Scatter Plot")
subtitle("Classic Weighting Schemes")

% Set Y-axis limits
ylim([0 0.1])
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Scatter Plot all Stratrgies %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure()
col = lines(num_strategy);
scatter(table2array(monthly_STD_fullSample), table2array(monthly_ExpectedRet_fullSample), [], col, "h", "filled")
hold on
scatter(table2array(monthly_STD_inSample), table2array(monthly_ExpectedRet_inSample), [], col, "o", "filled")
scatter(table2array(monthly_STD_Oos), table2array(monthly_ExpectedRet_Oos), [], col, "^","filled")

% To avoid overlapping between data and text 
delta_y = 0.002;
delta_x = 0.001/2;  
text(table2array(monthly_STD_fullSample)+delta_x, table2array(monthly_ExpectedRet_fullSample)+delta_y, shortName);
text(table2array(monthly_STD_inSample)+delta_x, table2array(monthly_ExpectedRet_inSample)+delta_y, shortName);
text(table2array(monthly_STD_Oos)+delta_x, table2array(monthly_ExpectedRet_Oos)+delta_y, shortName);
legend("Full Sample", "In Sample", "Out of Sample", "Location","northeast")
xlabel("\sigma", "FontSize", 20)
ylabel("\mu", "FontSize",20, Rotation=1)
title("Risk-Reward Scatter Plot")
subtitle("All Strategies")

% Set Y-axis limits
ylim([0 0.1])
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Scatter Plot Momentum & Adaptive Strategy %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure()
col = lines(num_strategy-5);
scatter(table2array(monthly_STD_fullSample(1,6:11)), table2array(monthly_ExpectedRet_fullSample(1,6:11)), [], col, "h", "filled")
hold on
scatter(table2array(monthly_STD_inSample(1,6:11)), table2array(monthly_ExpectedRet_inSample(1,6:11)), [], col, "o", "filled")
scatter(table2array(monthly_STD_Oos(1,6:11)), table2array(monthly_ExpectedRet_Oos(1,6:11)), [], col, "^","filled")

% To avoid overlapping between data and text 
delta_y = 0.002;
delta_x = 0.001/2;  
text(table2array(monthly_STD_fullSample(1,6:11))+delta_x, table2array(monthly_ExpectedRet_fullSample(1,6:11))+delta_y, shortName(1,6:11));
text(table2array(monthly_STD_inSample(1,6:11))+delta_x, table2array(monthly_ExpectedRet_inSample(1,6:11))+delta_y, shortName(1,6:11));
text(table2array(monthly_STD_Oos(1,6:11))+delta_x, table2array(monthly_ExpectedRet_Oos(1,6:11))+delta_y, shortName(1,6:11));
legend("Full Sample", "In Sample", "Out of Sample", "Location","northeast")
xlabel("\sigma", "FontSize", 20)
ylabel("\mu", "FontSize",20, Rotation=1)
title("Risk-Reward Scatter Plot")
subtitle("Momentum & Adaptive Strategies")

% Set Y-axis limits
ylim([0 0.1])
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Area graph of portfolios positions
date = convert2annual(backtester.Positions.Equally_Weighted, "Aggregation","lastvalue");
ew_yearly_positions = table2array(convert2annual(backtester.Positions.Equally_Weighted(:, 2:end)));
mwl_yearly_positions = table2array(convert2annual(backtester.Positions.Markowitz_Long(:, 2:end)));
mom_minv_yearly_positions = table2array(convert2annual(backtester.Positions.Momentum_Minimum_Variance(:, 2:end)));
mom_ew_yearly_positions = table2array(convert2annual(backtester.Positions.Momentum_Equally_Weighted(:, 2:end)));

mom_mwl_yearly_positions= table2array(convert2annual(backtester.Positions.Momentum_Markowitz_Long(:, 2:end)));
minv_yearly_positions= table2array(convert2annual(backtester.Positions.Minimum_Variance(:, 2:end)));
rp_yearly_positions= table2array(convert2annual(backtester.Positions.Risk_Parity(:, 2:end)));
vp_yearly_positions= table2array(convert2annual(backtester.Positions.Volatility_Parity(:, 2:end)));
mom_vp_yearly_positions= table2array(convert2annual(backtester.Positions.Momentum_Volatility_Parity(:, 2:end)));
mom_rp_yearly_positions= table2array(convert2annual(backtester.Positions.Momentum_Risk_Parity(:, 2:end)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Subplot
figure()
subplot(4, 2, 1)
area(date.Time, ew_yearly_positions./sum(ew_yearly_positions, 2));
title("Equally Weighted")
ylim([0,1])
hold off 

subplot(4, 2, 2)
area(date.Time, mom_ew_yearly_positions./sum(mom_ew_yearly_positions, 2));
title("Momentum Equally Weighted")
ylim([0,1])
hold off 

subplot(4, 2, 3)
area(date.Time, mwl_yearly_positions./sum(mwl_yearly_positions, 2));
title("Markowitz Long")
ylim([0,1])
hold off

subplot(4, 2, 4)
area(date.Time,mom_mwl_yearly_positions./sum(mom_mwl_yearly_positions, 2))
title("Momentum MWL")
ylim([0,1])
hold off

subplot(4, 2, 5)
area(date.Time, minv_yearly_positions./sum(minv_yearly_positions, 2));
title("Minimum Variance")
ylim([0,1])
hold off

subplot(4, 2, 6)
area(date.Time, mom_minv_yearly_positions./sum(mom_minv_yearly_positions, 2));
title("Momentum Minimum Variance")
ylim([0,1])
hold off

subplot(4, 2, 7)
area(date.Time,vp_yearly_positions./sum(vp_yearly_positions, 2))
title("Volatility Parity")
ylim([0,1])
hold off

subplot(4, 2, 8)
area(date.Time,mom_vp_yearly_positions./sum(mom_vp_yearly_positions, 2))
title(" Momentum Volatility Parity")
ylim([0,1])
hold off


%% Capitalized transaction costs
yearly_cost = convert2annual(backtester.BuyCost, "Aggregation","sum")+convert2annual(backtester.SellCost, "Aggregation","sum");
ew_tr_cost = cumsum(yearly_cost.Equally_Weighted);
mwl_tr_cost = cumsum(yearly_cost.Markowitz_Long);
mom_ew_tr_cost = cumsum(yearly_cost.Momentum_Equally_Weighted);
mom_minv_tr_cost = cumsum(yearly_cost.Momentum_Minimum_Variance);

vp_tr_cost = cumsum(yearly_cost.Volatility_Parity);
rp_tr_cost = cumsum(yearly_cost.Risk_Parity);
minv_tr_cost = cumsum(yearly_cost.Minimum_Variance);

mom_vp_tr_cost = cumsum(yearly_cost.Momentum_Volatility_Parity);
mom_rp_tr_cost = cumsum(yearly_cost.Momentum_Risk_Parity);
mom_mwl_tr_cost = cumsum(yearly_cost.Momentum_Markowitz_Long);
capitalized_ew_tr_cost= ew_tr_cost.*((1+CAGR(1)/2).^(14));
capitalized_mwl_tr_cost=mwl_tr_cost.*((1+CAGR(2)).^(14));
capitalized_minv_tr_cost= minv_tr_cost.*((1+CAGR(3)/2).^(14));
capitalized_rp_tr_cost= rp_tr_cost.*((1+CAGR(4)/2).^(14));
capitalized_vp_tr_cost=vp_tr_cost.*((1+CAGR(5)).^(14));
capitalized_mom_ew_tr_cost = mom_ew_tr_cost.*((1+CAGR(6)/2).^(14));
capitalized_mom_mwl_tr_cost = mom_mwl_tr_cost.*((1+CAGR(7)/2).^(14));
capitalized_mom_minv_tr_cost= mom_minv_tr_cost.*((1+CAGR(8)/2).^(14));
capitalized_mom_rp_tr_cost= mom_rp_tr_cost.*((1+CAGR(9)/2).^(14));
capitalized_mom_vp_tr_cost=mom_vp_tr_cost.*((1+CAGR(10)).^(14));




%% Equity Curve
eq_curve_MomEW = sum(mom_ew_yearly_positions, 2);
eq_curve_MomMinv = sum(mom_minv_yearly_positions, 2);
eq_curve_Mwl = sum(mwl_yearly_positions, 2);
eq_curve_EW = sum(ew_yearly_positions, 2);

eq_curve_MomMWL = sum(mom_mwl_yearly_positions, 2);
eq_curve_Minv = sum(minv_yearly_positions, 2);
eq_curve_rp = sum(rp_yearly_positions, 2);
eq_curve_vp = sum(rp_yearly_positions, 2);

eq_curve_MomVP = sum(mom_vp_yearly_positions, 2);
eq_curve_MomRP = sum(mom_rp_yearly_positions, 2);


%% Plot of the impact of transaction costs
figure()
colorOrder = lines(8); 
% Plot of capitalized transaction costs (linee trattegiate)
%yyaxis left
plot(date.Time, capitalized_ew_tr_cost, 'LineWidth', 1, 'LineStyle', '--', 'Color', colorOrder(1,:)); 
hold on
plot(date.Time, capitalized_mom_ew_tr_cost, 'LineWidth', 1, 'LineStyle', '--', 'Color', colorOrder(2,:)); 
plot(date.Time, capitalized_mwl_tr_cost, 'LineWidth', 1, 'LineStyle', '--', 'Color', colorOrder(3,:)); 
plot(date.Time, capitalized_mom_mwl_tr_cost, 'LineWidth', 1, 'LineStyle', '--', 'Color', colorOrder(4,:));
plot(date.Time, capitalized_minv_tr_cost, 'LineWidth', 1, 'LineStyle', '--', 'Color', colorOrder(5,:)); 
plot(date.Time, capitalized_mom_minv_tr_cost, 'LineWidth', 1, 'LineStyle', '--', 'Color', colorOrder(6,:)); 
plot(date.Time, capitalized_vp_tr_cost, 'LineWidth', 1, 'LineStyle', '--', 'Color', colorOrder(7,:)); 

plot(date.Time, capitalized_mom_rp_tr_cost, 'LineWidth', 1, 'LineStyle', '--', 'Color', colorOrder(8,:));

%ylabel('Capitalized tr. costs')

% Plot of equity curves (linee continue)
%yyaxis right
plot(date.Time, eq_curve_EW, 'LineWidth', 1.5, 'LineStyle', '-', 'Color', colorOrder(1,:)); 
plot(date.Time, eq_curve_MomEW, 'LineWidth', 1.5, 'LineStyle', '-', 'Color', colorOrder(2,:)); 
plot(date.Time, eq_curve_Mwl, 'LineWidth', 1.5, 'LineStyle', '-', 'Color', colorOrder(3,:)); 
plot(date.Time, eq_curve_MomMWL, 'LineWidth', 1.5, 'LineStyle', '-', 'Color', colorOrder(4,:));

plot(date.Time, eq_curve_Minv, 'LineWidth', 1.5, 'LineStyle', '-', 'Color', colorOrder(5,:)); 
plot(date.Time, eq_curve_MomMinv, 'LineWidth', 1.5, 'LineStyle', '-', 'Color', colorOrder(6,:)); 
plot(date.Time, eq_curve_vp, 'LineWidth', 1.5, 'LineStyle', '-', 'Color', colorOrder(7,:)); 
plot(date.Time, eq_curve_MomVP, 'LineWidth', 1.5, 'LineStyle', '-', 'Color', colorOrder(8,:));
%ylabel('Equity Curve')

title("Capitalized Transaction Cost vs Equity Curve")
xlabel('Time')
ylabel('Equity Curve')

% Combined legend for both transaction costs and equity curves
legend({'EW (tr)', 'Mom EW (tr)', 'MWL (tr)', 'Mom MWL (tr)', ...
        'MinVar (tr)', 'Mom MinVar (tr)', 'VolParity (tr)', 'Mom VolParity (tr)', ...
        'EW (eq)', 'Mom EW (eq)', 'MWL (eq)', 'Mom MWL (eq)', ...
        'MinVar (eq)', 'Mom MinVar (eq)', 'VolParity (eq)', 'Mom VolParity (eq)'}, ...
        'Location', 'northwest');

% Legend for strategies associated with colors (Southwest for clarity)
stratLegend = legend({'EW', 'Mom EW', 'MWL', 'Mom MWL', 'MinVar', 'Mom MinVar', 'VolParity', 'Mom VolParity'}, ...
                     'Location', 'southwest');
title(stratLegend, 'Strategy');
hold off

%% Plot of SR in different sampling
Returns_inSample= backtester_inSample.Returns;
% Monthly performance measures In sample
monthly_ret_IS = convert2monthly(Returns_inSample,"Aggregation", @compound);
monthly_stdev_IS = std(monthly_ret_IS);
monthly_ExpectedRet_IS = mean(monthly_ret_IS, 1);
SR_IS = monthly_ExpectedRet_IS./monthly_stdev_IS.*(12^0.5);


Returns_Oos= backtester_Oos.Returns;
% Monthly performance measures Out of sample
monthly_ret_Oos = convert2monthly(Returns_Oos,"Aggregation", @compound);
monthly_stdev_Oos = std(monthly_ret_Oos);
monthly_ExpectedRet_Oos = mean(monthly_ret_Oos, 1);
SR_Oos = monthly_ExpectedRet_Oos./monthly_stdev_Oos.*(12^0.5);

SR_IS = table2array(SR_IS);   
SR_Oos= table2array(SR_Oos); 
SR_FullSample = table2array(SR);         

strat = strrep(Strategy_name, '_', ' ');
barWidth = 0.25; 
%num_strategies = length(strat); % Numero totale di strategie
SR_matrix = [SR_IS; SR_Oos; SR_FullSample]; % Combina i dati in una matrice


figure();
% Creazione indice per le posizioni delle barre
x = 1:num_strategies; % Indici per le strategie
hold on; 

for i = 1:3 %In-Sample, Out-of-Sample, Full Sample
    bar(x + (i - 2) * barWidth, SR_matrix(i, :), 'BarWidth', barWidth); 
end

title('Sharpe Ratio in different sampling');
xlabel('Strategie');
ylabel('Sharpe Ratio');
xticks(x); 
xticklabels(strat); 
legend({'In-Sample', 'Out-of-Sample', 'Full Sample'}, 'Location', 'northeast');
grid on;
hold off; 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function new_weights = AaaStrategy(current_weights, prices)
    %input: Vettore contenente i pesi iniziali
    %input: Tabella di tipo timetable contenente i prezzi 
    
    %Strategia Adaptive basata sulla logica del Paper:
    % Adaptive Asset Allocation: A Primer(2013), l'idea sottostante è
    % quella di sfruttare il Momentum calcolato come rapporto di medie
    % mobili a lungo e breve periodo.
    %

    %numAssets =  size(prices.Variables,2);
    numAssets = width(prices);
    new_weights = ones(1,numAssets)/numAssets;
    pricesTable = timetable2table(prices);

    % utilizzo del trend decomposition per ottenere una
    % misura di trend di lungo periodo
    trend = trenddecomp(pricesTable(:,2:end));

    % indici_trend contiene  gli indici di trend di lungo periodo 
    indici_trend = 1:3:width(trend);

    % Determino il segno del trend di crescita
    bt_trend = table2array((trend(end, indici_trend)-trend(end-21, indici_trend))./trend(end-21, indici_trend));
    assetClass = AssetClass();
    %assetClass= category

    % If tied with 5 or more assets, EW all of them. If less than 5 assets, invest
    % with Minimum Variance logic
    if sum(bt_trend(:, assetClass~="Bond")>0) >= 5

        % Moving average (lt -> formation period, bt -> formation period/10.
        simple_ma_lt = movavg(prices, "simple", height(prices),"fill");
        simple_ma_bt = movavg(prices, "simple", round(height(prices)/10),"fill");

        % Indicatore Simple Moving Average
        moving_avg_diff = simple_ma_bt(end, :)./simple_ma_lt(end, :)-1;

        % Vettore logico dei 4 titoli migliori
        best_SMA = ismember(table2array(moving_avg_diff), maxk(table2array(moving_avg_diff), 5));
        
        % portafoglio Equally Weighted composto dai 4 migliori asset e
        % assegno 0 agli altri.
        current_weights = ones(size(best_SMA))/sum(ones(size(best_SMA)));
        temporary_weights = VolParity(current_weights, prices(:, best_SMA));
        
        new_weights(best_SMA) = temporary_weights;
        new_weights(not(best_SMA)) = 0;
    else

        current_weights = ones(width(prices))/sum(ones(width(prices)));
        new_weights = MinVariance(current_weights, prices);
    end
end


function new_weights = EquallyWeightedMomentum(~, prices)

    numAssets =  size(prices.Variables,2);
    new_weights = ones(1,numAssets)/numAssets;
    
    simple_ma_lt = movavg(prices, "simple", height(prices),"fill");
    simple_ma_bt = movavg(prices, "simple", round(height(prices)/10),"fill");
    moving_avg_diff = simple_ma_bt(end, :)./simple_ma_lt(end, :)-1;
    best_SMA = ismember(table2array(moving_avg_diff), maxk(table2array(moving_avg_diff), 5));

    current_weights = ones(size(best_SMA))/sum(ones(size(best_SMA)));
    temporary_weights = EW(current_weights, prices(:, best_SMA));

    new_weights(best_SMA) = temporary_weights;
    new_weights(not(best_SMA)) = 0;

end

function new_weights = MeanVariance_optLong(~, prices)

numAssets = width(prices);
numSample = height(prices);

sigma = robustcov(table2array(tick2ret(prices(end-90:end, :))));
mu = mean(table2array(tick2ret(prices)))';
mu0 = mean(mu);

unone = ones(numAssets, 1);

Aeq = [unone'; mu'];
beq = [1; mu0];
lower_bound = zeros(numAssets, 1);
upper_bound = ones(numAssets, 1);

new_weights = quadprog(2*sigma,[], ...
    [],[], ...
    Aeq, beq, ...
    lower_bound, upper_bound);

end

function new_weights = MeanVarianceMomLong(~, prices)
    %input: Tabella di tipo timetable contenente i prezzi 
    numAssets =  size(prices.Variables,2);
    new_weights = ones(1,numAssets)/numAssets;

    simple_ma_lt = movavg(prices, "simple", height(prices),"fill");
    simple_ma_bt = movavg(prices, "simple", round(height(prices)/10),"fill");
    moving_avg_diff = simple_ma_bt(end, :)./simple_ma_lt(end, :)-1;
    best_SMA = ismember(table2array(moving_avg_diff), maxk(table2array(moving_avg_diff), 5));

    current_weights = ones(size(best_SMA))/sum(ones(size(best_SMA)));
    temporary_weights = MeanVariance_optLong(current_weights, prices(:, best_SMA));

    new_weights(best_SMA) = temporary_weights;
    new_weights(not(best_SMA)) = 0;

end

function new_weights = MinVariance(~, prices)

    sigma = robustcov(table2array(tick2ret(prices(end-90:end, :))));
    numAssets = width(prices);
    ub = ones(numAssets, 1);
    lb = zeros(numAssets, 1);

    new_weights = fmincon(@(x)sum(x'*sigma*x) ...
        , ub/numAssets, [], [], ub', 1, lb, ub);

end

function new_weights = MinVarianceMomentum(~, prices)
    numAssets = size(prices.Variables,2);
    new_weights = ones(1,numAssets)/numAssets;

    simple_ma_lt = movavg(prices, "simple", height(prices),"fill");
    simple_ma_bt = movavg(prices, "simple", round(height(prices)/10),"fill");
    moving_avg_diff = simple_ma_bt(end, :)./simple_ma_lt(end, :)-1;
    best_SMA = ismember(table2array(moving_avg_diff), maxk(table2array(moving_avg_diff), 5));

    current_weights = ones(size(best_SMA))/sum(ones(size(best_SMA)));
    temporary_weights = MinVariance(current_weights, prices(:, best_SMA));

    new_weights(best_SMA) = temporary_weights;
    new_weights(not(best_SMA)) = 0;

end

function new_weights = RiskParity(~, prices)
    
    numAssets = size(prices.Variables,2);
    sigma = robustcov(table2array(tick2ret(prices(end-90:end, :))));

    ub = ones(numAssets, 1);
    
    %  risk budgeting portfolio è una strategia di allocazione che si
    %  concentra solo sul rischio senza considerare i rendimenti, con
    %  l'obiettivo di ottenere un'allocazione in cui ciascun asset contribuisce ad un determinato livello di rischio.
    %  the budget is set the same for all assets.
    
    new_weights = riskBudgetingPortfolio(sigma, ub/numAssets);
    new_weights = new_weights/sum(new_weights);

end

function new_weights = RiskParityMomentum(~, prices)

    numAssets =  size(prices.Variables,2);
    new_weights = ones(1,numAssets)/numAssets;
    
    simple_ma_lt = movavg(prices, "simple", height(prices),"fill");
    simple_ma_bt = movavg(prices, "simple", round(height(prices)/10),"fill");
    moving_avg_diff = simple_ma_bt(end, :)./simple_ma_lt(end, :)-1;
    best_SMA = ismember(table2array(moving_avg_diff), maxk(table2array(moving_avg_diff), 5));

    current_weights = ones(size(best_SMA))/sum(ones(size(best_SMA)));
    temporary_weights = RiskParity(current_weights, prices(:, best_SMA));


    new_weights(best_SMA) = temporary_weights;
    new_weights(not(best_SMA)) = 0;

end

function new_weights = VolParity(~, prices)
    
    returns = tick2ret(prices, "Method","continuous");
    volatility = sqrt(table2array(var(returns)));

    new_weights = (1./volatility)/sum(1./volatility);

end

function new_weights = VolParityMomentum(~, prices)

    numAssets =  size(prices.Variables,2);
    new_weights = ones(1,numAssets)/numAssets;
    
    simple_ma_lt = movavg(prices, "simple", height(prices),"fill");
    simple_ma_bt = movavg(prices, "simple", round(height(prices)/10),"fill");
    moving_avg_diff = simple_ma_bt(end, :)./simple_ma_lt(end, :)-1;
    best_SMA = ismember(table2array(moving_avg_diff), maxk(table2array(moving_avg_diff), 5));

    current_weights = ones(size(best_SMA))/sum(ones(size(best_SMA)));
    temporary_weights = VolParity(current_weights, prices(:, best_SMA));

    new_weights(best_SMA) = temporary_weights;
    new_weights(not(best_SMA)) = 0;


end

function new_weights = EW(~, prices)
    
    numAssets = size(prices.Variables,2);
    new_weights = ones(1,numAssets)/numAssets;

end
%% For strategies optimization
function new_weights = AaaStrategy2(current_weights, prices, strategy_func1, strategy_func2)
    %disp('Current Weights:');
    %disp(size(current_weights));
    %disp('Prices:');
    %disp(size(prices));

   
    
    numAssets = width(prices);
    new_weights = ones(1, numAssets) / numAssets;
    pricesTable = timetable2table(prices);

    % Trend decomposition
    trend = trenddecomp(pricesTable(:,2:end));
    
    %disp('Trend Size:');
    %disp(size(trend));  % Verifica la dimensione dell'output di trenddecomp

    indici_trend = 1:3:width(trend);
    bt_trend = table2array((trend(end, indici_trend) - trend(end - 21, indici_trend)) ./ trend(end - 21, indici_trend));

    assetClass = AssetClass();

    if sum(bt_trend(:, assetClass ~= "Bond") > 0) >= 5
        % Momentum calculation
        simple_ma_lt = movavg(prices, "simple", height(prices), "fill");
        simple_ma_bt = movavg(prices, "simple", round(height(prices) / 10), "fill");
        moving_avg_diff = simple_ma_bt(end, :) ./ simple_ma_lt(end, :) - 1;
        best_SMA = ismember(table2array(moving_avg_diff), maxk(table2array(moving_avg_diff), 5));

        %disp('Best SMA:');
        %disp(best_SMA);  % Verifica quali asset sono stati selezionati

        % Use the selected strategy for best SMA assets
        current_weights = ones(size(best_SMA))/sum(ones(size(best_SMA)));
        %disp("prova")
        temporary_weights = strategy_func1(current_weights, prices(:, best_SMA));
        %disp('prova1')
        new_weights(best_SMA) = temporary_weights;
        new_weights(~best_SMA) = 0;

        %disp('Temporary Weights after strategy_func1:');
        %disp(temporary_weights);  % Controlla i pesi temporanei generati
    else
        current_weights = ones(numAssets) / sum(ones(numAssets));
        new_weights = strategy_func2(current_weights, prices);
    end
end

function asset_class = AssetClass()
    % input: vector of char containing the type of asset class in the dataset
    % return a categorical object
    % alternativa si potrebbe impostare la funzione in modo che riceva in
    % input una stringa"sting" e restituisca come output:
    % asset_class = categorical(string);

    % Converte l'array di caratteri in una cella di stringhe
    categories = cellstr(["Stock" ,"Bond","Bond","Stock","Stock","Stock","Stock","Commodity","Bond","Stock"]);
    
    % Convert into categorical object
    asset_class = categorical(categories);

    % In alternativa passando  la variabile esternamente
    % asset_class = categorical(string);
end


function ret = compound(returns)
    %input:tabella o vettore di tipo timetable contenente i rendimenti
    ret = prod(returns+1, 1) - 1;
end
