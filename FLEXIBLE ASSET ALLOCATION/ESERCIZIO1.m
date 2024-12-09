clc
clear all

%% Full Sample
%load("FSDatasetD_price");
%load("FSDatasetD_delta");
%load("FSDataset_M_price");
%load("FSDataset_M_delta");

%m_delta = FSDataset_M_delta(:,2:end).Variables'./100;
%m_price = FSDataset_M_price(:,2:end).Variables';

%d_price = FSDatasetD_price(:,2:end).Variables';
%d_delta = FSDatasetD_delta(:,2:end).Variables'./100;

%% In-Sample load Book1 from sep-2005 to 31-dec-2012
load("Book1S3");
load("Book1S4");
load("Book1S5");
load("Book1S6");

m_delta = Book1S6(:,2:end).Variables'./100;
m_price = Book1S5(:,2:end).Variables';

d_price = Book1S3(:,2:end).Variables';
d_delta = Book1S4(:,2:end).Variables'./100;


%% Out of Sample

%load("OosDataset");
%load("OosDatasetD_delta");
%load("OosDataset_M_price");
%load("OosDataset_M_delta");

%m_delta = OosDataset_M_delta(:,2:end).Variables'./100;
%m_price = OosDataset_M_price(:,2:end).Variables';

%d_price = OosDataset(:,3:end).Variables';
%d_delta = OosDatasetD_delta(:,3:end).Variables'./100;

% Monthly variable
num = 3;
[num_rows, num_cols] = size(m_delta);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Relative momentum
relative_period = 4;
R_selected_assets = zeros(3, num_cols-relative_period);

% Otteniamo il numero di righe e colonne
[num_rows, num_cols] = size(m_delta);

% Creiamo una nuova matrice per salvare i risultati delle divisioni
result_matrix= zeros(num_rows, num_cols-relative_period);  % Usiamo NaN per le celle che non possono essere divise


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Itera su ogni riga
for i = 1:num_rows
    % Itera su ogni colonna per calcolare i rendimenti cumulati
    for j = relative_period + 1:num_cols
        % Estrai i rendimenti per il periodo
        returns = m_delta(i, j-relative_period:j-1);
        % Calcola i rendimenti cumulati
        cumulative_returns = ret2tick(returns');
        % Estrai l'ultimo rendimento cumulato
        last_cumulative_return = cumulative_returns(end); 
        
        result_matrix(i, j-relative_period) = last_cumulative_return;
    end
end

% Sorting dei rendimenti cumulati (per ciascuna colonna) in ordine decrescente
[rm, idx_rel_mom] = sort(result_matrix, 1, 'descend');  

% Ora selezioniamo i primi 3 asset per ciascun periodo
%
m_delta1 = m_delta(:, 5:end);  

% Selezioniamo i rendimenti dei primi 3 asset ordinati per ciascun mese (colonna)
ret = zeros(3, size(m_delta1, 2));  % Matrice per i rendimenti selezionati
for k = 1:size(m_delta1, 2)
    % Seleziona i rendimenti dei primi 3 asset (righe in idx_rel_mom) per ciascuna colonna k
    selected_assets = idx_rel_mom(1:3, k);  
    ret(:, k) = m_delta1(selected_assets, k);
end

% Calcola i rendimenti medi dei primi 3 asset selezionati per ciascun mese
strategy_ret = mean(ret, 1);

% Calcola i rendimenti cumulati della strategia
R_cumulative = ret2tick(strategy_ret');
%R_cumulative = ret2tick(strategy_ret(2:end)'); %Oos

%strategy_ret= mean(m_delta(R_selected_assets));
%R_cumulative = ret2tick(strategy_ret');
%R_cumulative2=cumsum(strategy_ret);

%Benchmark
bench= mean(m_delta(:, relative_period+1:end));
bench_cumulative = ret2tick(bench');


%dates = monthly_delta_price.Date(relative_period:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%date per sottoperiodo
%sep-2005 to 31-dec-2012
startDate = datetime('01-Dec-2005');endDate = datetime('31-Dec-2012'); % In Sample
%startDate = datetime(['01-Mar-2013']);endDate = datetime('28-Jun-2024'); % Out of Sample
%startDate = datetime('01-Apr-1997');endDate = datetime('28-Jun-2024'); % Full Sample

dates =  Book1S5.Date(Book1S5.Date >= startDate & Book1S5.Date <= endDate);% In Sample
%dates =  OosDataset_M_price.Date(OosDataset_M_price.Date >= startDate & OosDataset_M_price.Date <= endDate); % Out of Sample
%dates =  FSDataset_M_price.Date(FSDataset_M_price.Date >= startDate & FSDataset_M_price.Date <= endDate); % Full Sample

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot del grafico con scala logaritmica
figure;
hold on;
plot(dates, log(R_cumulative), 'LineWidth', 2);
plot(dates, log(bench_cumulative), '--', 'LineWidth', 1.5);
%plot((relative_period+1:num_mesi+1), log(cumulative_returns_strategy)', 'LineWidth', 2);
%plot((relative_period+1:num_mesi+1), log(cumulative_returns_benchmark)', '--', 'LineWidth', 1.5);
xlabel('Mesi');
ylabel('Rendimento Cumulato (log scale)');
legend('Strategia (Relative Momentum)', 'Benchmark (EW B&H)');
title('Confronto tra Strategia e Benchmark');
grid on;
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% STRATEGIA RA (Relative + Absolute Momentum)

% Prepara il proxy di liquidità (VFISX)
cash_proxy_m = m_delta(4,5:end);  % VFISX come proxy del risk free (rendimenti mensili)

% Rimuovi VFISX dal dataset per il calcolo del momentum relativo
m_delta_no_vfisx = m_delta;
m_delta_no_vfisx(4,:) = [];  % Rimuove la riga 4 (VFISX)% Rimuove la riga 4 (VFISX) %%%%%%%%%%  !!!!!!! %%%%%%%%%% devo ricordarmi di non rannare più volte uqesta linea , altrimenti mi riduce ogni volta la dimensione del dataset e il codice si interrompe.

% Variabili strategia RA
absolute_threshold = 1;  % Soglia per il momentum assoluto
[num_rows, num_cols] = size(m_delta_no_vfisx);  % Dimensioni senza VFISX

% Calcolo del momentum relativo (R) senza VFISX
relmom = zeros(num_rows, num_cols-relative_period);  


for i = 1:num_rows
    for j = relative_period + 1:num_cols
        returns2 = m_delta_no_vfisx(i, j-relative_period:j-1);  % Rendimenti su 4 mesi (lookback)
        cumulative_returns2 = ret2tick(returns2');  
        last_cumulative_return2 = cumulative_returns2(end);  % Prendi ultimo valore cumulato
        relmom(i, j-relative_period) = last_cumulative_return2;  % Salva il valore
    end
end


[a, idx_rel_mom2] = sort(relmom, 1, 'descend');  

m_delta2 = m_delta_no_vfisx(:, 5:end);  % Matrice dei rendimenti mensili escluso VFISX

% Preallocazione per i rendimenti della strategia
RA_ret = zeros(3, size(m_delta2, 2)); 


for k = 1:size(m_delta2, 2)
   
    RAselected_assets = idx_rel_mom2(1:3, k);  % Primi 3 asset selezionati per mese k
    
    for asset_idx = 1:3  
        % Momentum relativo è già calcolato nella matrice 'relmom'
        last_cumulative_return = relmom(RAselected_assets(asset_idx), k);  
        
        % Se il momentum assoluto è > 1, prendi l'asset selezionato, altrimenti --> con VFISX
        if last_cumulative_return > absolute_threshold
            RA_ret(asset_idx, k) = m_delta2(RAselected_assets(asset_idx), k);  
        else
            RA_ret(asset_idx, k) = cash_proxy_m(k);  % Sostituisci con VFISX
        end
    end
end


RA_returns = mean(RA_ret, 1);

% Calcolo dei rendimenti cumulati della strategia
RA_cumulative = ret2tick(RA_returns'); %In Sample
%dates = Book1S6.Date(relative_period:end);
 

% PLot delle 2 strategie R e RA e del benchmark
figure;
hold on;
plot(dates, log(R_cumulative)', 'LineWidth', 2);
plot(dates, log(RA_cumulative)','LineWidth',2)
plot(dates, log(bench_cumulative)', '--', 'LineWidth', 1.5);

xlabel('Mesi');
ylabel('Rendimento Cumulato (log scale)');
legend('Strategia (Relative Momentum)', 'Strategia (RA)','Benchmark (EW B&H)');
title('Confronto tra Strategia e Benchmark');
grid on;
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Strategia rav

wR = 1; % Peso del rendimento relativo
wV = 0.5; % Peso della volatilitàcash_proxy_index = 4; % VFISX come proxy del risk-free
num_days = size(d_delta, 2);
N=3;

cash_proxy_d= d_delta(4,:); %provare a spostare all inizio del codice per togliere il riskfree in tutte le strategie
d_delta(4,:) = []; % Rimuove la riga 4 (VFISX) %%%%%%%%%%  !!!!!!! %%%%%%%%%% devo ricordarmi di non rannare più volte uqesta linea , altrimenti mi riduce ogni volta la dimensione del dataset e il codice si interrompe.

num_assets= size(m_delta2,1);



% 1. Calcolo della volatilità con finestra mobile (84 giorni)
lookback_days = relative_period * 21;  % Circa 84 giorni per il lookback period
volatility = zeros(num_assets, num_cols-relative_period);  % Preallocazione
% movstd per calcolare la deviazione standard mobile su 84 giorni
for i = 1:num_assets
    % Calcola la volatilità su una finestra mobile di 84 giorni
    volatility_daily = movstd(d_delta(i, :), [lookback_days-1, 0], 0, 2);  % Volatilità mobile
    % Prendi solo i valori di volatilità corrispondenti alle date mensili (fine del mese)
    volatility(i, :) = volatility_daily(lookback_days:21:end);
end


ranking_function = wR * relmom + wV * volatility;  % Funzione di loss combinata
[~, sorted_indices] = sort(ranking_function,1, 'descend');  % Classifica in ordine decrescente

% Selezione dei top N asset
RAVselected_assets = sorted_indices(1:N, :);  % Seleziona i top N asset per ogni periodo


RAV_returns = zeros(1, num_cols-relative_period);  % Preallocazione per i rendimenti della strategia
m_delta3 = m_delta(:, relative_period+1:end); % dati mensilo senza cash proxy

for k = 1:(num_cols-relative_period)
    selected_assets = RAVselected_assets(:, k);  
    
    for asset_idx = 1:N
        asset = selected_assets(asset_idx);
        
        if relmom(asset, k) > absolute_threshold
            RAV_returns(asset_idx, k) = m_delta2(asset, k);  % Rendimento dell'asset selezionato
        else
            RAV_returns(asset_idx, k) = cash_proxy_m(k);  % Investi nel proxy di liquidità (VFISX)
        end
    end
end


RAV_strategy_returns = mean(RAV_returns, 1);

% 6. Calcolo dei rendimenti cumulati della strategia RAV
cumulative_RAV_returns = ret2tick(RAV_strategy_returns');

figure;
plot(log([bench_cumulative R_cumulative RA_cumulative cumulative_RAV_returns]))
grid on
xlabel('Time')
ylabel('Cumulative Returns')
legend ('Benchmark','R','RA', 'RAV')
title('Strategies vs Benchmark');

figure(2);
plot([bench_cumulative cumulative_RAV_returns])
grid on
xlabel('Time')
ylabel('Cumulative Returns')
legend ('Benchmark','RAV')
title('Strategies vs Benchmark');

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Strategia 4: RAVC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wR = 1;    % Peso del rendimento relativo
wV = 0.5;  % Peso della volatilità
wC = 0.5;  % Peso della correlazione



num_assets = size(m_delta_no_vfisx, 1);
num_cols = size(m_delta_no_vfisx, 2);
lookback_months = relative_period;  % 4 mesi

  
cash_proxy_index = 4;  % VFISX come proxy del risk-free

% Calcolo della correlazione media sui rendimenti giornalieri
correl = zeros(num_assets, num_cols-relative_period);  % Preallocazione per le correlazioni

for t = relative_period + 1:num_cols
    % Calcola l'indice di partenza e di fine per la finestra di 84 giorni
    start_idx = t * 21 - lookback_days + 1;  % Inizio finestra di 84 giorni
    end_idx = t * 21;                        % Fine finestra
    
    if end_idx <= size(d_delta, 2)
        % Estra i rendimenti giornalieri per il periodo di lookback
        period_returns = d_delta(:, start_idx:end_idx);  
        
        % corr sui rendimenti giornalieri
        corr_mat = corrcoef(period_returns');

        % Calcola la correlazione media per ogni asset rispetto agli altri (escludendo la diagonale)
        avg_corr = mean(corr_mat - diag(diag(corr_mat)), 2);
        correl(:, t-relative_period) = avg_corr;  % Salva le correlazioni medie
    end
end



L_i_RAVC = wR * relmom + wV * volatility + wC * correl;  % Funzione di ranking RAVC


[~, sorted] = sort(L_i_RAVC, 'descend'); 
selected_assets_RAVC = sorted(1:N, :);    


RAVC_returns = zeros(1, num_cols-relative_period);  
m_delta4 = m_delta(:, relative_period+1:end);       
cash_proxy_m = m_delta(cash_proxy_index, relative_period+1:end);  
for k = 1:(num_cols-relative_period)
    selected_assets = selected_assets_RAVC(:, k); 
    for asset_idx = 1:N
        asset = selected_assets(asset_idx);
        
       
        if relmom(asset, k) > absolute_threshold
            RAVC_returns(asset_idx, k) = m_delta2(asset, k);  % Rendimento dell'asset selezionato
        else
            RAVC_returns(asset_idx, k) = cash_proxy_m(k);  % Investi nel proxy di liquidità (VFISX)
        end
    end
end


RAVC_strategy_returns = mean(RAVC_returns, 1);

cumulative_RAVC_returns = ret2tick(RAVC_strategy_returns');

figure;
hold on;
plot(dates, log(R_cumulative), 'LineWidth', 2);
plot(dates, log(RA_cumulative),'LineWidth',2);
plot(dates, log(cumulative_RAV_returns),'LineWidth',2);
plot(dates, log(cumulative_RAVC_returns),'LineWidth',2);
plot(dates, log(bench_cumulative)', '--', 'LineWidth', 1.5);

xlabel('Mesi');
ylabel('Rendimento Cumulato (log scale)');
legend('Strategia (Relative Momentum)', 'Strategia (RA)','Strategia (RAV)','Strategia (RAVC)','Benchmark (EW B&H)');
title('Confronto tra Strategia e Benchmark');
grid on;
hold off;

%Test diagnostici delle strategie
mean_ret = mean([bench' strategy_ret' RA_returns' RAV_strategy_returns' RAVC_strategy_returns'],1)*12; 
dev_std_ret = std([bench' strategy_ret' RA_returns' RAV_strategy_returns' RAVC_strategy_returns']).*sqrt(12); %dev std annuale
efficiency = mean_ret./dev_std_ret; % 
worst_months = min([bench' strategy_ret' RA_returns' RAV_strategy_returns' RAVC_strategy_returns']);
best_month = max([bench' strategy_ret' RA_returns' RAV_strategy_returns' RAVC_strategy_returns']);
perc_sopra_benchmark= sum([bench' strategy_ret' RA_returns' RAV_strategy_returns' RAVC_strategy_returns']>0, 1)./size([bench' strategy_ret' RA_returns' RAV_strategy_returns' RAVC_strategy_returns'], 1);

Max_dd = maxdrawdown([bench_cumulative R_cumulative RA_cumulative cumulative_RAV_returns cumulative_RAVC_returns]);
VaR_5_pct = portvrisk(mean_ret, dev_std_ret)';






summary_stats = table()


summary_stats(1, :) = array2table(mean_ret);                % Mean Ret
summary_stats(2, :) = array2table(dev_std_ret);             % Ann. STD
summary_stats(3, :) =array2table(efficiency);              % Sharpe Ratio
summary_stats(4, :) =array2table(VaR_5_pct);               % VaR at 5%
summary_stats(5, :) = array2table(Max_dd);                  % Maximum Drawdown
summary_stats(6, :) = array2table(worst_months);            % Worst Month Return
summary_stats(7, :) = array2table(best_month);              % Best Month Return
summary_stats(8, :) = array2table(perc_sopra_benchmark);    % Profitable Months
% Impostazione dei nomi delle colonne
summary_stats.Properties.VariableNames = { ...
    'Benchmark', ...
    'R_momentum', ...
    'RA_momentum', ...
    'RAV_momentum', ...
    'RAVC_momentum'};

% Impostazione dei nomi delle righe
summary_stats.Properties.RowNames = { ...
    'Ann_MeanRet', ...
    'Ann_STD', ...
    'Sharpe_Ratio', ...
    'VaR_at_5%', ...
    'Maximum_Drawdown', ...
    'Worst_Month_Return', ...
    'Best_Month_Return', ...
    'Profitable_Months'};
