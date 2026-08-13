% general prediction moving window code with [S(0) I(0) R(0)] = [N - I0/rho - R0  I0/rho R0]
% Here SIR simulated data on the moving day started on day when cases is nonzero.
clear all; close all; clc;
[state, m, name_origin_list, time_start, lockdown_date, top_five_pop_flow_nodes, time_step_size , zeta_vec, zeta,  Nv, gamma, beta_1_vec, phi_0] = county_name_list();
lockdown_date   = 88;
statefips       = 40; % Define state FIPS code
str             = 'sir_model_fitting';
m               = 77;
N_o             = m;
v               = 1;
N_d             = m;
no_iteration    = 200;            % used the radom process to find the optimal value in the 1st iteration
num_sim         = 200;
% moving_window              = 5;
% no_ahead_prediction        = 10;
% ylimit        = 20000;
% xlimit        = 155;
covid_data = readmatrix('final_weekly_time_series_covid_data.csv');
original_covid_data = readmatrix('final_weekly_time_series_covid_data.csv');
all_counties_full_ts = covid_data(:,2:end);
%% here we can expand the time of the real data
covid_data = covid_data(7:76,:);
% Define dengue_div_data from covid_data (columns 4 to 11)
covid_county_data       = covid_data(:,2:end);  % Real data (time along rows, nodes along columns)
real_covid_data         = original_covid_data(7:81,:);
real_covid_county_data  = real_covid_data(:,2:end);  % Real data (time along rows, nodes along columns)
time_end_plot           = size(covid_data, 1);

origin_population                   = csvread('pop_vector_origin_only.csv',1,0);
origin_population_with_county_name  = readtable('pop_vector_origin.csv', 'Delimiter', ',');
vec_N                               =    origin_population(:,2);            % population of each nodes

% time_end_plot           = size(covid_data, 1);
% time_end_plot           = size(covid_county_data, 1);
% lockdown_date_0_indexed = lockdown_date;
% lockdown_week_0_indexed = 4;
time_end  = time_end_plot;
T         = time_end;
%plot_time_end = T;
time_step_size = 0.1; % can compute N_time_steps
% beta0         = 1;
%%
%rho    = 10^-1; % 1
tspan  = 1:time_step_size:time_end;  % Ensure this is within the bounds of `t`
t      = tspan;
node_vec = linspace(1,m,m);

% % Set beta and gamma vectors
gamma_lb                   = 0.5;
lb                         = [0.1 0.001  gamma_lb];                %[beta xi par_rho gamma]
ub                         = [3 0.2 1];                               %[beta xi rho gamma]
% initial parameter
beta                       = 1;
rho                        = 10^-1;
covid_cases_vec_all         = covid_county_data;
index_post_covid_cases = find(any(covid_cases_vec_all > 0, 2), 1, 'first');
int_index_post_covid_cases = index_post_covid_cases(:,1); % get the first nonzero covid-19 case
d                          = int_index_post_covid_cases;
covid_cases_vec_all        = covid_cases_vec_all(:,d:end);

optim_optim_beta_vec      = zeros(m,1);
optim_optim_rho_vec       = zeros(m,1);
optim_optim_gamma_vec     = zeros(m,1);
optim_optim_error_vec     = zeros(m,1);
n                         = size(covid_county_data,1);
sim_covid_county_data     = zeros(size(covid_county_data,1),m);
optim_optim_beta_vec_bootstrapping_sample  = zeros(num_sim,m);
optim_optim_rho_vec_bootstrapping_sample   = zeros(num_sim,m);
optim_optim_gamma_vec_bootstrapping_sample = zeros(num_sim,m);
%options = optimoptions('fmincon', 'Display', 'iter'); % Optional: set display options

%selected_counties = [1 3 9 14 16 55 72];   % example
%selected_counties = [77];   % example
selected_counties   = 1:m;
%selected_counties   = setdiff(1:m, [2,29,30, 65]);  % removes 65, 10, and 20;
% for county_index = 1 : m

for county_index = selected_counties

    %% select one county
    single_county_covid_cases  = covid_county_data(:,county_index);
    N                          = vec_N(county_index);
    index_post_covid_cases     = find(single_county_covid_cases> 0);       % find the nonzero covid-19 cases
    int_index_post_covid_cases = index_post_covid_cases(1); % get the first nonzero covid-19 case
    d                          = int_index_post_covid_cases;
    covid_cases                = single_county_covid_cases(d:end);

    Realdata                   = covid_cases;
    R0    = 0;                     % initial iteration R(0) =0
    % I0 is real rep. case, but practically it might not be available so update by X(2,2)
    I0    = (1/rho)*covid_cases(1);
    S0    = N - I0 - R0;
    % IC  = [N - I0 - R0  I0 R0 I0];   % [S0 I0 R0 I0]
    IC    = [S0 I0 R0 I0];
    % I0 is real rep. case, but practically it might not be available so update by X(2,2)
    tsol  = [1: length(Realdata)];
    %tsol = t;
    tdata = tsol;
    %t    = tsol;
    % we can change the pro9duct vector data [0 0 0 rho]
    %par =[par_beta_xi par_rho];
    %gamma = gamma_lb;
    par = [beta rho gamma];
    SIRparsol = @(par, t) SIRsol(par(1), par(2), par(3), IC, t);
    SumSquaresErrorSIRfunction = @(par)sum((SIRparsol(par, tdata)- Realdata').^2);
    %% start the optimization and error for different beta and gamma
    optim_SSE_error_matrix = zeros(no_iteration,1);
    optim_beta_matrix      = zeros(no_iteration,1);
    optim_gamma_matrix     = zeros(no_iteration,1);
    optim_rho_matrix       = zeros(no_iteration,1);


    %% we are optimizing beta, rho, and gamma using a random process

    for i = 1: no_iteration
        %for j = 1: no_iteration
        % lb                         = [0.1 1/(6*4) 0.001  gamma_lb]; [beta xi par_rho gamma]
        % ub                         = [3 1/10 0.2 1];
        X = rand(length(lb),1);
        new_randm_para = lb + (ub-lb).*X.';
        beta     = new_randm_para(1);
        par_beta = beta;
        par_rho  = new_randm_para(2);
        gamma    = new_randm_para(3);
        %             for i = 1: length(beta_vec)
        %                 beta = beta_vec(i);
        %                 for j = 1: length(gamma_vec)
        %                     %xi                               = 1/(6*4);
        %                     gamma                            = gamma_vec(j);
        [SIR_optm_theta, SIR_optm_error] = fmincon(SumSquaresErrorSIRfunction, [par_beta par_rho gamma], [], [], [], [], lb,ub);

        optim_SSE_error_matrix(i)        = SIR_optm_error;
        % Extract and store beta value
        optim_beta_matrix(i)             = SIR_optm_theta(1);
        % Extract and store par_rho value
        optim_rho_matrix(i)              = SIR_optm_theta(2);
        optim_gamma_matrix(i)            = SIR_optm_theta(3);
    end

    %% find the minimum error and its corresponding index and find corresponding beta and gamma
    min_error                 = min(min(optim_SSE_error_matrix));
    [row_index, column_index] = find(optim_SSE_error_matrix == min_error);
    opt_beta                  = optim_beta_matrix(row_index);
    opt_gamma                 = optim_gamma_matrix(row_index);
    opt_rho                   = optim_rho_matrix(row_index);  % Assuming optim_xi_matrix contains xi values
    opt_par                   = [opt_beta(1) opt_rho(1) opt_gamma(1)];
    opt_sol                   = SIRsol(opt_par(1), opt_par(2), opt_par(3), IC, tdata);
    % Store optimal parameters for each w
    optim_optim_beta_vec(county_index)             = opt_par(1);
    optim_optim_rho_vec(county_index)              = opt_par(2);
    optim_optim_gamma_vec(county_index)            = opt_par(3);
    optim_optim_error_vec(county_index)            = min_error;
    sim_covid_county_data(:, county_index) = [zeros(d-1,1); opt_sol'];
end


% ----- Prepare table of optimized parameters -----
% T_params = table((1:m)', optim_optim_beta_vec', optim_optim_rho_vec', optim_optim_gamma_vec', optim_optim_error_vec', ...
%     'VariableNames', {'CountyIndex','Beta','Rho','Gamma','MinError'});

T_params = table( ...
    (1:m)', ...
    optim_optim_beta_vec(:), ...
    optim_optim_rho_vec(:), ...
    optim_optim_gamma_vec(:), ...
    optim_optim_error_vec(:), ...
    'VariableNames', {'CountyIndex','Beta','Rho','Gamma','MinError'});
% ----- Ensure folder exists -----
folderPath = 'frac-cum y_k data/SIRmf03PFITboothstrapPIFullTimeSeriesFitting/';
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
end

% ----- Create filenames with state FIPS -----
paramsFilename = fullfile(folderPath, ['SIR_optim_parameters_state' num2str(statefips) '.csv']);
simFilename    = fullfile(folderPath, ['sim_covid_county_data_state' num2str(statefips) '.csv']);

% ----- Save CSV files -----
writetable(T_params, paramsFilename);          % parameters
writematrix(sim_covid_county_data, simFilename);  % simulation results

disp(['CSV files saved for state ' num2str(statefips)]);



%% ----------------- Figure -----------------

L = size(covid_county_data,1);

h03 = figure;
hold on;

for county_index = selected_counties % 1:m
    beta0 = optim_optim_beta_vec(county_index);
    % Determine common length for time and data
    t_plot = 1:size(covid_county_data,1);
    obs_node = covid_county_data(1:L, county_index);
    sim_node = sim_covid_county_data(1:L, county_index);
    
    % Plot other counties
    if ~ismember(county_index, top_five_pop_flow_nodes)
        plot(t_plot, obs_node, 'Color', [0.6 0.8 1], 'LineWidth', 1.5);  % observed
        plot(t_plot, sim_node, 'Color', [0.2 0.6 0.2], 'LineWidth', 1.5); % simulated
    end
    ylim([0, 6000]);
    xlim([0, 70]);
end

% Highlight top five counties
colors_top = {'c','b','g','k','m'};
for k = 1:numel(top_five_pop_flow_nodes)
    county_index = top_five_pop_flow_nodes(k);
    % Determine common length
    t_plot = 1:size(covid_county_data,1);
    obs_node = covid_county_data(1:L, county_index);
    sim_node = sim_covid_county_data(1:L, county_index);
    % Plot top counties with distinct colors
    plot(t_plot, obs_node, [colors_top{k} '-'], 'LineWidth', 3);  % observed solid
    plot(t_plot, sim_node, [colors_top{k} '--'], 'LineWidth', 3); % simulated dashed
    ylim([0, 6000]);
    xlim([0, 70]);
end

% % Optional: plot total simulated infections
% I_total = sum(sim_covid_county_data(1:L, :), 2);
% plot(t_plot, I_total, 'r--', 'LineWidth', 3);

% Appearance settings
xlabel('Time (days)', 'FontSize', 20);
ylabel('Infections', 'FontSize', 20);
title('Observed vs Simulated COVID Cases by County', 'FontSize', 22);
set(gca, 'LineWidth', 2, 'fontsize', 20, 'FontName', 'Times', 'box', 'on');
% legend({'Other counties (Observed)','Other counties (Simulated)', ...
%         'Top counties (Observed)','Top counties (Simulated)','Total Simulated'}, ...
%         'FontSize', 14, 'Location', 'best');

%% ----------------- SAVE FIGURES -----------------
folderPath = 'frac-cum y_k data/SIRmf03PFITboothstrapPIFullTimeSeriesFitting/result figures/fig/';
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
end

figFilename = fullfile(folderPath, ...
    ['Fig_All_Counties_TimeSeries_phi_' num2str(phi_0) '_m_' num2str(m) ...
    '_Net_' str '_inf_at_' num2str(v) '_Het_Node_Mult_' num2str(zeta) ...
    '_beta_i_' num2str(beta0) '_GT_' num2str(gamma) '.fig']);

folderPath = 'frac-cum y_k data/SIRmf03PFITboothstrapPIFullTimeSeriesFitting/result figures/png/';
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
end

pngFilename = fullfile(folderPath, ...
    ['Fig_All_Counties_TimeSeries_phi_' num2str(phi_0) '_m_' num2str(m) ...
    '_Net_' str '_inf_at_' num2str(v) '_Het_Node_Mult_' num2str(zeta) ...
    '_beta_i_' num2str(beta0) '_GT_' num2str(gamma) '.png']);

saveas(h03, figFilename);
saveas(h03, pngFilename);

%% ----------------- SAVE TIME SERIES DATA -----------------
% for county_index = 1:m
%     L = min([length(t), size(covid_county_data,1), size(sim_covid_county_data,1)]);
%     countyData = [covid_county_data(1:L, county_index), sim_covid_county_data(1:L, county_index)];
%     countyFilename = fullfile(folderPath, ...
%         ['County_' num2str(county_index) '_TimeSeries.csv']);
%     writematrix(countyData, countyFilename); % columns: [Observed, Simulated]
% end
%
% disp(['Figures and time series CSVs saved in: ' folderPath]);



% Initialize metric vectors
SSR_per_county   = zeros(1, m);
mse_per_county   = zeros(1, m);
rmse_per_county  = zeros(1, m);
mae_per_county   = zeros(1, m);
nrmse_per_county = zeros(1, m);
r2_per_county    = zeros(1, m);
mape_per_county  = zeros(1, m);
pae_per_county   = zeros(1, m);
pa_per_county    = zeros(1, m);
pe_per_county    = zeros(1, m); % percentage error

% Loop through each county (each column)
for i = 1:m
    % Predicted and real data
    opt_soln = sim_covid_county_data(:,i);
    Realdata = covid_county_data(:,i);
    T        = length(Realdata);

    % Residuals
    residuals = opt_soln - Realdata;

    % SSR (sum of squared residuals)
    SSR = sum(residuals.^2);

    % MSE
    mse = mean(residuals.^2);

    % RMSE
    rmse = sqrt(mse);

    % MAE
    mae = mean(abs(residuals));

    % Normalized RMSE
    nrmse = rmse / (max(Realdata) - min(Realdata));

    % R^2
    SST = sum((Realdata - mean(Realdata)).^2); % total sum of squares
    r2  = 1 - (SSR / SST);

    % MAPE (avoid divide by zero)
    mape = mean(abs(residuals) ./ max(Realdata, eps)) * 100;

    % PAE
    pae = sum(abs(residuals)) / sum(Realdata);

    % PA (percentage agreement)
    pa = mean(min(opt_soln, Realdata) ./ max(opt_soln, Realdata), 'omitnan');

    % PE (percentage error)
    pe = sum(opt_soln - Realdata) / sum(Realdata);

    % Store results
    SSR_per_county(i)   = SSR;
    mse_per_county(i)   = mse;
    rmse_per_county(i)  = rmse;
    mae_per_county(i)   = mae;
    nrmse_per_county(i) = nrmse;
    r2_per_county(i)    = r2;
    mape_per_county(i)  = mape;
    pae_per_county(i)   = pae;
    pa_per_county(i)    = pa;
    pe_per_county(i)    = pe;
end

% Add metrics as new columns in the table
origin_population_with_county_name.sir_baseline_SSR   = SSR_per_county';
origin_population_with_county_name.sir_baseline_mse   = mse_per_county';
origin_population_with_county_name.sir_baseline_mae   = mae_per_county';
origin_population_with_county_name.sir_baseline_rmse  = rmse_per_county';
origin_population_with_county_name.sir_baseline_nrmse = nrmse_per_county';
origin_population_with_county_name.sir_baseline_r2    = r2_per_county';
origin_population_with_county_name.sir_baseline_mape  = mape_per_county';
origin_population_with_county_name.sir_baseline_pae   = pae_per_county';
origin_population_with_county_name.sir_baseline_pa    = pa_per_county';
origin_population_with_county_name.sir_baseline_pe    = pe_per_county';

% Display the updated table
disp(origin_population_with_county_name);

% Save the updated table (CSV)
directory = 'frac-cum y_k data/SIRmf03PFITboothstrapPIFullTimeSeriesFitting';
file_path = fullfile(directory, sprintf( ...
    'ten_metric_sir_model_constant_beta_gamma_fitting_errors_with_all_counties_vari_selection.csv'));

if ~exist(directory, 'dir')
    mkdir(directory);
end

writetable(origin_population_with_county_name, file_path);



%% plot the reporting rate on each node


node_vec = linspace(1, m,m);

h2 = figure(2)
%I=linspace(1,11,11)
data = optim_optim_rho_vec' ;
%stem(data,'b','filled','LineWidth',2)
stem(node_vec, data,'r','filled','LineWidth',3)
%set(gca, 'YTick', 0:.1:ylimit); % Full form of gca is (Get Current Axes)
%ylim([0, ylimit])
% set(gca, 'XTick', 0:1: m); % Full form of gca is (Get Current Axes)
set(gca, 'XTick', [1 11 22 33 44 55 66 77]); % Full form of gca is (Get Current Axes)
xlim([1,m])
% ylim([0, ylimt])
% set(gca, 'XTick', 0:floor(tickSpacing):m);
xlabel('Node i')
ylabel('\rho_i')
% Adjust axis properties
set(gca, 'fontsize', 30, 'FontName', 'Times', 'LineWidth', 2, 'Box', 'on'); % Add square box and set line width
h2_figurename_fig = ['/frac-cum y_k data/SIRmf03PFITboothstrapPIFullTimeSeriesFitting/result figures/fig/Fig_baseline_0_modelm_node_vs_rho_vec.fig'];
saveas(h2,[pwd h2_figurename_fig])
h2_figurename_png = ['/frac-cum y_k data/SIRmf03PFITboothstrapPIFullTimeSeriesFitting/result figures/png/Fig_baseline_0_modelm_node_vs_rho_vec.png'];
saveas(h2,[pwd h2_figurename_png])


h3 = figure(3)
%I=linspace(1,11,11)
data = optim_optim_beta_vec' ;
%stem(data,'b','filled','LineWidth',2)
stem(node_vec, data,'r','filled','LineWidth',3)
%set(gca, 'YTick', 0:.1:ylimit); % Full form of gca is (Get Current Axes)
%ylim([0, ylimit])
% set(gca, 'XTick', 0:1: m); % Full form of gca is (Get Current Axes)
set(gca, 'XTick', [1 11 22 33 44 55 66 77]); % Full form of gca is (Get Current Axes)
xlim([1,m])
% ylim([0, ylimt])
% set(gca, 'XTick', 0:floor(tickSpacing):m);
xlabel('Node i')
ylabel('\beta_i')
% Adjust axis properties
set(gca, 'fontsize', 30, 'FontName', 'Times', 'LineWidth', 2, 'Box', 'on'); % Add square box and set line width
h3_figurename_fig = ['/frac-cum y_k data/SIRmf03PFITboothstrapPIFullTimeSeriesFitting/result figures/fig/Fig_baseline_0_modelm_node_vs_optim_optim_beta_vec.fig'];
saveas(h3,[pwd h3_figurename_fig])
h3_figurename_png = ['/frac-cum y_k data/SIRmf03PFITboothstrapPIFullTimeSeriesFitting/result figures/png/Fig_baseline_0_modelm_node_vs_optim_optim_beta_vec.png'];
saveas(h3,[pwd h3_figurename_png])


h4 = figure(4)
%I=linspace(1,11,11)
data = optim_optim_gamma_vec' ;
%stem(data,'b','filled','LineWidth',2)
stem(node_vec, data,'r','filled','LineWidth',3)
%set(gca, 'YTick', 0:.1:ylimit); % Full form of gca is (Get Current Axes)
%ylim([0, ylimit])
% set(gca, 'XTick', 0:1: m); % Full form of gca is (Get Current Axes)
set(gca, 'XTick', [1 11 22 33 44 55 66 77]); % Full form of gca is (Get Current Axes)
xlim([1,m])
% ylim([0, ylimt])
% set(gca, 'XTick', 0:floor(tickSpacing):m);
xlabel('Node i')
ylabel('\gamma_i')
% Adjust axis properties
set(gca, 'fontsize', 30, 'FontName', 'Times', 'LineWidth', 2, 'Box', 'on'); % Add square box and set line width
h4_figurename_fig = ['/frac-cum y_k data/SIRmf03PFITboothstrapPIFullTimeSeriesFitting/result figures/fig/Fig_baseline_0_modelm_node_vs_optim_optim_gamma_vec.fig'];
saveas(h4,[pwd h4_figurename_fig])
h4_figurename_png = ['/frac-cum y_k data/SIRmf03PFITboothstrapPIFullTimeSeriesFitting/result figures/png/Fig_baseline_0_modelm_node_vs_optim_optim_gamma_vec.png'];
saveas(h4,[pwd h4_figurename_png])


%% ============================================================
%   NEW SECTION: SAVE EACH COUNTY'S INDIVIDUAL TIME SERIES
% ===============================================================

figFolder = 'frac-cum y_k data/SIRmf03PFITboothstrapPIFullTimeSeriesFitting/indv_county_time_series_real_vs_fit/fig/';
pngFolder = 'frac-cum y_k data/SIRmf03PFITboothstrapPIFullTimeSeriesFitting/indv_county_time_series_real_vs_fit/png/';

if ~exist(figFolder, 'dir'), mkdir(figFolder); end
if ~exist(pngFolder, 'dir'), mkdir(pngFolder); end

colors_top = {'c','b','g','k','m'};   % same as before

for county_index = selected_counties
    optim_beta  = optim_optim_beta_vec(county_index);
    optim_rho   = optim_optim_rho_vec(county_index);
    optim_gamma = optim_optim_gamma_vec(county_index);
    r_squared   = r2_per_county(county_index);
    RMSE        = rmse_per_county(county_index);
    MAE         = mae_per_county(county_index);
    % Data
    t_plot   = 1:L;
    obs_node = covid_county_data(1:L, county_index);
    sim_node = sim_covid_county_data(1:L, county_index);
    y_max = max(max(obs_node), max(sim_node)) + 100;
    % Invisible figure
    hCounty = figure('Visible','off');
    hold on;

    % ------------------------------------------------------------
    %   COLOR LOGIC (matching exactly your previous plot)
    % ------------------------------------------------------------

    if ismember(county_index, top_five_pop_flow_nodes)
        % Determine color for this top county
        k = find(top_five_pop_flow_nodes == county_index);
        c_top = colors_top{k};

        % Top counties (solid / dashed)
        plot(t_plot, obs_node, [c_top '-'], 'LineWidth', 3);
        plot(t_plot, sim_node, [c_top '--'], 'LineWidth', 3);

    else
        % Other counties
        plot(t_plot, obs_node, 'Color', [0.6 0.8 1], 'LineWidth', 1.5);
        plot(t_plot, sim_node, 'Color', [0.2 0.6 0.2], 'LineWidth', 1.5);
    end

    % ------------------------------------------------------------
    %   VISUALS
    % ------------------------------------------------------------
    % xlabel('Time (days)', 'FontSize', 16);
    % ylabel('Infections', 'FontSize', 16);
    %title(['County Index = ' num2str(county_index)], 'FontSize', 30);
    set(gca, 'LineWidth', 1.5, 'FontSize', 30, 'FontName', 'Times');
    box on;
    xlim([0, 70]);
    ylim([0, y_max]);

    % ------------------------------------------------------------
    %   FILE NAMES
    % ------------------------------------------------------------
    figFile = fullfile(figFolder, ...
        ['County_' num2str(county_index) '_optim_beta_' num2str(optim_beta) ...
        '_optim_rho_' num2str(optim_rho) '_Net_' str '_inf_' num2str(v) ...
        '_optim_gamma_' num2str(optim_gamma) '_r_squared_' num2str(r_squared) '_RMSE_' num2str(RMSE) '_MAE_' num2str(MAE) '.fig']);

    pngFile = fullfile(pngFolder, ...
        ['County_' num2str(county_index) '_optim_beta_' num2str(optim_beta) ...
        '_optim_rho_' num2str(optim_rho) '_Net_' str '_inf_' num2str(v) ...
        '_optim_gamma_' num2str(optim_gamma) '_r_squared_' num2str(r_squared) '_RMSE_' num2str(RMSE) '_MAE_' num2str(MAE) '.png']);

    % Save files
    saveas(hCounty, figFile);
    saveas(hCounty, pngFile);

    close(hCounty);
end





%% 1. Fig for indv_county_time_series_real_vs_fit_bootstrap_pi
%% ==============================================================================
%  SAVE EACH COUNTY'S INDIVIDUAL TIME SERIES WITH BOOTHSTRAPPING PI
% ===============================================================================

figFolder = 'frac-cum y_k data/SIRmf03PFITboothstrapPIFullTimeSeriesFitting/indv_county_time_series_real_vs_fit_bootstrap_pi/fig/';
pngFolder = 'frac-cum y_k data/SIRmf03PFITboothstrapPIFullTimeSeriesFitting/indv_county_time_series_real_vs_fit_bootstrap_pi/png/';

if ~exist(figFolder, 'dir'), mkdir(figFolder); end
if ~exist(pngFolder, 'dir'), mkdir(pngFolder); end

colors_top              = {'r','b','g','k','m'};% same as before
color_rest              = [0.30, 0.70, 1.00];   % light blue (not gray)
ci_color                = [1.0, 0.6, 0.6];      % Soft Peach for confidence interval
color_mean_sim_samp_fit = [0.4, 0.0, 0.6];      % deep purple
%color_sim_samp_fit      = [0.75, 0.75, 0.25];   % olive-gold
color_sim_samp_fit       = [0.90, 0.88, 0.98];  % lighter


for county_index = selected_counties %1:m
    optim_beta   = optim_optim_beta_vec(county_index);
    optim_rho    = optim_optim_rho_vec(county_index);
    optim_gamma  = optim_optim_gamma_vec(county_index);
    r_squared    = r2_per_county(county_index);
    RMSE         = rmse_per_county(county_index);
    MAE          = mae_per_county(county_index);
    % Data
    t_plot   = 1:L;
    obs_node = covid_county_data(1:L, county_index);
    sim_node = sim_covid_county_data(1:L, county_index);
    y_max = max(max(obs_node), max(sim_node)) + 100;
    % Invisible figure
    hCounty = figure('Visible','off');
    hold on;


    %% Bootstraping prediction interval

    single_county_covid_cases  = covid_county_data(:,county_index);
    N                          = vec_N(county_index);
    index_post_covid_cases     = find(single_county_covid_cases> 0);       % find the nonzero covid-19 cases
    int_index_post_covid_cases = index_post_covid_cases(1); % get the first nonzero covid-19 case
    d                          = int_index_post_covid_cases;
    covid_cases                = single_county_covid_cases(d:end);

    Realdata                   = covid_cases;
    cases                      = Realdata;
    R0    = 0;                     % initial iteration R(0) =0
    % I0 is real rep. case, but practically it might not be available so update by X(2,2)
    I0    = (1/rho)*covid_cases(1);
    S0    = N - I0 - R0;
    % IC  = [N - I0 - R0  I0 R0 I0];   % [S0 I0 R0 I0]
    IC    = [S0 I0 R0 I0];
    %   I0 is real rep. case, but practically it might not be available so update by X(2,2)
    tsol  = [1: length(Realdata)];
    %tsol  = t;
    tdata = tsol;
    tspan = tdata;
    %t     = tsol;
    % we can change the pro9duct vector data [0 0 0 rho]
    %par =[par_beta_xi par_rho];
    %gamma = gamma_lb;
    par              = [optim_beta optim_rho optim_gamma];
    sir_sim_cases    = SIRsol(par(1), par(2), par(3), IC, tspan)
    %% ========================================================================
    %  95% Confidence interval for bestfit curve
    %  ========================================================================
    sim_cases = sir_sim_cases';
    %sim_cases = fun_fitting(x_optimal, tspan);
    F = cumsum(sim_cases); % Compute the cumulative sum of the array z1 and store the result in F.
    length_time = length(tspan);

    %% Additing negative binomial error distribution

    nbd_simulations = zeros(length_time,num_sim);

    % average squared deviation of observations from their mean
    variance_fit = var(cases);
    mean_fit = mean(cases);
    factor1 = variance_fit/mean_fit;


    % Initialize the first row of the simulation matrix
    nbd_simulations(1, :) = sim_cases(1);

    % Perform simulations for each column
    for matrix_col = 1:num_sim
        for matrix_row = 2:length_time
            % Compute new cases at time t
            newcases_t = F(matrix_row) - F(matrix_row - 1);

            % Compute mean and variance
            mean1 = newcases_t;
            var1 = mean1 * factor1;

            % Compute parameters for the negative binomial distribution
            param2 = mean1 / var1;
            param1 = mean1 * param2 / (1 - param2);

            % Generate random value from the negative binomial distribution
            nbd_simulations(matrix_row, matrix_col) = nbinrnd(param1, param2);
        end
    end


    %% model fitting for new data with negative binomial distribution

    variable_data = zeros(num_sim,length(lb));
    fval_data  = zeros(num_sim,1);

    for sim = 1:num_sim

        cases_sim = nbd_simulations(:,sim);
        SIRparsol = @(par, t) SIRsol(par(1), par(2), par(3), IC, t);
        SumSquaresErrorSIRfunction = @(par)sum((SIRparsol(par, tdata)- cases_sim').^2);
        [SIR_optm_theta, SIR_optm_error] = fmincon(SumSquaresErrorSIRfunction, [optim_beta optim_rho optim_gamma], [], [], [], [], lb,ub);

        % variable_guess = x_optimal;
        % [variable, fval] = lsqcurvefit(@fun_fitting,variable_guess,tspan,cases_sim,lb,ub,options);
        % variable_data(sim,:) = variable;
        % fval_data(sim)  = fval;
        variable_data(sim,:) = SIR_optm_theta;
        fval_data(sim)  = SIR_optm_error;
        % par              = [optim_beta optim_rho optim_gamma]
        optim_optim_beta_vec_bootstrapping_sample(sim,county_index)  = SIR_optm_theta(1);  
        optim_optim_rho_vec_bootstrapping_sample(sim,county_index)   = SIR_optm_theta(2);
        optim_optim_gamma_vec_bootstrapping_sample(sim,county_index) = SIR_optm_theta(3);
        fprintf('iteration number %d \n',sim);

    end

    %% Calculating 95% CI

    % Store all best fit curves in a matrix; each collumn represents one simulation

    matrix_best_fit = zeros(length(tspan),num_sim);
    for fig = 1:num_sim
        % matrix_best_fit(:,fig) = fun_fitting(variable_data(fig,:),tspan);
        matrix_best_fit(:,fig)    = SIRsol(variable_data(fig,1), variable_data(fig,2), variable_data(fig,3), IC, tspan)
    end


    %data = nbd_simulations;
    data = matrix_best_fit;


    % Desired quantiles to get the 95% CI (2.5% on the lower tail, 2.5% on the upper tail)
    q = [0.025, 0.975];

    % 90% CI
    % q = [0.05, 0.95]

    % 99% CI
    %q = [0.005, 0.995]

    % Compute quantiles along the second dimension (time series)
    quantiles = quantile(data, q, 2);

    % Extract lower and upper bounds
    orginal_lower_bound = quantiles(:, 1); % 2.5th percentile
    orginal_upper_bound = quantiles(:, 2); % 97.5th percentile

    lower_bound = [zeros(d-1,1); orginal_lower_bound];

    upper_bound = [zeros(d-1,1); orginal_upper_bound];
    %% Bootstraping prediction interval

    % ------------------------------------------------------------
    %   COLOR LOGIC (matching exactly your previous plot)
    % ------------------------------------------------------------

    if ismember(county_index, top_five_pop_flow_nodes)
        % Determine color for this top county
        k = find(top_five_pop_flow_nodes == county_index);
        c_top = colors_top{k};

        % Top counties (solid / dashed)
        % plot(t_plot, obs_node, [c_top '-'], 'LineWidth', 3);
        % plot(t_plot, sim_node, [c_top '--'], 'LineWidth', 3);
        fill([t_plot, fliplr(t_plot)], [lower_bound; flipud(upper_bound)], ...
            ci_color, 'FaceAlpha', 0.25, 'EdgeColor', 'none'); % Shaded area of CIs
        % Plot all trajectories lightly
        plot(t_plot, [zeros(d-1,num_sim); matrix_best_fit], 'Color', color_sim_samp_fit, 'LineWidth', 0.8)
        % Plot mean trajectory boldly
        plot(t_plot, [zeros(d-1,1); mean(matrix_best_fit,2)], 'Color', color_mean_sim_samp_fit, ...
            'LineWidth', 3)
        plot(t_plot, obs_node, [c_top '-'], 'LineWidth', 3);
        plot(t_plot, sim_node, [c_top '--'], 'LineWidth', 3);
    else
        % Other counties 'Color', [0.6 0.8 1]  (very light blue)
        %  [0.8 0.8 0.8]  light gray for the Prediction interval
        % plot(t_plot, obs_node, 'Color', [0.6 0.8 1], 'LineWidth', 1.5);
        % plot(t_plot, sim_node, 'Color', [0.2 0.6 0.2], 'LineWidth', 1.5);
        fill([t_plot, fliplr(t_plot)], [lower_bound; flipud(upper_bound)], ...
            ci_color, 'FaceAlpha', 0.25, 'EdgeColor', 'none'); % Shaded area of CIs
        % Plot all trajectories lightly
        plot(t_plot, [zeros(d-1,num_sim); matrix_best_fit], 'Color', color_sim_samp_fit, 'LineWidth', 0.8)
        % Plot mean trajectory boldly
        plot(t_plot, [zeros(d-1,1); mean(matrix_best_fit,2)], 'Color', color_mean_sim_samp_fit, ...
            'LineWidth', 3)
        plot(t_plot, obs_node, 'Color', color_rest,'LineStyle', '-', 'LineWidth', 3)
        plot(t_plot, sim_node, 'Color', color_rest, 'LineStyle', '--', 'LineWidth', 3)

        % plot(1: length(real_covid_county_data(:,county_index)), real_covid_county_data(:,county_index), 'Color', color_rest,'LineStyle', '-', 'LineWidth', 3)

    end

    % ------------------------------------------------------------
    %   VISUALS
    % ------------------------------------------------------------
    % xlabel('Time (days)', 'FontSize', 16);
    % ylabel('Infections', 'FontSize', 16);
    %title(['County Index = ' num2str(county_index)], 'FontSize', 30);
    set(gca, 'LineWidth', 1.5, 'FontSize', 30, 'FontName', 'Times');
    box on;
    xlim([0, 80]);
    ylim([0, y_max]);

    % ------------------------------------------------------------
    %   FILE NAMES
    % ------------------------------------------------------------
    figFile = fullfile(figFolder, ...
        ['County_' num2str(county_index) '_optim_beta_' num2str(optim_beta) ...
        '_optim_rho_' num2str(optim_rho) '_Net_' str '_inf_' num2str(v) ...
        '_optim_gamma_' num2str(optim_gamma) '_r_squared_' num2str(r_squared) '_RMSE_' num2str(RMSE) '_MAE_' num2str(MAE) '.fig']);

    pngFile = fullfile(pngFolder, ...
        ['County_' num2str(county_index) '_optim_beta_' num2str(optim_beta) ...
        '_optim_rho_' num2str(optim_rho) '_Net_' str '_inf_' num2str(v) ...
        '_optim_gamma_' num2str(optim_gamma) '_r_squared_' num2str(r_squared) '_RMSE_' num2str(RMSE) '_MAE_' num2str(MAE) '.png']);

    % Save files
    saveas(hCounty, figFile);
    saveas(hCounty, pngFile);

    close(hCounty);
end





%% ==============================================================================================================
%  Fig 2. with forecast SAVE EACH COUNTY'S INDIVIDUAL TIME SERIES WITH BOOTHSTRAPPING PI and FORECASTS PI USING BOOTHSTRAPPING
% ===============================================================================================================
% Preallocate table copunty level data
all_county_results      = [];   % Initialize once before the loop

no_week_ahead_forecasts = 5;

forecast_horizon = no_week_ahead_forecasts;  % number of extra weeks to forecast


figFolder = 'frac-cum y_k data/SIRmf03PFITboothstrapPIFullTimeSeriesFitting/indv_county_time_series_real_vs_fit_bootstrap_pi_forecasts_pi/fig/';
pngFolder = 'frac-cum y_k data/SIRmf03PFITboothstrapPIFullTimeSeriesFitting/indv_county_time_series_real_vs_fit_bootstrap_pi_forecasts_pi/png/';

if ~exist(figFolder, 'dir'), mkdir(figFolder); end
if ~exist(pngFolder, 'dir'), mkdir(pngFolder); end

% === Baseline SIR model ===
colors_top              = {[0.35, 0.35, 0.35],[0.35, 0.35, 0.35], ...
                                [0.35, 0.35, 0.35],[0.35, 0.35, 0.35], ...
                                [0.35, 0.35, 0.35]};
color_rest              = [0.35, 0.35, 0.35];   % charcoal gray
% Shared
real_data_color         = [0 0 1];              % blue
ci_color                = [0.85, 0.85, 0.85];   % gray CI
color_mean_sim_samp_fit = [0.15, 0.15, 0.15];   % near-black
color_sim_samp_fit      = [0.75, 0.75, 0.75];   % light gray ensemble
forecast_scatter_color  = [0.9,0.75,0.3];  % muted golden yellow [0.6,0.4,0.7]; % Muted purple [0.85, 0.33, 0.10];   % muted burnt orange
color_fit               = [0.25, 0.25, 0.25];
marker_size             = 30;
train_ci_color          = [0 1 0];
%% CYAN  [0.4,0.8,0.8]

% Pre acllocate mean of the calibration model + forecast for plotting
mean_across_sim          = zeros(size(lower_bound,1), m);
% Preallocate the mean forecast matrix (time × counties)
mean_forecast_across_sim = zeros(forecast_horizon, m);


for county_index = selected_counties
    optim_beta  = optim_optim_beta_vec(county_index);
    optim_rho   = optim_optim_rho_vec(county_index);
    optim_gamma = optim_optim_gamma_vec(county_index);
    r_squared   = r2_per_county(county_index);
    RMSE        = rmse_per_county(county_index);
    MAE         = mae_per_county(county_index);
    % Data
    t_plot   = 1:L;
    obs_node = covid_county_data(1:L, county_index);
    sim_node = sim_covid_county_data(1:L, county_index);
    y_max = max(max(obs_node), max(sim_node)) + 100;
    % Invisible figure
    hCounty = figure('Visible','off');
    hold on;


    %% Bootstraping prediction interval

    single_county_covid_cases  = covid_county_data(:,county_index);
    N                          = vec_N(county_index);
    index_post_covid_cases     = find(single_county_covid_cases> 0);       % find the nonzero covid-19 cases
    int_index_post_covid_cases = index_post_covid_cases(1); % get the first nonzero covid-19 case
    d                          = int_index_post_covid_cases;
    covid_cases                = single_county_covid_cases(d:end);

    Realdata                   = covid_cases;
    cases                      = Realdata;
    R0    = 0;                     % initial iteration R(0) =0
    % I0 is real rep. case, but practically it might not be available so update by X(2,2)
    I0    = (1/rho)*covid_cases(1);
    S0    = N - I0 - R0;
    % IC  = [N - I0 - R0  I0 R0 I0];   % [S0 I0 R0 I0]
    IC    = [S0 I0 R0 I0];
    %   I0 is real rep. case, but practically it might not be available so update by X(2,2)
    tsol  = [1: length(Realdata)];
    %tsol  = t;
    tdata = tsol;
    tspan = tdata;
    %t     = tsol;
    % we can change the pro9duct vector data [0 0 0 rho]
    %par =[par_beta_xi par_rho];
    %gamma = gamma_lb;
    par              = [optim_beta optim_rho optim_gamma];
    sir_sim_cases    = SIRsol(par(1), par(2), par(3), IC, tspan)
    %% ========================================================================
    %  95% Confidence interval for bestfit curve
    %  ========================================================================
    sim_cases = sir_sim_cases';
    %sim_cases = fun_fitting(x_optimal, tspan);
    F = cumsum(sim_cases); % Compute the cumulative sum of the array z1 and store the result in F.
    length_time = length(tspan);

    %% Additing negative binomial error distribution


    nbd_simulations = zeros(length_time,num_sim);

    % average squared deviation of observations from their mean
    variance_fit = var(cases);
    mean_fit = mean(cases);
    factor1 = variance_fit/mean_fit;


    % Initialize the first row of the simulation matrix
    nbd_simulations(1, :) = sim_cases(1);

    % Perform simulations for each column
    for matrix_col = 1:num_sim
        for matrix_row = 2:length_time
            % Compute new cases at time t
            newcases_t = F(matrix_row) - F(matrix_row - 1);

            % Compute mean and variance
            mean1 = newcases_t;
            var1 = mean1 * factor1;

            % Compute parameters for the negative binomial distribution
            param2 = mean1 / var1;
            param1 = mean1 * param2 / (1 - param2);

            % Generate random value from the negative binomial distribution
            nbd_simulations(matrix_row, matrix_col) = nbinrnd(param1, param2);
        end
    end


    %% model fitting for new data with negative binomial distribution

    variable_data = zeros(num_sim,length(lb));
    fval_data  = zeros(num_sim,1);

    for sim = 1:num_sim

        cases_sim = nbd_simulations(:,sim);
        SIRparsol = @(par, t) SIRsol(par(1), par(2), par(3), IC, t);
        SumSquaresErrorSIRfunction = @(par)sum((SIRparsol(par, tdata)- cases_sim').^2);
        [SIR_optm_theta, SIR_optm_error] = fmincon(SumSquaresErrorSIRfunction, [optim_beta optim_rho optim_gamma], [], [], [], [], lb,ub);

        % variable_guess = x_optimal;
        % [variable, fval] = lsqcurvefit(@fun_fitting,variable_guess,tspan,cases_sim,lb,ub,options);
        % variable_data(sim,:) = variable;
        % fval_data(sim)  = fval;
        variable_data(sim,:) = SIR_optm_theta;
        fval_data(sim)  = SIR_optm_error;

        fprintf('iteration number %d \n',sim);

    end

    %% Calculating 95% CI

    % Store all best fit curves in a matrix; each collumn represents one simulation

    matrix_best_fit = zeros(length(tspan),num_sim);
    for fig = 1:num_sim
        % matrix_best_fit(:,fig) = fun_fitting(variable_data(fig,:),tspan);
        matrix_best_fit(:,fig)    = SIRsol(variable_data(fig,1), variable_data(fig,2), variable_data(fig,3), IC, tspan)
    end



    %data = nbd_simulations;
    data = matrix_best_fit;


    % Desired quantiles to get the 95% CI (2.5% on the lower tail, 2.5% on the upper tail)
    q = [0.025, 0.975];

    % 90% CI
    % q = [0.05, 0.95]

    % 99% CI
    %q = [0.005, 0.995]

    % Compute quantiles along the second dimension (time series)
    quantiles = quantile(data, q, 2);

    % Extract lower and upper bounds
    orginal_lower_bound = quantiles(:, 1); % 2.5th percentile
    orginal_upper_bound = quantiles(:, 2); % 97.5th percentile

    lower_bound = [zeros(d-1,1); orginal_lower_bound];

    upper_bound = [zeros(d-1,1); orginal_upper_bound];

    %% After computing matrix_best_fit and original CIs, add the forecast logic

    %% ========================================================================
    % Extend dynamics for no_week_ahead_forecasts
    % ========================================================================

    t_forecast = (length(tspan)+1) : (length(tspan)+forecast_horizon);

    % Preallocate forecast matrix
    forecast_all_traj = zeros(forecast_horizon, num_sim);

    % Last IC from calibration period
    % sol_all_class_opt_data_fit   = SIRsolAllClass(par(1), par(2), par(3), IC, tspan);
    %
    % S_end = sol_all_class_opt_data_fit(1,end);
    % I_end = sol_all_class_opt_data_fit(2,end);
    % R_end = sol_all_class_opt_data_fit(3,end);
    % P_end = sol_all_class_opt_data_fit(4,end);
    % IC_forecast_base_opt_data_fit = [S_end I_end R_end P_end];

    % Last IC from calibration period
    t_span_and_t_forecast = [tspan t_forecast];
    SIR_final_forecast   = SIRsol(par(1), par(2), par(3), IC, t_span_and_t_forecast);
    final_model_forecast = SIR_final_forecast(t_forecast);
    % Forecast all bootstrap trajectories
    for j = 1:num_sim
        beta_j  = variable_data(j,1);
        rho_j   = variable_data(j,2);
        gamma_j = variable_data(j,3);
        t_span_and_t_forecast = [tspan t_forecast];
        sum_traj_t_span_and_t_forecast = SIRsol(beta_j, rho_j, gamma_j, IC, t_span_and_t_forecast);
        sum_traj_ext                   = sum_traj_t_span_and_t_forecast(t_forecast);
        forecast_all_traj(:,j) = sum_traj_ext(1,:)';  % forecasted infected
    end

    % Compute 95% CI for forecast
    forecast_quantiles = quantile(forecast_all_traj, q, 2);
    forecast_lower = forecast_quantiles(:,1);
    forecast_upper = forecast_quantiles(:,2);

    % Combine calibration + forecast for plotting
    combined_mean = [zeros(d-1,1); mean(matrix_best_fit,2); mean(forecast_all_traj,2)];
    combined_low  = [zeros(d-1,1); orginal_lower_bound; forecast_lower];
    combined_up   = [zeros(d-1,1); orginal_upper_bound; forecast_upper];

    t_full = [t_plot, t_plot(end) + (1:forecast_horizon)];

    mean_across_sim(:,county_index)           = combined_mean(t_plot);
    mean_forecast_across_sim(:, county_index) = combined_mean([t_plot(end) + (1:forecast_horizon)]) ;


    %% Bootstraping prediction interval

    % ------------------------------------------------------------
    %   COLOR LOGIC (matching exactly your previous plot)
    % ------------------------------------------------------------

    if ismember(county_index, top_five_pop_flow_nodes)
        % Determine color for this top county
        k = find(top_five_pop_flow_nodes == county_index);
        c_top = colors_top{k};

        % % Time points for fitting part
        % t_fit = t_plot;
        % Top counties (solid / dashed)
        % plot(t_plot, obs_node, [c_top '-'], 'LineWidth', 3);
        % plot(t_plot, sim_node, [c_top '--'], 'LineWidth', 3);

        %% CI plot
        %  fill([t_plot, fliplr(t_plot)], [lower_bound; flipud(upper_bound)], ...
        %  ci_color, 'FaceAlpha', 0.25, 'EdgeColor', 'none'); % Shaded area of CIs
        % % Plot all trajectories lightly
        % plot(t_plot, [zeros(d-1,num_sim); matrix_best_fit], 'Color', color_sim_samp_fit, 'LineWidth', 0.8)
        % % Plot mean trajectory boldly
        % plot(t_plot, [zeros(d-1,1); mean(matrix_best_fit,2)], 'Color', color_mean_sim_samp_fit, ...
        %   'LineWidth', 3)

        % Light bootstrap trajectories (calibration + forecast)
        plot([t_plot, t_full(end-forecast_horizon+1:end)], [zeros(d-1,num_sim); [matrix_best_fit; forecast_all_traj]], ...
            'Color', color_sim_samp_fit, 'LineWidth', 0.8);

        % ===== Extract Corresponding CI Values fore forecasting domain only =====
        ci_low_forecast = combined_low(end-forecast_horizon:end);
        ci_up_forecast  = combined_up(end-forecast_horizon:end);

        t_new = [t_forecast(1) + d - 2, t_forecast + d - 1];

        % ===== Soft Salmon CI for FORECASTING DOMAIN ONLY =====
        fill([t_new, fliplr(t_new)], ...
            [ci_low_forecast; flipud(ci_up_forecast)], ...
            ci_color, 'FaceAlpha', 0.25, 'EdgeColor', 'none');  % Forecast CI shading

        % ===== CYAN SHADING FOR FITTING DOMAIN ONLY =====
        fill([t_plot, fliplr(t_plot)], [combined_low(1:length(t_plot)); flipud(combined_up(1:length(t_plot)))], ...
            train_ci_color, 'FaceAlpha', 0.25, 'EdgeColor', 'none'); % Shaded CI

        % --- Lower CI boundary (dashed) ---
        plot(t_full, combined_low, ...
            '--', 'Color', color_mean_sim_samp_fit, 'LineWidth', 2);

        % --- Upper CI boundary (dashed) ---
        plot(t_full, combined_up, ...
            '--', 'Color', color_mean_sim_samp_fit, 'LineWidth', 2);

        % Mean trajectory (calibration + forecast)
        plot(t_full, combined_mean, 'Color', color_mean_sim_samp_fit, 'LineWidth', 4);

        scatter(t_forecast_points, combined_mean(t_forecast_points), marker_size, 'x', ... % diamond shape
            'MarkerEdgeColor', forecast_scatter_color, 'LineWidth', 2);

        %plot(1: length(real_covid_county_data(:,county_index)), real_covid_county_data(:,county_index), 'Color', c_top,'LineStyle', '-', 'LineWidth', 3)

        %plot(t_plot, obs_node, [c_top '-'], 'LineWidth', 3);
        %plot(t_plot, sim_node, [c_top '--'], 'LineWidth', 3);
        % Combine existing fitted SIR and forecast
        sim_node_extended = [sim_node; final_model_forecast'];

        % Extend t_plot to include forecast days
        %t_plot_extended = [t_plot t_forecast];

        % Plot fitted + forecast SIR
        %plot(t_full, sim_node_extended, 'Color', color_fit, 'LineStyle', '--', 'LineWidth', 3);

        % Time points for fitting part
        t_fit = t_plot;

        % Time points for forecast part
        t_forecast_points = t_full(length(t_plot)+1:end);
        forecast_vals = final_model_forecast';

        %% --- REAL DATA (obs_node) SCATTER ---
        scatter(1: length(real_covid_county_data(:,county_index)), real_covid_county_data(:,county_index), marker_size, 'o', ...
            'MarkerEdgeColor', real_data_color, 'LineWidth', 2);
        % 'MarkerFaceColor', c_top, ...
        %% --- FITTED MODEL (sim_node) SCATTER ---
        scatter(t_fit, sim_node, marker_size, 's', ...   % squares to distinguish from circles
                     'MarkerEdgeColor', color_fit, 'LineWidth', 0.8);

        %% --- FORECAST (final_model_forecast) SCATTER ---
        % scatter(t_forecast_points, forecast_vals, marker_size, 'x', ... % diamond shape
        %               'MarkerEdgeColor', forecast_scatter_color, 'LineWidth', 2);

    else

        % Light bootstrap trajectories (calibration + forecast)
        plot([t_plot, t_full(end-forecast_horizon+1:end)], [zeros(d-1,num_sim); [matrix_best_fit; forecast_all_traj]], ...
            'Color', color_sim_samp_fit, 'LineWidth', 0.8);

        % ===== Extract Corresponding CI Values fore forecasting domain only =====
        ci_low_forecast = combined_low(end-forecast_horizon:end);
        ci_up_forecast  = combined_up(end-forecast_horizon:end);

        t_new = [t_forecast(1) + d - 2, t_forecast + d - 1];
        % ===== Soft Salmon CI for FORECASTING DOMAIN ONLY =====
        fill([t_new, fliplr(t_new)], ...
            [ci_low_forecast; flipud(ci_up_forecast)], ...
            ci_color, 'FaceAlpha', 0.25, 'EdgeColor', 'none');  % Forecast CI shading

        % ===== CYAN SHADING FOR FITTING DOMAIN ONLY =====
        fill([t_plot, fliplr(t_plot)], [combined_low(1:length(t_plot)); flipud(combined_up(1:length(t_plot)))], ...
            train_ci_color, 'FaceAlpha', 0.25, 'EdgeColor', 'none'); % Shaded CI

        % --- Lower CI boundary (dashed) ---
        plot(t_full, combined_low, ...
            '--', 'Color', color_mean_sim_samp_fit, 'LineWidth', 2);

        % --- Upper CI boundary (dashed) ---
        plot(t_full, combined_up, ...
            '--', 'Color', color_mean_sim_samp_fit, 'LineWidth', 2);

        % Mean trajectory (calibration + forecast)
        plot(t_full, combined_mean, 'Color', color_mean_sim_samp_fit, 'LineWidth', 4);
        % Time points for forecast part

        t_forecast_points = t_full(length(t_plot)+1:end);

        forecast_vals = final_model_forecast';

        scatter(t_forecast_points, combined_mean(t_forecast_points), marker_size, 'x', ... % diamond shape
            'MarkerEdgeColor', forecast_scatter_color, 'LineWidth', 2);

        %plot(1: length(real_covid_county_data(:,county_index)), real_covid_county_data(:,county_index), 'Color', c_top,'LineStyle', '-', 'LineWidth', 3)

        %plot(t_plot, obs_node, [c_top '-'], 'LineWidth', 3);
        %plot(t_plot, sim_node, [c_top '--'], 'LineWidth', 3);
        % Combine existing fitted SIR and forecast
        sim_node_extended = [sim_node; final_model_forecast'];

        % Extend t_plot to include forecast days
        %t_plot_extended = [t_plot t_forecast];

        % Plot fitted + forecast SIR
        %plot(t_full, sim_node_extended, 'Color', color_fit, 'LineStyle', '--', 'LineWidth', 3);

        % Time points for fitting part
        t_fit = t_plot;

        % Time points for forecast part
        t_forecast_points = t_full(length(t_plot)+1:end);
        forecast_vals = final_model_forecast';

        %% --- REAL DATA (obs_node) SCATTER ---
        scatter(1: length(real_covid_county_data(:,county_index)), real_covid_county_data(:,county_index), marker_size, 'o', ...
            'MarkerEdgeColor', real_data_color, 'LineWidth', 2);
        % 'MarkerFaceColor', c_top, ...
        %% --- FITTED MODEL (sim_node) SCATTER ---
         scatter(t_fit, sim_node, marker_size, 's', ...   % squares to distinguish from circles
                     'MarkerEdgeColor', color_fit, 'LineWidth', 0.8);

        %% --- FORECAST (final_model_forecast) SCATTER ---
        % scatter(t_forecast_points, forecast_vals, marker_size, 'x', ... % diamond shape
        %              'MarkerEdgeColor', forecast_scatter_color, 'LineWidth', 2);
    end

    % ------------------------------------------------------------
    %   VISUALS
    % ------------------------------------------------------------
    % xlabel('Time (days)', 'FontSize', 16);
    % ylabel('Infections', 'FontSize', 16);
    %title(['County Index = ' num2str(county_index)], 'FontSize', 30);
    set(gca, 'LineWidth', 1.5, 'FontSize', 30, 'FontName', 'Times');
    box on;
    xline(t_fit(end), 'k--', 'LineWidth', 1);
    xlim([0, 70]);
    xticks(0:10:70);              % Set tick marks every 10 units
    ylim([0, ceil(max(real_covid_county_data(:, county_index)) / 200) * 200]);

    % ------------------------------------------------------------
    %   FILE NAMES
    % ------------------------------------------------------------

    final_model_forecast_data = mean_forecast_across_sim(:,county_index);

    final_real_fcst_window_covid_county_data = real_covid_county_data(t_forecast+d-1);
    % Residuals
    residuals = final_model_forecast_data - final_real_fcst_window_covid_county_data';
    % SSR (sum of squared residuals)
    SSR = sum(residuals.^2);
    % MSE
    mse = mean(residuals.^2);
    % RMSE
    RMSE = sqrt(mse);
    % MAE
    MAE = mean(abs(residuals));
    % Normalized RMSE
    nrmse = rmse / (max(final_real_fcst_window_covid_county_data) - min(final_real_fcst_window_covid_county_data));
    % R^2
    SST = sum((final_real_fcst_window_covid_county_data - mean(final_real_fcst_window_covid_county_data)).^2); % total sum of squares
    r_squared  = 1 - (SSR / SST);

    %     figFile = fullfile(figFolder, ...
    %         ['County_' num2str(county_index) '_opt_beta_' num2str(optim_beta) ...
    %         '_optim_rho_' num2str(optim_rho) '_Net_' str '_inf_' num2str(v) ...
    %         '_optim_gamma_' num2str(optim_gamma) '_r_squared_' num2str(r_squared) '_RMSE_' num2str(RMSE) '_MAE_' num2str(MAE) '.fig']);
    %
    %     pngFile = fullfile(pngFolder, ...
    %         ['County_' num2str(county_index) '_opt_beta_' num2str(optim_beta) ...
    %         '_optim_rho_' num2str(optim_rho) '_Net_' str '_inf_' num2str(v) ...
    %         '_optim_gamma_' num2str(optim_gamma) '_r_squared_' num2str(r_squared) '_RMSE_' num2str(RMSE) '_MAE_' num2str(MAE) '.png']);
    %
    %     figName = sprintf('County_%d_b%.4f_r%.4f_g%.4f_R2_%.4f_RMSE_%.3f_MAE_%.3f.fig', ...
    %     county_index, optim_beta, optim_rho, optim_gamma, r_squared, RMSE, MAE);
    %
    % pngName = sprintf('County_%d_b%.4f_r%.4f_g%.4f_R2_%.4f_RMSE_%.3f_MAE_%.3f.png', ...
    %     county_index, optim_beta, optim_rho, optim_gamma, r_squared, RMSE, MAE);


    figName = sprintf('County_%d_b%.4f_r%.4f_%s_inf_%d_g%.4f_R2_%.4f_RMSE_%.3f_MAE_%.3f.fig', ...
        county_index, optim_beta, optim_rho, str, v, optim_gamma, r_squared, RMSE, MAE);

    pngName = sprintf('County_%d_b%.4f_r%.4f_%s_inf_%d_g%.4f_R2_%.4f_RMSE_%.3f_MAE_%.3f.png', ...
        county_index, optim_beta, optim_rho, str, v, optim_gamma, r_squared, RMSE, MAE);

    figFile = fullfile(figFolder, figName);
    pngFile = fullfile(pngFolder, pngName);

    % Save files
    saveas(hCounty, figFile);
    saveas(hCounty, pngFile);
    close(hCounty);

    %% ================================================================
    %  COLLECT RESULTS FOR THIS COUNTY (5 columns)
    % ================================================================

    real_data_full = real_covid_county_data(:, county_index);

    len_real = length(real_data_full);
    len_comb = length(combined_mean);

    % Pad with NaNs at the end
    % if len_comb < len_real
    %     combined_mean = [combined_mean; NaN(len_real - len_comb, 1)];
    %     combined_low  = [combined_low;  NaN(len_real - len_comb, 1)];
    %     combined_up   = [combined_up;   NaN(len_real - len_comb, 1)];
    %     t_full        = [t_full(:);     NaN(len_real - length(t_full), 1)];
    % end

    transpose_t_full = t_full';
    county_block = [ ...
        transpose_t_full, ...
        real_data_full(transpose_t_full), ...
        combined_mean(:), ...
        combined_low(:), ...
        combined_up(:) ...
        ];

    % Append for this county (5 columns added each loop)
    all_county_results = [all_county_results, county_block];

    hold off
end


% =====================================================================
% After the for-loop ends: Build the master CSV table
% =====================================================================

csvgFolder = 'frac-cum y_k data/SIRmf03PFITboothstrapPIFullTimeSeriesFitting/csv/';
if ~exist(csvgFolder,'dir')
    mkdir(csvgFolder);
end

% all_county_results is numeric: rows × (5*m)
[nrows, ncols]      = size(all_county_results);
n_selected_counties = length(selected_counties);

if ncols ~= 5*n_selected_counties
    error('Column mismatch: expected 5*m columns in all_county_results.');
end

T = table();

for i = 1:n_selected_counties
    county_id = selected_counties(i);  % actual county index (unique county code)

    col_start = (i-1)*5 + 1;
    block = all_county_results(:, col_start:col_start+4);

    T.(['t_full_C' num2str(county_id)]) = block(:,1);
    T.(['real_C'   num2str(county_id)]) = block(:,2);
    T.(['mean_C'   num2str(county_id)]) = block(:,3);
    T.(['low_C'    num2str(county_id)]) = block(:,4);
    T.(['up_C'     num2str(county_id)]) = block(:,5);
end

csvFile = fullfile(csvgFolder, 'All_County_SIR_Results.csv');
writetable(T, csvFile);

fprintf('Saved file: %s\n', csvFile);




%% last part: forecasting window's metric collections
origin_population_with_county_name = readtable('pop_vector_origin.csv', 'Delimiter', ',');

% Initialize metric vectors
forecasting_SSR_per_county   = zeros(1, m);
forecasting_mse_per_county   = zeros(1, m);
forecasting_rmse_per_county  = zeros(1, m);
forecasting_mae_per_county   = zeros(1, m);
forecasting_nrmse_per_county = zeros(1, m);
forecasting_r2_per_county    = zeros(1, m);
forecasting_mape_per_county  = zeros(1, m);
forecasting_pae_per_county   = zeros(1, m);
forecasting_pa_per_county    = zeros(1, m);
forecasting_pe_per_county    = zeros(1, m); % percentage error

% Loop through each county (each column)
for i = 1:m
    % Predicted and real data
    final_model_forecast_data = mean_forecast_across_sim(:,i);
    final_real_fcst_window_covid_county_data = real_covid_county_data(t_forecast+d-1,i);
    % Residuals
    residuals = final_model_forecast_data - final_real_fcst_window_covid_county_data;

    % SSR (sum of squared residuals)
    SSR = sum(residuals.^2);

    % MSE
    mse = mean(residuals.^2);

    % RMSE
    rmse = sqrt(mse);

    % MAE
    mae = mean(abs(residuals));

    % Normalized RMSE
    nrmse = rmse / (max(final_real_fcst_window_covid_county_data) - min(final_real_fcst_window_covid_county_data));

    % R^2
    SST = sum((final_real_fcst_window_covid_county_data - mean(final_real_fcst_window_covid_county_data)).^2); % total sum of squares
    r2  = 1 - (SSR / SST);

    % MAPE (avoid divide by zero)
    mape = mean(abs(residuals) ./ max(final_real_fcst_window_covid_county_data, eps)) * 100;

    % PAE
    pae = sum(abs(residuals)) / sum(final_real_fcst_window_covid_county_data);

    % PA (percentage agreement)
    pa = mean(min(final_model_forecast_data, final_real_fcst_window_covid_county_data) ./ max(final_model_forecast_data, final_real_fcst_window_covid_county_data), 'omitnan');

    % PE (percentage error)
    pe = sum(final_model_forecast_data - final_real_fcst_window_covid_county_data) / sum(final_real_fcst_window_covid_county_data);

    % Store results
    forecasting_SSR_per_county(i)   = SSR;
    forecasting_mse_per_county(i)   = mse;
    forecasting_rmse_per_county(i)  = rmse;
    forecasting_mae_per_county(i)   = mae;
    forecasting_nrmse_per_county(i) = nrmse;
    forecasting_r2_per_county(i)    = r2;
    forecasting_mape_per_county(i)  = mape;
    forecasting_pae_per_county(i)   = pae;
    forecasting_pa_per_county(i)    = pa;
    forecasting_pe_per_county(i)    = pe;
end

% Add metrics as new columns in the table
origin_population_with_county_name.sir_fcst_SSR   = forecasting_SSR_per_county';
origin_population_with_county_name.sir_fcst_mse   = forecasting_mse_per_county';
origin_population_with_county_name.sir_fcst_mae   = forecasting_mae_per_county';
origin_population_with_county_name.sir_fcst_rmse  = forecasting_rmse_per_county';
origin_population_with_county_name.sir_fcst_nrmse = forecasting_nrmse_per_county';
origin_population_with_county_name.sir_fcst_r2    = forecasting_r2_per_county';
origin_population_with_county_name.sir_fcst_mape  = forecasting_mape_per_county';
origin_population_with_county_name.sir_fcst_pae   = forecasting_pae_per_county';
origin_population_with_county_name.sir_fcst_pa    = forecasting_pa_per_county';
origin_population_with_county_name.sir_fcst_pe    = forecasting_pe_per_county';

% Display the updated table
disp(origin_population_with_county_name);

% Save the updated table (CSV)
directory = 'frac-cum y_k data/SIRmf03PFITboothstrapPIFullTimeSeriesFitting';
file_path = fullfile(directory, sprintf( ...
    'ten_metric_sirm_for_sir_forecasting_errors_with_all_counties_opt_beta_%opt_gamma_%opt_rho_%g.csv', ...
    opt_beta, opt_gamma, opt_rho));

if ~exist(directory, 'dir')
    mkdir(directory);
end

writetable(origin_population_with_county_name, file_path);


%% all bootstrapping fitting SIR model parameters saving pipelines 
% par              = [optim_beta optim_rho optim_gamma]

% Number of columns in your matrix
m = size(optim_optim_beta_vec_bootstrapping_sample, 2);

% Create column names: "1", "2", ..., "m"
varNames = string(1:m);

% Convert matrix to table with column names
T1 = array2table(optim_optim_beta_vec_bootstrapping_sample, ...
                'VariableNames', varNames);

% File path (your original format)
file_path = fullfile(directory, sprintf( ...
    'optim_optim_beta_vec_with_col_names_bootstrapping_sample_%g_opt_gamma_%g_opt_rho_%g.csv', ...
    opt_beta, opt_gamma, opt_rho));

% Write table (writes column names)
writetable(T1, file_path);


%% rho 

% Number of columns in your matrix
m = size(optim_optim_rho_vec_bootstrapping_sample, 2);

% Create column names: "1", "2", ..., "m"
varNames = string(1:m);

% Convert matrix to table with column names
T2 = array2table(optim_optim_rho_vec_bootstrapping_sample, ...
                'VariableNames', varNames);

% File path (your original format)
file_path = fullfile(directory, sprintf( ...
    'optim_optim_rho_vec_with_col_names_bootstrapping_sample_%g_opt_gamma_%g_opt_rho_%g.csv', ...
    opt_beta, opt_gamma, opt_rho));

% Write table (writes column names)
writetable(T2, file_path);


%% gamma 
% Number of columns in your matrix
m = size(optim_optim_gamma_vec_bootstrapping_sample, 2);

% Create column names: "1", "2", ..., "m"
varNames = string(1:m);

% Convert matrix to table with column names
T3 = array2table(optim_optim_gamma_vec_bootstrapping_sample, ...
                'VariableNames', varNames);

% File path (your original format)
file_path = fullfile(directory, sprintf( ...
    'optim_optim_gamma_vec_with_col_names_bootstrapping_sample_%g_opt_gamma_%g_opt_rho_%g.csv', ...
    opt_beta, opt_gamma, opt_rho));

% Write table (writes column names)
writetable(T3, file_path);