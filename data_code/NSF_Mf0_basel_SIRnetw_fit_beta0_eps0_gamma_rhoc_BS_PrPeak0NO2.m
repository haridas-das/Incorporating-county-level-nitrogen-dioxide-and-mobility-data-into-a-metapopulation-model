clear all; close all; clc;

%% Model fitting with 80 parameters: bline_beta0; gamma, and rhoc
[state, m, name_origin_list, time_start, lockdown_date, top_five_pop_flow_nodes, time_step_size , zeta_vec, zeta,  Nv, gamma, beta_1_vec, phi_0] = county_name_list();
zeta = 1;
%% select model based on NO_2 and high populous
vari_selection = 'NO2';

%% Sort by variability (Std_NO2) in descending order
%top5_var_index = 'Std_NO2';
%% Sort by variability (Std_NO2) in descending order
top5_var_index = 'Mean_NO2';
% top5_var_index = 'vec_N';
% Load county level NO_2 data
data = readtable('NO2_county_data_2020.csv');

%% Origin population from the data
origin_population                  = csvread('pop_vector_origin_only.csv',1,0);
origin_population_with_county_name = readtable('pop_vector_origin.csv', 'Delimiter', ',');
vec_N                              = origin_population(:,2);            % population of each nodes

% Define prepandemic cutoff date
prepandemic_cutoff                 = datetime('2020-03-01');
% Filter for Oklahoma only
oklahoma_data = data(strcmp(data.StateName, 'Oklahoma'), :);
% Get unique county names
counties = unique(oklahoma_data.CountyName);
num_counties = numel(counties);

% Initialize vectors for mean and variability
mean_NO2 = zeros(num_counties, 1);
std_NO2  = zeros(num_counties, 1);  % Standard deviation

% Ensure Date is datetime
if ~isdatetime(oklahoma_data.Date)
    oklahoma_data.Date = datetime(oklahoma_data.Date);
end

% Loop through each county
for i = 1:num_counties
    county_name = counties{i};

    % Filter data for the current county
    county_data = oklahoma_data(strcmp(oklahoma_data.CountyName, county_name), :);

    % Select only pre-pandemic data
    county_prepandemic = county_data(county_data.Date <= prepandemic_cutoff, :);

    % Compute mean and std deviation, ignoring NaNs
    mean_NO2(i) = mean(county_prepandemic.NO2_column_number_density, 'omitnan');
    std_NO2(i)  = std(county_prepandemic.NO2_column_number_density, 'omitnan');
end
county_idx = linspace(1,num_counties,num_counties)';
% Create table with county name, mean, and variability
NO2_stats_table = table(county_idx, counties, mean_NO2, std_NO2,vec_N, ...
    'VariableNames', {'county_idx', 'CountyName', 'Mean_NO2', 'Std_NO2','vec_N'});

[~, sorted_idx] = sort(NO2_stats_table.(top5_var_index), 'descend');
% Extract top 5 counties by NO2 variability
top5_variability = sorted_idx(1:5, :)

if strcmpi(vari_selection, 'NO2')  % or any condition that checks if variability depends on NO2
    top_five_pop_flow_nodes = top5_variability;
else
    top_five_pop_flow_nodes = top_five_pop_flow_nodes';  % keeps previous value
end


% Plot bar chart for top 5 counties with three possible variability

h = figure;
bar(NO2_stats_table.(top5_var_index)(top_five_pop_flow_nodes), 'FaceColor', [0.8 0.4 0.2]);
% set(gca, 'XTickLabel', NO2_stats_table.CountyName(top5_variability), 'XTick', 1:5);
% Combine index and county name
xt_labels = strcat(string(top_five_pop_flow_nodes), ': ', NO2_stats_table.CountyName(top_five_pop_flow_nodes));
% Set XTick and labels
set(gca, 'XTick', 1:5, 'XTickLabel', xt_labels);
xtickangle(45);
xlabel('County Name');
% Force both into char
% ylabel(sprintf('%s (%s)', char(top5_var_index), char(vari_selection)), 'Interpreter', 'latex');
ylabel(sprintf('%s (%s)', char(top5_var_index), char(vari_selection)), 'Interpreter', 'latex');
title(sprintf('%s (%s)', char(top5_var_index), char(vari_selection)), 'Interpreter', 'latex');
grid on;
% % Generate file paths for each figure
% % Convert numeric array to strings and join with underscores
node_str = strjoin(string(top_five_pop_flow_nodes(1:5)), '_');
%% Build PNG filename safely

fig_filename = sprintf(['frac-cum y_k data/', ...
    'SIRnetwMf03PFITboothstrapPIPrePeak/county_data_fiiting_constant_beta_function_all_counties/fig/', ...
    '0_Fig_Bar_Chart_Beta_0_%s_%s_%s.fig'], ...
    vari_selection, top5_var_index, node_str);

png_filename = sprintf(['frac-cum y_k data/', ...
    'SIRnetwMf03PFITboothstrapPIPrePeak/county_data_fiiting_constant_beta_function_all_counties/png/', ...
    '0_Fig_Bar_Chart_Beta_0_%s_%s_%s.png'], ...
    vari_selection, top5_var_index, node_str);


saveas(h, fig_filename);
saveas(h, png_filename);
% Close the figure to free memory
close(h);


% Plot histogram of Mean_NO2 values
h2 = figure(2);
histogram(NO2_stats_table.Mean_NO2, 10); % You can change the number of bins (e.g., 10)
xlabel('Mean NO_2 Column Density');
ylabel('Number of Counties');
%title('Histogram of Mean NO_2 Column Density across Oklahoma Counties');
% Convert numeric array to strings and join with underscores
node_str = strjoin(string(top_five_pop_flow_nodes(1:5)), '_');

fig_filename = sprintf(['frac-cum y_k data/', ...
    'SIRnetwMf03PFITboothstrapPIPrePeak/county_data_fiiting_constant_beta_function_all_counties/fig/', ...
    '0_Fig_histogram_Beta_0_%s_%s_%s.fig'], ...
    vari_selection, top5_var_index, node_str);

png_filename = sprintf(['frac-cum y_k data/', ...
    'SIRnetwMf03PFITboothstrapPIPrePeak/county_data_fiiting_constant_beta_function_all_counties/png/', ...
    '0_Fig_histogram_Beta_0_%s_%s_%s.png'], ...
    vari_selection, top5_var_index, node_str);


saveas(h2, fig_filename);
saveas(h2, png_filename);
% Close the figure to free memory
close(h2);



%% Cellphone mobility data driven SIRnetwork model input
%f = 6; % Put f=0 for Gen Flux, f=1 for star shape; f=2 for fully connected,
M = m; % Number of nodes
N_o = m;
N_d = m;
nodes_vec = 1:m;
str = 'ok_county_wise_network';
%covid_data = readmatrix('final_weekly_division_col4_col11_wise_covid_data_2023_apr1_to_2024_apr6.csv');
covid_data = readmatrix('final_weekly_time_series_covid_data.csv');
% first covid cases in Oklahoma State is in Tulsa county on March 6
% On weekly data it is basically on March 7.
% since the first state lockdown started on March 28 so first 4 weeks
% simulation used the average mobility from the Jan 1 to March 28.
% rest of the time we used avg mobility from March 29 to Decmeber 31st, 2020

%% Fig shows all the counties tiem series together
all_counties_full_ts = covid_data(:,2:end);
% List of counties to highlight
highlight_nodes    = [top_five_pop_flow_nodes', 60, 42, 43];
% Corresponding county names
highlight_counties = {'Oklahoma', 'Tulsa', 'Cleveland', 'Canadian', 'Comanche', 'Payne', 'Logan', 'Love'};

% Define highlight colors (make sure they pop in print)

highlight_colors = [
    1, 0, 0;                  % 'r'  red
    0, 0, 1;                  % 'b'  blue
    0, 1, 0;                  % 'g'  green
    0, 0, 0;                  % 'k'  black
    1, 0, 1;                  % 'm'  magenta
    0.30, 0.70, 1.00;         % color_rest (light blue)
    1.0, 0.6, 0.6;            % ci_color (soft peach)
    0.4, 0.0, 0.6;            % color_mean_sim_samp_fit (deep purple)
    0.90, 0.88, 0.98          % color_sim_samp_fit (light lavender)
    ];

T = size(all_counties_full_ts, 1);          % number of time points

hft = figure('Position', [100, 100, 1800, 800]); % [left, bottom, width, height]
hold on;


% Plot background counties and save handle for legend
hGray = [];  % To store the handle for "Others"
for c = 1:size(all_counties_full_ts, 2)
    if ~ismember(c, highlight_nodes)
        h = plot(1:T, all_counties_full_ts(:, c), 'Color', [0.85, 0.85, 0.85], 'LineWidth', 3);
        if isempty(hGray)
            hGray = h; % Save the handle of the first gray line for legend
        end
    end
end

% Plot highlighted counties and store their handles
hHighlights = gobjects(length(highlight_nodes),1);
for k = 1:length(highlight_nodes)
    idx = highlight_nodes(k);
    hHighlights(k) = plot(1:T, all_counties_full_ts(:, idx), ...
        'Color', highlight_colors(k,:), ...
        'LineWidth', 3, ...
        'Marker', 'o', ...
        'MarkerSize', 4, ...
        'MarkerFaceColor', highlight_colors(k,:));
end

xlabel('Time (days)');
ylabel('COVID cases');
set(gca, 'FontSize', 30, 'FontName', 'Times', 'LineWidth', 2, 'Box', 'on');
ylim([0,16000])
xlim([1, 160])
%legend(arrayfun(@(x) sprintf('County %d', x), highlight_nodes, 'UniformOutput', false), ...
% 'Location', 'best');
% Add legend with custom county names
%legend(highlight_counties, 'Location', 'best', 'FontSize', 14);
% Pass only the handles of highlighted plots to the legend
% legend(hHighlights, highlight_counties, 'Location', 'best', 'FontSize', 14, 'Interpreter', 'none');
% Build legend entries
legend_handles = [hHighlights; hGray];
legend_labels = [highlight_counties, {'Other counties'}];
legend(legend_handles, legend_labels, 'Location', 'best', 'FontSize', 14, 'Interpreter', 'none');
hold off;

% Generate file paths for each figure

fig_filename = sprintf(['frac-cum y_k data/', ...
    'SIRnetwMf03PFITboothstrapPIPrePeak/county_data_fiiting_constant_beta_function_all_counties/fig/', ...
    'Fig_all_counties_real_cases_together_time_series_matrix_size_%d.fig'], ...
    m);

png_filename = sprintf(['frac-cum y_k data/', ...
    'SIRnetwMf03PFITboothstrapPIPrePeak/county_data_fiiting_constant_beta_function_all_counties/png/', ...
    'Fig_all_counties_real_cases_together_time_series_matrix_size_%d.png'], ...
    m);

% Save the figure in .fig and .png formats
saveas(hft, fig_filename);
saveas(hft, png_filename);
% Close the figure to free memory
close(hft);

%% Model calibration and loading of observed data
covid_data = readmatrix('final_weekly_time_series_covid_data.csv');
original_covid_data = readmatrix('final_weekly_time_series_covid_data.csv');
all_counties_full_ts = covid_data(:,2:end);
%% here we can expand the time of the real data
covid_data = covid_data(7:36,:);
% Define dengue_div_data from covid_data (columns 4 to 11)
covid_county_data = covid_data(:,2:end);  % Real data (time along rows, nodes along columns)
real_covid_data        = original_covid_data(7:76,:);
real_covid_county_data = real_covid_data(:,2:end);  % Real data (time along rows, nodes along columns)
time_end_plot           = size(covid_data, 1);
% Convert lockdown_date to 0-indexed for MATLAB slicing (adjusting for 1-based indexing)
lockdown_date_0_indexed = lockdown_date;
lockdown_week_0_indexed = 4;   % this is in week

u = 1;
v = 1;
% time_start = 1;
time_end  = time_end_plot;
T         = time_end;
ylimit        = 0.3;
bline_beta0   = 0.5;
baseline_rho  = 0.1; %10^-1; % 1
time_step_size = 1; % can compute N_time_steps
tspan         = 1:time_step_size:time_end;  % Ensure this is within the bounds of `t`
t             = tspan; %linspace(1, time_end_plot, time_end_plot); % Time vector (0 to 20)
%f_vec   = 1.3*ones(1,m); % beta_vec values
%rho_vec = rho*ones(1,m);
rhoc          = baseline_rho;
% epsilon       = 10^-1;
%% Range and Infected nodes defined for the simulations
%simulation by the mean mobility flux value 0.00651 from 2020 mobility data
%phi_0 = 0.00651;
phi_0      = 1;  % 1 means gravity model flux values
beta_1_vec = linspace(0,0.90,1); % simulate the infection rates
gamma = 1/1; % Recovery rate 1/1 means that infected people recover in 1 week if we consider weekly data
% psi   = 1/10;
% mean NO_2 data
NO2_data = NO2_stats_table.Mean_NO2;

% I0    = (1/rho)*covid_cases(1);
vec_I0 = (1/baseline_rho)*covid_county_data(1,:);
% vec_N =N*ones(m,1);
%% Initial value of the Problem by the following subroutine function
[vec_S0,vec_I0,vec_N,y0] =  DataDrivenInitialValueMobilitySIRNetwork(m,vec_N,vec_I0);
%% computing flux matrix for pre and post pandemic:
time_step_size_pre_post_covid_flux = 1; % can compute N_time_steps
% tspan = 1:time_step_size_pre_post_covid_flux:366;
time_end_mobility_flux   =  366;
orig_matrix_phi          = TDFluxMatrix(N_o,N_d,time_step_size_pre_post_covid_flux,time_start,time_end_mobility_flux,name_origin_list);

% Calculate the average along the third dimension (time)
avg_matrix               = mean(orig_matrix_phi, 1:lockdown_date);

% Slice the matrix from time steps 1 to lockdown_date
pre_lockdown_flux_matrix = orig_matrix_phi(:, :, 1:lockdown_date_0_indexed);

% Compute the average along the third dimension (time steps)
avg_pre_lockdown_flux_matrix = mean(pre_lockdown_flux_matrix, 3);

% Slice the matrix from time steps 1 to lockdown_date
post_lockdown_flux_matrix = orig_matrix_phi(:, :, lockdown_date_0_indexed+1:366);

% Compute the average along the third dimension (time steps)
avg_post_lockdown_flux_matrix = mean(post_lockdown_flux_matrix, 3);

%% checking that the row sum is 1
for i =1: m
    sum(avg_post_lockdown_flux_matrix(i,:));
end

for i =1: m
    sum(avg_pre_lockdown_flux_matrix(i,:));
end

%% Precompute the flux matrix function
% for the simulated time steps 32 to later is the post lockdown avg mobility data
phi_matrix_func = @(t) get_matrix_phi(t, lockdown_week_0_indexed, time_start, time_step_size, avg_pre_lockdown_flux_matrix, avg_post_lockdown_flux_matrix);

[t, y_k] = ode45(@(t, y) SIR_Network_WithoutNO2_baseB0(t, m, phi_matrix_func(t), y, bline_beta0, vec_S0, vec_I0, vec_N, gamma, zeta), tspan, y0);

h00 = figure(1)
plot(1: size(covid_county_data, 1),covid_county_data(1: size(covid_county_data, 1),4),'r-','MarkerFaceColor','r','LineWidth',3);
hold on
set(gca,'LineWidth',2);
hold on
plot(1: size(covid_county_data, 1),covid_county_data(1: size(covid_county_data, 1),5),'b-','MarkerFaceColor','b','LineWidth',3);
plot(1: size(covid_county_data, 1),covid_county_data(1: size(covid_county_data, 1),6),'g-','MarkerFaceColor','g','LineWidth',3);
plot(1: size(covid_county_data, 1),covid_county_data(1: size(covid_county_data, 1),7),'k-','MarkerFaceColor','k','LineWidth',3);
plot(1: size(covid_county_data, 1),covid_county_data(1: size(covid_county_data, 1),8),'m-','MarkerFaceColor','m','LineWidth',3);
plot(1: size(covid_county_data, 1),covid_county_data(1: size(covid_county_data, 1),9),'Color','#145A32','LineStyle', '-','MarkerFaceColor','m','LineWidth',3);
plot(1: size(covid_county_data, 1),covid_county_data(1: size(covid_county_data, 1),10),'Color','#EDB120','LineStyle', '-','MarkerFaceColor','m','LineWidth',3);
plot(1: size(covid_county_data, 1),covid_county_data(1: size(covid_county_data, 1),11),'Color','#7E2F8E','LineStyle', '-','MarkerFaceColor','m','LineWidth',3);
%plot(1: size(covid_county_data, 1),covid_county_data(1: size(covid_county_data, 1),13),'r--','MarkerFaceColor','r','LineWidth',3);
%plot(t,y_k(:,m+top_five_pop_flow_nodes(6))/vec_N(m),'k-','MarkerFaceColor','k','LineWidth',3);
%grid on;
% set(gca, 'fontsize', 20, 'fontweight','bold'); % gca: get current axis
%set(gca, 'fontsize', 30); % gca: get current axis
%set(gca, 'YTick', 0:.1:ylimit); % Full form of gca is (Get Current Axes)
% Set X-axis ticks and labels
% tick_positions1 = 1:88:time_end;
% tick_positions11 = [tick_positions1, time_end];  % Add end time to the tick positions
% tick_labels = num2cell(tick_positions11);  % Convert tick positions to cell array of strings
% set(gca, 'XTick', tick_positions11, 'XTickLabel', tick_labels); % Set the tick positions and labels
% set(gca, 'YTick', 0:split_y_ylimit_non_norm:ylimit_non_norm); % Full form of gca is (Get Current Axes)
% set(gca, 'XTick', [0:11:time_end, time_end])%set(gca, 'XTick', [0:split_x_xlimit_non_norm:time_end_plot, time_end_plot])

%ylim([0,ylimit]);
% ylim([0,ylimit]);
xlim([1,time_end_plot])
set(gca, 'fontsize', 30,'LineWidth',2,'FontName', 'Times','box','on'); % gca: get current axis
% legend only sum of total cases in each day node sum
% legend(['$I_{' num2str(top_five_pop_flow_nodes(1)) '}$'], ...
%        ['$I_{' num2str(top_five_pop_flow_nodes(2)) '}$'], ...
%        ['$I_{' num2str(top_five_pop_flow_nodes(3)) '}$'], ...
%        ['$I_{' num2str(top_five_pop_flow_nodes(4)) '}$'], ...
%        ['$$I_{' num2str(top_five_pop_flow_nodes(5)) '}$$'], ...
%        '$$\sum_{i=1}^{m} I_i$$', 'Interpreter', 'latex', ...
%        'Location', 'Best', 'NumColumns', 2);
% %       '$$\sum_{i=1}^{m} I_i/\sum_{i=1}^{m} N_i$$', 'Interpreter', 'latex', ...
% %       'Location', 'Best', 'NumColumns', 2);
hold off
h00_figurename_fig = ['/frac-cum y_k data/SIRnetwMf03PFITboothstrapPIPrePeak/fig/Fig_real_data_dynmics_cons_beta_TD_flux_' num2str(phi_0) '_matrix_size_' num2str(m) '_Network_' str '_Size_' num2str(T) '_infected_strat_node_' num2str(v) '_Het_Node_Mult_' num2str(zeta) '_beta_i_' num2str(bline_beta0) '_GT_' num2str(gamma) '.fig'];
h00_figurename_png = ['/frac-cum y_k data/SIRnetwMf03PFITboothstrapPIPrePeak/png/Fig_real_data_dynmics_cons_beta_TD_flux_' num2str(phi_0) '_matrix_size' num2str(m) '_Network_' str '_Size_' num2str(T) '_infected_start_node_' num2str(v) '_Het_Node_Mult_' num2str(zeta) '_beta_i_' num2str(bline_beta0) '_GT_' num2str(gamma) '.png'];
saveas(h00,[pwd h00_figurename_fig])
saveas(h00,[pwd h00_figurename_png])


h2 = figure(2)
plot(t,y_k(:,m+top_five_pop_flow_nodes(1)),'r-','MarkerFaceColor','r','LineWidth',3);
hold on
set(gca,'LineWidth',2);
hold on
plot(t,y_k(:,m+top_five_pop_flow_nodes(2)),'b-','MarkerFaceColor','b','LineWidth',3);
plot(t,y_k(:,m+top_five_pop_flow_nodes(3)),'g-','MarkerFaceColor','g','LineWidth',3);
plot(t,y_k(:,m+top_five_pop_flow_nodes(4)),'k-','MarkerFaceColor','k','LineWidth',3);
plot(t,y_k(:,m+top_five_pop_flow_nodes(5)),'m-','MarkerFaceColor','m','LineWidth',3);
% plot(t,y_k(:,m+top_five_pop_flow_nodes(6)),'Color','#145A32','LineStyle', '-','MarkerFaceColor','m','LineWidth',3);
% plot(t,y_k(:,m+top_five_pop_flow_nodes(7)),'Color','#EDB120','LineStyle', '-','MarkerFaceColor','m','LineWidth',3);
% plot(t,y_k(:,m+top_five_pop_flow_nodes(8)),'Color','#7E2F8E','LineStyle', '-','MarkerFaceColor','m','LineWidth',3);
plot(t,sum(y_k(:,m+1:2*m),2),'r--','MarkerFaceColor','r','LineWidth',3);
%plot(t,y_k(:,m+top_five_pop_flow_nodes(6))/vec_N(m),'k-','MarkerFaceColor','k','LineWidth',3);
%grid on;
% set(gca, 'fontsize', 20, 'fontweight','bold'); % gca: get current axis
%set(gca, 'fontsize', 30); % gca: get current axis
%set(gca, 'YTick', 0:.1:ylimit); % Full form of gca is (Get Current Axes)
% Set X-axis ticks and labels
% tick_positions1 = 1:88:time_end;
% tick_positions11 = [tick_positions1, time_end];  % Add end time to the tick positions
% tick_labels = num2cell(tick_positions11);  % Convert tick positions to cell array of strings
% set(gca, 'XTick', tick_positions11, 'XTickLabel', tick_labels); % Set the tick positions and labels
%set(gca, 'YTick', 0:1000000:9000000); % Full form of gca is (Get Current Axes)
% set(gca, 'XTick', 1:11: time_end); % Full form of gca is (Get Current Axes)
%set(gca, 'XTick', [0:11:time_end, time_end])
%ylim([0,ylimit]);
% ylim([0,ylimit]);
xlim([1,time_end])
set(gca, 'fontsize', 30,'LineWidth',2,'FontName', 'Times','box','on'); % gca: get current axis
% legend only sum of total cases in each day node sum
% legend(['$I_{' num2str(top_five_pop_flow_nodes(1)) '}$'], ...
%        ['$I_{' num2str(top_five_pop_flow_nodes(2)) '}$'], ...
%        ['$I_{' num2str(top_five_pop_flow_nodes(3)) '}$'], ...
%        ['$I_{' num2str(top_five_pop_flow_nodes(4)) '}$'], ...
%        ['$$I_{' num2str(top_five_pop_flow_nodes(5)) '}$$'], ...
%        '$$\sum_{i=1}^{m} I_i$$', 'Interpreter', 'latex', ...
%        'Location', 'Best', 'NumColumns', 2);
% %       '$$\sum_{i=1}^{m} I_i/\sum_{i=1}^{m} N_i$$', 'Interpreter', 'latex', ...
% %       'Location', 'Best', 'NumColumns', 2);
hold off
h2_figurename_fig = ['/frac-cum y_k data/SIRnetwMf03PFITboothstrapPIPrePeak/fig/Fig_A_B_TD_non_frac_dynmics_Time_dependent_beta_TD_flux_' num2str(phi_0) '_matrix_size_' num2str(m) '_Network_' str '_Size_' num2str(T) '_infected_strat_node_' num2str(v) '_Het_Node_Mult_' num2str(zeta) '_beta_i_' num2str(bline_beta0) '_GT_' num2str(gamma) '.fig'];
h2_figurename_png = ['/frac-cum y_k data/SIRnetwMf03PFITboothstrapPIPrePeak/png/Fig_A_B_TD_non_frac_dynmics_Time_dependent_beta_TD_flux_' num2str(phi_0) '_matrix_size' num2str(m) '_Network_' str '_Size_' num2str(T) '_infected_start_node_' num2str(v) '_Het_Node_Mult_' num2str(zeta) '_beta_i_' num2str(bline_beta0) '_GT_' num2str(gamma) '.png'];
saveas(h2,[pwd h2_figurename_fig])
saveas(h2,[pwd h2_figurename_png])
% baseline_rho is rhoc which is also optimizing with the

%% Model caliberation data
P_init          = [bline_beta0, gamma, rhoc];
% Define the length of P_init
n               = length(P_init);  % This should be m + 3
%Lower Bound≤Initial Value≤Upper Bound
%% fitting bline_beta0, gamma works great

lb        = 0.0001*ones(1, n);     % Set the lower bounds for each parameter (0 or a reasonable value)
ub        = 2.5*ones(1, n);
lb(1)     = 0.01;                 % Lower bounds for bline_beta0
ub(1)     = 0.25 + max(NO2_stats_table.Mean_NO2);  % Upper bounds for bline_beta0
lb(2)     = 0.0001;               % Lower bounds for gamma
ub(2)     = 2;                    % Upper bounds for gamma
lb(3)     = 0.0001;               % Lower bounds for rhoc
ub(3)     = 0.4;                  % Upper bounds for rhoc
% reporting rate upper bound
% K         = 0.1;
% rho_T     = 0.4;
% ub(4)     = K;
% Compute epsilon bound for each county
% Take the global upper bound (minimum across counties)
% NO2_scaled_for_ub = (rho_T - K) ./ NO2_data;
% epsilon_ub = min(NO2_scaled_for_ub);
% % epsilon_ub = max(epsilon_bounds);
% lb(5)     = 0.01;              % Lower bounds for rho_vec
% ub(5)     = epsilon_ub;


% Number of time points in the real data
tt                       = linspace(1, time_end, time_end);  % Adjust to match `dengue_div_data`
covid_county_data_interp = interp1(tt, covid_county_data, tspan, 'linear');
objectiveSSR             = @(P) Compute_SSR_with_SIR_Network_NO2_baseB0_WithoutNO2_rhoc(P, covid_county_data_interp, t, m, lockdown_week_0_indexed, time_start, time_step_size, avg_pre_lockdown_flux_matrix, avg_post_lockdown_flux_matrix, vec_S0, vec_I0, vec_N,zeta, y0);
% Optimization to minimize SSR with bounds
options                  = optimset('Display', 'iter', 'TolFun', 1e-6, 'TolX', 1e-6);

%% optimizing over i , j and k
% beta_vec  = 0.5:.5:2.5;     % vector of beta values
% gamma_vec = 0.5:.5:2.5;     % vector of gamma values

%% optimizing over no. of iterations
no_iteration = 200;
% allocating space for initial values and the optimum variables
optim_SSE_error_matrix  = zeros(no_iteration,1);
optim_beta0_matrix      = zeros(no_iteration,1);
optim_gamma_matrix      = zeros(no_iteration,1);
optim_rho_vec_matrix    = zeros(no_iteration,1);

%% we are optimizing beta_vec, gamma, and rho_vec using a random process with 25 times 5 = 125  iterations
for i = 1: no_iteration
    i
    % rand funct generates pseudorandom numbers from a uniform distribution
    X = rand(n,1);
    new_randm_para  = lb + (ub-lb).*X.';
    bline_beta0     = new_randm_para(1);
    gamma           = new_randm_para(2);
    rhoc            = new_randm_para(3);
    P_init          = [bline_beta0, gamma, rhoc];

    [P_opt, SIR_optm_error] = fmincon(objectiveSSR, P_init, [], [], [], [], lb, ub, [], options);

    optim_SSE_error_matrix(i)   = SIR_optm_error;
    optim_beta0_matrix(i)       = P_opt(1);
    optim_gamma_matrix(i)       = P_opt(2);
    optim_rho_vec_matrix(i)     = P_opt(3);
end

% Define your target directory
directory = 'frac-cum y_k data/SIRnetwMf03PFITboothstrapPIPrePeak';

% Create directory if it doesn't exist
if ~exist(directory, 'dir')
    mkdir(directory);
end

% Save each matrix in that directory
writematrix(optim_SSE_error_matrix, fullfile(directory, 'optim_SSE_error_matrix.csv'));
writematrix(optim_beta0_matrix,      fullfile(directory, 'optim_beta0_matrix.csv'));
writematrix(optim_gamma_matrix,     fullfile(directory, 'optim_gamma_matrix.csv'));
writematrix(optim_rho_vec_matrix,     fullfile(directory, 'optim_rho_vec_matrix.csv'));


%% find the minimum error and its corresponding index and find corresponding beta and gamma
min_error                 = min(min(optim_SSE_error_matrix));
[row_index, column_index] = find(optim_SSE_error_matrix == min_error);
row_idx                   = min(row_index);

bline_beta0 = optim_beta0_matrix(row_idx);
gamma       = optim_gamma_matrix(row_idx);
% rhoc        = squeeze(optim_rho_vec_matrix(row_idx));
rhoc        = optim_rho_vec_matrix(row_idx);
P_opt       = [bline_beta0, gamma, rhoc];
optim_beta  = bline_beta0;
optim_rho   = rhoc;
optim_gamma = gamma;

%% Display optimized parameters
disp('Optimized Parameters:');
disp(P_opt);
disp(['Minimum SSE: ', num2str(min_error)]);

%% optimal data fitting curve
bline_beta0 = P_opt(1);
gamma       = P_opt(2);
rhoc        = P_opt(3);

[t, y_k] = ode45(@(t, y) SIR_Network_WithoutNO2_baseB0(t, m, phi_matrix_func(t), y, bline_beta0, vec_S0, vec_I0, vec_N, gamma, zeta), tspan, y0);

rho_optc = rhoc*ones(m,1);
vec_beta = bline_beta0*ones(m,1);



h3 = figure(3)
plot(t,y_k(:,m+top_five_pop_flow_nodes(1)),'r-','MarkerFaceColor','r','LineWidth',3);
hold on
set(gca,'LineWidth',2);
hold on
plot(t,y_k(:,m+top_five_pop_flow_nodes(2)),'b-','MarkerFaceColor','b','LineWidth',3);
plot(t,y_k(:,m+top_five_pop_flow_nodes(3)),'g-','MarkerFaceColor','g','LineWidth',3);
plot(t,y_k(:,m+top_five_pop_flow_nodes(4)),'k-','MarkerFaceColor','k','LineWidth',3);
plot(t,y_k(:,m+top_five_pop_flow_nodes(5)),'m-','MarkerFaceColor','m','LineWidth',3);
% plot(t,y_k(:,m+top_five_pop_flow_nodes(6)),'Color','#145A32','LineStyle', '-','MarkerFaceColor','m','LineWidth',3);
% plot(t,y_k(:,m+top_five_pop_flow_nodes(7)),'Color','#EDB120','LineStyle', '-','MarkerFaceColor','m','LineWidth',3);
% plot(t,y_k(:,m+top_five_pop_flow_nodes(8)),'Color','#7E2F8E','LineStyle', '-','MarkerFaceColor','m','LineWidth',3);
% %plot(t,sum(y_k(:,m+1:2*m),2),'r--','MarkerFaceColor','r','LineWidth',3);
%plot(t,y_k(:,m+top_five_pop_flow_nodes(6))/vec_N(m),'k-','MarkerFaceColor','k','LineWidth',3);
%grid on;
% set(gca, 'fontsize', 20, 'fontweight','bold'); % gca: get current axis
%set(gca, 'fontsize', 30); % gca: get current axis
%set(gca, 'YTick', 0:.1:ylimit); % Full form of gca is (Get Current Axes)
% Set X-axis ticks and labels
% tick_positions1 = 1:88:time_end;
% tick_positions11 = [tick_positions1, time_end];  % Add end time to the tick positions
% tick_labels = num2cell(tick_positions11);  % Convert tick positions to cell array of strings
% set(gca, 'XTick', tick_positions11, 'XTickLabel', tick_labels); % Set the tick positions and labels
set(gca, 'YTick', 0:1000000:9000000); % Full form of gca is (Get Current Axes)
% set(gca, 'XTick', 1:11: time_end); % Full form of gca is (Get Current Axes)
set(gca, 'XTick', [0:11:time_end, time_end]);
%ylim([0,ylimit]);
% ylim([0,ylimit]);
xlim([1,time_end])
set(gca, 'fontsize', 30,'LineWidth',2,'FontName', 'Times','box','on'); % gca: get current axis
h3_figurename_fig = ['/frac-cum y_k data/SIRnetwMf03PFITboothstrapPIPrePeak/fig/Fig_real_data_fitting_dynmics_Time_dependent_beta_TD_flux_' num2str(phi_0) '_matrix_size_' num2str(m) '_Network_' str '_Size_' num2str(T) '_infected_strat_node_' num2str(v) '_Het_Node_Mult_' num2str(zeta) '_beta_i_' num2str(bline_beta0) '_GT_' num2str(gamma) '.fig'];
h3_figurename_png = ['/frac-cum y_k data/SIRnetwMf03PFITboothstrapPIPrePeak/png/Fig_real_data_fitting_dynmics_Time_dependent_beta_TD_flux_' num2str(phi_0) '_matrix_size' num2str(m) '_Network_' str '_Size_' num2str(T) '_infected_start_node_' num2str(v) '_Het_Node_Mult_' num2str(zeta) '_beta_i_' num2str(bline_beta0) '_GT_' num2str(gamma) '.png'];
saveas(h3,[pwd h3_figurename_fig])
saveas(h3,[pwd h3_figurename_png])
hold off


h4 = figure(4);
hold on;
% Real data as scatter points with circles
scatter(1:size(covid_county_data, 1), covid_county_data(:, top_five_pop_flow_nodes(1)), 80, 'r', 'o', 'filled', 'MarkerEdgeColor', 'r');
scatter(1:size(covid_county_data, 1), covid_county_data(:, top_five_pop_flow_nodes(1)), 80, 'b', 'o', 'filled', 'MarkerEdgeColor', 'b');
scatter(1:size(covid_county_data, 1), covid_county_data(:, top_five_pop_flow_nodes(1)), 80, 'g', 'o', 'filled', 'MarkerEdgeColor', 'g');
scatter(1:size(covid_county_data, 1), covid_county_data(:, top_five_pop_flow_nodes(1)), 80, 'k', 'o', 'filled', 'MarkerEdgeColor', 'k');
scatter(1:size(covid_county_data, 1), covid_county_data(:, top_five_pop_flow_nodes(1)), 80, 'm', 'o', 'filled', 'MarkerEdgeColor', 'm');
% Fitting results as lines
plot(t, rho_optc(top_five_pop_flow_nodes(1))*y_k(:, m + top_five_pop_flow_nodes(1)), 'r-', 'LineWidth', 3);
plot(t, rho_optc(top_five_pop_flow_nodes(2))*y_k(:, m + top_five_pop_flow_nodes(2)), 'b-', 'LineWidth', 3);
plot(t, rho_optc(top_five_pop_flow_nodes(3))*y_k(:, m + top_five_pop_flow_nodes(3)), 'g-', 'LineWidth', 3);
plot(t, rho_optc(top_five_pop_flow_nodes(4))*y_k(:, m + top_five_pop_flow_nodes(4)), 'k-', 'LineWidth', 3);
plot(t, rho_optc(top_five_pop_flow_nodes(5))*y_k(:, m + top_five_pop_flow_nodes(5)), 'm-', 'LineWidth', 3);
set(gca, 'XTick', [0:11:time_end, time_end]);
%ylim([0,ylimit]);
% ylim([0,ylimit]);
xlim([1,time_end])
% Formatting the axes
set(gca, 'fontsize', 30,'LineWidth',2,'FontName', 'Times','box','on'); % gca: get current axis
xlabel('Time');
ylabel('Cases');
h4_figurename_fig = ['/frac-cum y_k data/SIRnetwMf03PFITboothstrapPIPrePeak/county_data_fiiting_constant_beta_function_8_counties/fig/Fig_real_data_fitting_dynmics_with_scatter_Time_dependent_beta_TD_flux_' num2str(phi_0) '_matrix_size_' num2str(m) '_Network_' str '_Size_' num2str(T) '_infected_strat_node_' num2str(v) '_Het_Node_Mult_' num2str(zeta) '_beta_i_' num2str(bline_beta0) '_GT_' num2str(gamma) '.fig'];
h4_figurename_png = ['/frac-cum y_k data/SIRnetwMf03PFITboothstrapPIPrePeak/county_data_fiiting_constant_beta_function_8_counties/png/Fig_real_data_fitting_dynmics_with_scatter_Time_dependent_beta_TD_flux_' num2str(phi_0) '_matrix_size' num2str(m) '_Network_' str '__Size_' num2str(T) '_infected_start_node_' num2str(v) '_Het_Node_Mult_' num2str(zeta) '_beta_i_' num2str(bline_beta0) '_GT_' num2str(gamma) '.png'];
saveas(h4,[pwd h4_figurename_fig])
saveas(h4,[pwd h4_figurename_png])
hold off


%% Create figure
h5 = figure(5);
hold on;

% Define the colors for scatter and fits
colors = {'r', 'b', 'g', 'k', 'm', [20/255, 90/255, 50/255], [237/255, 178/255, 32/255], [126/255, 47/255, 142/255]};

% Create 8 subplots (4 rows, 2 columns)
numSubplots = 8;
rows = 4;
cols = 2;
%top_five_pop_flow_nodes_plus_60_42_43 = [top_five_pop_flow_nodes 60 42 43]; % 60 - payne , 42-logan, 43- love,
top_five_pop_flow_nodes_plus_60_42_43 = [55 72 14 9 16 60 42 43]; % 60 - payne , 42-logan, 43- love,

for i = 1:numSubplots
    subplot(rows, cols, i); % Create subplot for each pair

    % Real data scatter points for the pair with corresponding color
    scatter(1:size(covid_county_data, 1), covid_county_data(:, top_five_pop_flow_nodes_plus_60_42_43(i)), 80, 'o', 'MarkerEdgeColor', colors{i});
    hold on;

    % Fitting results for the pair with corresponding color
    plot(t, rho_optc(top_five_pop_flow_nodes_plus_60_42_43(i))* y_k(:, m + top_five_pop_flow_nodes_plus_60_42_43(i)), 'Color', colors{i}, 'LineWidth', 3);

    % Adjust axis and formatting for each subplot
    %     set(gca, 'XTick', [1:11:time_end, time_end]);
    set(gca, 'XTick', [1    12    23    34    45  58  65]);
    xlim([1, time_end]);
    %     xlabel('Time');
    %     ylabel('Cases');
    set(gca, 'fontsize', 15, 'LineWidth', 2, 'FontName', 'Times', 'box', 'on');
end

h5_figurename_fig = ['/frac-cum y_k data/SIRnetwMf03PFITboothstrapPIPrePeak/county_data_fiiting_constant_beta_function_8_counties/fig/Fig_subplots_real_data_fitting_dynmics_with_scatter_Time_dependent_beta_TD_flux_' num2str(phi_0) '_matrix_size_' num2str(m) '_Network_' str '_Size_' num2str(T) '_infected_strat_node_' num2str(v) '_Het_Node_Mult_' num2str(zeta) '_beta_i_' num2str(bline_beta0) '_GT_' num2str(gamma) '.fig'];
h5_figurename_png = ['/frac-cum y_k data/SIRnetwMf03PFITboothstrapPIPrePeak/county_data_fiiting_constant_beta_function_8_counties/png/Fig_subplots_real_data_fitting_dynmics_with_scatter_Time_dependent_beta_TD_flux_' num2str(phi_0) '_matrix_size' num2str(m) '_Network_' str '__Size_' num2str(T) '_infected_start_node_' num2str(v) '_Het_Node_Mult_' num2str(zeta) '_beta_i_' num2str(bline_beta0) '_GT_' num2str(gamma) '.png'];
saveas(h5,[pwd h5_figurename_fig])
saveas(h5,[pwd h5_figurename_png])
hold off


%% all counties
% Nodes to highlight for eight counties
highlight_nodes = [55 72 14 9 16 60 42 43];

% Define base colors for highlighted nodes (same length as highlight_nodes)
highlight_colors = [...
    1, 0, 0;                  % Red
    0, 0, 1;                  % Blue
    0, 1, 0;                  % Green
    0, 0, 0;                  % Black
    1, 0, 1;                  % Magenta
    20/255, 90/255, 50/255;   % Dark Green
    237/255, 178/255, 32/255; % Golden Yellow
    126/255, 47/255, 142/255  % Purple
    ];

% Initialize color cell
m_colors = cell(1, m);

% Assign base colors to highlighted nodes
for k = 1:length(highlight_nodes)
    m_colors{highlight_nodes(k)} = highlight_colors(k,:);
end

% Assign remaining colors from parula
remaining_nodes = setdiff(1:m, highlight_nodes);
num_extra = length(remaining_nodes);
extra_colors = parula(num_extra);

for k = 1:num_extra
    m_colors{remaining_nodes(k)} = extra_colors(k,:);
end

% Now m_colors{i} is a valid RGB triplet for all i = 1..m
disp(m_colors);

% Z-scores for prediction intervals
z_scores = struct('CI50', 0.674, 'CI80', 1.282, 'CI90', 1.645, 'CI95', 1.960);

for i = 1:m
    h = figure('Visible', 'off'); % batch processing
    hold on;

    % Scatter plot actual data
    %scatter(1:size(covid_county_data, 1), covid_county_data(:, nodes_vec(i)), 150, 'o', 'MarkerEdgeColor', m_colors{i}, 'LineWidth', 3);

    % Fitted curve
    fitted_curve = rho_optc(nodes_vec(i))*y_k(:, m + nodes_vec(i));
    %plot(1:size(covid_county_data, 1), fitted_curve, 'Color', m_colors{i}, 'LineWidth', 3);
    % Calculate residuals and std deviation
    errors = covid_county_data(:, nodes_vec(i)) - fitted_curve;
    std_err = std(errors);

    x_vals = 1:size(covid_county_data, 1);

    % Compute intervals for each CI level
    lower_95 = fitted_curve - z_scores.CI95 * std_err;
    upper_95 = fitted_curve + z_scores.CI95 * std_err;

    lower_90 = fitted_curve - z_scores.CI90 * std_err;
    upper_90 = fitted_curve + z_scores.CI90 * std_err;

    lower_80 = fitted_curve - z_scores.CI80 * std_err;
    upper_80 = fitted_curve + z_scores.CI80 * std_err;

    lower_50 = fitted_curve - z_scores.CI50 * std_err;
    upper_50 = fitted_curve + z_scores.CI50 * std_err;

    % Plot filled areas from widest CI to narrowest to avoid overlap issues
    % Updated colors for confidence intervals (brighter and publication-friendly)
    ci_colors = struct(...
        'CI95', [1.0, 0.95, 0.7], ...   % Light warm yellow
        'CI90', [1.0, 0.8, 0.5], ...  % Warm orange
        'CI80', [0.5, 0.9, 1.0], ...  % Soft sky blue
        'CI50', [0.3, 0.75, 1.0]);      % Brighter blue

    % And in the plotting loop, replace fill commands with:

    fill([x_vals, fliplr(x_vals)], [lower_95', fliplr(upper_95')], ci_colors.CI95, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    fill([x_vals, fliplr(x_vals)], [lower_90', fliplr(upper_90')], ci_colors.CI90, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    fill([x_vals, fliplr(x_vals)], [lower_80', fliplr(upper_80')], ci_colors.CI80, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    fill([x_vals, fliplr(x_vals)], [lower_50', fliplr(upper_50')], ci_colors.CI50, 'FaceAlpha', 0.3, 'EdgeColor', 'none');

    % Plot fitted curve **after CI fills** (to bring on top)
    plot(x_vals, fitted_curve, 'Color', m_colors{i}, 'LineWidth', 3);

    % Plot real data scatter **last** (on top)
    scatter(x_vals, covid_county_data(:, nodes_vec(i)), 150, 'o', 'MarkerEdgeColor', m_colors{i}, 'LineWidth', 3);

    % Axis formatting
    set(gca, 'XTick', [1 12 23 34 45 58 65]);
    xlim([1, time_end]);
    ylim([0, max(upper_95)*1.1]); % Slightly above upper 95% CI max

    set(gca, 'fontsize', 30, 'LineWidth', 2, 'FontName', 'Times', 'box', 'on');
    hold off;

    % Save figures
    fig_filename = fullfile(pwd, 'frac-cum y_k data', 'SIRnetwMf03PFITboothstrapPIPrePeak', 'county_data_fiiting_constant_beta_function_all_counties', 'fig', ...
        sprintf('Fig_individual_dynamic_with_CI_for_Node_%d_opt_rho_%g_Network_%s_rho_opt_vec_%g_gamma_%g.fig', i, rho_optc(i), str, rho_optc(i), gamma));
    png_filename = fullfile(pwd, 'frac-cum y_k data', 'SIRnetwMf03PFITboothstrapPIPrePeak', 'county_data_fiiting_constant_beta_function_all_counties', 'png', ...
        sprintf('Fig_individual_dynamic_with_CI_for_Node_%d_opt_rho_%g_Network_%s_rho_opt_vec_%g_gamma_%g.png', i, rho_optc(i), str, rho_optc(i), gamma));

    saveas(h, fig_filename);
    saveas(h, png_filename);

    close(h);
end

%% compute the errors for different metrices
% Extract infected populations from ODE solution
simulated_optimal_infected_data_all_infected_class = y_k(:, m+1:2*m); % Infected data for all nodes

% Ensure rho_vec is a row vector
rho_opt_vec_row = rho_optc'; % Convert rho_vec to row vector if not already

% Multiply each column of simulated_data by the corresponding scalar in rho_vec
rho_vec_times_sim_opt_all_infected_class = simulated_optimal_infected_data_all_infected_class.* rho_opt_vec_row;

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
    opt_soln = rho_vec_times_sim_opt_all_infected_class(:,i);
    Realdata = covid_county_data(:,i);
    T = length(Realdata);

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
origin_population_with_county_name.sir_net_fit_SSR   = SSR_per_county';
origin_population_with_county_name.sir_net_fit_mse   = mse_per_county';
origin_population_with_county_name.sir_net_fit_mae   = mae_per_county';
origin_population_with_county_name.sir_net_fit_rmse  = rmse_per_county';
origin_population_with_county_name.sir_net_fit_nrmse = nrmse_per_county';
origin_population_with_county_name.sir_net_fit_r2    = r2_per_county';
origin_population_with_county_name.sir_net_fit_mape  = mape_per_county';
origin_population_with_county_name.sir_net_fit_pae   = pae_per_county';
origin_population_with_county_name.sir_net_fit_pa    = pa_per_county';
origin_population_with_county_name.sir_net_fit_pe    = pe_per_county';

% Display the updated table
disp(origin_population_with_county_name);

% Save the updated table (CSV)
directory = 'frac-cum y_k data/SIRnetwMf03PFITboothstrapPIPrePeak';
file_path = fullfile(directory, sprintf( ...
    'ten_metric_sir_network_for_constant_beta_01_function_fitting_errors_with_all_counties_vari_selection_%s_top5_var_index_%s_beta0_%gamma_%g.csv', ...
    vari_selection, top5_var_index, bline_beta0, gamma));

if ~exist(directory, 'dir')
    mkdir(directory);
end

writetable(origin_population_with_county_name, file_path);




%% plot the reporting rate on each node
node_vec = linspace(1, m,m)
h11 = figure(7)
%I=linspace(1,11,11)
data = rho_optc' ;
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
h11_figurename_fig = ['/frac-cum y_k data/SIRnetwMf03PFITboothstrapPIPrePeak/county_data_fiiting_constant_beta_function_8_counties/fig/Fig_2_C_D_rho_opt_vec_vs_Node_Fig2_m_Time_dependent_beta_' num2str(m) '_Network_' str '_Size_' num2str(time_end) '_infected_start_node_' num2str(v) '_Het_Node' num2str(u) '_RR_' num2str(gamma) '.fig'];
saveas(h11,[pwd h11_figurename_fig])
h11_figurename_png = ['/frac-cum y_k data/SIRnetwMf03PFITboothstrapPIPrePeak/county_data_fiiting_constant_beta_function_8_counties/png/Fig_2_C_D_rho_opt_vec_vs_Node_Fig2_m_Time_dependent_beta_' num2str(m) '_Network_' str '_Size_' num2str(time_end) '_infected_start_node_' num2str(v) '_Het_Node' num2str(u) '_RR_' num2str(gamma) '.png'];
saveas(h11,[pwd h11_figurename_png])


%% plot the opt_vec_beta on each node
node_vec = linspace(1, m,m)
h12 = figure(8)
%I=linspace(1,11,11)
%stem(data,'b','filled','LineWidth',2)
stem(node_vec, vec_beta,'r','filled','LineWidth',3)
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
h12_figurename_fig = ['/frac-cum y_k data/SIRnetwMf03PFITboothstrapPIPrePeak/county_data_fiiting_constant_beta_function_8_counties/fig/Fig_opt_vec_beta_vs_Node_Fig2_m_Time_dependent_beta_' num2str(m) '_Network_' str '_Size_' num2str(time_end) '_infected_start_node_' num2str(v) '_Het_Node' num2str(u) '_RR_' num2str(gamma) '.fig'];
saveas(h12,[pwd h12_figurename_fig])
h12_figurename_png = ['/frac-cum y_k data/SIRnetwMf03PFITboothstrapPIPrePeak/county_data_fiiting_constant_beta_function_8_counties/png/Fig_opt_vec_beta_vs_Node_Fig2_m_Time_dependent_beta_' num2str(m) '_Network_' str '_Size_' num2str(time_end) '_infected_start_node_' num2str(v) '_Het_Node' num2str(u) '_RR_' num2str(gamma) '.png'];
saveas(h12,[pwd h12_figurename_png])


%% save all the optimum variables bline_beta0, gamma, and rhoc are scalars
% rho_vec is a vector of length m
% Create a cell array of parameter names
param_names = cell(m+3, 1);
param_names{1} = 'bline_beta0';
param_names{2} = 'gamma';
param_names{3} = 'rhoc';

for i = 1:m
    param_names{3+i} = sprintf('vec_beta_%d', i)
end


% Combine parameter values into a column vector
param_values = [bline_beta0; gamma; rhoc; vec_beta];
% Create table with two columns: names and values
opt_fit_table = table(param_names, param_values, 'VariableNames', {'Parameter', 'opt_fit_pars'});
% Define directory and file name
directory = 'frac-cum y_k data/SIRnetwMf03PFITboothstrapPIPrePeak';
file_path = fullfile(directory, sprintf('opt_fit_params_gamma_%g_rhoCount_%d.csv', gamma, m));

% Ensure the directory exists
if ~exist(directory, 'dir')
    mkdir(directory);
end

% Write the table to CSV file
writetable(opt_fit_table, file_path);


% Preallocate arrays for opt_soln and Realdata
nT = size(rho_vec_times_sim_opt_all_infected_class, 1); % number of time steps
opt_soln_all = zeros(nT, m);
Realdata_all = zeros(nT, m);

% Collect time series for each county
for i = 1:m
    opt_soln_all(:,i) = rho_vec_times_sim_opt_all_infected_class(:,i);
    Realdata_all(:,i) = covid_county_data(:,i);
end

% Build table with time column
Time = (1:nT)'; % or use your actual time vector if available
all_results = table(Time);

% Add each county’s data as separate columns
for i = 1:m
    all_results.(sprintf('OptSoln_County%d', i))  = opt_soln_all(:,i);
    all_results.(sprintf('RealData_County%d', i)) = Realdata_all(:,i);
end

% Define directory and file path
directory = 'frac-cum y_k data/SIRnetwMf03PFITboothstrapPIPrePeak';
file_path = fullfile(directory, sprintf( ...
    'opt_soln_vs_realdata_side_by_side_for_all_counties_vari_selection_%s_top5_var_index_%s_beta0_%g_gamma_1_%g.csv', ...
    vari_selection, top5_var_index, bline_beta0, gamma));

% Make directory if not existing
if ~exist(directory, 'dir')
    mkdir(directory);
end

% Save wide-format CSV
writetable(all_results, file_path);


%% ==============================================================================================================
%  Fig 2. with forecast SAVE EACH COUNTY'S INDIVIDUAL TIME SERIES WITH BOOTHSTRAPPING PI and FORECASTS PI USING BOOTHSTRAPPING
% ===============================================================================================================
% Preallocate table copunty level data
all_county_results      = [];   % Initialize once before the loop
num_sim                 = 200;    % will try with 100
no_week_ahead_forecasts = 5;
length_time             = nT;
forecast_horizon        = no_week_ahead_forecasts;  % number of extra weeks to forecast


figFolder = 'frac-cum y_k data/SIRnetwMf03PFITboothstrapPIPrePeak/indv_county_time_series_real_vs_fit_bootstrap_pi_forecasts_pi/fig/';
pngFolder = 'frac-cum y_k data/SIRnetwMf03PFITboothstrapPIPrePeak/indv_county_time_series_real_vs_fit_bootstrap_pi_forecasts_pi/png/';

if ~exist(figFolder, 'dir'), mkdir(figFolder); end
if ~exist(pngFolder, 'dir'), mkdir(pngFolder); end

colors_top              = {[0.30, 0.70, 1.00],[0.30, 0.70, 1.00],[0.30, 0.70, 1.00],[0.30, 0.70, 1.00],[0.30, 0.70, 1.00]};% same as before
color_rest              = [0.30, 0.70, 1.00];   % light blue (not gray) for line curves
real_data_color         = [0 0 1];
ci_color                = [1.0, 0.6, 0.6];      % Soft Peach for confidence interval
color_mean_sim_samp_fit = [1, 0, 0]; % [0.4, 0.0, 0.6];      % deep purple
%color_sim_samp_fit      = [0.75, 0.75, 0.25];   % olive-gold
color_sim_samp_fit       = [0.90, 0.88, 0.98];  % lighter
forecast_scatter_color   = [0 1 0]; % green [0.95, 0.75, 0.1];    % Vivid Gold (contrasts with everything)
color_fit                = [0,0.6,0.6]; % teal
marker_size              = 30;

%% change start here

%% Additing negative binomial distribution for num_sim and each coiunties

nbd_simulations = zeros(length_time,num_sim, m);

selected_counties = 1:m;

for county_index =  selected_counties

    sim_cases = rho_vec_times_sim_opt_all_infected_class(:,county_index) % sir_sim_cases';
    cases     = covid_county_data(:,county_index);
    %sim_cases = fun_fitting(x_optimal, tspan);
    F = cumsum(sim_cases); % Compute the cumulative sum of the array z1 and store the result in F.
    % length_time = length(tspan);

    %% Additing negative binomial error distribution for ecah counties
    length_time = length(sim_cases);
    % num_sim = 10;
    % average squared deviation of observations from their mean
    variance_fit = var(cases);
    mean_fit     = mean(cases);
    factor1      = variance_fit/mean_fit;
    % Initialize the first row of the simulation matrix
    % nbd_simulations(1, :) = sim_cases(1);
    nbd_simulations(1, :,county_index) = sim_cases(1);

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
            nbd_simulations(matrix_row, matrix_col,county_index) = nbinrnd(param1, param2);
        end
    end
end
% nbd_simulations(:,:,77)

%% model fitting for new data with negative binomial distribution

variable_data = zeros(num_sim,length(lb));
fval_data  = zeros(num_sim,1);

for sim = 1:num_sim
    cases_sim                = nbd_simulations(:,sim,:); % becomes nT × 1 × m
    % covid_cases_sim          =  squeeze(covid_county_data_interp);   % becomes 65 × 77
    covid_cases_sim          =  squeeze(cases_sim);   % becomes nT × m
    % Optimization to minimize SSR with bounds
    options                  = optimset('Display', 'iter', 'TolFun', 1e-6, 'TolX', 1e-6);
    P_init                   = [bline_beta0, gamma, rhoc];
    objectiveSSR             = @(P) Compute_SSR_with_SIR_Network_NO2_baseB0_WithoutNO2_rhoc(P, covid_cases_sim, t, m, lockdown_week_0_indexed, time_start, time_step_size, avg_pre_lockdown_flux_matrix, avg_post_lockdown_flux_matrix, vec_S0, vec_I0, vec_N,zeta, y0);
    % Optimization to minimize SSR with bounds
    [P_opt, SIR_optm_error]  = fmincon(objectiveSSR, P_init, [], [], [], [], lb, ub, [], options);
    variable_data(sim,:)     = P_opt;
    fval_data(sim)           = SIR_optm_error;
    fprintf('sim steps or iteration number %d \n',sim);
end

%% Calculating 95% CI

% Store all best fit curves in a matrix; 2nd index represents simulation, county in the third index

matrix_best_fit              = zeros(length(t), num_sim, m);
for fig = 1:num_sim
    % matrix_best_fit(:,fig)    = SIRsol(variable_data(fig,1), variable_data(fig,2), variable_data(fig,3), IC, tspan)
    [t, y_k] = ode45(@(t, y) SIR_Network_WithoutNO2_baseB0(t, m, phi_matrix_func(t), y, variable_data(fig,1), vec_S0, vec_I0, vec_N, variable_data(fig,2), zeta), t, y0);
    % compute the errors for different metrices
    % Extract infected populations from ODE solution
    sum_simulated_optimal_infected_data_all_infected_class = y_k(:, m+1:2*m); % Infected data for all nodes
    sumulated_rho_optc  = variable_data(fig,3)*ones(m,1);
    % Ensure rho_vec is a row vector
    sim_rho_opt_vec_row = sumulated_rho_optc'; % Convert rho_vec to row vector if not already
    % Multiply each column of simulated_data by the corresponding scalar in rho_vec
    sim_rho_vec_times_sim_opt_all_infected_class = sum_simulated_optimal_infected_data_all_infected_class.* sim_rho_opt_vec_row;
    matrix_best_fit(:,fig, :)                    = sim_rho_vec_times_sim_opt_all_infected_class;
end

data = matrix_best_fit;

% Desired quantiles to get the 95% CI (2.5% on the lower tail, 2.5% on the upper tail)
q = [0.025, 0.975];
% 90% CI
% q = [0.05, 0.95]
% 99% CI
%q = [0.005, 0.995]


%% Extract lower and upper bounds
% collecting quantiles for each county: Compute quantiles along the second dimension (time series)
quantiles = quantile(data, q, 2);
% Extract lower and upper bounds for all counties

orginal_lower_bound = quantiles(:, 1, :); % 2.5th percentile

orginal_upper_bound = quantiles(:, 2, :); % 97.5th percentile

lower_bound = squeeze(orginal_lower_bound);  % time × counties (remove dimension of length 1)

upper_bound = squeeze(orginal_upper_bound);  % time × counties

%% After computing matrix_best_fit and original CIs, add the forecast logic

%% ========================================================================
% Extend dynamics for no_week_ahead_forecasts
% ========================================================================

t_forecast = (length(t)+1) : (length(t)+forecast_horizon);

forecast_all_traj = zeros(forecast_horizon, num_sim, m);
% Last IC from calibration period
t_span_and_t_forecast = [t' t_forecast];


for fig = 1:num_sim
    [t, y_k] = ode45(@(t, y) SIR_Network_WithoutNO2_baseB0(t, m, phi_matrix_func(t), y, variable_data(fig,1), vec_S0, vec_I0, vec_N, variable_data(fig,2), zeta), t_span_and_t_forecast, y0);

    % Extract infected populations from ODE solution
    sum_simulated_optimal_infected_data_all_infected_class = y_k(:, m+1:2*m); % Infected data for all nodes

    sumulated_rho_optc  = variable_data(fig,3)*ones(m,1);
    % Ensure rho_vec is a row vector
    sim_rho_opt_vec_row = sumulated_rho_optc'; % Convert rho_vec to row vector if not already

    % Multiply each column of simulated_data by the corresponding scalar in rho_vec
    sim_rho_vec_times_sim_opt_all_infected_class = sum_simulated_optimal_infected_data_all_infected_class.* sim_rho_opt_vec_row;
    % matrix_best_fit(:,fig, :)                    = sim_rho_vec_times_sim_opt_all_infected_class;
    sum_traj_ext               = sim_rho_vec_times_sim_opt_all_infected_class(t_forecast, :);
    forecast_all_traj(:,fig,:) = sum_traj_ext % (1,:)';  % forecasted infected

end


% Compute 95% CI for forecast
forecast_quantiles = quantile(forecast_all_traj, q, 2);

forecast_lower_bound = forecast_quantiles(:,1, :);

forecast_upper_bound = forecast_quantiles(:,2, :);

forecast_lower = squeeze(forecast_lower_bound); % 2.5th percentile  time × counties
forecast_upper = squeeze(forecast_upper_bound); % 97.5th percentile time × counties


% Combine calibration model + forecast for plotting

mean_across_sim = zeros(size(matrix_best_fit,1), size(matrix_best_fit,3));

for county = 1:size(matrix_best_fit,3)
    mean_across_sim(:,county) = mean(matrix_best_fit(:, :, county), 2);
end

% Preallocate the mean matrix (time × counties)
mean_forecast_across_sim = zeros(size(forecast_all_traj, 1), size(forecast_all_traj, 3));

% Loop through each county and compute mean across simulations
for county = 1:size(forecast_all_traj, 3)
    mean_forecast_across_sim(:, county) = mean(forecast_all_traj(:, :, county), 2);
end


combined_mean = [mean_across_sim; mean_forecast_across_sim];               % Vertical concatenation
combined_low  = [lower_bound; forecast_lower];
combined_up   = [upper_bound; forecast_upper];

t_plot = 1: time_end_plot;
t_full = [t_plot, t_plot(end) + (1:forecast_horizon)];


%% With Bootstraping prediction interval

% ------------------------------------------------------------
%  Final Fig 2:  COLOR LOGIC (matching exactly your previous plot)
% ------------------------------------------------------------

for county_index = selected_counties
    %% Invisible figure
    hCounty = figure('Visible','off');   % ✅ NEW FIGURE PER COUNTY
    hold on;
    if ismember(county_index, top_five_pop_flow_nodes)
        % Determine color for this top county
        k = find(top_five_pop_flow_nodes == county_index);
        %c_top = colors_top{k};

        % Time points for fitting part
        t_fit = t_plot;

        % ===== Soft SalmonG FOR FITTING DOMAIN =====
        % fill([t_full, fliplr(t_full)], [combined_low(:,county_index); flipud(combined_up(:,county_index))], ...
        % ci_color, 'FaceAlpha', 0.25, 'EdgeColor', 'none'); % Shaded CI

        % Light bootstrapping all trajectories (calibration + forecast)
        plot([t_plot, t_full(end-forecast_horizon+1:end)], [matrix_best_fit(:,:,county_index); forecast_all_traj(:,:,county_index)], ...
            'Color', color_sim_samp_fit, 'LineWidth', 0.8);

        % ===== Extract Corresponding CI Values fore forecasting domain only =====
        ci_low_forecast = combined_low(end-forecast_horizon:end, county_index);
        ci_up_forecast  = combined_up(end-forecast_horizon:end, county_index);

        t_new = [t_forecast(1) - 1, t_forecast];

        % ===== Soft Salmon CI for FORECASTING DOMAIN ONLY =====
        fill([t_new, fliplr(t_new)], ...
            [ci_low_forecast; flipud(ci_up_forecast)], ...
            ci_color, 'FaceAlpha', 0.25, 'EdgeColor', 'none');  % Forecast CI shading

        % ===== CYAN SHADING FOR FITTING DOMAIN =====
        fill([t_plot, fliplr(t_plot)], ...
            [combined_low(1:length(t_plot), county_index); ...
            flipud(combined_up(1:length(t_plot), county_index))], ...
            'c', 'FaceAlpha', 0.25, 'EdgeColor', 'none');

        % --- Lower CI boundary (dashed) ---
        plot(t_full, combined_low(:,county_index), ...
            '--', 'Color', color_mean_sim_samp_fit, 'LineWidth', 2);

        % --- Upper CI boundary (dashed) ---
        plot(t_full, combined_up(:,county_index), ...
            '--', 'Color', color_mean_sim_samp_fit, 'LineWidth', 2);


        % Mean trajectory (calibration + forecast)
        plot(t_full, combined_mean(:,county_index), 'Color', color_mean_sim_samp_fit, 'LineWidth', 4);
        % plot(1: length(real_covid_county_data(:,county_index)), real_covid_county_data(:,county_index), 'Color', real_data_color,'LineStyle', '-', 'LineWidth', 1)

        % Combine existing fitted SIR-network for each county sepeartely and forecast
        sim_node             = rho_vec_times_sim_opt_all_infected_class(:,county_index);
        final_model_forecast = mean_forecast_across_sim(:,county_index);
        sim_node_plus_forecast_extended    = [sim_node; final_model_forecast];

        % Extend t_plot to include forecast days
        % Plot fitted + forecast SIR]
        % plot(t_full, sim_node_plus_forecast_extended, 'Color', color_fit, 'LineStyle', '-', 'LineWidth', 3);

        % Time points for fitting part
        t_fit = t_plot;

        % Time points for forecast part
        t_forecast_points = t_full(length(t_plot)+1:end);
        forecast_vals = final_model_forecast;

        %% --- REAL DATA (obs_node) SCATTER ---
        scatter(1: length(real_covid_county_data(:,county_index)), real_covid_county_data(:,county_index), marker_size, 'o', ...
            'MarkerEdgeColor', real_data_color, 'LineWidth', 2);
        % 'MarkerFaceColor', c_top, ...
        %% --- FITTED MODEL (sim_node) SCATTER ---
        % scatter(t_fit, sim_node, marker_size, 's', ...   % squares to distinguish from circles
        %            'MarkerEdgeColor', color_fit, 'LineWidth', 1.2);

        %% --- FORECAST (final_model_forecast) SCATTER ---
        % scatter(t_forecast_points, forecast_vals, marker_size, 'x', ... % diamond shape
        %               'MarkerEdgeColor', forecast_scatter_color, 'LineWidth', 2);
        scatter(t_forecast_points(:), forecast_vals(:), marker_size, 'x', ...
            'MarkerEdgeColor', forecast_scatter_color, 'LineWidth', 2);

    else

        % ===== Soft SalmonG FOR FITTING DOMAIN =====
        % fill([t_full, fliplr(t_full)], [combined_low(:,county_index); flipud(combined_up(:,county_index))], ...
        % ci_color, 'FaceAlpha', 0.25, 'EdgeColor', 'none'); % Shaded CI

        % Time points for fitting part
        t_fit = t_plot;

        sim_node             = rho_vec_times_sim_opt_all_infected_class(:,county_index);
        final_model_forecast = mean_forecast_across_sim(:,county_index);
        sim_node_plus_forecast_extended    = [sim_node; final_model_forecast];


        % Light bootstrap trajectories (calibration + forecast)
        plot([t_plot, t_full(end-forecast_horizon+1:end)], [matrix_best_fit(:,:,county_index); forecast_all_traj(:,:,county_index)], ...
            'Color', color_sim_samp_fit, 'LineWidth', 0.8);

        % ===== Extract Corresponding CI Values fore forecasting domain only =====
        ci_low_forecast = combined_low(end-forecast_horizon:end, county_index);
        ci_up_forecast  = combined_up(end-forecast_horizon:end, county_index);

        t_new = [t_forecast(1) - 1, t_forecast];

        % ===== Soft Salmon CI for FORECASTING DOMAIN ONLY =====
        fill([t_new, fliplr(t_new)], ...
            [ci_low_forecast; flipud(ci_up_forecast)], ...
            ci_color, 'FaceAlpha', 0.25, 'EdgeColor', 'none');  % Forecast CI shading


        % ===== CYAN SHADING FOR FITTING DOMAIN =====
        fill([t_plot, fliplr(t_plot)], ...
            [combined_low(1:length(t_plot), county_index); ...
            flipud(combined_up(1:length(t_plot), county_index))], ...
            'c', 'FaceAlpha', 0.25, 'EdgeColor', 'none');

        % --- Lower CI boundary (dashed) ---
        plot(t_full, combined_low(:,county_index), ...
            '--', 'Color', color_mean_sim_samp_fit, 'LineWidth', 2);

        % --- Upper CI boundary (dashed) ---
        plot(t_full, combined_up(:,county_index), ...
            '--', 'Color', color_mean_sim_samp_fit, 'LineWidth', 2);


        % Mean trajectory (calibration + forecast)
        plot(t_full, combined_mean(:,county_index), 'Color', color_mean_sim_samp_fit, 'LineWidth', 4);
        % plot(1: length(real_covid_county_data(:,county_index)), real_covid_county_data(:,county_index), 'Color', real_data_color,'LineStyle', '-', 'LineWidth', 1)

        % Extend t_plot to include forecast days
        % Plot fitted + forecast SIR]
        % plot(t_full, sim_node_plus_forecast_extended, 'Color', color_fit, 'LineStyle', '--', 'LineWidth', 3)
        % Time points for forecast part
        t_forecast_points = t_full(length(t_plot)+1:end);
        forecast_vals = final_model_forecast;

        %% --- REAL DATA (obs_node) SCATTER ---
        scatter(1: length(real_covid_county_data(:,county_index)), real_covid_county_data(:,county_index), marker_size, 'o', ...
            'MarkerEdgeColor', real_data_color, 'LineWidth', 2);
        % 'MarkerFaceColor', c_top, ...
        %% --- FITTED MODEL (sim_node) SCATTER ---
        %scatter(t_fit, sim_node, marker_size, 'o', ...   % squares to distinguish from circles
        %            'MarkerEdgeColor', color_fit, 'LineWidth', 2);

        %% --- FORECAST (final_model_forecast) SCATTER ---
        % scatter(t_forecast_points, forecast_vals, marker_size, 'x', ... % diamond shape
        %              'MarkerEdgeColor', forecast_scatter_color, 'LineWidth', 2);
        scatter(t_forecast_points(:), forecast_vals(:), marker_size, 'x', ...
            'MarkerEdgeColor', forecast_scatter_color, 'LineWidth', 2);


    end
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
    final_real_fcst_window_covid_county_data = real_covid_county_data(t_forecast,county_index);
    % Residuals
    residuals = final_model_forecast_data - final_real_fcst_window_covid_county_data;
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

    %% ================================================================
    %  COLLECT RESULTS FOR THIS COUNTY (5 columns)
    % ================================================================

    real_data_full = real_covid_county_data(:, county_index);

    len_real = length(real_data_full);
    len_comb = length(combined_mean(:,county_index));

    % Pad with NaNs at the end
    % if len_comb < len_real
    %     combined_mean = combined_mean(:,county_index); % [combined_mean; NaN(len_real - len_comb, 1)];
    %     combined_low  = [combined_low;  NaN(len_real - len_comb, 1)];
    %     combined_up   = [combined_up;   NaN(len_real - len_comb, 1)];
    %     t_full        = [t_full(:);     NaN(len_real - length(t_full), 1)];
    % end


    county_block = [ ...
        t_full(:), ...
        real_data_full(t_full(:)), ...
        combined_mean(:,county_index), ...
        combined_low(:,county_index), ...
        combined_up(:,county_index) ...
        ];

    % Append for this county (5 columns added each loop)
    all_county_results = [all_county_results, county_block];

    hold off
end


% =====================================================================
% After the for-loop ends: Build the master CSV table
% =====================================================================

csvgFolder = 'frac-cum y_k data/SIRnetwMf03PFITboothstrapPIPrePeak/csv/';
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

csvFile = fullfile(csvgFolder, 'All_County_SIRnetwork_Results.csv');
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
    final_model_forecast_data                = mean_forecast_across_sim(:,i);
    final_real_fcst_window_covid_county_data = real_covid_county_data(t_forecast,i);
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
origin_population_with_county_name.sir_net_fcst_SSR   = forecasting_SSR_per_county';
origin_population_with_county_name.sir_net_fcst_mse   = forecasting_mse_per_county';
origin_population_with_county_name.sir_net_fcst_mae   = forecasting_mae_per_county';
origin_population_with_county_name.sir_net_fcst_rmse  = forecasting_rmse_per_county';
origin_population_with_county_name.sir_net_fcst_nrmse = forecasting_nrmse_per_county';
origin_population_with_county_name.sir_net_fcst_r2    = forecasting_r2_per_county';
origin_population_with_county_name.sir_net_fcst_mape  = forecasting_mape_per_county';
origin_population_with_county_name.sir_net_fcst_pae   = forecasting_pae_per_county';
origin_population_with_county_name.sir_net_fcst_pa    = forecasting_pa_per_county';
origin_population_with_county_name.sir_net_fcst_pe    = forecasting_pe_per_county';

% Display the updated table
disp(origin_population_with_county_name);

% Save the updated table (CSV)
directory = 'frac-cum y_k data/SIRnetwMf03PFITboothstrapPIPrePeak';
file_path = fullfile(directory, sprintf( ...
    'ten_metric_sir_network_for_constant_beta_01_function_forecasting_errors_with_all_counties_vari_selection_%s_top5_var_index_%s_beta0_%gamma_%g.csv', ...
    vari_selection, top5_var_index, bline_beta0, gamma));

if ~exist(directory, 'dir')
    mkdir(directory);
end

writetable(origin_population_with_county_name, file_path);
