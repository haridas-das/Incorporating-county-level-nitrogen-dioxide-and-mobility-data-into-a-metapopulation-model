clear all; close all; clc;

%% Model fitting with 80 parameters: bline_beta0; omega, gamma, and rho_vec
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
% Define prepandemic cutoff date at the begining
%prepandemic_cutoff = datetime('2020-03-01');
prepandemic_cutoff = datetime(2020,3,28);  % first day of the official lockdown
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

h = figure(1);
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
    'Fig1TimeSeriesDataDynamics/county_data_fiiting_constant_beta_function_all_counties/fig/', ...
    '0_Fig_Bar_Chart_Beta_0_%s_%s_%s.fig'], ...
    vari_selection, top5_var_index, node_str);

png_filename = sprintf(['frac-cum y_k data/', ...
    'Fig1TimeSeriesDataDynamics/county_data_fiiting_constant_beta_function_all_counties/png/', ...
    '0_Fig_Bar_Chart_Beta_0_%s_%s_%s.png'], ...
    vari_selection, top5_var_index, node_str);


saveas(h, fig_filename);
saveas(h, png_filename);
% % Close the figure to free memory
% close(h);


% Plot histogram of Mean_NO2 values
h2 = figure(2);
histogram(NO2_stats_table.Mean_NO2, 10); % You can change the number of bins (e.g., 10)
xlabel('Mean NO_2 Column Density');
ylabel('Number of Counties');
%title('Histogram of Mean NO_2 Column Density across Oklahoma Counties');
% Convert numeric array to strings and join with underscores
node_str = strjoin(string(top_five_pop_flow_nodes(1:5)), '_');

fig_filename = sprintf(['frac-cum y_k data/', ...
    'Fig1TimeSeriesDataDynamics/county_data_fiiting_constant_beta_function_all_counties/fig/', ...
    '0_Fig_histogram_Beta_0_%s_%s_%s.fig'], ...
    vari_selection, top5_var_index, node_str);

png_filename = sprintf(['frac-cum y_k data/', ...
    'Fig1TimeSeriesDataDynamics/county_data_fiiting_constant_beta_function_all_counties/png/', ...
    '0_Fig_histogram_Beta_0_%s_%s_%s.png'], ...
    vari_selection, top5_var_index, node_str);


saveas(h2, fig_filename);
saveas(h2, png_filename);
% Close the figure to free memory
% close(h2);



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

%%% Old Fig legend 
% %% Fig shows all the counties tiem series together
% all_counties_full_ts = covid_data(:,2:end);
% % List of counties to highlight
% highlight_nodes    = [top_five_pop_flow_nodes', 60, 42, 43];
% % Corresponding county names
% highlight_counties = {'Oklahoma', 'Tulsa', 'Cleveland', 'Canadian', 'Comanche', 'Payne', 'Logan', 'Love'};
% 
% % Define highlight colors (make sure they pop in print)
% 
% highlight_colors = [
%     1, 0, 0;                  % 'r'  red
%     0, 0, 1;                  % 'b'  blue
%     0, 1, 0;                  % 'g'  green
%     0, 0, 0;                  % 'k'  black
%     1, 0, 1;                  % 'm'  magenta
%     0.30, 0.70, 1.00;         % color_rest (light blue)
%     1.0, 0.6, 0.6;            % ci_color (soft peach)
%     0.4, 0.0, 0.6;            % color_mean_sim_samp_fit (deep purple)
%     0.90, 0.88, 0.98          % color_sim_samp_fit (light lavender)
%     ];

%% New Fig legend 
%% Fig shows all the counties tiem series together
all_counties_full_ts = covid_data(:,2:end);
% List of counties to highlight
% highlight_nodes    = [13, 38, 24, 37, 55, 72, 14, 47, 25, 60];
% Corresponding county names
% highlight_counties = {'Oklahoma', 'Tulsa', 'Cleveland', 'Canadian', 'Comanche', 'Payne', 'Logan', 'Love'};

% highlight_counties = {'Cimarron', 'Kiowa', 'Garfield', 'Kingfisher', 'Oklahoma', 'Tulsa', 'Cleveland', 'Canadian', 'Comanche', 'Payne'};

highlight_nodes    = [13, 38, 24, 37, 55, 72, 14, 47, 60];
highlight_counties = {'Cimarron', 'Kiowa', 'Garfield', 'Kingfisher', 'Oklahoma', 'Tulsa', 'Cleveland', 'Canadian', 'Payne'};

% Define highlight colors (make sure they pop in print)

% highlight_colors = [
%     1, 0, 0;                  % 'r'  red
%     0, 0, 1;                  % 'b'  blue
%     0, 1, 0;                  % 'g'  green
%     0, 0, 0;                  % 'k'  black
%     1, 0, 1;                  % 'm'  magenta
%     0.30, 0.70, 1.00;         % color_rest (light blue)
%     1.0, 0.6, 0.6;            % ci_color (soft peach)
%     0.4, 0.0, 0.6;            % color_mean_sim_samp_fit (deep purple)
%     0.90, 0.88, 0.98          % color_sim_samp_fit (light lavender)
%     ];

% highlight_colors = [
%     1, 0, 1;                   % Very Low - 'm'  magenta
%     0, 0, 1;                  %  Low -'b'  blue
%     0, 1, 0;                  % Medium -'g'  green
%     0, 0, 0;                  % High -'k'  black
%     1, 0, 0;                  % Very High -'r'  red 
%     0.30, 0.70, 1.00;         % color_rest (light blue)
%     0, 1, 1;            % ci_color (soft peach)
%     0.4, 0.0, 0.6;            % color_mean_sim_samp_fit (deep purple)
%     1, 1, 0                   % color_sim_samp_fit (light lavender)
%     1, 0.5, 0;                % Low - Orange
% 
% ];

highlight_colors = [
    1, 0, 1;                   % Very Low - 'm'  magenta
    0, 0, 1;                  %  Low -'b'  blue
    0, 1, 0;                  % Medium -'g'  green
    0, 0, 0;                  % High -'k'  black
    1, 0, 0;                  % Very High -'r'  red 
    0.30, 0.70, 1.00;         % color_rest (light blue)
    0, 1, 1;            % ci_color (soft peach)
    1, 1, 0                   % color_sim_samp_fit (light lavender)
    1, 0.5, 0;                % Low - Orange

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
box on;
xline(7, 'k--', 'LineWidth', 2);  % first covid cases date in week 2020/03/07
xline(77, 'k--', 'LineWidth', 2); % end of first wave of covid cases date in week 2021/07/03

% xlabel('Time (days)');
% ylabel('COVID cases');
set(gca, 'FontSize', 50, 'FontName', 'Times', 'LineWidth', 2, 'Box', 'on');
ylim([0,16000])
xlim([1, 160])
%legend(arrayfun(@(x) sprintf('County %d', x), highlight_nodes, 'UniformOutput', false), ...
% 'Location', 'best');
% Add legend with custom county names
%legend(highlight_counties, 'Location', 'best', 'FontSize', 14);
% Pass only the handles of highlighted plots to the legend
% legend(hHighlights, highlight_counties, 'Location', 'best', 'FontSize', 14, 'Interpreter', 'none');
% Build legend entries
% legend_handles = [hHighlights; hGray];
% legend_labels = [highlight_counties, {'Other counties'}];
% legend(legend_handles, legend_labels, 'Location', 'best', 'FontSize', 14, 'Interpreter', 'none');
hold off;

% Generate file paths for each figure

fig_filename = sprintf(['frac-cum y_k data/', ...
    'Fig1TimeSeriesDataDynamics/county_data_fiiting_constant_beta_function_all_counties/fig/', ...
    'Fig_all_counties_real_cases_together_time_series_matrix_size_%d.fig'], ...
    m);

png_filename = sprintf(['frac-cum y_k data/', ...
    'Fig1TimeSeriesDataDynamics/county_data_fiiting_constant_beta_function_all_counties/png/', ...
    'Fig_all_counties_real_cases_together_time_series_matrix_size_%d.png'], ...
    m);

% Save the figure in .fig and .png formats
saveas(hft, fig_filename);
saveas(hft, png_filename);





hft = figure('Position', [100, 100, 1800, 800]); % [left, bottom, width, height]
% hft = figure()
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
box on;
xline(7, 'k--', 'LineWidth', 2);  % first covid cases date in week 2020/03/07
xline(77, 'k--', 'LineWidth', 2); % end of first wave of covid cases date in week 2021/07/03

% xlabel('Time (days)');
% ylabel('COVID cases');
set(gca, 'FontSize', 50, 'FontName', 'Times', 'LineWidth', 2, 'Box', 'on');
ylim([0,6000])
xlim([1, 80])
%legend(arrayfun(@(x) sprintf('County %d', x), highlight_nodes, 'UniformOutput', false), ...
% 'Location', 'best');
% Add legend with custom county names
%legend(highlight_counties, 'Location', 'best', 'FontSize', 14);
% Pass only the handles of highlighted plots to the legend
% legend(hHighlights, highlight_counties, 'Location', 'best', 'FontSize', 14, 'Interpreter', 'none');
% Build legend entries
% legend_handles = [hHighlights; hGray];
% legend_labels = [highlight_counties, {'Other counties'}];
% legend(legend_handles, legend_labels, 'Location', 'best', 'FontSize', 14, 'Interpreter', 'none');
% hold off;

% Generate file paths for each figure

fig_filename = sprintf(['frac-cum y_k data/', ...
    'Fig1TimeSeriesDataDynamics/county_data_fiiting_constant_beta_function_all_counties/fig/', ...
    'Fig_Zoom_all_counties_real_cases_together_time_series_matrix_size_%d.fig'], ...
    m);

png_filename = sprintf(['frac-cum y_k data/', ...
    'Fig1TimeSeriesDataDynamics/county_data_fiiting_constant_beta_function_all_counties/png/', ...
    'Fig_Zoom_all_counties_real_cases_together_time_series_matrix_size_%d.png'], ...
    m);

% Save the figure in .fig and .png formats
saveas(hft, fig_filename);
saveas(hft, png_filename);



%% Expanding Training Windows with 5-Day Forecast for County 72

% --- Select county 72 ---
county_idx = 72;
ts_full = all_counties_full_ts(:, county_idx);
marker_size              = 30;
% --- Start from day 7 (epidemic day 0) ---
start_day = 7;
ts = ts_full(start_day:end);
t_rel = 0:length(ts)-1;   % epidemic-relative time

% --- Define training windows and forecast horizon ---
training_windows = [20, 30, 40, 50, 60];
forecast_horizon = 5;  % 5-day forecast beyond training

% --- Define save directories ---
fig_dir = 'frac-cum y_k data/Fig1TimeSeriesDataDynamics/county_data_fiiting_constant_beta_function_all_counties/fig/';
png_dir = 'frac-cum y_k data/Fig1TimeSeriesDataDynamics/county_data_fiiting_constant_beta_function_all_counties/png/';

% --- Create directories if they don't exist ---
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end
if ~exist(png_dir, 'dir')
    mkdir(png_dir);
end

% --- Loop over training windows ---
for w = 1:length(training_windows)
    
    window_end = training_windows(w);
    forecast_end = window_end + forecast_horizon;

    % Ensure forecast_end does not exceed data length
    if forecast_end > length(ts)
        forecast_end = length(ts);
    end

    % --- Create figure ---
    %hft = figure('Position', [100, 100, 1800, 1200]);
    hft = figure()
    hold on;


    % --- Line: training + forecast window ---
    idx_window = t_rel <= forecast_end;
    plot(t_rel(idx_window), ts(idx_window), 'b', 'LineWidth', 2);

        % --- Scatter: full epidemic data ---
    % scatter(t_rel, ts, 130, 'b', 'LineWidth', 3); % 'filled'
    scatter(t_rel, ts, marker_size, ...           % marker size
        'o', ...                              % marker type
        'MarkerEdgeColor', 'b', ...          % marker edge color (blue)
        'LineWidth', 2);                      % line width

    % --- Vertical lines marking training end and forecast end ---
    % xline(window_end, 'b--', 'LineWidth', 3, 'Label', 'Training End');
    % xline(forecast_end, 'g--', 'LineWidth', 3, 'Label', 'Forecast End');
    % --- Vertical lines marking training end and forecast end ---
    xline(window_end, 'k--', 'LineWidth', 1);
    xline(forecast_end, 'g--', 'LineWidth', 1);
    % Set x-axis ticks every 10 days
    xticks(0:10:max(t_rel));
    % --- Formatting ---
    box on;
    set(gca, 'FontSize', 30, 'FontName', 'Times', 'LineWidth', 2);
    ylim([0, 6000])
    xlim([0, 70])
    % xlabel('Days since first reported case')
    % ylabel('COVID cases')
    % title(sprintf('County 72: Training Window = %d days + %d-day Forecast', window_end, forecast_horizon));

    % --- Save figure ---
    fig_filename = sprintf([fig_dir, 'County72_training_%d_forecast_%d.fig'], window_end, forecast_horizon);
    png_filename = sprintf([png_dir, 'County72_training_%d_forecast_%d.png'], window_end, forecast_horizon);

    saveas(hft, fig_filename);
    saveas(hft, png_filename);

    hold off;
end



%% --- Prepare full NO2 time series matrix ---


%% --- Parameters ---
num_counties = length(counties);
vari_selection = 'NO2';                     % select variability based on NO2
% highlight_colors = lines(5);               % color map for top 5 counties

%% --- Compute mean and std NO2 for each county ---
mean_NO2 = zeros(num_counties,1);
std_NO2  = zeros(num_counties,1);

% Find range of dates
first_date = min(oklahoma_data.Date);
last_date  = max(oklahoma_data.Date);

% Create numeric time vector
T = days(last_date - first_date) + 1;  % number of days
all_counties_full_NO2_ts = nan(T, num_counties);

% Create a date vector for plotting
time_vector = (first_date:last_date)';


for i = 1:num_counties
    county_name = counties{i};

    % Filter data for current county
    county_data = oklahoma_data(strcmp(oklahoma_data.CountyName, county_name), :);

    % Pre-pandemic subset for statistics
    county_prepandemic = county_data(county_data.Date <= prepandemic_cutoff, :);
    mean_NO2(i) = mean(county_prepandemic.NO2_column_number_density, 'omitnan');
    std_NO2(i)  = std(county_prepandemic.NO2_column_number_density, 'omitnan');

    % Store full time series in matrix aligned by date
    for j = 1:height(county_data)
        idx = find(time_vector == county_data.Date(j));
        all_counties_full_NO2_ts(idx, i) = county_data.NO2_column_number_density(j);
    end
end

% %% --- Identify top 5 counties by NO2 variability ---
% [~, sorted_idx] = sort(std_NO2, 'descend');
% top5_variability = sorted_idx(1:5);
% 
% if strcmpi(vari_selection, 'NO2')
%     highlight_nodes = top5_variability;
% else
%     highlight_nodes = [];  % keep previous value or empty
% end
% highlight_counties = counties(highlight_nodes);

%% --- Plot time series ---
T = length(time_vector);
hft = figure('Position', [100, 100, 1800, 800]);
hold on;

% Plot all other counties in gray
hGray = [];
for c = 1:num_counties
    if ~ismember(c, highlight_nodes)
        h = plot(1:T, all_counties_full_NO2_ts(:, c), 'Color', [0.85 0.85 0.85], 'LineWidth', 2);
        if isempty(hGray)
            hGray = h; % for legend
        end
    end
end

% Plot top 9 counties with colors and markers
num_highlight = length(highlight_nodes);  % e.g., 5
if num_highlight > size(highlight_colors,1)
    error('You have more highlighted counties than colors defined.');
end

hHighlights = gobjects(num_highlight,1);
for k = 1:num_highlight
    idx = highlight_nodes(k);
    hHighlights(k) = plot(1:T, all_counties_full_NO2_ts(:, idx), ...
        'Color', highlight_colors(k,:), ...
        'LineWidth', 3, ...
        'Marker', 'o', ...
        'MarkerSize', 4, ...
        'MarkerFaceColor', highlight_colors(k,:));
end

% Optional: mark key dates
%xline(find(time_vector == datetime(2021,7,3)), 'k--', 'LineWidth', 2);  % end of first wave
xline(find(time_vector == datetime(prepandemic_cutoff)), 'k--', 'LineWidth', 2);

% Figure formatting
xlabel('Time');
ylabel('NO_2 Column Number Density');
set(gca, 'FontSize', 20, 'FontName', 'Times', 'LineWidth', 1.5, 'Box', 'on');
ylim([0 12]);
xlim([1 T]);

% Legend
legend_handles = [hHighlights; hGray];
legend_labels  = [highlight_counties, {'Other counties'}];
legend(legend_handles, legend_labels, 'Location', 'best', 'FontSize', 14, 'Interpreter', 'none');

hold off;

% %% --- Save figure ---
% fig_filename = sprintf('Figures/NO2_all_counties_timeseries.fig');
% png_filename = sprintf('Figures/NO2_all_counties_timeseries.png');
% saveas(h5, fig_filename);
% saveas(h5, png_filename);

% Generate file paths for each figure

fig_filename = sprintf(['frac-cum y_k data/', ...
    'Fig1TimeSeriesDataDynamics/county_data_fiiting_constant_beta_function_all_counties/fig/', ...
    'Fig_NO2_all_counties_timeserie_matrix_size_%d.fig'], ...
    m);

png_filename = sprintf(['frac-cum y_k data/', ...
    'Fig1TimeSeriesDataDynamics/county_data_fiiting_constant_beta_function_all_counties/png/', ...
    'Fig_NO2_all_counties_timeserie_matrix_size_%d.png'], ...
    m);

% Save the figure in .fig and .png formats
saveas(hft, fig_filename);
saveas(hft, png_filename);

%% --- Plot NO2 time series for Oklahoma County (index 55) ---
%% --- Plot NO2 time series for Oklahoma County (index 55) ---
county_idx = 55;  % Oklahoma County
T = length(time_vector);

% Extract Oklahoma County data
county_data = oklahoma_data(strcmp(oklahoma_data.CountyName, counties{county_idx}), :);

% Pre-pandemic subset
county_prepandemic = county_data(county_data.Date <= prepandemic_cutoff, :);

% Compute pre-pandemic mean NO2
mean_NO2_prepandemic = mean(county_prepandemic.NO2_column_number_density, 'omitnan');

% Prepare full NO2 time series aligned to time_vector
all_counties_full_NO2_ts = nan(length(time_vector), 1);
for j = 1:height(county_data)
    idx = find(time_vector == county_data.Date(j));
    all_counties_full_NO2_ts(idx) = county_data.NO2_column_number_density(j);
end

%% --- Plot ---
hft = figure('Position', [100, 100, 1200, 600]);
hold on;

% Plot Oklahoma County NO2 time series
plot(1:T, all_counties_full_NO2_ts, ...
    'Color', [1 0 0], ...          % red
    'LineWidth', 3, ...
    'Marker', 'o', ...
    'MarkerSize', 4, ...
    'MarkerFaceColor', [1 0 0]);

% Plot horizontal line for pre-pandemic mean
% yline(mean_NO2_prepandemic, 'b--', 'LineWidth', 2, ...
%     'Label', sprintf('Pre-pandemic Avg = %.2f', mean_NO2_prepandemic), ...
%     'LabelHorizontalAlignment', 'left', ...
%     'LabelVerticalAlignment', 'bottom', ...
%     'FontSize', 14);
yline(mean_NO2_prepandemic, 'b--', 'LineWidth', 4);
% Optional: mark first COVID case date
xline(find(time_vector == datetime(prepandemic_cutoff)), 'k--', 'LineWidth', 2);

% Figure formatting
% xlabel('Time');
% ylabel('NO_2 Column Number Density');
% title('NO_2 Time Series - Oklahoma County');
set(gca, 'FontSize', 40, 'FontName', 'Times', 'LineWidth', 1.5, 'Box', 'on');
ylim([0 12]);  % adjust as needed
xlim([1 T]);

% Legend
legend('Oklahoma County', 'Pre-pandemic Avg', 'Location', 'best', 'FontSize', 14);

hold off;

% Generate file paths for each figure

fig_filename = sprintf(['frac-cum y_k data/', ...
    'Fig1TimeSeriesDataDynamics/county_data_fiiting_constant_beta_function_all_counties/fig/', ...
    'Fig_NO2_Oklahoma_county_timeserie_matrix_size_%d.fig'], ...
    m);

png_filename = sprintf(['frac-cum y_k data/', ...
    'Fig1TimeSeriesDataDynamics/county_data_fiiting_constant_beta_function_all_counties/png/', ...
    'Fig_NO2_Oklahoma_county_timeserie_matrix_size_%d.png'], ...
    m);

% Save the figure in .fig and .png formats
saveas(hft, fig_filename);
saveas(hft, png_filename);

% %% data helps me to plot the following staffs
% highlight_nodes    = [13, 38, 24, 37, 55, 72, 14, 47, 60];
% highlight_counties = {'Cimarron', 'Kiowa', 'Garfield', 'Kingfisher', 'Oklahoma', 'Tulsa', 'Cleveland', 'Canadian', 'Payne'};
% 
% highlight_colors = [
%  1, 0, 1;                   % Very Low - 'm'  magenta
%  0, 0, 1;                  %  Low -'b'  blue
%  0, 1, 0;                  % Medium -'g'  green
%  0, 0, 0;                  % High -'k'  black
%  1, 0, 0;                  % Very High -'r'  red 
%  0.30, 0.70, 1.00;         % color_rest (light blue)
%  0, 1, 1;            % ci_color (soft peach)
%  1, 1, 0                   % color_sim_samp_fit (light lavender)
%  1, 0.5, 0;                % Low - Orange
% ];


%% --- Plot NO2 time series for Oklahoma County (55) and Payne County (60) ---
county_idx_list = [55, 13];  % 55 = Oklahoma, 60 = Payne 13 = 'Cimarron'
T = length(time_vector);

colors = [1 0 0;    % red for Oklahoma
          1 0 1];   % magenta for Cimarron

hft = figure('Position', [100, 100, 1200, 600]);
hold on;

for c = 1:length(county_idx_list)
    county_idx = county_idx_list(c);
    
    % Extract county data
    county_data = oklahoma_data(strcmp(oklahoma_data.CountyName, counties{county_idx}), :);
    
    % Pre-pandemic subset
    county_prepandemic = county_data(county_data.Date <= prepandemic_cutoff, :);
    
    % Compute pre-pandemic mean NO2
    mean_NO2_prepandemic = mean(county_prepandemic.NO2_column_number_density, 'omitnan');
    
    % Prepare full NO2 time series aligned to time_vector
   % all_counties_full_NO2_ts = nan(length(time_vector), 1);
    for j = 1:height(county_data)
        idx = find(time_vector == county_data.Date(j));
        all_counties_full_NO2_ts(idx) = county_data.NO2_column_number_density(j);
    end
    
    % Plot county NO2 time series
    plot(1:T, all_counties_full_NO2_ts, ...
        'Color', colors(c,:), ...
        'LineWidth', 2, ...
        'Marker', 'o', ...
        'MarkerSize', 2, ...
        'MarkerFaceColor', colors(c,:));
    
    % Plot horizontal line for pre-pandemic mean
    yline(mean_NO2_prepandemic, '--', 'LineWidth', 3, 'Color', colors(c,:));
end

% Optional: mark pre-pandemic cutoff
xline(find(time_vector == prepandemic_cutoff), 'k--', 'LineWidth', 2);

% % Figure formatting
% xlabel('Time');
% ylabel('NO_2 Column Number Density');
% title('NO_2 Time Series - Oklahoma and Payne Counties');
set(gca, 'FontSize', 40, 'FontName', 'Times', 'LineWidth', 1.5, 'Box', 'on');
ylim([0 12]);  % adjust as needed
%xlim([1 T]);
xlim([1, 80])

% % Legend
% legend({'Oklahoma County','Oklahoma Pre-pandemic Avg', ...
%         'Payne County','Payne Pre-pandemic Avg'}, ...
%         'Location','best', 'FontSize', 14);
% 
% hold off;

%% --- Save the figure ---
fig_filename = sprintf('frac-cum y_k data/Fig1TimeSeriesDataDynamics/county_data_fiiting_constant_beta_function_all_counties/fig/Fig_NO2_OK_Payne_counties_timeserie_matrix_size_%d.fig', m);
png_filename = sprintf('frac-cum y_k data/Fig1TimeSeriesDataDynamics/county_data_fiiting_constant_beta_function_all_counties/png/Fig_NO2_OK_Payne_counties_timeserie_matrix_size_%d.png', m);

saveas(hft, fig_filename);
saveas(hft, png_filename);
