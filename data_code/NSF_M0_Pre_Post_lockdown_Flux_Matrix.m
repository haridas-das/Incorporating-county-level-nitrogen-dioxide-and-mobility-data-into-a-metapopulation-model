% Driver File to run the Mobility SIR Network model Dynamics
clear all; clc; close all;
% Data for the time dependent flux Matrix
[state, m, name_origin_list, time_start, time_end, top_five_pop_flow_nodes, time_step_size , zeta_vec, zeta,  Nv, gamma, beta_1_vec, phi_0] = county_name_list()
Nv     = 1;
% Number of nodes or size of the flux matrix
N_o    = m;                                                                  % number of origin
N_d    = m;                                                                  % number of destination
ylimit = 0.6;
% top_five_pop_flow_nodes = [55, 72, 14, 9, 16];
zeta          = 1;
gamma         = 0.1429;  %1/7                wekly recovered
beta_1_vec    = linspace(0,0.1428,1); % simulate the infection rates% Recovery rate
epsilon       = 0.2;
bline_beta0   = 0.1428; % 0.5;
omega         = 0.01;

str = sprintf('TD_Realistic_Human_Mobility_NO2_data_driven_beta_epsilon_%.2f', epsilon);

%% select model based on NO_2 and high populous
vari_selection = 'NO2';

%% Sort by variability (Std_NO2) in descending order
top5_var_index = 'Mean_NO2';
% Load county level NO_2 data
data = readtable('NO2_county_data_2020.csv');

%% Origin population from the data
origin_population                  = csvread('pop_vector_origin_only.csv',1,0);
origin_population_with_county_name = readtable('pop_vector_origin.csv', 'Delimiter', ',');
vec_N = origin_population(:,2);            % population of each nodes

% Define prepandemic cutoff date
prepandemic_cutoff = datetime('2020-03-01');
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

%% Plot bar chart for top 5 NO_2 counties with three possible variability
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
% Generate file paths for each figure
% Convert numeric array to strings and join with underscores
node_str = strjoin(string(top_five_pop_flow_nodes(1:5)), '_');
% Example: save figure with indices of top 5 counties in the filename
% Build PNG filename safely
fig_filename = sprintf('%s/result figures/NO2SimulationStaticPreLockdownMobility/fig/0_Fig_Bar_Chart_Beta_0_%s_%s_%s.png', ...
    pwd, char(vari_selection),char(top5_var_index), node_str);;% Example: save figure with indices of top 5 counties in the filename
% Build PNG filename safely
png_filename = sprintf('%s/result figures/NO2SimulationStaticPreLockdownMobility/png/0_Fig_Bar_Chart_Beta_0_%s_%s_%s.png', ...
    pwd, char(vari_selection),char(top5_var_index), node_str);

% Save the figure in .fig and .png formats
saveas(h, fig_filename);
saveas(h, png_filename);
% Close the figure to free memory
close(h);

%% Plot histogram of Mean_NO2 values
h2 = figure(2);
histogram(NO2_stats_table.Mean_NO2, 10); % You can change the number of bins (e.g., 10)
xlabel('Mean NO_2 Column Density');
ylabel('Number of Counties');
%title('Histogram of Mean NO_2 Column Density across Oklahoma Counties');
% Convert numeric array to strings and join with underscores
node_str = strjoin(string(top_five_pop_flow_nodes(1:5)), '_');
% Example: save figure with indices of top 5 counties in the filename
% Build PNG filename safely
fig_filename = sprintf('%s/result figures/NO2SimulationStaticPreLockdownMobility/fig/0_Fig_histogram_Beta_0_%s_%s_%s.png', ...
    pwd, char(vari_selection),char(top5_var_index), node_str);;% Example: save figure with indices of top 5 counties in the filename
% Build PNG filename safely
png_filename = sprintf('%s/result figures/NO2SimulationStaticPreLockdownMobility/png/0_Fig_histogram_Beta_0_%s_%s_%s.png', ...
    pwd, char(vari_selection),char(top5_var_index), node_str);
% Save the figure in .fig and .png formats
saveas(h2, fig_filename);
saveas(h2, png_filename);
% Close the figure to free memory
close(h2);


%% Cellphone mobility data driven SIRnetwork model input
%f = 6; % Put f=0 for Gen Flux, f=1 for star shape; f=2 for fully connected,
M = m; % Number of nodes
N_o = m;
N_d = m;
% nodes_vec = 1:m;
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
% highlight_nodes    = [top_five_pop_flow_nodes, 60, 42, 43];
highlight_nodes    = [55, 72, 14, 9, 16, 60, 42, 43];
% Corresponding county names
highlight_counties = {'Oklahoma', 'Tulsa', 'Cleveland', 'Canadian', 'Comanche', 'Payne', 'Logan', 'Love'};

% Define highlight colors (make sure they pop in print)
highlight_colors = [
    1, 0, 0;                  % red
    0, 0, 1;                  % blue
    0, 0.6, 0;                % green (brighter than pure [0,1,0])
    0, 0, 0;                  % black
    1, 0, 1;                  % magenta
    20/255, 90/255, 50/255;   % dark green
    237/255, 178/255, 32/255; % golden yellow
    126/255, 47/255, 142/255  % purple
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
fig_filename = [pwd '/result figures/NO2SimulationStaticPreLockdownMobility/fig/Fig_all_counties_real_cases_together_time_series_matrix_size_' num2str(m) '.fig'];
png_filename = [pwd '/result figures/NO2SimulationStaticPreLockdownMobility/png/Fig_all_counties_real_cases_together_time_series_matrix_size_' num2str(m) '.png'];

% Save the figure in .fig and .png formats
saveas(hft, fig_filename);
saveas(hft, png_filename);
% Close the figure to free memory
close(hft);

%covid_data = readmatrix('final_weekly_time_series_covid_data.csv');
% covid_data = covid_data(7:64,:); covid_data(7:66,:);
covid_data = covid_data(7:71,:);
% Define dengue_div_data from covid_data (columns 4 to 11)
covid_county_data = covid_data(:,2:end);  % Real data (time along rows, nodes along columns)
time_end_plot           = size(covid_data, 1);
% Convert lockdown_date to 0-indexed for MATLAB slicing (adjusting for 1-based indexing)
lockdown_date           = time_end;
lockdown_date_0_indexed = lockdown_date;
lockdown_week_0_indexed = 4;  
u = 1;
v = 1;
time_end  = lockdown_date; %time_end_plot;
T         = time_end;
time_step_size = 0.1; % can compute N_time_steps
tspan         = 1:time_step_size:time_end;  % Ensure this is within the bounds of `t`
t             = tspan; %linspace(1, time_end_plot, time_end_plot); % Time vector (0 to 20)
phi_0         = 1;  % 1 means gravity model flux values
NO2_data      = NO2_stats_table.Mean_NO2;
%% Data for the Initial value of the Network Model
u = top_five_pop_flow_nodes(1);                                                                     % heteroginous node
v = top_five_pop_flow_nodes(1);                                                                     % Node where the infectious started , put Oklahoma County  55
T = time_end;
plot_time_end = T;
% time span for the simulation in week
% Range and Infected nodes defined for the simulations
tspan = [1 : time_step_size : T];
%% Origin population from the data
origin_population   = csvread('pop_vector_origin_only.csv',1,0);             % this has no origin name
% origin_population = csvread('pop_vector_origin.csv',1,0)                 % this includes name of the origin nad rank
vec_N               = origin_population(:,2);                                            % population of each nodes
%% Initial value of the Problem by the following subroutine function
[vec_S0,vec_I0,vec_N,y0] = InitialValueMobilitySIRNetwork(m,vec_N,u,Nv);

time_step_size_pre_post_covid_flux = 1; % can compute N_time_steps

% tspan = 1:time_step_size_pre_post_covid_flux:366;
time_end_mobility_flux   =  366;

orig_matrix_phi          = TDFluxMatrix(N_o,N_d,time_step_size_pre_post_covid_flux,time_start,time_end_mobility_flux,name_origin_list);

fprintf('The mean flux value is %10.3f \n',mean(mean(mean(orig_matrix_phi))));

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

% avg_post_lockdown_flux_matrix = avg_pre_lockdown_flux_matrix;

%% checking that the row sum is 1
for i =1: m
    sum(avg_post_lockdown_flux_matrix(i,:));
end

for i =1: m
    sum(avg_pre_lockdown_flux_matrix(i,:));
end


% avg_post_lockdown_flux_matrix
% 
% 
% avg_pre_lockdown_flux_matrix


% ==========================================
% Folder to save flux matrices
% ==========================================
save_folder = fullfile(pwd, 'NSF_M0_Numerical_Simulation');

if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end

% ==========================================
% Safety checks
% ==========================================
assert(exist('avg_pre_lockdown_flux_matrix','var') == 1, ...
       'avg_pre_lockdown_flux_matrix does not exist');

assert(exist('avg_post_lockdown_flux_matrix','var') == 1, ...
       'avg_post_lockdown_flux_matrix does not exist');

% ==========================================
% File names (clear + publication-ready)
% ==========================================
filename_pre = fullfile(save_folder, ...
    sprintf('avg_pre_lockdown_flux_matrix_network_f%d.csv', m));

filename_post = fullfile(save_folder, ...
    sprintf('avg_post_lockdown_flux_matrix_network_f%d.csv', m));

% ==========================================
% Save as CSV
% ==========================================
writematrix(avg_pre_lockdown_flux_matrix, filename_pre);
writematrix(avg_post_lockdown_flux_matrix, filename_post);

% ==========================================
% Confirm save
% ==========================================
fprintf('Saved:\n%s\n%s\n', filename_pre, filename_post);
