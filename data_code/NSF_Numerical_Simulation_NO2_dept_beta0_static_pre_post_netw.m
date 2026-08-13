% Driver File to run the Mobility SIR Network model Dynamics
clear all; clc; close all;
% Data for the time dependent flux Matrix
[state, m, name_origin_list, time_start, time_end, top_five_pop_flow_nodes, time_step_size , zeta_vec, zeta,  Nv, gamma, beta_1_vec, phi_0] = county_name_list()
%[state, m, name_origin_list, time_start, lockdown_date, top_five_pop_flow_nodes, time_step_size , zeta_vec, zeta,  Nv, gamma, beta_1_vec, phi_0] = county_name_list();
Nv     = 1; 
% Number of nodes or size of the flux matrix
N_o    = m;                                                                  % number of origin
N_d    = m;                                                                  % number of destination
ylimit = 0.6;
% top_five_pop_flow_nodes = [55, 72, 14, 9, 16];
zeta       = 1;
gamma      = 7/5; %0.1429;  %1/7                wekly recovered
% beta_1_vec =linspace(0,0.1428,1); % simulate the infection rates% Recovery rate

% %% Case 1
% beta_1_vec = 0.5; %linspace(0,0.1428,1); % simulate the infection rates% Recovery rate
% epsilon    = 0.5;
% 
% %% Case 2
% beta_1_vec = 1; %linspace(0,0.1428,1); % simulate the infection rates% Recovery rate
% epsilon    = 0.5;

% %% Case 3
% beta_1_vec = 0.5; %linspace(0,0.1428,1); % simulate the infection rates% Recovery rate
% epsilon    = 1;

%% Case 4
beta_1_vec = 1.0; %linspace(0,0.1428,1); % simulate the infection rates% Recovery rate
epsilon = 1;
beta0 = beta_1_vec;
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
fig_filename = sprintf('%s/result figures/NO2SimulationPeakVsTimeVariesOverEpsilon/fig/0_Fig_Bar_Chart_Beta_0_%s_%s_%s.png', ...
    pwd, char(vari_selection),char(top5_var_index), node_str);;% Example: save figure with indices of top 5 counties in the filename
% Build PNG filename safely
png_filename = sprintf('%s/result figures/NO2SimulationPeakVsTimeVariesOverEpsilon/png/0_Fig_Bar_Chart_Beta_0_%s_%s_%s.png', ...
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
fig_filename = sprintf('%s/result figures/NO2SimulationPeakVsTimeVariesOverEpsilon/fig/0_Fig_histogram_Beta_0_%s_%s_%s.png', ...
    pwd, char(vari_selection),char(top5_var_index), node_str);;% Example: save figure with indices of top 5 counties in the filename
% Build PNG filename safely
png_filename = sprintf('%s/result figures/NO2SimulationPeakVsTimeVariesOverEpsilon/png/0_Fig_histogram_Beta_0_%s_%s_%s.png', ...
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

%T = size(all_counties_full_ts, 1);          % number of time points

%covid_data = readmatrix('final_weekly_time_series_covid_data.csv');
% covid_data = covid_data(7:64,:); covid_data(7:66,:);
covid_data = covid_data(7:71,:);
% Define dengue_div_data from covid_data (columns 4 to 11)
covid_county_data = covid_data(:,2:end);  % Real data (time along rows, nodes along columns)
time_end_plot           = size(covid_data, 1);
% Convert lockdown_date to 0-indexed for MATLAB slicing (adjusting for 1-based indexing)
lockdown_date           = time_end;
lockdown_date_0_indexed = lockdown_date;
lockdown_week_0_indexed = 4;   % this is in week
% top_five_pop_flow_nodes = [1 2 3 4 5 6 7 8];
u = 72;
v = 72;
% time_start = 1;
time_end  = lockdown_date; %time_end_plot;
T         = time_end;
%plot_time_end = T;
time_step_size = 0.1; % can compute N_time_steps
% ylimit         = 0.6;
% bline_beta0   = 1;
% omega_0       = 1;
% omega_1       = 1;
% rho           = 10^-1; % 1
% optim_rho0    = 10^-1; 
% optim_epsilon0= 10^-1;
% optim_epsilon1= 10^-1;

tspan         = 1:time_step_size:time_end;  % Ensure this is within the bounds of `t`
t             = tspan; %linspace(1, time_end_plot, time_end_plot); % Time vector (0 to 20)
%f_vec   = 1.3*ones(1,m); % beta_vec values
%rho_vec = rho*ones(1,m);
%top_five_pop_flow_nodes = [37,45,49,51,2];
%% Range and Infected nodes defined for the simulations
%simulation by the mean mobility flux value 0.00651 from 2020 mobility data
%phi_0 = 0.00651;
phi_0      = 1;  % 1 means gravity model flux values

%% Extract and save NO2 data
NO2_data   = NO2_stats_table.Mean_NO2;
min(NO2_data);
max(NO2_data);
% Define folder
save_folder = "result figures/NO2SimulationPeakVsTimeVariesOverEpsilon"

NO2_table = table(NO2_data, 'VariableNames', {'Mean_NO2'});

filename = fullfile(save_folder, 'NO2_data.csv');
writetable(NO2_table, filename);

%str = 'Fully_Connected';
% List of the origin name
% name_origin_list = {'Adair County', 'Alfalfa County', 'Atoka County', 'Beaver County', 'Beckham County', 'Blaine County', 'Bryan County', 'Caddo County', 'Canadian County', 'Carter County', 'Cherokee County', 'Choctaw County', 'Cimarron County', 'Cleveland County', 'Coal County', 'Comanche County', 'Cotton County', 'Craig County', 'Creek County', 'Custer County', 'Delaware County', 'Dewey County', 'Ellis County', 'Garfield County', 'Garvin County', 'Grady County', 'Grant County', 'Greer County', 'Harmon County', 'Harper County', 'Haskell County', 'Hughes County', 'Jackson County', 'Jefferson County', 'Johnston County', 'Kay County', 'Kingfisher County', 'Kiowa County', 'Latimer County', 'Le Flore County', 'Lincoln County', 'Logan County', 'Love County', 'Major County', 'Marshall County', 'Mayes County', 'McClain County', 'McCurtain County', 'McIntosh County', 'Murray County', 'Muskogee County', 'Noble County', 'Nowata County', 'Okfuskee County', 'Oklahoma County', 'Okmulgee County', 'Osage County', 'Ottawa County', 'Pawnee County', 'Payne County', 'Pittsburg County', 'Pontotoc County', 'Pottawatomie County', 'Pushmataha County', 'Roger Mills County', 'Rogers County', 'Seminole County', 'Sequoyah County', 'Stephens County', 'Texas County', 'Tillman County', 'Tulsa County', 'Wagoner County', 'Washington County', 'Washita County', 'Woods County', 'Woodward County'};
%name_origin_list = county_name_list(m);
% [m, name_origin_list, state, time_end, top_five_pop_flow_nodes] = county_name_list(m)
%% Data for the Initial value of the Network Model
u = top_five_pop_flow_nodes(2);                                                                     % heteroginous node
v = top_five_pop_flow_nodes(2);                                                                     % Node where the infectious started , put Oklahoma County  55
T = time_end;
plot_time_end = T;
% time span for the simulation in week
% Range and Infected nodes defined for the simulations
tspan = [1 : time_step_size : T];
%% Origin population from the data
origin_population = csvread('pop_vector_origin_only.csv',1,0)              % this has no origin name
% origin_population = csvread('pop_vector_origin.csv',1,0)                 % this includes name of the origin nad rank
vec_N = origin_population(:,2);                                            % population of each nodes
%% Constant population for the simulation
%% Initial value of the Problem by the following subroutine function
[vec_S0,vec_I0,vec_N,y0] = InitialValueMobilitySIRNetwork(m,vec_N,u,Nv);
%vec_I0(72) = 3;

% 
% %% this is time dependent flux from the real mobility data
% matrix_phi = TDFluxMatrix(N_o,N_d,time_step_size,time_start,time_end,name_origin_list)
% 
% fprintf('The mean flux value is %10.3f \n',mean(mean(mean(matrix_phi))))
% 
% %Np = TD_PresentPopFlow(m,matrix_phi,vec_N,tspan);
% 
% Np = zeros(m,length(tspan)); % initial present population at i
% Npsum = 0;
% for t = 1 : length(tspan)
%     for j = 1: m
%         for k = 1: m
%             %total population at node j
%             Npsum = Npsum+matrix_phi(k,j,t)*vec_N(k); % phi(k,j)*vec_N(k,1)is the total population at vertex k that travels every day to vertex j
%         end
%         Np(j,t) = Npsum;
%         %break;
%         Npsum = 0;
%     end
% end

%% Data for heterogeneous transmission node

size(NO2_stats_table.Mean_NO2);
vec_beta = beta0*ones(m,1)+ epsilon*NO2_stats_table.Mean_NO2;
vec_beta(u,1) = zeta*beta0; % heterogeneous node u which will be different infection rate
bline_beta0 = beta0; 
omega = epsilon;
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

%% Data for heterogeneous transmission node
% static_matrix_phi = avg_pre_lockdown_flux_matrix;
% beta0             = beta_1_vec;
% vec_beta          = beta0*ones(m,1)+ epsilon*NO2_stats_table.Mean_NO2;
% [t,y_k] = ode45(@(t,y) SirGenNetworkModel(m,static_matrix_phi,y,vec_beta,vec_S0,vec_I0,vec_N,gamma,zeta),tspan,y0);



%% Precompute the flux matrix function
% pre and  post lockdown avg mobility data
phi_matrix_func = @(t) get_matrix_phi(t, lockdown_week_0_indexed, time_start, time_step_size, avg_pre_lockdown_flux_matrix, avg_post_lockdown_flux_matrix);
[t, y_k] = ode45(@(t, y) SIR_Network_NO2_baseB0(t, m, phi_matrix_func(t), y, NO2_data, bline_beta0, omega, vec_S0, vec_I0, vec_N, gamma, zeta), tspan, y0);


%% PART 2: Solve the time-dependent Ordinary differential equation
%[t, y_k] = ode45(@(t,y) TDSirGenNetworkModel(t,y,m,Np, matrix_phi,vec_beta, vec_S0, vec_I0, vec_N, gamma, zeta,time_start,time_step_size), tspan, y0)

%% Some particulr nodes infectious dynamics
%figure('DefaultAxesFontSize',20) % this will change the Font Size on the both Axis\
h2 = figure(2)
plot(t,y_k(:,m+top_five_pop_flow_nodes(1))/vec_N(top_five_pop_flow_nodes(1)),'r-','MarkerFaceColor','r','LineWidth',3);
hold on
set(gca,'LineWidth',2);
set(gca,'FontWeight','Bold');
hold on
plot(t,y_k(:,m+top_five_pop_flow_nodes(2))/vec_N(top_five_pop_flow_nodes(2)),'g-','MarkerFaceColor','g','LineWidth',3);
plot(t,y_k(:,m+top_five_pop_flow_nodes(3))/vec_N(top_five_pop_flow_nodes(3)),'b-','MarkerFaceColor','b','LineWidth',3);
plot(t,y_k(:,m+top_five_pop_flow_nodes(4))/vec_N(top_five_pop_flow_nodes(4)),'c-','MarkerFaceColor','c','LineWidth',3);
plot(t,y_k(:,m+top_five_pop_flow_nodes(5))/vec_N(top_five_pop_flow_nodes(5)),'m-','MarkerFaceColor','m','LineWidth',3);
plot(t,sum(y_k(:,m+1:2*m),2)/sum(vec_N),'r--','MarkerFaceColor','b','LineWidth',3);
%plot(t,y_k(:,m+top_five_pop_flow_nodes(6))/vec_N(m),'k-','MarkerFaceColor','k','LineWidth',3);
%grid on;
set(gca, 'fontsize', 20, 'fontweight','bold'); % gca: get current axis
% set(gca, 'YTick', 0:.1:ylimit); % Full form of gca is (Get Current Axes)
set(gca, 'XTick', 0:8:time_end); % Full form of gca is (Get Current Axes)
% ylim([0,ylimit]);
xlim([1,plot_time_end])
set(gca, 'fontsize', 20, 'fontweight','bold','LineWidth',2); % gca: get current axis
% Legend for individual nodes
legend(['$I_{' num2str(top_five_pop_flow_nodes(1)) '}$'], ...
    ['$I_{' num2str(top_five_pop_flow_nodes(2)) '}$'], ...
    ['$I_{' num2str(top_five_pop_flow_nodes(3)) '}$'], ...
    ['$I_{' num2str(top_five_pop_flow_nodes(4)) '}$'], ...
    ['$$I_{' num2str(top_five_pop_flow_nodes(5)) '}$$'], ...
    '$$\sum_{i=1}^{m} I_i$$', 'Interpreter', 'latex', 'Location', 'Best','FontSize', 20, 'FontWeight', 'bold');
%       '$$\sum_{i=1}^{m} I_i/\sum_{i=1}^{m} N_i$$', 'Interpreter', 'latex', 'Location', 'Best','FontSize', 20, 'FontWeight', 'bold');

h2_figurename_fig = ['/result figures/NO2SimulationPeakVsTimeVariesOverEpsilon/fig/Fig_2_B_cons_beta_TD_Flux_ind_Dynmics_v_value_' num2str(v) '_matrix_size_' num2str(m) '_Network_' str '_Size_' num2str(T) '_infected_strat_node_' num2str(v) '_Het_Node_' num2str(u) '_RR_' num2str(gamma) '.fig'];
h2_figurename_png = ['/result figures/NO2SimulationPeakVsTimeVariesOverEpsilon/png/Fig_2_B_cons_beta_TD_Flux_ind_Dynmics_v_value_' num2str(v) '_matrix_size' num2str(m) '_Network_' str '__Size_' num2str(T) '_infected_start_node_' num2str(v) '_Het_Node_' num2str(u) '_RR_' num2str(gamma) '.png'];
saveas(h2,[pwd h2_figurename_fig])
saveas(h2,[pwd h2_figurename_png])

filename=['Indidual_Dynamics_1_77_TD_flux_het_u' num2str(u) '_Origin_Data_Pop_zeta_' num2str(zeta) '.csv'];
cum_y_k_fraction = sum(y_k(:,m+1:2*m),2)/sum(vec_N);
csvwrite(filename,cum_y_k_fraction)

top_five_pop_flow_nodes = [55 72 14 9 16];

%% Some particulr nodes infectious dynamics
%figure('DefaultAxesFontSize',20) % this will change the Font Size on the both Axis\
h3 = figure(3)
plot(t,y_k(:,m+top_five_pop_flow_nodes(1))/vec_N(top_five_pop_flow_nodes(1)),'r-','MarkerFaceColor','r','LineWidth',3);
hold on
set(gca,'LineWidth',2);
set(gca,'FontWeight','Bold');
hold on
plot(t,y_k(:,m+top_five_pop_flow_nodes(2))/vec_N(top_five_pop_flow_nodes(2)),'b-','MarkerFaceColor','b','LineWidth',3);
plot(t,y_k(:,m+top_five_pop_flow_nodes(3))/vec_N(top_five_pop_flow_nodes(3)),'g-','MarkerFaceColor','g','LineWidth',3);
plot(t,y_k(:,m+top_five_pop_flow_nodes(4))/vec_N(top_five_pop_flow_nodes(4)),'k-','MarkerFaceColor','k','LineWidth',3);
plot(t,y_k(:,m+top_five_pop_flow_nodes(5))/vec_N(top_five_pop_flow_nodes(5)),'m-','MarkerFaceColor','m','LineWidth',3);

plot(t,y_k(:,m+top_five_pop_flow_nodes(2))/vec_N(top_five_pop_flow_nodes(2)),'b-','MarkerFaceColor','b','LineWidth',3);
plot(t,y_k(:,m+top_five_pop_flow_nodes(3))/vec_N(top_five_pop_flow_nodes(3)),'g-','MarkerFaceColor','g','LineWidth',3);
plot(t,y_k(:,m+top_five_pop_flow_nodes(4))/vec_N(top_five_pop_flow_nodes(4)),'k-','MarkerFaceColor','k','LineWidth',3);
plot(t,y_k(:,m+top_five_pop_flow_nodes(5))/vec_N(top_five_pop_flow_nodes(5)),'m-','MarkerFaceColor','m','LineWidth',3);

plot(t,sum(y_k(:,m+1:2*m),2)/sum(vec_N),'r--','MarkerFaceColor','r','LineWidth',3);
%plot(t,y_k(:,m+top_five_pop_flow_nodes(6))/vec_N(m),'k-','MarkerFaceColor','k','LineWidth',3);
%grid on;
set(gca, 'fontsize', 20, 'fontweight','bold'); % gca: get current axis
% set(gca, 'YTick', 0:.1:ylimit); % Full form of gca is (Get Current Axes)
set(gca, 'XTick', 0:88: time_end); % Full form of gca is (Get Current Axes)
% ylim([0,ylimit]);
xlim([1,plot_time_end])
set(gca, 'fontsize', 25, 'fontweight','bold','LineWidth',2); % gca: get current axis
% legend only sum of total cases in each day node sum
legend(['$I_{' num2str(top_five_pop_flow_nodes(1)) '}$'], ...
    ['$I_{' num2str(top_five_pop_flow_nodes(2)) '}$'], ...
    ['$I_{' num2str(top_five_pop_flow_nodes(3)) '}$'], ...
    ['$I_{' num2str(top_five_pop_flow_nodes(4)) '}$'], ...
    ['$$I_{' num2str(top_five_pop_flow_nodes(5)) '}$$'], ...
    '$$\sum_{i=1}^{m} I_i$$', 'Interpreter', 'latex', ...
    'Location', 'Best', 'NumColumns', 2);
%        '$$\sum_{i=1}^{m} I_i/\sum_{i=1}^{m} N_i$$', 'Interpreter', 'latex', ...
%        'Location', 'Best', 'NumColumns', 2);
h3_figurename_fig = ['/result figures/NO2SimulationPeakVsTimeVariesOverEpsilon/fig/Fig_2_B_cons_beta_TD_Flux_ind_and_cum_Dynmics_v_value_' num2str(v) '_matrix_size_' num2str(m) '_Network_' str '_Size_' num2str(T) '_infected_strat_node_' num2str(v) '_Het_Node_' num2str(u) '_RR_' num2str(gamma) '.fig'];
h3_figurename_png = ['/result figures/NO2SimulationPeakVsTimeVariesOverEpsilon/png/Fig_2_B_cons_beta_TD_Flux_ind_and_cum_Dynmics_v_value_' num2str(v) '_matrix_size' num2str(m) '_Network_' str '__Size_' num2str(T) '_infected_start_node_' num2str(v) '_Het_Node_' num2str(u) '_RR_' num2str(gamma) '.png'];
saveas(h3,[pwd h3_figurename_fig])
saveas(h3,[pwd h3_figurename_png])
hold off



%% Some particulr nodes infectious dynamics
%figure('DefaultAxesFontSize',20) % this will change the Font Size on the both Axis\
colors = {'#145A32', 'r', 'b', 'g', 'k', 'm', '#EDB120'}; 
top5_idx = [42; 55; 72; 14 ; 9; 16; 41] 

h4 = figure(4)

plot(t,y_k(:,m+top_five_pop_flow_nodes(1))/vec_N(top_five_pop_flow_nodes(1)),'r-','MarkerFaceColor','r','LineWidth',3);
hold on
set(gca,'LineWidth',2);
% set(gca,'FontWeight','Bold');
plot(t,y_k(:,m+top_five_pop_flow_nodes(2))/vec_N(top_five_pop_flow_nodes(2)),'b-','MarkerFaceColor','b','LineWidth',3);
plot(t,y_k(:,m+top_five_pop_flow_nodes(3))/vec_N(top_five_pop_flow_nodes(3)),'g-','MarkerFaceColor','g','LineWidth',3);
plot(t,y_k(:,m+top_five_pop_flow_nodes(4))/vec_N(top_five_pop_flow_nodes(4)),'k-','MarkerFaceColor','k','LineWidth',3);
plot(t,y_k(:,m+top_five_pop_flow_nodes(5))/vec_N(top_five_pop_flow_nodes(5)),'m-','MarkerFaceColor','m','LineWidth',3);
plot(t,sum(y_k(:,m+1:2*m),2)/sum(vec_N),'r--','MarkerFaceColor','r','LineWidth',3);

% plot(t, y_k(:, m + top5_idx(1)) / vec_N(top5_idx(1)), '-', 'Color', '#145A32', 'MarkerFaceColor', '#145A32', 'LineWidth', 3);
% plot(t, y_k(:, m + top5_idx(7)) / vec_N(top5_idx(7)), '-', 'Color', '#EDB120', 'MarkerFaceColor', '#EDB120', 'LineWidth', 3);

% plot(t,y_k(:,m+1)/vec_N(1),'r-.','MarkerFaceColor','r','LineWidth',3);
% plot(t,y_k(:,m+2)/vec_N(2),'b-.','MarkerFaceColor','b','LineWidth',3);
% plot(t,y_k(:,m+3)/vec_N(3),'g-.','MarkerFaceColor','g','LineWidth',3);
% plot(t,y_k(:,m+4)/vec_N(4),'k-.','MarkerFaceColor','k','LineWidth',3);
% plot(t,y_k(:,m+5)/vec_N(5),'m-.','MarkerFaceColor','m','LineWidth',3);

%grid on;
% set(gca, 'fontsize', 20); % gca: get current axis
%set(gca, 'fontsize', 20, 'fontweight','bold'); % gca: get current axis
%set(gca, 'YTick', 0:.1:ylimit); % Full form of gca is (Get Current Axes)
% set(gca, 'XTick', 1:11: time_end); % Full form of gca is (Get Current Axes)
% set(gca, 'XTick', [1:11:time_end, time_end])
set(gca, 'XTick', [0:10:70])
% ylim([0,ylimit]);
xlim([1,70])
set(gca, 'fontsize', 30,'LineWidth',2,'FontName', 'Times'); % gca: get current axis
%set(gca, 'fontsize', 22, 'fontweight','bold','LineWidth',2); % gca: get current axis
% legend only sum of total cases in each day node sum
% legend(['$I_{' num2str(top_five_pop_flow_nodes(1)) '}$'], ...
%        ['$I_{' num2str(top_five_pop_flow_nodes(2)) '}$'], ...
%        ['$I_{' num2str(top_five_pop_flow_nodes(3)) '}$'], ...
%        ['$I_{' num2str(top_five_pop_flow_nodes(4)) '}$'], ...
%        ['$$I_{' num2str(top_five_pop_flow_nodes(5)) '}$$'], ...
%        '$$\sum_{i=1}^{m} I_i$$', 'Interpreter', 'latex', ...
%        'Location', 'Best', 'NumColumns', 2);
% %        '$$\sum_{i=1}^{m} I_i/\sum_{i=1}^{m} N_i$$', 'Interpreter', 'latex', ...
%        'Location', 'Best', 'NumColumns', 2);
h4_figurename_fig = ['/result figures/NO2SimulationPeakVsTimeVariesOverEpsilon/fig/Fig_2_B4_cons_beta_TD_Flux_ind_and_cum_Dynmics_without_legend_v_value_' num2str(v) '_matrix_size_' num2str(m) '_Network_' str '_Size_' num2str(T) '_infected_strat_node_' num2str(v) '_Het_Node_' num2str(u) '_beta_i_' num2str(beta0) '_GT_' num2str(gamma) '.fig'];
h4_figurename_png = ['/result figures/NO2SimulationPeakVsTimeVariesOverEpsilon/png/Fig_2_B4_cons_beta_TD_Flux_ind_and_cum_Dynmics_without_legend_v_value_' num2str(v) '_matrix_size' num2str(m) '_Network_' str '__Size_' num2str(T) '_infected_start_node_' num2str(v) '_Het_Node_' num2str(u) '_beta_i_' num2str(beta0) '_GT_' num2str(gamma) '.png'];
saveas(h4,[pwd h4_figurename_fig])
saveas(h4,[pwd h4_figurename_png])
hold off

% 
% %% 
% highlight_nodes    = [55, 72, 14, 9, 16, 60, 42, 43];
% % Corresponding county names
% highlight_counties = {'Oklahoma', 'Tulsa', 'Cleveland', 'Canadian', 'Comanche', 'Payne', 'Logan', 'Love'};
% 
% % Define highlight colors (make sure they pop in print)
% highlight_colors = [
%     1, 0, 0;                  % red
%     0, 0, 1;                  % blue
%     0, 0.6, 0;                % green (brighter than pure [0,1,0])
%     0, 0, 0;                  % black
%     1, 0, 1;                  % magenta
%     20/255, 90/255, 50/255;   % dark green
%     237/255, 178/255, 32/255; % golden yellow
%     126/255, 47/255, 142/255  % purple
%     ];
% 
% h4 = figure(4); clf;
% hold on
% legend_entries = cell(length(highlight_nodes) + 1, 1);
% 
% for k = 1:length(highlight_nodes)
%     idx = highlight_nodes(k);
% 
%     plot(t, ...
%          y_k(:, m + idx) / vec_N(idx), ...
%          '-', ...
%          'Color', highlight_colors(k, :), ...
%          'LineWidth', 3);
% 
%     legend_entries{k} = highlight_counties{k};
% end
% 
% % ---- Total normalized infected population
% plot(t, sum(y_k(:, m+1:2*m), 2) / sum(vec_N), ...
%      'k--', 'LineWidth', 3);
% 
% legend_entries{end} = 'Total (normalized)';
% set(gca, 'FontSize', 30, 'LineWidth', 2, 'FontName', 'Times');
% ylim([0 ylimit]);
% xlim([1 plot_time_end]);
% set(gca, 'YTick', 0:0.2:ylimit);
% set(gca, 'XTick', [1:11:time_end, time_end]);
% legend(legend_entries, ...
%        'Location', 'northwest', ...
%        'NumColumns', 2, ...
%        'Box', 'on');
% 
% %% Highlighted nodes
% highlight_nodes    = [55, 72, 14, 9, 16, 60, 42, 43];
% highlight_counties = {'Oklahoma', 'Tulsa', 'Cleveland', ...
%                       'Canadian', 'Comanche', 'Payne', ...
%                       'Logan', 'Love'};
% 
% % Highlight colors (print-safe)
% highlight_colors = [
%     1, 0, 0;                  
%     0, 0, 1;                  
%     0, 0.6, 0;                
%     0, 0, 0;                  
%     1, 0, 1;                  
%     20/255, 90/255, 50/255;   
%     237/255, 178/255, 32/255; 
%     126/255, 47/255, 142/255  
% ];
% 
% h44 = figure(5); clf;
% hold on
% 
% %% --------------------------------------------------
% % 1) Plot ALL non-highlighted nodes in light gray
% %% --------------------------------------------------
% all_nodes = 1:m;
% other_nodes = setdiff(all_nodes, highlight_nodes);
% % Plot background counties and save handle for legend
% hGray = [];  % To store the handle for "Others"
% for i = other_nodes
%    h =  plot(t, ...
%          y_k(:, m + i) / vec_N(i), ...
%          '-', ...
%          'Color', [0.85, 0.85, 0.85], ...
%          'LineWidth', 3);
%             hGray = h; % Save the handle of the first gray line for legend
% 
% end
% 
% %% --------------------------------------------------
% % 2) Plot highlighted nodes ON TOP
% %% --------------------------------------------------
% legend_entries = cell(length(highlight_nodes) + 1, 1);
% 
% for k = 1:length(highlight_nodes)
%     idx = highlight_nodes(k);
% 
%     plot(t, ...
%          y_k(:, m + idx) / vec_N(idx), ...
%          '-', ...
%          'Color', highlight_colors(k, :), ...
%          'LineWidth', 3);
% 
%     legend_entries{k} = highlight_counties{k};
% end
% 
% %% --------------------------------------------------
% % 3) Total normalized infected population
% %% --------------------------------------------------
% plot(t, sum(y_k(:, m+1:2*m), 2) / sum(vec_N), ...
%      'k--', 'LineWidth', 3);
% 
% legend_entries{end} = 'Total (normalized)';
% 
% %% --------------------------------------------------
% % 4) Axis & legend formatting
% %% --------------------------------------------------
% set(gca, 'FontSize', 30, 'LineWidth', 2, 'FontName', 'Times');
% ylim([0 ylimit]);
% xlim([1 plot_time_end]);
% set(gca, 'YTick', 0:0.2:ylimit);
% set(gca, 'XTick', [1:11:time_end, time_end]);
% % 
% % legend(legend_entries, ...
% %        'Location', 'northwest', ...
% %        'NumColumns', 2, ...
% %        'Box', 'on');
% legend_handles = [hHighlights; hGray];
% legend_labels = [highlight_counties, {'Other counties'}];
% legend(legend_handles, legend_labels, 'Location', 'best', 'FontSize', 14, 'Interpreter', 'none');
% hold off;
% box on;
% h44_figurename_fig = ['/result figures/NO2SimulationPeakVsTimeVariesOverEpsilon/fig/All_Dynmics_with1_normalz_legend_v_value_' num2str(v) '_matrix_size_' num2str(m) '_Network_' str '_Size_' num2str(T) '_infected_strat_node_' num2str(v) '_Het_Node_' num2str(u) '_TR_' num2str(beta0) '_GT_' num2str(gamma) '.fig'];
% h44_figurename_png = ['/result figures/NO2SimulationPeakVsTimeVariesOverEpsilon/png/All_Dynmics_with1_normalz_legend_v_value_' num2str(v) '_matrix_size' num2str(m) '_Network_' str '__Size_' num2str(T) '_infected_start_node_' num2str(v) '_Het_Node_' num2str(u) '_TR_' num2str(beta0) '_GT_' num2str(gamma) '.png'];
% saveas(h44,[pwd h44_figurename_fig])
% saveas(h44,[pwd h44_figurename_png])
% hold off

%% Highlighted nodes
highlight_nodes    = [55, 72, 14, 9, 16, 60, 42, 43];
highlight_counties = {'Oklahoma', 'Tulsa', 'Cleveland', ...
                      'Canadian', 'Comanche', 'Payne', ...
                      'Logan', 'Love'};

% Highlight colors (print-safe)
highlight_colors = [
    1, 0, 0;                  
    0, 0, 1;                  
    0, 0.6, 0;                
    0, 0, 0;                  
    1, 0, 1;                  
    20/255, 90/255, 50/255;   
    237/255, 178/255, 32/255; 
    126/255, 47/255, 142/255  
];

h44 = figure(5); clf;
hold on

%% --------------------------------------------------
% 1) Plot ALL non-highlighted nodes (light gray)
%% --------------------------------------------------
all_nodes   = 1:m;
other_nodes = setdiff(all_nodes, highlight_nodes);

hGray = [];   % handle for "Other counties"

for i = other_nodes
    h = plot(t, ...
             y_k(:, m + i) / vec_N(i), ...
             '-', ...
             'Color', [0.85, 0.85, 0.85], ...
             'LineWidth', 3);
    if isempty(hGray)
        hGray = h;   % save only once for legend
    end
end

%% --------------------------------------------------
% 2) Plot highlighted nodes ON TOP
%% --------------------------------------------------
hHighlights = gobjects(length(highlight_nodes),1);

for k = 1:length(highlight_nodes)
    idx = highlight_nodes(k);

    hHighlights(k) = plot(t, ...
        y_k(:, m + idx) / vec_N(idx), ...
        '-', ...
        'Color', highlight_colors(k,:), ...
        'LineWidth', 3);
end

%% --------------------------------------------------
% 3) Total normalized infected population
%% --------------------------------------------------
hTotal = plot(t, ...
    sum(y_k(:, m+1:2*m), 2) / sum(vec_N), ...
    'k--', 'LineWidth', 3);

%% --------------------------------------------------
% 4) Axis formatting
%% --------------------------------------------------
set(gca, 'FontSize', 40, 'LineWidth', 2, 'FontName', 'Times');
% ylim([0 0.3]);
% xlim([1 plot_time_end]);
xlim([1 40]);
% set(gca, 'YTick', 0:0.1:ylimit);
% set(gca, 'XTick', [1:11:time_end, time_end]);
set(gca, 'XTick', [0:10:70]);
box on;

%% --------------------------------------------------
% 5) Legend (clean & minimal)
%% --------------------------------------------------
% legend_handles = [hHighlights; hGray; hTotal];
% legend_labels  = [highlight_counties, ...
%                   {'Other counties'}, ...
%                   {'Total (normalized)'}];

% legend(legend_handles, legend_labels, ...
%        'Location', 'northwest', ...
%        'NumColumns', 2, ...
%        'Box', 'on', ...
%        'Interpreter', 'none');

hold off
h44_figurename_fig = ['/result figures/NO2SimulationPeakVsTimeVariesOverEpsilon/fig/All_Dynmics_with1_normalz_legend_v_value_' num2str(v) '_matrix_size_' num2str(m) '_Network_' str '_Size_' num2str(T) '_infected_strat_node_' num2str(v) '_Het_Node_' num2str(u) '_TR_' num2str(beta0) '_GT_' num2str(gamma) '.fig'];
h44_figurename_png = ['/result figures/NO2SimulationPeakVsTimeVariesOverEpsilon/png/All_Dynmics_with1_normalz_legend_v_value_' num2str(v) '_matrix_size' num2str(m) '_Network_' str '__Size_' num2str(T) '_infected_start_node_' num2str(v) '_Het_Node_' num2str(u) '_TR_' num2str(beta0) '_GT_' num2str(gamma) '.png'];
saveas(h44,[pwd h44_figurename_fig])
saveas(h44,[pwd h44_figurename_png])
hold off

% %% Some particulr nodes infectious dynamics
% %figure('DefaultAxesFontSize',20) % this will change the Font Size on the both Axis\
% h44 = figure(5)
% plot(t,y_k(:,m+top_five_pop_flow_nodes(1)),'r-','MarkerFaceColor','r','LineWidth',3);
% hold on
% set(gca,'LineWidth',2);
% % set(gca,'FontWeight','Bold');
% hold on
% plot(t,y_k(:,m+top_five_pop_flow_nodes(2)),'b-','MarkerFaceColor','b','LineWidth',3);
% plot(t,y_k(:,m+top_five_pop_flow_nodes(3)),'g-','MarkerFaceColor','g','LineWidth',3);
% plot(t,y_k(:,m+top_five_pop_flow_nodes(4)),'k-','MarkerFaceColor','k','LineWidth',3);
% plot(t,y_k(:,m+top_five_pop_flow_nodes(5)),'m-','MarkerFaceColor','m','LineWidth',3);
% %plot(t,sum(y_k(:,m+1:2*m),2),'r--','MarkerFaceColor','r','LineWidth',3);
% %plot(t,y_k(:,m+top_five_pop_flow_nodes(6))/vec_N(m),'k-','MarkerFaceColor','k','LineWidth',3);
% %grid on;
% % set(gca, 'fontsize', 20); % gca: get current axis
% %set(gca, 'fontsize', 20, 'fontweight','bold'); % gca: get current axis
% %set(gca, 'YTick', 0:.1:ylimit); % Full form of gca is (Get Current Axes)
% set(gca, 'XTick', 0:8: time_end); % Full form of gca is (Get Current Axes)
% % ylim([0,ylimit]);
% xlim([1,plot_time_end])
% set(gca, 'fontsize', 30,'LineWidth',2,'FontName', 'Times'); % gca: get current axis
% %set(gca, 'fontsize', 22, 'fontweight','bold','LineWidth',2); % gca: get current axis
% % legend only sum of total cases in each day node sum
% % legend(['$I_{' num2str(top_five_pop_flow_nodes(1)) '}$'], ...
% %        ['$I_{' num2str(top_five_pop_flow_nodes(2)) '}$'], ...
% %        ['$I_{' num2str(top_five_pop_flow_nodes(3)) '}$'], ...
% %        ['$I_{' num2str(top_five_pop_flow_nodes(4)) '}$'], ...
% %        ['$$I_{' num2str(top_five_pop_flow_nodes(5)) '}$$'], ...
% %        '$$\sum_{i=1}^{m} I_i$$', 'Interpreter', 'latex', ...
% %        'Location', 'Best', 'NumColumns', 2);
% % %        '$$\sum_{i=1}^{m} I_i/\sum_{i=1}^{m} N_i$$', 'Interpreter', 'latex', ...
% %        'Location', 'Best', 'NumColumns', 2);
% h44_figurename_fig = ['/result figures/NO2SimulationPeakVsTimeVariesOverEpsilon/fig/Fig_2_B_cons_beta_TD_Flux_ind_and_top_five_Dynmics_without_normalization_legend_v_value_' num2str(v) '_matrix_size_' num2str(m) '_Network_' str '_Size_' num2str(T) '_infected_strat_node_' num2str(v) '_Het_Node_' num2str(u) '_TR_' num2str(beta0) '_GT_' num2str(gamma) '.fig'];
% h44_figurename_png = ['/result figures/NO2SimulationPeakVsTimeVariesOverEpsilon/png/Fig_2_B_cons_beta_TD_Flux_ind_and_top_five_Dynmics_without_normalization_legend_v_value_' num2str(v) '_matrix_size' num2str(m) '_Network_' str '__Size_' num2str(T) '_infected_start_node_' num2str(v) '_Het_Node_' num2str(u) '_TR_' num2str(beta0) '_GT_' num2str(gamma) '.png'];
% saveas(h44,[pwd h44_figurename_fig])
% saveas(h44,[pwd h44_figurename_png])
% hold off

% 
% %% Some particulr nodes infectious dynamics
% %figure('DefaultAxesFontSize',20) % this will change the Font Size on the both Axis\
% h5 = figure(6)
% plot(t,y_k(:,m+1)/vec_N(1),'r-','MarkerFaceColor','r','LineWidth',3);
% hold on
% set(gca,'LineWidth',2);
% % set(gca,'FontWeight','Bold');
% hold on
% plot(t,y_k(:,m+2)/vec_N(2),'b-','MarkerFaceColor','b','LineWidth',3);
% plot(t,y_k(:,m+3)/vec_N(3),'g-','MarkerFaceColor','g','LineWidth',3);
% plot(t,y_k(:,m+4)/vec_N(4),'k-','MarkerFaceColor','k','LineWidth',3);
% plot(t,y_k(:,m+5)/vec_N(5),'m-','MarkerFaceColor','m','LineWidth',3);
% plot(t,sum(y_k(:,m+1:2*m),2)/sum(vec_N),'r--','MarkerFaceColor','r','LineWidth',3);
% %plot(t,y_k(:,m+top_five_pop_flow_nodes(6))/vec_N(m),'k-','MarkerFaceColor','k','LineWidth',3);
% %grid on;
% % set(gca, 'fontsize', 20); % gca: get current axis
% %set(gca, 'fontsize', 20, 'fontweight','bold'); % gca: get current axis
% % set(gca, 'YTick', 0:.1:ylimit); % Full form of gca is (Get Current Axes)
% set(gca, 'XTick', 0:88: time_end); % Full form of gca is (Get Current Axes)
% % ylim([0,ylimit]);
% xlim([1,plot_time_end])
% set(gca, 'fontsize', 30,'LineWidth',2,'FontName', 'Times'); % gca: get current axis
% %set(gca, 'fontsize', 22, 'fontweight','bold','LineWidth',2); % gca: get current axis
% % legend only sum of total cases in each day node sum
% % legend(['$I_{' num2str(top_five_pop_flow_nodes(1)) '}$'], ...
% %        ['$I_{' num2str(top_five_pop_flow_nodes(2)) '}$'], ...
% %        ['$I_{' num2str(top_five_pop_flow_nodes(3)) '}$'], ...
% %        ['$I_{' num2str(top_five_pop_flow_nodes(4)) '}$'], ...
% %        ['$$I_{' num2str(top_five_pop_flow_nodes(5)) '}$$'], ...
% %        '$$\sum_{i=1}^{m} I_i$$', 'Interpreter', 'latex', ...
% %        'Location', 'Best', 'NumColumns', 2);
% % %        '$$\sum_{i=1}^{m} I_i/\sum_{i=1}^{m} N_i$$', 'Interpreter', 'latex', ...
% %        'Location', 'Best', 'NumColumns', 2);
% h5_figurename_fig = ['/result figures/NO2SimulationPeakVsTimeVariesOverEpsilon/fig/Fig_2_B6_cons_beta_TD_Flux_ind_and_cum_Dynmics_without_legend_v_value_' num2str(v) '_matrix_size_' num2str(m) '_Network_' str '_Size_' num2str(T) '_infected_strat_node_' num2str(v) '_Het_Node_' num2str(u) '_TR_' num2str(beta0) '_GT_' num2str(gamma) '.fig'];
% h5_figurename_png = ['/result figures/NO2SimulationPeakVsTimeVariesOverEpsilon/png/Fig_2_B6_cons_beta_TD_Flux_ind_and_cum_Dynmics_without_legend_v_value_' num2str(v) '_matrix_size' num2str(m) '_Network_' str '__Size_' num2str(T) '_infected_start_node_' num2str(v) '_Het_Node_' num2str(u) '_TR_' num2str(beta0) '_GT_' num2str(gamma) '.png'];
% saveas(h5,[pwd h5_figurename_fig])
% saveas(h5,[pwd h5_figurename_png])
% hold off
% 
% 
% 
% %% Some particulr nodes infectious dynamics
% %figure('DefaultAxesFontSize',20) % this will change the Font Size on the both Axis\
% h6 = figure(7)
% plot(t,y_k(:,m+1),'r-','MarkerFaceColor','r','LineWidth',3);
% hold on
% set(gca,'LineWidth',2);
% % set(gca,'FontWeight','Bold');
% hold on
% plot(t,y_k(:,m+2),'b-','MarkerFaceColor','b','LineWidth',3);
% plot(t,y_k(:,m+3),'g-','MarkerFaceColor','g','LineWidth',3);
% plot(t,y_k(:,m+4),'k-','MarkerFaceColor','k','LineWidth',3);
% plot(t,y_k(:,m+5),'m-','MarkerFaceColor','m','LineWidth',3);
% %plot(t,sum(y_k(:,m+1:2*m),2),'r--','MarkerFaceColor','r','LineWidth',3);
% %plot(t,y_k(:,m+top_five_pop_flow_nodes(6))/vec_N(m),'k-','MarkerFaceColor','k','LineWidth',3);
% %grid on;
% % set(gca, 'fontsize', 20); % gca: get current axis
% %set(gca, 'fontsize', 20, 'fontweight','bold'); % gca: get current axis
% %set(gca, 'YTick', 0:.1:ylimit); % Full form of gca is (Get Current Axes)
% set(gca, 'XTick', 0:8: time_end); % Full form of gca is (Get Current Axes)
% % ylim([0,ylimit]);
% xlim([1,plot_time_end])
% set(gca, 'fontsize', 30,'LineWidth',2,'FontName', 'Times'); % gca: get current axis
% %set(gca, 'fontsize', 22, 'fontweight','bold','LineWidth',2); % gca: get current axis
% % legend only sum of total cases in each day node sum
% % legend(['$I_{' num2str(top_five_pop_flow_nodes(1)) '}$'], ...
% %        ['$I_{' num2str(top_five_pop_flow_nodes(2)) '}$'], ...
% %        ['$I_{' num2str(top_five_pop_flow_nodes(3)) '}$'], ...
% %        ['$I_{' num2str(top_five_pop_flow_nodes(4)) '}$'], ...
% %        ['$$I_{' num2str(top_five_pop_flow_nodes(5)) '}$$'], ...
% %        '$$\sum_{i=1}^{m} I_i$$', 'Interpreter', 'latex', ...
% %        'Location', 'Best', 'NumColumns', 2);
% % %        '$$\sum_{i=1}^{m} I_i/\sum_{i=1}^{m} N_i$$', 'Interpreter', 'latex', ...
% %        'Location', 'Best', 'NumColumns', 2);
% h6_figurename_fig = ['/result figures/NO2SimulationPeakVsTimeVariesOverEpsilon/fig/Fig_2_B_cons_beta_TD_Flux_ind_and_cum_Dynmics_without_normalization_legend_v_value_' num2str(v) '_matrix_size_' num2str(m) '_Network_' str '_Size_' num2str(T) '_infected_strat_node_' num2str(v) '_Het_Node_' num2str(u) '_beta_i_' num2str(beta0) '_GT_' num2str(gamma) '.fig'];
% h6_figurename_png = ['/result figures/NO2SimulationPeakVsTimeVariesOverEpsilon/png/Fig_2_B_cons_beta_TD_Flux_ind_and_cum_Dynmics_without_normalization_legend_v_value_' num2str(v) '_matrix_size' num2str(m) '_Network_' str '__Size_' num2str(T) '_infected_start_node_' num2str(v) '_Het_Node_' num2str(u) '_beta_i_' num2str(beta0) '_GT_' num2str(gamma) '.png'];
% saveas(h6,[pwd h6_figurename_fig])
% saveas(h6,[pwd h6_figurename_png])
% hold off
% 
% 
% 
% %% new 
% 
% 
% h00 = figure(8)
% plot(t,y_k(:,m+top_five_pop_flow_nodes(1))/vec_N(top_five_pop_flow_nodes(1)),'r-','MarkerFaceColor','r','LineWidth',3);
% hold on
% set(gca,'LineWidth',2);
% hold on
% plot(t,y_k(:,m+top_five_pop_flow_nodes(2))/vec_N(top_five_pop_flow_nodes(2)),'b-','MarkerFaceColor','b','LineWidth',3);
% plot(t,y_k(:,m+top_five_pop_flow_nodes(3))/vec_N(top_five_pop_flow_nodes(3)),'g-','MarkerFaceColor','g','LineWidth',3);
% plot(t,y_k(:,m+top_five_pop_flow_nodes(4))/vec_N(top_five_pop_flow_nodes(4)),'k-','MarkerFaceColor','k','LineWidth',3);
% plot(t,y_k(:,m+top_five_pop_flow_nodes(5))/vec_N(top_five_pop_flow_nodes(5)),'m-','MarkerFaceColor','m','LineWidth',3);
% plot(t,sum(y_k(:,m+1:2*m),2)/sum(vec_N),'r--','MarkerFaceColor','r','LineWidth',3);
% %plot(t,y_k(:,m+top_five_pop_flow_nodes(6))/vec_N(m),'k-','MarkerFaceColor','k','LineWidth',3);
% %grid on;
% % set(gca, 'fontsize', 20, 'fontweight','bold'); % gca: get current axis
% %set(gca, 'fontsize', 30); % gca: get current axis
% %set(gca, 'YTick', 0:.1:ylimit); % Full form of gca is (Get Current Axes)
% % Set X-axis ticks and labels
% % tick_positions1 = 1:88:time_end;
% % tick_positions11 = [tick_positions1, time_end];  % Add end time to the tick positions
% % tick_labels = num2cell(tick_positions11);  % Convert tick positions to cell array of strings
% % set(gca, 'XTick', tick_positions11, 'XTickLabel', tick_labels); % Set the tick positions and labels
% set(gca, 'YTick', 0:.2:ylimit); % Full form of gca is (Get Current Axes)
% % set(gca, 'XTick', 1:11: time_end); % Full form of gca is (Get Current Axes)
% set(gca, 'XTick', [1:11:time_end, time_end])
% ylim([0,ylimit]);
% % ylim([0,ylimit]);
% xlim([0,time_end])
% set(gca, 'fontsize', 30,'LineWidth',2,'FontName', 'Times'); % gca: get current axis
% % legend only sum of total cases in each day node sum
% % legend(['$I_{' num2str(top_five_pop_flow_nodes(1)) '}$'], ...
% %        ['$I_{' num2str(top_five_pop_flow_nodes(2)) '}$'], ...
% %        ['$I_{' num2str(top_five_pop_flow_nodes(3)) '}$'], ...
% %        ['$I_{' num2str(top_five_pop_flow_nodes(4)) '}$'], ...
% %        ['$$I_{' num2str(top_five_pop_flow_nodes(5)) '}$$'], ...
% %        '$$\sum_{i=1}^{m} I_i$$', 'Interpreter', 'latex', ...
% %        'Location', 'Best', 'NumColumns', 2);
% % %       '$$\sum_{i=1}^{m} I_i/\sum_{i=1}^{m} N_i$$', 'Interpreter', 'latex', ...
% % %       'Location', 'Best', 'NumColumns', 2);
% hold off
% h00_figurename_fig = ['/result figures/NO2SimulationPeakVsTimeVariesOverEpsilon/fig/Fig_A_B_TD_dynmics_cons_beta_TD_flux_' num2str(phi_0) '_matrix_size_' num2str(m) '_Network_' str '_Size_' num2str(T) '_infected_strat_node_' num2str(v) '_Het_Node_Mult_' num2str(zeta) '_beta_i_' num2str(beta0) '_GT_' num2str(gamma) '.fig'];
% h00_figurename_png = ['/result figures/NO2SimulationPeakVsTimeVariesOverEpsilon/png/Fig_A_B_TD_dynmics_cons_beta_TD_flux_' num2str(phi_0) '_matrix_size' num2str(m) '_Network_' str '__Size_' num2str(T) '_infected_start_node_' num2str(v) '_Het_Node_Mult_' num2str(zeta) '_beta_i_' num2str(beta0) '_GT_' num2str(gamma) '.png'];
% saveas(h00,[pwd h00_figurename_fig])
% saveas(h00,[pwd h00_figurename_png])
% 
% 
% 
% 
% 
% %% find I_max and time of I_max using command
% for j = 1:m
%     B = y_k(:, m + j);
%     [max_x, idx] = max(B);
%     maxt = t(idx);
%     D(j) = max_x;
%     T(j) = maxt;
% end
% 
% x = T;
% %% Maximum diffusion number is the peak/population effect
% I_diff = zeros(m,1); 
% for i =1:m
%    I_diff(i) = D(i)/vec_N(i) 
% end
% y=I_diff;
% 
% %% here b=[b(1) b(2) b(3)] is the unknown parameter to be estimated by fminsearch
% % this function for FC data 
% f = @(b,x) b(1).*exp(b(2).*x)+b(3);                                       % Objective Function
% % this function for TD mobility data
% %f = @(b,x) b(1).*exp(b(2).*x)+b(3);                                      % Objective Function
% %Bmin = fminsearch(@(b) norm(y - f(b,x)), [-200; -1; 100]);                 % Estimate Parameters
% %Real cases exp fit 
% %BmiWhatsApp Image 2023-08-25 at 11.29.05 (1)n = fminsearch(@(b) norm(y - f(b,x)), [.75; .0005; 1.5]);                 % Estimate Parameters
% %TD orig pop fit
% %Bmin = fminsearch(@(b) norm(y - f(b,x)), [1; .9; 0]);                 % Estimate Parameters
% %g= Bm(1).*exp(Bm(2).*x)+Bm(3);
% 
% 
% 
% mdl = fitnlm(x(:), y(:), f,  [2; .009; .5]) 
% B = mdl.Coefficients.Estimate;
% Bci = coefCI(mdl);
% CoefficientMatrix = array2table([Bci(:,1) B, Bci(:,2)], 'VariableNames',{'-95% CI','Beta','+95% CI'})
% xv = linspace(min(x), max(x));
% [yv,yci] = predict(mdl,xv(:));
% 
% 
% h2 = figure(9);
% %top5_idx = [42; 55; 72; 14 ; 9; 16; 41] 
% 
% % Define colors as a cell array of strings
% colors = {'#145A32', 'r', 'b', 'g', 'k', 'm', '#EDB120'}; 
% % forest green, Red, Blue, Green, Black, Magenta, snd "gold" or "amber". 
% % Logan, Oklahoma, 
% 
% % Normalize y-values for use as marker sizes (or scale them)
% min_size = 5;  % Minimum marker size
% max_size = 40; % Maximum marker size
% scaled_sizes = min_size + (y - min(y)) / (max(y) - min(y)) * (max_size - min_size);
% 
% % Identify the top 5 peaks by finding the indices of the top 5 y-values
% [~, idx] = sort(y, 'descend'); % Sort y in descending order
% %top5_idx = idx(1:5); % Indices of top 5 peaks
% % here we have manually define the peak nodes from beta_0 = .01428
% top5_idx = [42; 55; 72; 14 ; 9; 16; 41] 
% % Define different colors for the top 5 peaks
% % colors = ['#7E2F8E', 'r', 'b', 'g', 'k' 'm', '#EDB120']; % Red, Blue, Green, Black, Magenta
% % colors = ['#7E2F8E', 'r', 'b', 'g', 'k', 'm', '#EDB120'];
% 
% 
% % Plot all data points with different circle sizes
% hp{1} = scatter(x, y, scaled_sizes * 10, 'o', 'MarkerEdgeColor', '#A2142F', 'MarkerFaceColor', 'w');
% hold on;
% 
% % Plot the top 5 peak points with different colors and scaled sizes
% for i = 1:7
%     scatter(x(top5_idx(i)), y(top5_idx(i)), scaled_sizes(top5_idx(i)) * 10, 'o', ... % Scale sizes for top 5
%         'MarkerEdgeColor', colors{i}, 'MarkerFaceColor', 'w', 'LineWidth', 2); % Assign different colors to top 5
%     % Print the value of y at the top 5 indices
%     fprintf('Value of y at index %d: %.4f\n', top5_idx(i), y(top5_idx(i)));
% end
% 
% % Plot the line connecting the points
% hp{7} = plot(xv, yv, '-r', 'LineWidth', 3);
% set(gca, 'YTick', 0:.2:ylimit); % Full form of gca is (Get Current Axes)
% % Set plot limits and labels
% ylim([0, 0.6]);
% % xlim([18, 55]);
% set(gca, 'XTick', 0:5:time_end);
% 
% % Adjust axis properties
% set(gca, 'fontsize', 30, 'FontName', 'Times', 'LineWidth', 2, 'Box', 'on'); % Add square box and set line width
% 
% % Add text to the plot
% text(27, 205, sprintf('f(x) = %.1f\\cdote^{%.3f\\cdotx}%+.1f', B));
% 
% hold off;
% 
% % Save figures
% h2_figurename_fig = ['/result figures/NO2SimulationPeakVsTimeVariesOverEpsilon/fig/Fig_2_C_D_diff_top_7_colors_Exp_fit_I_diff_' num2str(m) '_Network_' str '_Size_' num2str(time_end) '_infected_start_node_' num2str(v) '_Het_Node' num2str(u) '_beta_i_' num2str(beta0) '_GT_' num2str(gamma) '.fig'];
% saveas(h2, [pwd h2_figurename_fig]);
% 
% h2_figurename_png = ['/result figures/NO2SimulationPeakVsTimeVariesOverEpsilon/png/Fig_2_C_D_diff_top_7_colors_Exp_fit_I_diff_' num2str(m) '_Network_' str '_Size_' num2str(time_end) '_infected_start_node_' num2str(v) '_Het_Node' num2str(u) '_beta_i_' num2str(beta0) '_GT_' num2str(gamma) '.png'];
% saveas(h2, [pwd h2_figurename_png]);
% 
% 
% 
% node_vec = linspace(1,m,m);
% 
% 
% %node_vec = linspace(1, m,m)
% h3 = figure(10)
% %I=linspace(1,11,11)
% data = I_diff
% %stem(data,'b','filled','LineWidth',2)
% stem(node_vec, data,'r','filled','LineWidth',3)
% xlim([0,m])
% ylim([0, 0.6])
% set(gca, 'XTick', 0:11: m); % Full form of gca is (Get Current Axes)
% % xlim([0,m])
% % ylim([0, ylimt])
% % set(gca, 'XTick', 0:floor(tickSpacing):m);
% % xlabel('Node i')
% % ylabel('Time of peak')
% % Adjust axis properties
% set(gca, 'fontsize', 30, 'FontName', 'Times', 'LineWidth', 2, 'Box', 'on'); % Add square box and set line width
% h3_figurename_fig = ['/result figures/NO2SimulationPeakVsTimeVariesOverEpsilon/fig/Fig_2_C_D_I_diff_vs_Node_Fig2_m_' num2str(m) '_Network_' str '_Size_' num2str(time_end) '_infected_start_node_' num2str(v) '_Het_Node' num2str(u) '_RR_' num2str(gamma) '.fig'];
% saveas(h3,[pwd h3_figurename_fig])
% h3_figurename_png = ['/result figures/NO2SimulationPeakVsTimeVariesOverEpsilon/png/Fig_2_C_D_I_diff_vs_Node_Fig2_m_' num2str(m) '_Network_' str '_Size_' num2str(time_end) '_infected_start_node_' num2str(v) '_Het_Node' num2str(u) '_RR_' num2str(gamma) '.png'];
% saveas(h3,[pwd h3_figurename_png])
% 
% 
% 
% 
% 
% 
% % Define the data
% 
% % Create a table with the data
% cons_beta_SIR_network_simulation_node_peak_fraction = table(node_vec', name_origin_list', I_diff, x','VariableNames', {'node_vec', 'name_origin_list', 'peak_fraction', 'time_peak_fraction'});
% 
% % Define the directory and file name
% directory = 'frac-cum y_k data/NO2SimulationPeakVsTimeVariesOverEpsilon';
% % file_path = fullfile(directory, 'cons_beta_TD_mobility_SIR_network_simulation_node_peak_fraction_values.csv');
% file_path = fullfile(directory, sprintf('cons_beta_%s_sim_node_peak_frac_gamma_values_%g_beta0_values_%g.csv', str, gamma, beta0));
% 
% % Ensure the directory exists
% if ~exist(directory, 'dir')
%     mkdir(directory);
% end
% 
% % Write the table to a CSV file
% writetable(cons_beta_SIR_network_simulation_node_peak_fraction, file_path);
% 
