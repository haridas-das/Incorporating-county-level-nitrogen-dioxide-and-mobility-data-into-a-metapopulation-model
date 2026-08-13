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
% gamma         = 1.4; %0.1429;  %1/7                wekly recovered
gamma         = 1; %0.1429;  %1/7   
beta_1_vec    = linspace(0,0.1428,1); % simulate the infection rates% Recovery rate
epsilon       = 0.1;
bline_beta0   = 0.8; %0.1428; % 0.5;
omega         = epsilon;

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

% T = size(all_counties_full_ts, 1);          % number of time points


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
u = top_five_pop_flow_nodes(2);                                                                     % heteroginous node
v = top_five_pop_flow_nodes(2);                                                                     % Node where the infectious started , put Oklahoma County  55
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

%% Data for heterogeneous transmission node
% static_matrix_phi = avg_pre_lockdown_flux_matrix;
% beta0             = beta_1_vec;
% vec_beta          = beta0*ones(m,1)+ epsilon*NO2_stats_table.Mean_NO2;
% [t,y_k] = ode45(@(t,y) SirGenNetworkModel(m,static_matrix_phi,y,vec_beta,vec_S0,vec_I0,vec_N,gamma,zeta),tspan,y0);

bline_beta0_vec = [0:0.01:1.5];            % simuulate the flux only for the Unidirectional cycle networks
omega_vec       = [0:0.01:1.5];            % simulate the infection rates

% phi_0      = [0.01:0.01:0.99]         % simuulate the flux only for the Unidirectional cycle networks
% beta_1_vec =[0:0.01:1];               % simulate the infection rates
% Peak of the cumulative cases on the whole time domain
D = zeros(m,1);                     %finding the maximum infected
T = zeros(m,1);                     %time when the maximum infected occured

%% R_0 estimates 

% Inputs:
% beta_vec: M x 1 vector of beta_j
% matrix_phi: M x M flux matrix
% vec_N: M x 1 population vector N_i
% Np: M x 1 present population vector at nodes
% gamma: scalar

% Prelockdown_mobility_matrix_phi  = avg_post_lockdown_flux_matrix;
% M = length(vec_N);

% bline_beta0 = 0.5;
% omega       = 0;
% vec_beta   = bline_beta0*ones(m,1) + omega*NO2_data;
% Postlockdown_mobility_matrix_phi = avg_pre_lockdown_flux_matrix;
% 
% kappa = zeros(m,m);
% 
% %%Start: First define the present population at node i, which is also constant in time come from rest of the nodes
% Np = zeros(m,1); % initial present population at i
% 
% Npsum=0;
% for i=1:m
%     for k=1:m
%         Npsum=Npsum+matrix_phi(k,i).*vec_N(k);                             % phi(k,i)*vec_N(k,1)is the total population at vertex k that travels every day to vertex i
%     end
%     Np(i)=Npsum;
%     %break;
%     Npsum=0;
% end
% 
% for i = 1:M
%     for l = 1:M
%         s = 0;
%         for j = 1:M
%             s = s + vec_beta(j) * matrix_phi(i,j) * matrix_phi(l,j) * (vec_N(i)/Np(j));
%         end
%         kappa(i,l) = s / gamma;
%     end
% end
% 
% Prelockdown_mobility_R0 = max(eig(kappa));



% Preallocate empty arrays
bline_beta0_selected = [];
omega_selected       = [];
R0_selected          = [];

tolerance = 0.001;  % how close to 1 you consider "near 1"


for a = 1: length(bline_beta0_vec)
    for b = 1: length(omega_vec)
        omega         = omega_vec(b);
        
        bline_beta0   = bline_beta0_vec(a); %

        vec_beta   = bline_beta0*ones(m,1) + omega*NO2_data;
        
        avg_pre_lockdown_matrix_phi = avg_pre_lockdown_flux_matrix;
        
        kappa = zeros(m,m);
        
        %%Start: First define the present population at node i, which is also constant in time come from rest of the nodes
        Np = zeros(m,1); % initial present population at i
        
        Npsum = 0;
        for r = 1:m
            for s = 1:m
                Npsum = Npsum + avg_pre_lockdown_matrix_phi(s,r).*vec_N(s);                             % phi(k,i)*vec_N(k,1)is the total population at vertex k that travels every day to vertex i
            end
            Np(r) = Npsum;
            %break;
            Npsum = 0;
        end
        
        for i = 1:M
            for l = 1:M
                s = 0;
                for j = 1:M
                    s = s + vec_beta(j) * avg_pre_lockdown_matrix_phi(i,j) * avg_pre_lockdown_matrix_phi(l,j) * (vec_N(i)/Np(j));
                end
                kappa(i,l) = s / gamma;
            end
        end
        
        % Compute R0
        Prelockdown_mobility_R0 = max(eig(kappa));

        % Store all R0 values if needed
        R0_collection(b,a) = Prelockdown_mobility_R0;

        % --- Check condition ---
        if abs(Prelockdown_mobility_R0 - 1) <= tolerance
            % append this combination
            bline_beta0_selected(end+1,1) = bline_beta0;
            omega_selected(end+1,1)       = omega;
            R0_selected(end+1,1)          = Prelockdown_mobility_R0;
        end

        %% Precompute the flux matrix function
        % pre and  post lockdown avg mobility data
        phi_matrix_func = @(t) get_matrix_phi(t, lockdown_week_0_indexed, time_start, time_step_size, avg_pre_lockdown_flux_matrix, avg_post_lockdown_flux_matrix);
        %% PART 2: Solve the differential equation
        % [t,y_k] = ode45(@(t,y) SirGenNetworkModel(m,matrix_phi,y,vec_beta,vec_S0,vec_I0,vec_N,gamma,zeta),tspan,y0);
        [t, y_k] = ode45(@(t, y) SIR_Network_NO2_baseB0(t, m, phi_matrix_func(t), y, NO2_data, bline_beta0, omega, vec_S0, vec_I0, vec_N, gamma, zeta), tspan, y0);

        %% Epidemic size color maps
        % values: last row of the simulation i.e. total recovered or epidemic size
        %D = y_k((1+(time_end-1)/time_step_size),2*m+1:3*m)
        %    I_i^max  computation
        %epidemic_size_each_node = y_k((1+(time_end-1)/time_step_size),2*m+1:3*m)
        S0 = y0(1:m)';
        % t_n is the last index of the time
        t_n = ((time_end - time_start)/time_step_size)+1;
        ST = y_k(t_n,1:m)'; % remaining susceptible population in each node
        epidemic_size_each_node               = S0 - ST;
        total_netw_epidemic_size              = sum(epidemic_size_each_node);
        NetEpiSize_beta0_epsilon_vec(b,a)     = total_netw_epidemic_size;

        %% Infectious total peak color maps
        NetItotPeak_beta0_epsilon_vec(b,a)    = max(sum(y_k(:,m+1:2*m),2)); % 2 means it is the column sum, and 1 means row sum

        %% epidemic vs. non epidemic color maps
        matrix_max_1(b,a)    = max(sum(y_k(:,m+1:2*m),2)); % 2 means it is the column sum, and 1 means row sum
        if matrix_max_1(b,a)>= 2   % maximum number of infected I1(0)=1, I2(0)=0
            matrix_binary(b,a)=1
        else
            matrix_binary(b,a)=0
        end

    end
end

% ======================================
% Add highlighted points as circles when gamma = 1.4
% (x, y) = (beta0, omega)
% ======================================
% highlight_points = [
%     0.5 0.5
%     0.5 1
%     1 0.5
%     1 1
%     ];

% ======================================
% Add highlighted points as circles when gamma = 1.0
% (x, y) = (beta0, omega)
% ======================================
% highlight_points = [
%     0.2 0.5
%     0.2 1
%     1 0.5
%     1 1
%     ];

highlight_points = [
    0.2 0.5
    0.2 0.75
    1 0.25
    1 0.65
    ];

%% Save the simulated binary color map
figure('DefaultAxesFontSize',18) % this will change the Font Size on the both Axis
h = figure(1);
clf;
set(h, 'Color', 'w', 'Position', [120 120 900 700]);
pcolor(bline_beta0_vec, omega_vec, matrix_binary) % pcolor(X axis, Yaxis, ploting data)
%yline(0.76,'-r','LineWidth',2); %'Epidemic Control Threshold'
%yline((m*gamma)/(zeta+m-1),'-r','LineWidth',2); %'Epidemic Control Threshold'
% xlabel('\phi_1');
% ylabel('\beta_0');
hold on

plot(bline_beta0_selected, omega_selected, ...
    '--w', ...          % dashed black line
    'LineWidth', 3);

plot(highlight_points(:,1), highlight_points(:,2), ...
    'o', ...
    'MarkerSize', 14, ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', [1 1 1], ... % white-filled circles
    'LineWidth', 2);


hold off
set(gca, 'fontsize', 30,'LineWidth',2,'FontName', 'Times');
shading flat % make the color flat
figurename_fig = ['/NSF_M0_Numerical_Simulation/fig/Figure_epidemic_non_epidemic_binary_color_maps_bline_beta0_vec_end_' num2str(bline_beta0) '_Network_' str '_test_infected_node_' num2str(v) '_omega_vec_end_' num2str(omega) '_Recovery_rate_' num2str(gamma) '.fig'];
figurename_png = ['/NSF_M0_Numerical_Simulation/png/Figure_epidemic_non_epidemic_binary_color_maps_bline_beta0_vec_end_' num2str(bline_beta0) '_Network_' str '_test_infected_node_' num2str(v) '_omega_vec_end_' num2str(omega) '_Recovery_rate_' num2str(gamma) '.png'];
saveas(h,[pwd figurename_fig])
saveas(h,[pwd figurename_png])

%% Save the simulated matrix binary
save_folder = fullfile(pwd, 'NSF_M0_Numerical_Simulation');
if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end
filename_matrix = fullfile(save_folder, ...
    sprintf('binary_matrix_Pre_post_lockdown_mobility_network_bline_beta0_vec_end_%d_omega_vec_end_%d_gamma_%d.csv', bline_beta0, omega, gamma));
assert(exist('matrix_binary','var')==1, 'matrix_binary does not exist');
writematrix(matrix_binary, filename_matrix);
disp(['Saved to: ' filename_matrix])


%% NetEpiSize_beta0_epsilon_vec 
h = figure(2);
clf;
set(h, 'Color', 'w', 'Position', [100 100 900 700]);
% Main color map
pcolor(bline_beta0_vec, omega_vec, NetEpiSize_beta0_epsilon_vec);
% pcolor(omega_vec, bline_beta0_vec, NetEpiSize_beta0_epsilon_vec);
shading flat;
hold on;

% Nature-style perceptually uniform colormap

% colormap(turbo);    % Smooth and strong
colormap(parula);   % Smooth but slightly muted
% colormap(jet);      % Vibrant but not perceptually uniform

% colorbar;
plot(bline_beta0_selected, omega_selected, ...
    '--w', ...          % dashed black line
    'LineWidth', 3);

% Highlight selected parameter points
plot(highlight_points(:,1), highlight_points(:,2), ...
    'o', ...
    'MarkerSize', 14, ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', [1 1 1], ...
    'LineWidth', 2);

% Axis formatting (Nature style)
set(gca, ...
    'FontSize', 30, ...
    'LineWidth', 2, ...
    'FontName', 'Times', ...
    'Layer', 'top', ...
    'Box', 'on');

% xlabel('$\beta_0$', 'Interpreter', 'latex', 'FontSize', 34);
% ylabel('$\omega$',  'Interpreter', 'latex', 'FontSize', 34);

% Colorbar formatting
cb = colorbar;
cb.LineWidth = 1.5;
cb.FontSize  = 26;
cb.FontName  = 'Times';
%ylabel(cb, 'Final epidemic size', 'FontSize', 30);

% Fix color scaling (important for reproducibility)
caxis([ ...
    min(NetEpiSize_beta0_epsilon_vec(:)), ...
    max(NetEpiSize_beta0_epsilon_vec(:)) ]);

axis tight;

%% ============================================================
%  Save figures
% ============================================================

fig_folder = fullfile(pwd, 'NSF_M0_Numerical_Simulation', 'fig');
png_folder = fullfile(pwd, 'NSF_M0_Numerical_Simulation', 'png');

if ~exist(fig_folder, 'dir'); mkdir(fig_folder); end
if ~exist(png_folder, 'dir'); mkdir(png_folder); end

figurename_fig = fullfile(fig_folder, ...
    sprintf(['Figure_NetEpiSize_beta0_epsilon_vec_' ...
             'bline_beta0_%d_Network_%s_test_infected_node_%d_' ...
             'omega_%d_gamma_%d.fig'], ...
             bline_beta0, str, v, omega, gamma));

figurename_png = fullfile(png_folder, ...
    sprintf(['Figure_NetEpiSize_beta0_epsilon_vec_' ...
             'bline_beta0_%d_Network_%s_test_infected_node_%d_' ...
             'omega_%d_gamma_%d.png'], ...
             bline_beta0, str, v, omega, gamma));

saveas(h, figurename_fig);
saveas(h, figurename_png);

disp(['Figure saved as: ' figurename_png]);

%% ============================================================
%  Save epidemic size matrix (CSV)
% ============================================================

save_folder = fullfile(pwd, 'NSF_M0_Numerical_Simulation');
if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end

filename_matrix = fullfile(save_folder, ...
    sprintf(['NetEpiSize_beta0_epsilon_vec_' ...
             'Pre_post_lockdown_mobility_' ...
             'bline_beta0_%d_omega_%d_gamma_%d.csv'], ...
             bline_beta0, omega, gamma));

assert(exist('NetEpiSize_beta0_epsilon_vec','var')==1, ...
       'NetEpiSize_beta0_epsilon_vec does not exist');

writematrix(NetEpiSize_beta0_epsilon_vec, filename_matrix);

disp(['CSV saved to: ' filename_matrix]);

hold off;

%% Fig 3: NetItotPeak_beta0_epsilon_vec
h = figure(3);
clf;
set(h, 'Color', 'w', 'Position', [120 120 900 700]);

% Main heatmap
pcolor(bline_beta0_vec, omega_vec, NetItotPeak_beta0_epsilon_vec);
% pcolor(omega_vec, bline_beta0_vec, NetItotPeak_beta0_epsilon_vec);

shading flat;
hold on;

% Nature-style colormap for peak intensity
% colormap(magma);
colormap(jet);  % or parula, hot, cool, etc.

plot(bline_beta0_selected, omega_selected, ...
    '--w', ...          % dashed black line
    'LineWidth', 3);

% Highlight selected parameter points
plot(highlight_points(:,1), highlight_points(:,2), ...
    'o', ...
    'MarkerSize', 14, ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', [1 1 1], ...
    'LineWidth', 2);

% Axis formatting
set(gca, ...
    'FontSize', 30, ...
    'LineWidth', 2, ...
    'FontName', 'Times', ...
    'Layer', 'top', ...
    'Box', 'on');

% xlabel('$\beta_0$', 'Interpreter', 'latex', 'FontSize', 34);
% ylabel('$\omega$',  'Interpreter', 'latex', 'FontSize', 34);

% Colorbar
cb = colorbar;
cb.LineWidth = 1.5;
cb.FontSize  = 26;
cb.FontName  = 'Times';
% ylabel(cb, 'Peak infection prevalence', 'FontSize', 30);

% Fix color scaling
caxis([min(NetItotPeak_beta0_epsilon_vec(:)), ...
    max(NetItotPeak_beta0_epsilon_vec(:)) ]);

axis tight;

%% ============================================================
%  Save figures
% ============================================================

fig_folder = fullfile(pwd, 'NSF_M0_Numerical_Simulation', 'fig');
png_folder = fullfile(pwd, 'NSF_M0_Numerical_Simulation', 'png');

if ~exist(fig_folder, 'dir'); mkdir(fig_folder); end
if ~exist(png_folder, 'dir'); mkdir(png_folder); end

figurename_fig = fullfile(fig_folder, ...
    sprintf(['Figure_NetItotPeak_beta0_epsilon_vec_' ...
             'bline_beta0_%d_Network_%s_test_infected_node_%d_' ...
             'omega_%d_gamma_%d.fig'], ...
             bline_beta0, str, v, omega, gamma));

figurename_png = fullfile(png_folder, ...
    sprintf(['Figure_NetItotPeak_beta0_epsilon_vec_' ...
             'bline_beta0_%d_Network_%s_test_infected_node_%d_' ...
             'omega_%d_gamma_%d.png'], ...
             bline_beta0, str, v, omega, gamma));

saveas(h, figurename_fig);
saveas(h, figurename_png);

disp(['Figure saved as: ' figurename_png]);

%% ============================================================
%  Save peak infection matrix (CSV)
% ============================================================

save_folder = fullfile(pwd, 'NSF_M0_Numerical_Simulation');
if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end

filename_matrix = fullfile(save_folder, ...
    sprintf(['NetItotPeak_beta0_epsilon_vec_' ...
             'Pre_post_lockdown_mobility_' ...
             'bline_beta0_%d_omega_%d_gamma_%d.csv'], ...
             bline_beta0, omega, gamma));

assert(exist('NetItotPeak_beta0_epsilon_vec','var')==1, ...
       'NetItotPeak_beta0_epsilon_vec does not exist');

writematrix(NetItotPeak_beta0_epsilon_vec, filename_matrix);

disp(['CSV saved to: ' filename_matrix]);


total_netw_population = sum(vec_N);
max(max(NetEpiSize_beta0_epsilon_vec./total_netw_population))


hold off 


highlight_points = [
    0.2 0.75
    1 0.65
    ];

%% Fig 4 NetEpiSizefractionbeta0_epsilon_vec 
h = figure(4);
clf;

set(h, 'Color', 'w', 'Position', [100 100 900 700]);
% Main color map
pcolor(bline_beta0_vec, omega_vec, NetEpiSize_beta0_epsilon_vec./total_netw_population)
% pcolor(omega_vec, bline_beta0_vec, NetEpiSize_beta0_epsilon_vec);
shading flat;
hold on;

% Nature-style perceptually uniform colormap

%colormap(turbo);    % Smooth and strong
colormap(parula);   % Smooth but slightly muted
%colormap(jet);      % Vibrant but not perceptually uniform

%colormap(magma);
% colormap(mycmap);

colorbar

% colorbar;

plot(bline_beta0_selected, omega_selected, ...
    '--w', ...          % dashed white line
    'LineWidth', 3);

% Highlight selected parameter points
plot(highlight_points(:,1), highlight_points(:,2), ...
    'o', ...
    'MarkerSize', 14, ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', [1 1 1], ...
    'LineWidth', 2);

% Axis formatting (Nature style)
set(gca, ...
    'FontSize', 30, ...
    'LineWidth', 2, ...
    'FontName', 'Times', ...
    'Layer', 'top', ...
    'Box', 'on');

% xlabel('$\beta_0$', 'Interpreter', 'latex', 'FontSize', 34);
% ylabel('$\omega$',  'Interpreter', 'latex', 'FontSize', 34);

% Colorbar formatting
cb = colorbar;
cb.LineWidth = 1.5;
cb.FontSize  = 26;
cb.FontName  = 'Times';
%ylabel(cb, 'Final epidemic size', 'FontSize', 30);

% Fix color scaling (important for reproducibility)
caxis([ ...
    min(min(NetEpiSize_beta0_epsilon_vec./total_netw_population)), ...
    1]);

axis tight;

%% ============================================================
%  Save figures
% ============================================================

fig_folder = fullfile(pwd, 'NSF_M0_Numerical_Simulation', 'fig');
png_folder = fullfile(pwd, 'NSF_M0_Numerical_Simulation', 'png');

if ~exist(fig_folder, 'dir'); mkdir(fig_folder); end
if ~exist(png_folder, 'dir'); mkdir(png_folder); end

figurename_fig = fullfile(fig_folder, ...
    sprintf(['Figure_NetEpiSizefractionbeta0_epsilon_vec_' ...
             'bline_beta0_%d_Network_%s_test_infected_node_%d_' ...
             'omega_%d_gamma_%d.fig'], ...
             bline_beta0, str, v, omega, gamma));

figurename_png = fullfile(png_folder, ...
    sprintf(['Figure_NetEpiSizefractionbeta0_epsilon_vec_' ...
             'bline_beta0_%d_Network_%s_test_infected_node_%d_' ...
             'omega_%d_gamma_%d.png'], ...
             bline_beta0, str, v, omega, gamma));

saveas(h, figurename_fig);
saveas(h, figurename_png);

disp(['Figure saved as: ' figurename_png]);

%% ============================================================
%  Save epidemic size matrix (CSV)
% ============================================================

save_folder = fullfile(pwd, 'NSF_M0_Numerical_Simulation');
if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end

filename_matrix = fullfile(save_folder, ...
    sprintf(['NetEpiSize_beta0_epsilon_vec_' ...
             'Pre_post_lockdown_mobility_' ...
             'bline_beta0_%d_omega_%d_gamma_%d.csv'], ...
             bline_beta0, omega, gamma));

assert(exist('NetEpiSize_beta0_epsilon_vec','var')==1, ...
       'NetEpiSize_beta0_epsilon_vec does not exist');

writematrix(NetEpiSize_beta0_epsilon_vec, filename_matrix);

disp(['CSV saved to: ' filename_matrix]);


hold off ;

highlight_points = [
    0.2 0.5
    1 0.25
    ];


%% Fig 5: NetItotPeak_fraction_beta0_epsilon_vec
h = figure(5);
clf;
set(h, 'Color', 'w', 'Position', [120 120 900 700]);

% Main heatmap
pcolor(bline_beta0_vec, omega_vec, NetItotPeak_beta0_epsilon_vec./total_netw_population);
% pcolor(omega_vec, bline_beta0_vec, NetItotPeak_beta0_epsilon_vec);

shading flat;
hold on;

% Nature-style colormap for peak intensity
% colormap(magma);
colormap(jet);  % or parula, hot, cool, etc.

% Highlight selected parameter points
plot(highlight_points(:,1), highlight_points(:,2), ...
    'o', ...
    'MarkerSize', 14, ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', [1 1 1], ...
    'LineWidth', 2);

plot(bline_beta0_selected, omega_selected, ...
    '--w', ...          % dashed black line
    'LineWidth', 3);

hold off;

% Axis formatting
set(gca, ...
    'FontSize', 30, ...
    'LineWidth', 2, ...
    'FontName', 'Times', ...
    'Layer', 'top', ...
    'Box', 'on');

% xlabel('$\beta_0$', 'Interpreter', 'latex', 'FontSize', 34);
% ylabel('$\omega$',  'Interpreter', 'latex', 'FontSize', 34);

% Colorbar
cb = colorbar;
cb.LineWidth = 1.5;
cb.FontSize  = 26;
cb.FontName  = 'Times';
% ylabel(cb, 'Peak infection prevalence', 'FontSize', 30);

% Fix color scaling
% caxis([min(min(NetItotPeak_beta0_epsilon_vec./total_netw_population)), ...
%     max(max(NetItotPeak_beta0_epsilon_vec./total_netw_population)) ]);

caxis([0, 0.4]);

axis tight;

%% ============================================================
%  Save figures
% ============================================================

fig_folder = fullfile(pwd, 'NSF_M0_Numerical_Simulation', 'fig');
png_folder = fullfile(pwd, 'NSF_M0_Numerical_Simulation', 'png');

if ~exist(fig_folder, 'dir'); mkdir(fig_folder); end
if ~exist(png_folder, 'dir'); mkdir(png_folder); end

figurename_fig = fullfile(fig_folder, ...
    sprintf(['Figure_NetItotPeak_fraction_beta0_epsilon_vec_' ...
             'bline_beta0_%d_Network_%s_test_infected_node_%d_' ...
             'omega_%d_gamma_%d.fig'], ...
             bline_beta0, str, v, omega, gamma));

figurename_png = fullfile(png_folder, ...
    sprintf(['Figure_NetItotPeak_fraction_beta0_epsilon_vec_' ...
             'bline_beta0_%d_Network_%s_test_infected_node_%d_' ...
             'omega_%d_gamma_%d.png'], ...
             bline_beta0, str, v, omega, gamma));

saveas(h, figurename_fig);
saveas(h, figurename_png);

disp(['Figure saved as: ' figurename_png]);

%% ============================================================
%  Save peak infection matrix (CSV)
% ============================================================

save_folder = fullfile(pwd, 'NSF_M0_Numerical_Simulation');
if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end

filename_matrix = fullfile(save_folder, ...
    sprintf(['NetItotPeak_beta0_epsilon_vec_' ...
             'Pre_post_lockdown_mobility_' ...
             'bline_beta0_%d_omega_%d_gamma_%d.csv'], ...
             bline_beta0, omega, gamma));

assert(exist('NetItotPeak_beta0_epsilon_vec','var')==1, ...
       'NetItotPeak_beta0_epsilon_vec does not exist');

writematrix(NetItotPeak_beta0_epsilon_vec, filename_matrix);

disp(['CSV saved to: ' filename_matrix]);
