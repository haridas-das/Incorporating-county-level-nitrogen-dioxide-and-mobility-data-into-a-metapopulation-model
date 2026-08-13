clear all; clc; close all; 
%% This is FullTimeSeries_70 days fitting Booth strapping fitting parameters
%% Load error metrices data
SIRmf03_vs_SIRnetmf03_error_metrices_filename = 'Mf0_SIRmf03_vs_SIRnetmf03_PFITboothstrapPIFullTimeSeriesFittingErrorMetricesR2.csv';
SIRmf03_vs_SIRnetmf03_error_metrices = readtable(SIRmf03_vs_SIRnetmf03_error_metrices_filename);

SIRmf03_vs_SIRnetmf03_error_metrices.Properties.VariableNames

%% Directory to save figures
directory = 'frac-cum y_k data/Mf0_Full_TS_Error_metrices_with_R2/result figures/png';
if ~exist(directory, 'dir')
    mkdir(directory);  % Create folder if it doesn't exist
end

%% Define colors
color_sir  = [1 0 0];      % red 
color_netw = [0 0 1];      % blue

%% List of metrics exactly as in the table
metrics = {'SSR','mse','mae','rmse','nrmse','r2','mape','pae','pa','pe'};

baseline_prefix = 'sir_baseline_';
net_prefix      = 'sir_net_fit_';

%% Loop through metrics and plot
for i = 1:length(metrics)
    metric = metrics{i};
    
    baseline_var = strcat(baseline_prefix, metric);  % keep case
    net_var      = strcat(net_prefix, metric);       % keep case
    
    y_baseline = SIRmf03_vs_SIRnetmf03_error_metrices.(baseline_var);
    y_net      = SIRmf03_vs_SIRnetmf03_error_metrices.(net_var);
    labels     = SIRmf03_vs_SIRnetmf03_error_metrices.o_name;
    
    % Create figure
    figure('Visible','off');
    hold on;
    
    % Plot each bar individually with defined colors
    x = 1:numel(labels);
    bar(x-0.15, y_baseline, 1, 'FaceColor', color_sir);   % SIR
    bar(x+0.15, y_net, 1, 'FaceColor', color_netw);       % SIRNet
    % set(gca, 'XTick', x, 'XTickLabel', labels);
    xtickangle(45);
    %ylabel(metric);
   % title(['Comparison of ', metric, ' between SIR and SIRNet']);
    %legend('SIR', 'SIRNet', 'Location','best');
    % grid on;
    box on 
    hold off;
    set(gca, 'LineWidth', 1.5, 'FontSize', 30, 'FontName', 'Times');

    % Save figure
    save_filename = fullfile(directory, ['Error_', metric, '_SIR_vs_SIRNet.png']);
    saveas(gcf, save_filename);
    close(gcf);
end

disp('All error metric plots with colors saved successfully!');

%% Loop through metrics and plot boxplots for SIR vs SIRNet
for i = 1:length(metrics)
    metric = metrics{i};
    
    baseline_var = strcat(baseline_prefix, metric);  % keep case
    net_var      = strcat(net_prefix, metric);       % keep case
    
    y_baseline = SIRmf03_vs_SIRnetmf03_error_metrices.(baseline_var);
    y_net      = SIRmf03_vs_SIRnetmf03_error_metrices.(net_var);
    
    % Combine data into one column
    all_data  = [y_baseline; y_net];
    % Create grouping variable (1=SIR, 2=SIRNet)
    all_group = [ones(size(y_baseline)); 2*ones(size(y_net))];
    
    % Create figure
    figure('Visible','off');
    
    % Boxplot with colors
    h = boxplot(all_data, all_group, 'Whisker', 1.5);
    % h = boxplot(all_data, all_group, 'Labels', {'SIR','SIRNet'}, 'Whisker', 1.5);
    % Set colors for boxes
    set(h, {'linew'}, {2});  % make lines thicker
    boxColors = [color_sir; color_netw];
    edgeColors =[color_sir; color_netw];% [0 0 0; 0 0 0] %[0 1 0; 1 0 1];               % Example: black edges; you can change

    % Change the box colors
    boxes = findobj(gca, 'Tag', 'Box');
    for j = 1:length(boxes)
        % patch(get(boxes(j), 'XData'), get(boxes(j), 'YData'), boxColors(j,:), 'FaceAlpha',0.5);
        patch(get(boxes(j), 'XData'), get(boxes(j), 'YData'), boxColors(j,:), ...
              'FaceAlpha',1, 'EdgeColor', edgeColors(j,:), 'LineWidth', 1.5);
    end
    set(gca, 'LineWidth', 1.5, 'FontSize', 30, 'FontName', 'Times');
     switch metric
        case {'r2'}
            ylim([0 1])
    
       case {'pae'}
            ylim([0 1.2])

       case {'nrmse'}
            ylim([0 0.3])

        case {'SSR','mse','rmse'}
            ylim([0 max(all_data)*1.1])
    
        case {'mae','mape','pa','pe'}
            ylim([0 max(all_data)*1.2])
    
        otherwise
            padding = 0.05 * (max(all_data) - min(all_data));
            ylim([min(all_data)-padding, max(all_data)+padding])
    end

    % ylabel(metric);
    % title(['Comparison of ', metric, ' between SIR and SIRNet']);
    % grid on;
    
    % Save figure
    save_filename = fullfile(directory, ['Boxplot_Error_', metric, '_SIR_vs_SIRNet.png']);
    saveas(gcf, save_filename);
    close(gcf);
end

%disp('All error metric box plots saved successfully!');


%% 

%% Loop through metrics and plot histograms for SIR vs SIRNet

for i = 1:length(metrics)
    metric = metrics{i};
    
    baseline_var = strcat(baseline_prefix, metric);
    net_var      = strcat(net_prefix, metric);
    
    y_baseline = SIRmf03_vs_SIRnetmf03_error_metrices.(baseline_var);
    y_net      = SIRmf03_vs_SIRnetmf03_error_metrices.(net_var);
    
    figure('Visible','off'); hold on;

    % Histograms (PDF normalized)
    histogram(y_baseline, 'Normalization','pdf', ...
        'FaceColor', color_sir, 'FaceAlpha', 1, 'EdgeColor','none');
    histogram(y_net, 'Normalization','pdf', ...
        'FaceColor', color_netw, 'FaceAlpha', 1, 'EdgeColor','none');

    box on
    set(gca, 'LineWidth', 1.5, 'FontSize', 30, 'FontName', 'Times');

    % Axis limits (same logic as boxplots)
    all_data = [y_baseline; y_net];

    switch metric
        case {'r2'}
            xlim([0 1])

        case {'pae'}
            xlim([0 1.2])

        case {'nrmse'}
            xlim([0 0.3])

        case {'SSR','mse','rmse'}
            xlim([0 max(all_data)*1.1])

        case {'mae','mape','pa','pe'}
            xlim([0 max(all_data)*1.2])

        otherwise
            padding = 0.05 * (max(all_data) - min(all_data));
            xlim([min(all_data)-padding, max(all_data)+padding])
    end

    % ylabel('Probability Density')
    % xlabel(metric)

    % legend({'SIR','SIRNet'}, 'Location','best')  % optional

    % Save figure
    save_filename = fullfile(directory, ...
        ['Hist_Error_', metric, '_SIR_vs_SIRNet.png']);
    saveas(gcf, save_filename);
    close(gcf);
end



%% ============================================================
% Summary Statistics for SIR vs SIRNet Error Metrics
%% ============================================================

metrics = {'SSR','mse','mae','rmse','nrmse','r2','mape','pae','pa','pe'};

baseline_prefix = 'sir_baseline_';
net_prefix      = 'sir_net_fit_';

summary_results = table();

for i = 1:length(metrics)

    metric = metrics{i};

    baseline_var = [baseline_prefix metric];
    net_var      = [net_prefix metric];

    % Extract data
    y_baseline = SIRmf03_vs_SIRnetmf03_error_metrices.(baseline_var);
    y_net      = SIRmf03_vs_SIRnetmf03_error_metrices.(net_var);

    % Remove Inf and NaN
    y_baseline = y_baseline(isfinite(y_baseline));
    y_net      = y_net(isfinite(y_net));

    %% Baseline statistics
    baseline_mean   = mean(y_baseline);
    baseline_median = median(y_baseline);
    baseline_std    = std(y_baseline);
    baseline_min    = min(y_baseline);
    baseline_max    = max(y_baseline);
    baseline_iqr    = iqr(y_baseline);

    %% SIRNet statistics
    net_mean   = mean(y_net);
    net_median = median(y_net);
    net_std    = std(y_net);
    net_min    = min(y_net);
    net_max    = max(y_net);
    net_iqr    = iqr(y_net);

    %% Create temporary table
    temp_table = table( ...
        {metric}, ...
        baseline_mean, baseline_median, baseline_std, ...
        baseline_min, baseline_max, baseline_iqr, ...
        net_mean, net_median, net_std, ...
        net_min, net_max, net_iqr, ...
        'VariableNames', ...
        {'Metric', ...
        'SIR_Mean', 'SIR_Median', 'SIR_STD', ...
        'SIR_Min', 'SIR_Max', 'SIR_IQR', ...
        'SIRNet_Mean', 'SIRNet_Median', 'SIRNet_STD', ...
        'SIRNet_Min', 'SIRNet_Max', 'SIRNet_IQR'} );

    summary_results = [summary_results; temp_table];

end

%% Display
disp(summary_results)

%% Save
summary_filename = fullfile(directory, ...
    'Summary_Statistics_SIR_vs_SIRNet.csv');

writetable(summary_results, summary_filename);

disp('Summary statistics saved successfully!');