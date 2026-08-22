clear all; clc; close all; 
%% BS forecasting performance across different window sizes, evaluated using nRMSE and PAE (BS)
%% Load error metrices data
nRMSE_filename = 'Mf0_sir_sirnetw_diff_windows_forecasting_mobility_vs_sirnetw_nRMSE_data.csv';
nRMSE_diff_win_BS_fcst = readtable(nRMSE_filename);

nRMSE_diff_win_BS_fcst.Properties.VariableNames

%% Directory
directory = 'frac-cum y_k data/Mf0_diff_windows_mobility_vs_NO2sirnetw_nRMSE_PAE/result figures/png';
if ~exist(directory,'dir')
    mkdir(directory);
end

%% Colors
color_sir  = [1 0 0];
color_netw = [0 0 1];

%% Windows
windows = {'early','PrePeak','Peak','PostPeak','Decline'}; % ,'Late'


%% Figure A — histogram (distribution only)


for i = 1:length(windows)
    win = windows{i};

    sir_var     = ['netm_fcst_nRMSE_no_NO2_' win];
    sirnetw_var = ['sirnetw_nRMSE_' win];

    y_sir  = nRMSE_diff_win_BS_fcst.(sir_var);
    y_netw = nRMSE_diff_win_BS_fcst.(sirnetw_var);

    figure('Visible','off'); hold on;

    histogram(y_sir,  'Normalization','pdf', ...
        'FaceColor',color_sir,  'FaceAlpha',1, 'EdgeColor','none');
    histogram(y_netw, 'Normalization','pdf', ...
        'FaceColor',color_netw, 'FaceAlpha',1, 'EdgeColor','none');

    box on
    set(gca,'LineWidth',1.5,'FontSize',30,'FontName','Times')
    % title(['nRMSE distribution (', win, ')'])
    % legend({'SIR','SIRNet'},'Location','best')

    xlim([0 max([y_sir; y_netw])*1.05])

    saveas(gcf, fullfile(directory, ...
        ['Hist_nRMSE_', win, '_SIR_vs_SIRNet.png']));
    close(gcf);
end

%% Option B: True histogram / mirrored density (recommended for paper)

for i = 1:length(windows)
    win = windows{i};

    sir_var     = ['netm_fcst_nRMSE_no_NO2_' win];
    sirnetw_var = ['sirnetw_nRMSE_' win];

    y_sir  = nRMSE_diff_win_BS_fcst.(sir_var);
    y_netw = nRMSE_diff_win_BS_fcst.(sirnetw_var);

    figure('Visible','off'); hold on;

    [f1,x1] = ksdensity(y_sir);
    [f2,x2] = ksdensity(y_netw);

    area(x1,  f1,  'FaceColor',color_sir,  'FaceAlpha',1,'EdgeColor','none');
    area(x2, -f2, 'FaceColor',color_netw, 'FaceAlpha',1,'EdgeColor','none');

    yline(0,'k','LineWidth',1)

    box on
    set(gca,'LineWidth',1.5,'FontSize',30,'FontName','Times')
    title(['nRMSE density (', win, ')'])
    legend({'SIR','SIRNet'},'Location','best')

    % xlim([0 max([y_sir; y_netw])*1.05])

    saveas(gcf, fullfile(directory, ...
        ['JHist_nRMSE_', win, '_SIR_vs_SIRNet.png']));
    close(gcf);
end


%% Figure B — Boxplots only

for i = 1:length(windows)
    win = windows{i};

    sir_var     = ['netm_fcst_nRMSE_no_NO2_' win];
    sirnetw_var = ['sirnetw_nRMSE_' win];

    y_sir  = nRMSE_diff_win_BS_fcst.(sir_var);
    y_netw = nRMSE_diff_win_BS_fcst.(sirnetw_var);

    all_data  = [y_sir; y_netw];
    all_group = [ones(size(y_sir)); 2*ones(size(y_netw))];

    figure('Visible','off');

    % Create boxplot
    h = boxplot(all_data, all_group, 'Whisker',1.5);
    set(h,{'LineWidth'},{2});

    % Colors: first box = SIR, second box = SIRNet
    boxColors = [color_sir; color_netw];

    % Get box objects and sort left-to-right
    boxes = findobj(gca,'Tag','Box');
    boxes = sortrows([boxes(:).XData].');   % sort by X position
    boxes = flipud(boxes);                  % MATLAB returns top-to-bottom, so flip

    % Apply colors
    boxes_handles = findobj(gca,'Tag','Box');
    boxes_handles = flipud(boxes_handles);  % ensure left-to-right
    for j = 1:length(boxes_handles)
        patch(get(boxes_handles(j),'XData'), get(boxes_handles(j),'YData'), ...
              boxColors(j,:), 'FaceAlpha',1, ...
              'EdgeColor',boxColors(j,:), 'LineWidth',1.5);
    end

    % ylim([0 30])
    box on
    set(gca,'LineWidth',1.5,'FontSize',30,'FontName','Times')

    saveas(gcf, fullfile(directory, ...
        ['Boxplot_nRMSE_', win, '_SIR_vs_SIRNet.png']));
    close(gcf);
end






%% Load error metrices data
PAE_filename = 'Mf0_sir_sirnetw_diff_windows_forecasting_mobility_vs_sirnetw_PAE_data.csv';
PAE_diff_win_BS_fcst = readtable(PAE_filename);

PAE_diff_win_BS_fcst.Properties.VariableNames


%% Figure A — Histogram (PAE distribution only)

for i = 1:length(windows)
    win = windows{i};

    sir_var     = ['netm_fcst_pae_no_NO2_' win];
    sirnetw_var = ['sirnetw_PAE_' win];

    y_sir  = PAE_diff_win_BS_fcst.(sir_var);
    y_netw = PAE_diff_win_BS_fcst.(sirnetw_var);

    figure('Visible','off'); hold on;

    histogram(y_sir,  'Normalization','pdf', ...
        'FaceColor',color_sir,  'FaceAlpha',1, 'EdgeColor','none');
    histogram(y_netw, 'Normalization','pdf', ...
        'FaceColor',color_netw, 'FaceAlpha',1, 'EdgeColor','none');

    box on
    set(gca,'LineWidth',1.5,'FontSize',30,'FontName','Times')

    xlim([0 max([y_sir; y_netw])*1.05])

    saveas(gcf, fullfile(directory, ...
        ['Hist_PAE_', win, '_SIR_vs_SIRNet.png']));
    close(gcf);
end



%% Figure A (Option B) — True histogram for PAE

for i = 1:length(windows)
    win = windows{i};

    sir_var     = ['netm_fcst_pae_no_NO2_' win];
    sirnetw_var = ['sirnetw_PAE_' win];

    y_sir  = PAE_diff_win_BS_fcst.(sir_var);
    y_netw = PAE_diff_win_BS_fcst.(sirnetw_var);

    figure('Visible','off'); hold on;

    [f1,x1] = ksdensity(y_sir);
    [f2,x2] = ksdensity(y_netw);

    area(x1,  f1,  'FaceColor',color_sir,  ...
        'FaceAlpha',1,'EdgeColor','none');
    area(x2, -f2, 'FaceColor',color_netw, ...
        'FaceAlpha',1,'EdgeColor','none');

    yline(0,'k','LineWidth',1)

    box on
    set(gca,'LineWidth',1.5,'FontSize',30,'FontName','Times')

    xlim([0 max([y_sir; y_netw])*1.05])

    saveas(gcf, fullfile(directory, ...
        ['JHist_PAE_', win, '_SIR_vs_SIRNet.png']));
    close(gcf);
end


%% Figure B — Boxplots (PAE)
%% Figure B — Boxplots (PAE) corrected

for i = 1:length(windows)
    win = windows{i};

    sir_var     = ['netm_fcst_pae_no_NO2_' win];
    sirnetw_var = ['sirnetw_PAE_' win];

    y_sir  = PAE_diff_win_BS_fcst.(sir_var);
    y_netw = PAE_diff_win_BS_fcst.(sirnetw_var);

    all_data  = [y_sir; y_netw];
    all_group = [ones(size(y_sir)); 2*ones(size(y_netw))];

    figure('Visible','off');

    h = boxplot(all_data, all_group, 'Whisker',1.5);
    set(h,{'LineWidth'},{2})

    % Get boxes in left-to-right order
    boxes = flipud(findobj(gca,'Tag','Box'));  

    boxColors = [color_sir; color_netw];  % first box = SIR, second = SIRNet

    for j = 1:length(boxes)
        patch(get(boxes(j),'XData'), get(boxes(j),'YData'), ...
              boxColors(j,:), 'FaceAlpha',1, ...
              'EdgeColor',boxColors(j,:), 'LineWidth',1.5);
    end

    %ylim([0, 15])
    box on
    set(gca,'LineWidth',1.5,'FontSize',30,'FontName','Times')

    saveas(gcf, fullfile(directory, ...
        ['Boxplot_PAE_', win, '_SIR_vs_SIRNet.png']));
    close(gcf);
end



%% =========================================================
% Quantitative summary statistics for each epidemic window
%% =========================================================

summary_results = table();

for i = 1:length(windows)

    win = windows{i};

    %% ---------- nRMSE ----------
    sir_var     = ['netm_fcst_nRMSE_no_NO2_' win];
    sirnetw_var = ['sirnetw_nRMSE_' win];

    y_sir_nrmse = nRMSE_diff_win_BS_fcst.(sir_var);
    y_net_nrmse = nRMSE_diff_win_BS_fcst.(sirnetw_var);

    % Remove Inf and NaN
    y_sir_nrmse = y_sir_nrmse(isfinite(y_sir_nrmse));
    y_net_nrmse = y_net_nrmse(isfinite(y_net_nrmse));

    %% ---------- PAE ----------
    sir_var_pae     = ['netm_fcst_pae_no_NO2_' win];
    sirnetw_var_pae = ['sirnetw_PAE_' win];

    y_sir_pae = PAE_diff_win_BS_fcst.(sir_var_pae);
    y_net_pae = PAE_diff_win_BS_fcst.(sirnetw_var_pae);

    % Remove Inf and NaN
    y_sir_pae = y_sir_pae(isfinite(y_sir_pae));
    y_net_pae = y_net_pae(isfinite(y_net_pae));

    %% ---------- Create rows ----------

    temp = table();

    temp.Window = {win};

    %% nRMSE statistics
    temp.SIR_nRMSE_mean   = mean(y_sir_nrmse);
    temp.SIR_nRMSE_median = median(y_sir_nrmse);
    temp.SIR_nRMSE_std    = std(y_sir_nrmse);

    temp.NET_nRMSE_mean   = mean(y_net_nrmse);
    temp.NET_nRMSE_median = median(y_net_nrmse);
    temp.NET_nRMSE_std    = std(y_net_nrmse);

    %% PAE statistics
    temp.SIR_PAE_mean   = mean(y_sir_pae);
    temp.SIR_PAE_median = median(y_sir_pae);
    temp.SIR_PAE_std    = std(y_sir_pae);

    temp.NET_PAE_mean   = mean(y_net_pae);
    temp.NET_PAE_median = median(y_net_pae);
    temp.NET_PAE_std    = std(y_net_pae);

    %% Append
    summary_results = [summary_results; temp];

end

%% Display
disp(summary_results)

%% Save
writetable(summary_results, ...
    fullfile(directory,'Summary_Statistics_nRMSE_PAE.csv'));