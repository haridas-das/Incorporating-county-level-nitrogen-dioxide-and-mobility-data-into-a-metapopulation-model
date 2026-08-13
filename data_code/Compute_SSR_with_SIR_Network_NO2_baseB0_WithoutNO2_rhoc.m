function SSR = Compute_SSR_with_SIR_Network_NO2_baseB0_WithoutNO2_rhoc(P, covid_county_data_interp, t, m, lockdown_week_0_indexed, time_start, time_step_size, avg_pre_lockdown_flux_matrix, avg_post_lockdown_flux_matrix, vec_S0, vec_I0, vec_N,zeta, y0)

    phi_matrix_func = @(t) get_matrix_phi(t, lockdown_week_0_indexed, time_start, time_step_size, avg_pre_lockdown_flux_matrix, avg_post_lockdown_flux_matrix);

    % Extract parameters for each node
    % P_init          P = [bline_beta0, gamma, rhoc];
    bline_beta0          = P(1);
    gamma                = P(2);
    baseline_rho_vec     = P(3)*ones(m,1);                                 % baseline_rho
    % Solve the SIR model with time-dependent infection rates
    [~, y_k] = ode45(@(t, y) SIR_Network_WithoutNO2_baseB0(t, m, phi_matrix_func(t), y, bline_beta0, vec_S0, vec_I0, vec_N, gamma, zeta), t, y0);
    
    % Extract infected populations from ODE solution
    simulated_data = y_k(:, m+1:2*m);                                      % Infected data for all nodes
    rho_epsilon_vec = baseline_rho_vec;

    % Ensure rho_vec is a row vector
    rho_vec_row = rho_epsilon_vec(:)'; % Convert rho_vec to row vector if not already

    % Multiply each column of simulated_data by the corresponding scalar in rho_vec
    rho_vec_times_simulated_data = simulated_data.* rho_vec_row;
    
    % Compute residuals
    residuals = covid_county_data_interp - rho_vec_times_simulated_data;
     
    % residuals(:): This converts the residuals matrix (size 531×8) into a column vector (size 4248×1) by stacking its columns vertically.
    % Compute SSR
    % SSR = sum(residuals(:).^2);
    % Ignore NaNs counties in residuals
    valid_idx = ~isnan(residuals);          % logical array of valid entries
    SSR = sum(residuals(valid_idx).^2);     % sum only over valid entries
    
    % Optionally, display
    disp('Sum of Squared Residuals (SSR):');
    disp(SSR);
end