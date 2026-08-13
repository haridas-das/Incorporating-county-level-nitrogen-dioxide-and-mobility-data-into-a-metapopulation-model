function matrix_phi = get_matrix_phi(t,lockdown_week_0_indexed ,time_start,time_step_size,avg_pre_lockdown_flux_matrix, avg_post_lockdown_flux_matrix)
    if t <= floor((lockdown_week_0_indexed  - time_start) / time_step_size) + 1;
        matrix_phi = avg_pre_lockdown_flux_matrix;
    else
        matrix_phi = avg_post_lockdown_flux_matrix;  % Change as needed
    end
end

