function dydt =        TDSirGenNetworkModel(t,y,m,Np, matrix_phi,vec_beta, vec_S0, vec_I0, vec_N, gamma, zeta,time_start,time_step_size)

% initialize the set of equations for susceptible and infected and recovered ODEs at times t
dydt = zeros(3 * m, 1);

% define the set of ODEs (dSdt) at times t for susceptible class
dsdtsum = 0;
for i = 1:m
    for j = 1:m
        for k = 1:m
            % Evaluate ODE at times t
            dsdtsum = dsdtsum - vec_beta(j) * matrix_phi(i, j, floor((t - time_start) / time_step_size) + 1) * y(i) * matrix_phi(k, j, floor((t - time_start) / time_step_size) + 1) * y(m + k) / Np(j,floor((t - time_start) / time_step_size) + 1);
            % dsdtsum = dsdtsum - vec_beta(j) * matrix_phi(i, j, ceil((t-time_start)/(time_step_size) +1)) * y(i) * matrix_phi(k, j, ceil((t-time_start)/(time_step_size) +1)) * y(m + k) / Np(j,ceil((t-time_start)/(time_step_size) +1));
            % dsdtsum = dsdtsum - vec_beta(j) * matrix_phi(i, j, int8((t-time_start)/(time_step_size) +1)) * y(i) * matrix_phi(k, j, int8((t-time_start)/(time_step_size) +1)) * y(m + k) / Np(j,int8((t-time_start)/(time_step_size) +1));
            %dsdtsum = dsdtsum - vec_beta(j)*matrix_phi(i, j, t) * y(i)*matrix_phi(k, j, t)*y(m + k)/Np(j);
        end
    end
    dydt(i) = dsdtsum;
    dsdtsum = 0;
end

%% define the set of ODEs (dIdt) at times t for infected class

didtsum = 0;
for i = (m+1):2*m
    for j = 1:m
        for k = 1:m
            % Evaluate ODE at times t
            didtsum = didtsum + vec_beta(j) * matrix_phi(i - m, j, floor((t - time_start) / time_step_size) + 1) * y(i - m) * matrix_phi(k, j, floor((t - time_start) / time_step_size) + 1) * y(m + k) / Np(j,floor((t - time_start) / time_step_size) + 1);
            %didtsum = didtsum + vec_beta(j) * matrix_phi(i - m, j, ceil((t-time_start)/(time_step_size) +1)) * y(i - m) * matrix_phi(k, j, ceil((t-time_start)/(time_step_size) +1)) * y(m + k) / Np(j,ceil((t-time_start)/(time_step_size) +1));
            %didtsum = didtsum + vec_beta(j) * matrix_phi(i - m, j, int8((t-time_start)/(time_step_size) +1)) * y(i - m) * matrix_phi(k, j, int8((t-time_start)/(time_step_size) +1)) * y(m + k) / Np(j,int8((t-time_start)/(time_step_size) +1));
            %didtsum = didtsum + vec_beta(j)*matrix_phi(i-m, j, t) * y(i-m)*matrix_phi(k, j, t) * y(m + k) / Np(j);
        end
    end
    dydt(i) = didtsum - gamma*y(i);
    didtsum = 0;
end

%% define the set of ODEs (dIdt) at times t for recovered class
%drdtsum = 0;
for i = (2*m+1):3*m
    dydt(i) = gamma*y(i-m);
end

