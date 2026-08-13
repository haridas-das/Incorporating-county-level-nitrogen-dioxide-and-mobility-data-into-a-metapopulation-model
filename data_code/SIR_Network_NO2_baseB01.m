function dydt = SIR_Network_NO2_baseB01(t, m, matrix_phi, y, NO2_data, bline_beta0, bline_beta1, top_five_pop_flow_nodes, omega_0, omega_1, vec_S0, vec_I0, vec_N, gamma, zeta)

% Evaluate time-dependent phi matrix
% matrix_phi = phi_matrix_func(t);

% Compute NO2-driven infection rates when NO2 data is a fixed m × 1 vector
%vec_beta = bline_beta0 + omega*NO2_data;  % vector of size m × 1

% % Initialize vec_beta
% vec_beta = zeros(m, 1);
% 
% % Complement nodes
% all_nodes = 1:m;
% complement_nodes = setdiff(all_nodes, top_five_pop_flow_nodes);
% 
% % Assign values
% vec_beta(top_five_pop_flow_nodes) = bline_beta0 + omega * NO2_data(top_five_pop_flow_nodes);                  % top five nodes
% vec_beta(complement_nodes) = bline_beta1 + omega * NO2_data(complement_nodes);  % rest


% m = total number of nodes
% top_five_pop_flow_nodes = indices of top five nodes
% bline_beta0 = baseline for top five nodes (scalar or vector of length 5)
% bline_beta1 = baseline for other nodes (scalar or vector of length m-5)
% NO2_data = vector of length m
% omega = scalar

% % Initialize vec_beta
% vec_beta = zeros(m, 1);
% 
% % All node indices
% all_nodes = 1:m;
% 
% % Complement nodes
% complement_nodes = setdiff(all_nodes, top_five_pop_flow_nodes);
% 
% % Assign values safely
% vec_beta(top_five_pop_flow_nodes) = bline_beta0 + omega_0 * NO2_data(top_five_pop_flow_nodes);
% vec_beta(complement_nodes) = bline_beta1 + omega_1.* NO2_data(complement_nodes);

% NO2_data = NO2_data(:);   % make it column vector, 77×1
% vec_beta = zeros(m,1);    % also column vector
% top_five_pop_flow_nodes = top_five_pop_flow_nodes(:);  % 5×1 column
% % all_nodes = 1:m;
% %vec_beta(top_five_pop_flow_nodes) = bline_beta0 + omega_0* NO2_data(top_five_pop_flow_nodes);
% %complement_nodes = setdiff(all_nodes, top_five_pop_flow_nodes);
% complement_nodes = setdiff((1:m)', top_five_pop_flow_nodes);  % column (m-5)×1
% 
% %vec_beta(complement_nodes)        = bline_beta1 + omega_1* NO2_data(complement_nodes);
% 
% 
% vec_beta(top_five_pop_flow_nodes) = bline_beta0 + omega_0 .* NO2_data(top_five_pop_flow_nodes);
% vec_beta(complement_nodes)        = bline_beta1 + omega_1 .* NO2_data(complement_nodes);


vec_beta = zeros(m,1);                   % vector of size m × 1
top_five_pop_flow_nodes = top_five_pop_flow_nodes(:); % column
complement_nodes = setdiff((1:m)', top_five_pop_flow_nodes); % column

vec_beta(top_five_pop_flow_nodes) = bline_beta0 + omega_0 * NO2_data(top_five_pop_flow_nodes);
vec_beta(complement_nodes)        = bline_beta1 + omega_1 * NO2_data(complement_nodes);



%% here matrix_phi could be replaced bty phi_matrix_func(t)
%% Initialize time-dependent beta vector
% %vec_beta                       = zeros(m, 1);
% %vec_beta = ComputeGaussianTimeDependentBeta(t, eta_vec, sigma_vec, mu_vec, f_vec, A, T, tau, m)
% vec_beta = zeros(m, 1); % m rows (nodes), length(t) columns (time steps)
%
% % Loop over each node
% for i = 1:m
% %     % Compute the first term
% %     term1 = (eta_vec(i) / (sigma_vec(i) * sqrt(2 * pi))) * ...
% %             exp(-((t - mu_vec(i)).^2) / (2 * sigma_vec(i)^2));
% %
% %     % Add the second term
% %     term2 = f_vec(i);
% %
% %     % Combine with the sinusoidal modulation
% % %     vec_beta(i, :) = (term1 + term2) .* (1 + A * sin(2 * pi * (t / T + tau)));
% %  %    vec_beta(i, 1) = (term1 + term2).*(1 + A * sin(2 * pi * (t / T + tau)));
% %  vec_beta(i, 1) = (term1 + term2);
% vec_beta(i,1) = f_vec(i);
% end

%%Start: First define the present population at node i, which is also constant in time come from rest of the nodes
Np = zeros(m,1); % initial present population at i
Npsum = 0;
for i= 1:m
    for k= 1:m
        Npsum = Npsum+matrix_phi(k,i).*vec_N(k); % phi(k,i)*vec_N(k,1)is the total population at vertex k that travels every day to vertex i
    end
    Np(i)= Npsum;
    %break;
    Npsum = 0;
end
%% end to define the present population at node i.

% initialize the set of equations for susceptible and infected and recovered ODEs at times t
dydt = zeros(3*m,1); % this is the susceptible and infected differential equations
%% Constract the susceptibl class differential equations dSdt
dsdtsum=0;
for i=1:m

    for j=1:m
        for k=1:m
            dsdtsum=dsdtsum-vec_beta(j)*matrix_phi(i,j)*y(i)*matrix_phi(k,j)*y(m+k)/Np(j);
        end
    end
    dydt(i)=dsdtsum;
    dsdtsum=0;
end


%% Constract the infected class differential equations dIdt
didtsum=0;
%didt = zeros(2*m,1); % this is the susceptible
for i=m+1:2*m
    for j=1:m
        for k=1:m
            didtsum=didtsum+vec_beta(j)*matrix_phi(i-m,j)*y(i-m)*matrix_phi(k,j)*y(m+k)/Np(j);
        end
    end
    dydt(i)=didtsum-gamma*y(i);
    didtsum=0;
end

for i = (2*m+1):3*m
    dydt(i) = gamma*y(i-m);
end


end





