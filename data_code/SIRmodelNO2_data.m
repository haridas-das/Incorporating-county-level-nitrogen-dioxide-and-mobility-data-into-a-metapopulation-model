function f = SIRmodelNO2_data(t, y, bline_beta0, omega, NO2_data, gamma)
   % Label the parameters and variables 
%    beta       = par(1);
% %    gamma      = par(2);
%    xi         = par(2);
   beta = bline_beta0 + omega*NO2_data;

   %beta       = par_beta_xi(1);
   %gamma      = par(2);
   % xi         = par_beta_xi(2);
   S          = y(1);
   I          = y(2);
   R          = y(3);
   %C          = y(4);
   % Conservation law
   N          = S + I + R;
   % Input the differential equations 
   dS         = -beta*S*I/N ;              % susceptible
   dI         =  beta*S*I/N - gamma*I;            % prevalence
   dR         =  gamma*I;                   % recovered individuals
   dis_inc    =  beta*S*I/N;                      % cumulative incidence computation
   f          = [dS dI dR dis_inc]';              % return the [S I R C]
end
