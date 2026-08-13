function f = SIRSmodel(t, y, par_beta_xi, gamma)
   % Label the parameters and variables 
%    beta       = par(1);
% %    gamma      = par(2);
%    xi         = par(2);
   beta       = par_beta_xi(1);
%    gamma      = par(2);
   xi         = par_beta_xi(2);
   S          = y(1);
   I          = y(2);
   R          = y(3);
   %C          = y(4);
   % Conservation law
   N          = S + I + R;
   % Input the differential equations 
   dS         = -beta*S*I/N + xi*R;              % susceptible
   dI         =  beta*S*I/N - gamma*I;            % prevalence
   dR         =  gamma*I-xi*R;                   % recovered individuals
   dis_inc    =  beta*S*I/N;                      % cumulative incidence computation
   f          = [dS dI dR dis_inc]';              % return the [S I R C]
end
