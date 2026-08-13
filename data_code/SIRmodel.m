function f = SIRmodel(t, y, par_beta, gamma)
   % Label the parameters and variables 
   beta       = par_beta;
   S          = y(1);
   I          = y(2);
   R          = y(3);
   %C          = y(4);
   % Conservation law
   N          = S + I + R;
   % Input the differential equations 
   dS         = -beta*S*I/N;                      % susceptible
   dI         =  beta*S*I/N - gamma*I;            % prevalence
   dR         =  gamma*I;                          % recovered individuals
   dis_inc    =  beta*S*I/N;                      % cumulative incidence computation
   f          = [dS dI dR dis_inc]';              % return the [S I R C]
end
