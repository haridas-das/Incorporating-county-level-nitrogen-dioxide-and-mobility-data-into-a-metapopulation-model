function sol          = SIRsol(par_beta, par_rho, gamma, IC, t)
% Return [susceptible prevalence recovered incidence] of the SIR model 
         %disp(num2str(par));                   
         DeHandle     = @(t,y) SIRmodel(t, y, par_beta, gamma);
         [~, Y]       = ode45(DeHandle, t, IC);      % returns the sol in column vector and the time in the row vector
       %  sol          = [0 par_rho 0 0]*Y';         % Return transformed solutio
         sol = sum([0 par_rho 0 0] * Y', 1);         % Sum along the first dimension
end


