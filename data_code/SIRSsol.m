function sol          = SIRSsol(par_beta_xi, par_rho, gamma, IC,t)
% Return [susceptible prevalence recovered incidence] of the SIR model 
         %disp(num2str(par));                   
         DeHandle     = @(t,y) SIRSmodel(t, y, par_beta_xi, gamma);
         [~, Y]       = ode45(DeHandle, t, IC);      % returns the sol in column vector and the time in the row vector
         % sol is the solution for 4 classes along the row vectors 
%          sol          = Y';                          % S, I, R, C are in the row vectors
%          solo         = Y;                           % S, I, R, C's are in the colum vector
%          solution     = solo(:,4);                   % cumulative incidence
%          % size of the solution  
%          st = size(solution);
%          % initialize the incidence by zero vector
%          dis_inc = zeros(st(1),1);    
%          % compute the incidence
%          for q = 1: (st(1)-1)
% %        dis_inc(q+1,1) = int64(solution(q+1)) - int64(solution(q));
% %        dis_inc(q+1,1) = int32(solution(q+1)) - int32(solution(q));
%          dis_inc(q+1,1) = solution(q+1) - solution(q);
%          end
%          % Update the first incidence from the cumulative incidence vector 
%          dis_inc(1,:) = solution(1);
%          % update the cum sol row vec by the incid sol row vector 
%          sol(4,:) = dis_inc';
%          
        new_sol = sum([0 par_rho 0 0] * Y', 1); % Sum along the first dimension
        sol     = new_sol';
end

% checking solution
%t   =   [1  2]
% t   =   tdata
% par_beta_xi = [par(1), par(2)]
% par_rho     = par(3)
% gamma       = par(4); 
% IC  = [3928051        7109       24195        7109]

