function dydt = SIR_Network_WithoutNO2_baseB0(t, m, matrix_phi, y, bline_beta0, vec_S0, vec_I0, vec_N, gamma, zeta)

% Compute NO2-driven infection rates wit mean NO2 data is a fixed m × 1 vector                                  % vector of size m × 1
vec_beta = bline_beta0*ones(m,1);

%%Start: First define the present population at node i, which is also constant in time come from rest of the nodes
Np=zeros(m,1); % initial present population at i

Npsum=0;
for i=1:m
    for k=1:m
        Npsum=Npsum+matrix_phi(k,i).*vec_N(k);                             % phi(k,i)*vec_N(k,1)is the total population at vertex k that travels every day to vertex i
    end
    Np(i)=Npsum;
    %break;
    Npsum=0;
end

%% end to define the present population at node i.
% initialize the set of equations for susceptible and infected and recovered ODEs at times t
dydt = zeros(3*m,1);                                                       % this is the susceptible and infected differential equations
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





