function Np = TD_PresentPopFlow(m,matrix_phi,vec_N,tspan)
Np = zeros(m,length(tspan)); % initial present population at i
Npsum = 0;
for t = 1 : length(tspan)
    for j = 1: m
        for k = 1: m
            %total population at node j
            Npsum = Npsum+matrix_phi(k,j,t)*vec_N(k); % phi(k,j)*vec_N(k,1)is the total population at vertex k that travels every day to vertex j
        end
        Np(j,t) = Npsum;
        %break;
        Npsum = 0;
    end
end

end