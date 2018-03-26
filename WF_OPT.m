function P_opt = WF_OPT(N,N_E,Pi,Ps,h_sq_sorted,delta_sq)
P_opt = zeros(1,N-N_E);
% mu = min(delta_sq./h_sq_sorted(N_E+1:N));% Find the water level  
% sum_P = sum(max(mu-delta_sq./h_sq_sorted(N_E+1:N),0));
% tol = 0.01;
% 
% while abs(Pi-sum_P)>tol
%     mu = mu + (Pi-sum_P)/(N-N_E);
%     sum_P = sum(max(mu-delta_sq./h_sq_sorted(N_E+1:N),0));
% end
% P_opt = (max(mu-delta_sq./h_sq_sorted(N_E+1:N),0));
length = N-N_E;
sum_p = sum(delta_sq./h_sq_sorted(N_E+1:N));
mu = (Pi+sum_p)/length;
cp = h_sq_sorted(N_E+1:N);

for i = N:-1:N_E+1
    if mu<max(delta_sq./cp)
        [~,min_index] = min(cp);
        cp(min_index) = [];
        P_opt(i-N_E) = 0;
        length = length - 1;
        sum_p = sum(delta_sq./cp);
        mu = (Pi+sum_p)/length;
    else
        P_opt(i-N_E) = mu-delta_sq/h_sq_sorted(i);
    end
end


%% The WF_with_PSD_constraint algorithm
for i = 1:N-N_E
    if Pi <= 0
        break
    end
    P_opt(i) = min(P_opt(i),Ps);
    Pi = Pi-P_opt(i);   
end
P_opt = P_opt';