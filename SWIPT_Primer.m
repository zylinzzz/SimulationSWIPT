clear, close all
Ps = 5e-3;  % The PSD constraint is 1 w
Pt = 1;
% Bs = 1*10^5; %The bandwidth is 100 KHz
delta_sq = 1; % The noise covariance
beta = 0.7; % Variance
% M = 1; % Number of transmitter antennas
N = 40; % Number of channels
Pc = linspace(0,60,7); % Power needed for circuit


for M = [1 2 4 8]
    H = normrnd(0,beta,M,N)+1j*normrnd(0,beta,M,N); % Modeling the channel
    h_sq = zeros(1,N);
    Q = zeros(1,13);
    R = zeros(1,13);
    N_E = zeros(size(Pc));
    for k = 1:13
        for i = 1:N
            h_sq(i) = norm(H(:,i),2).^2;
        end
        [h_sq_sorted, sorted_index] = sort(h_sq,'descend');
        h_sorted = H(:,sorted_index);

        %% Find the optimal number of sub-bands used for energy transfer
        for i = 1:N
            if Pc(k) == 0
                N_E(k) = 0;
                break;
            else
                Q(k) = Q(k) + Ps * h_sq_sorted(i);
                if Q(k) > Pc(k)
                    Q(k) = Q(k) - Ps * h_sq_sorted(i);
                    N_E(k) = i;
                    break
                end
            end
             %% All channel used for energy transmission
            if i == N
                N_E(k) = N;
                break;
            end
        end
        
        %% Calculate the v_opt[n]
        v = zeros(M,N);
        Pe = zeros(1,13);
        for i = 1:N_E(k)
            if i < N_E(k)
                v(:,i) = sqrt(Ps).* h_sorted(:,i)/norm(h_sorted(:,i),2);
            else
                v(:,i) = (Pc(k) - Q(k)).* h_sorted(:,i)/ ...
                    (norm(h_sorted(:,i),2))^3;
            end
        end
        
         %% Calcute the energy used for info. and energy
        for i = 1:N_E(k)
            Pe(k) = Pe(k)+norm(v(:,i))^2;
        end
%         Pe(k) = sum(norm(v(:,1:N_E(k)))^2);
        Pi = max(Pt - Pe(k),0);
        P_opt = WF_OPT(N,N_E(k),Pi,Ps,h_sq_sorted,delta_sq);
        for i = N_E(k)+1:N
            v(:,i) = sqrt(P_opt(i-N_E(k)))* ...
                conj(h_sorted(:,i))/norm(h_sorted(:,i),2);
        end

        %% Calculate the sum rate
        for i = N_E(k)+1:N
            R(k) = R(k) + log2(1+abs(h_sorted(:,i)' * v(:,i))^2/delta_sq);
        end     
    end
    % plot(Pc,N_E,'-*')
    % grid on
    % figure
    plot(Pc,R,'-*')
    hold on
    grid on
end
legend('Pt = 140,M = 1,N = 40', ...
    'Pt = 40,M = 2,N = 40', ...
    'Pt = 40,M = 4,N = 40', ...
    'Pt = 40,M = 8,N = 40')
xlabel('Circuit Power Consumption (watta)')
ylabel('Sum-Rate (bits per hertz)')
