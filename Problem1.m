clear, close all
Ps = 5e-3;  % The PSD constraint is 1 w
Pt = 1;
% Bs = 1*10^5; %The bandwidth is 100 KHz
delta_sq = 1e-5; % The noise covariance
beta = 0.5; % Variance
M = 4; % Number of transmitter antennas
N = 200; % Number of channels
points = 13;
Pc = linspace(0,60e-3,points); % Power needed for circuit
d1 = 6;
H = normrnd(0,beta,M,N)+1j*normrnd(0,beta,M,N); % Modeling the channel

for d2 =  [7,8]
    H = normrnd(0,beta,M,N)+1j*normrnd(0,beta,M,N); % Modeling the channel
    h_sq = zeros(1,N);
    Q1 = zeros(1,points);
    Q2 = zeros(1,points);
    R = zeros(1,points);
    RN = zeros(1,points);
    RF = zeros(1,points);
    N_E1 = zeros(size(Pc));
    N_E2 = N*ones(size(Pc));
    N_I1 = zeros(size(Pc));
    
    for k = 1:points
        for i = 1:N
            h_sq(i) = norm(H(:,i),2).^2;
        end
        [h_sq_sorted, sorted_index] = sort(h_sq,'descend');
        h_sorted = H(:,sorted_index);

        %% Find the optimal number of sub-bands used for energy transfer USER1
        for i = 1:N
            if Pc(k) == 0
                N_E1(k) = 0;
                break;
            else
                Q1(k) = Q1(k) + Ps * h_sq_sorted(i)/d1^2;
                if Q1(k) > Pc(k)
                    Q1(k) = Q1(k) - Ps * h_sq_sorted(i)/d1^2;
                    N_E1(k) = i;
                    break
                end
            end
             %% All channel used for energy transmission
            if i == N
                N_E1(k) = N;
                break;
            end
        end
        
        %% Find the optimal number of sub-bands used for energy transfer USER2
        if N_E1(k)~= N 
            for i = N_E1(k):N
                if Pc(k) == 0
                    N_E2(k) = 0;
                    break
                else
                    Q2(k) = Q2(k) + Ps * h_sq_sorted(i)/d2^2;
                    if Q2(k) > Pc(k)
                        Q2(k) = Q2(k) - Ps * h_sq_sorted(i)/d2^2;
                        N_E2(k) = i;
                        break
                    end
                end
                %% All channel used for energy transmission
                if i == N
                    N_E2(k) = N;
                    break;
                end
            end
        end
        
        %% Calculate the v_opt[n]
        v = zeros(M,N);
        Pe = zeros(1,points);
        
        %% Pe1
        for i = 1:N_E1(k)
            if i < N_E1(k)
                v(:,i) = sqrt(Ps).* h_sorted(:,i)/(d1* ...
                    norm(h_sorted(:,i),2));
            else
                v(:,i) = (Pc(k) - Q1(k)).* d1*h_sorted(:,i)/ ...
                    (norm(h_sorted(:,i),2))^2;
            end
        end
        
        %% Pe2
        for i = N_E1(k)+1:N_E2(k)
            if i < N_E2(k)
                v(:,i) = sqrt(Ps).* h_sorted(:,i)/(d2* ...
                    norm(h_sorted(:,i),2));
            else
                v(:,i) = (Pc(k) - Q1(k)).* d2*h_sorted(:,i)/ ...
                    (norm(h_sorted(:,i),2))^2;
            end
        end
        
         %% Calcute the energy used for info. and energy
        for i = 1:N_E2(k)
            Pe(k) = Pe(k)+norm(v(:,i))^2;
        end
        Pi = max(Pt - Pe(k),0);
        P_opt = WF_OPT(N,N_E1(k),Pi,Ps,h_sq_sorted,delta_sq);
        
        N_I1(k) = N_E2(k) + ceil((N-N_E2(k))/2);
        
        for i = N_E2(k)+1:N_I1(k)
            v(:,i) = sqrt(P_opt(i-N_E2(k)))* ...
                (h_sorted(:,i))/ ...
                (d1*norm(h_sorted(:,i),2));
        end
        
        for i = N_I1(k)+1:N
            v(:,i) = sqrt(P_opt(i-N_I1(k)))* ...
                (h_sorted(:,i))/ ...
                (d2*norm(h_sorted(:,i),2));
        end
        

        %% Calculate the sum rate
        for i = N_E1(k)+1:N_I1
            RN(k) = RN(k) + log2(1+abs(h_sorted(:,i)' * v(:,i))^2/ ...
                (d1^2*delta_sq));
        end   
        for i = N_I1(k)+1:N
            RF(k) = RF(k) + log2(1+abs(h_sorted(:,i)' * v(:,i))^2/ ...
                (d2^2*delta_sq));
        end
        
        R(k) = RN(k)+RF(k);
    end
    
    plot(Pc,R,'-*')
    hold on
    plot(Pc,RN,'-*')
    plot(Pc,RF,'-*')
    grid on
end
legend('R(d1 = 6,d2 = 7)','R-NU(d1 = 6,d2 = 7)','R-FU(d1 = 6,d2 = 7)', ...
    'R(d1 = 6,d2 = 8)','R-NU(d1 = 6,d2 = 8)','R-FU(d1 = 6,d2 = 8)')
xlabel('Circuit Power Consumption (watta)')
ylabel('Sum-Rate (bits per hertz)')
