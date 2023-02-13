clc
clear
clc
sl = 10000;
I = 2*randi([0 1], 1, sl) - 1;
L = 2;
x = [1/12 1/2 5/6 1/2 1/12];

SNR = 13;
y = filter(x, 1, I);
noise = randn(1, length(y)) + 1j*randn(1, length(y));

c = 0.5533;
F = [c 1/2/(c + 1/12/c) 1/12/c ];

filtered_noise = filter(F, 1, noise);
Pe = zeros(1, 13);
for SNR = 0:12
    
    filtered_noise  = filtered_noise /sqrt(mean(abs(filtered_noise ).^2)) * ...
        sqrt((mean(abs(y).^2)/(10^(SNR/10))));
    y_n = y + filtered_noise;
    states = zeros(4,6*L);
    prev_states = zeros(4,6*L);
    costs = zeros(4,1);
    prev_costs = zeros(4,1);
    
    z = zeros(1,length(I) - 6*L);
    I_vec = [-1; 1];
    ind1 = [1 3];
    ind2 = [2 4];
    
    for i=1:length(y_n)
        for j=1:2
            tmp_cost = real(I_vec(j) * (2*y_n(i) - I_vec(j)*x(3) - 2* I_vec(1) * x(4) - 2*I_vec * x(5)));
            [costs(j) , index] =  max(tmp_cost+[prev_costs(1);prev_costs(3)]);
            states(j,:)= [prev_states(ind1(index),2:end) I_vec(j)];
        end
        for j=3:4
            tmp_cost = real(I_vec(j-2) * (2*y_n(i) - I_vec(j-2)*x(3) - 2* I_vec(2) * x(4) - 2*I_vec * x(5)));
            [costs(j) , index] = max(tmp_cost+[prev_costs(2);prev_costs(4)]);
            states(j,:)= [prev_states(ind2(index),2:end) I_vec(j-2)];
        end
        prev_states = states;
        prev_costs = costs;
        [~,index ] = max(costs);
        z(i) = states(index,1);
    end
    error = sum(z(6*L+2:end) ~= I(1:end-(6*L+1)));
    Pe(SNR+1) = error / (length(z)-13);
    
    SNR;
end
SNR = 0:12;
semilogy(SNR, Pe,'mx-','Linewidth',2)
axis([0 12 10^-5 1])
grid on

xlabel('SNR-db')
ylabel('Pe')