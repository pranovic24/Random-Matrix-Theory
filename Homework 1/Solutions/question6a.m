%III. Numerical Solutions
%Question 6 Part a.
close all;
clear all;


n = 2:2:8; %equivalent to nr and nt 
K = 1;
Pi = 10^(5/10);
P_dB = 0:40;
P = 10.^(P_dB/10);

ergodic_capacity = zeros(1,length(P));
mean_capacity = zeros(1,length(P));

for i = 1:length(n)
    for j = 1:length(P)
        for k = 1:20000
            H = 1/sqrt(2)*(randn(n(i),n(i)) + 1i*randn(n(i),n(i)));
            A = H*(H');
            Hi = 1/sqrt(2)*(randn(n(i),K*n(i)) + 1i*randn(n(i),K*n(i)));
            B = Hi*(Hi');
            C = log(det(eye(n(i))+(P(j)/(Pi/K))*A*(inv(B))));
            ergodic_capacity(k) = C;
        end
        mean_capacity(j) = mean(ergodic_capacity);
    end
    plot(P_dB,mean_capacity,'LineWidth',1);
    hold on; 
end

set(gcf,'color','w');
legend('(n_r,n_t)=(2,2)','(n_r,n_t)=(4,4)','(n_r,n_t)=(6,6)','(n_r,n_t)=(8,8)','Location','Northwest');
name = ['Ergodic Capacity of the System'];
xlabel('Transmit Power P(dB)');
ylabel('Average Ergodic Capacity');
title(name);
grid on;