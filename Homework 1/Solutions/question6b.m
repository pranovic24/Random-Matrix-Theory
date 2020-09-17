%III. Numerical Solutions
%Question 6 Part b.
close all;
clear all;


n = 2:2:8; %equivalent to nr and nt 
K = 1;
Pi = 10^(5/10);
P_dB = 0:40;
P = 10.^(P_dB/10);

ergodic_capacity = zeros(1,length(P));
mean_capacity_per_antenna = zeros(1,length(P));

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
        mean_capacity_per_antenna(j) = (mean(ergodic_capacity))/(n(i));
        
    end
    plot(P_dB,mean_capacity_per_antenna,'LineWidth',1);
    hold on;
end

N = n;
m1 = n;
m2 = K*n;
steps = 0.001;
%lamda = (0:steps:1) + rand*1e-5;
lamda = 0.001:steps:0.999; %Change due to numerical issues


c1 = N/m1;
c2 = N/m2;
b0 = (c1-c2)*lamda-c1+2;
b1 = (2*c1-2*c2)*(lamda.^2)+(2-3*c1+c2)*lamda+c1-1;
b2 = (c1-c2)*(lamda.^3)+((-1)*2*c1+c2)*(lamda.^2)+c1*lamda;
p_jac = real((sqrt(4.*b0.*b2-(b1.^2)))./(2.*pi.*b2));

theoretical_capacity = [];
for i = 1:length(P)
    capacity = N*sum(log(1+(P(i)/(Pi/K))*(lamda./(1-lamda))).*p_jac.*steps);
    theoretical_capacity = [theoretical_capacity; capacity];
end
theoretical_capacity_per_antenna = theoretical_capacity/N;

plot(P_dB,theoretical_capacity_per_antenna,'o-','LineWidth',1);
hold off;
set(gcf,'color','w');
legend('(n_r,n_t)=(2,2)','(n_r,n_t)=(4,4)','(n_r,n_t)=(6,6)','(n_r,n_t)=(8,8)','Approximation','Location','Northwest');
name = ['Ergodic Capacity per Receive Antenna'];
xlabel('Transmit Power P(dB)');
ylabel('Average Ergodic Capacity');
title(name);
grid on;

