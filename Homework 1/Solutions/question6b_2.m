%III. Numerical Solutions
%Question 6 Part b.
close all;
clear all;

n = 8; %2,4,6,8
%equivalent to nr and nt 
K = 1;
Pi = 10^(5/10);
P_dB = 0:40;
P = 10.^(P_dB/10);

ergodic_capacity = zeros(1,length(P));
mean_capacity_per_antenna = zeros(1,length(P));

for i = 1:length(P)
    for j = 1:20000
        H = randn(n,n) + 1i*randn(n,n);
        A = H*(H');
        Hi = randn(n,K*n) + 1i*randn(n,K*n);
        B = Hi*(Hi');
        C = log(det(eye(n)+(P(i)/(Pi/K))*A*(inv(B))));
        ergodic_capacity(j) = C;
    end
    mean_capacity_per_antenna(i) = (mean(ergodic_capacity))/(n);
end

N = n;
m1 = n;
m2 = K*n;
steps = 0.001;
lamda = (0:steps:1)+ rand*1e-5;


c1 = N/m1;
c2 = N/m2;
b0 = (c1-c2)*lamda-c1+2;
b1 = (2*c1-2*c2)*(lamda.^2)+(2-3*c1+c2)*lamda+c1-1;
b2 = (c1-c2)*(lamda.^3)+((-1)*2*c1+c2)*(lamda.^2)+c1*lamda;
p_jac = (sqrt(4.*b0.*b2-(b1.^2)))./(2.*pi.*b2);

theoretical_capacity = zeros(1,length(P));
for i = 1:length(P)
    capacity = N*sum(log(1+(P(i)/(Pi/K))*(lamda./(1-lamda))).*p_jac.*steps);
    theoretical_capacity(i) = capacity;
end
theoretical_capacity_per_antenna = theoretical_capacity./N;

set(gcf,'color','w');
plot(P_dB,mean_capacity_per_antenna,'LineWidth',1)
hold on;
plot(P_dB,theoretical_capacity_per_antenna,'LineWidth',1);
hold off;
legend('Simulated','Approximation');
name = ['Ergodic Capacity per Receive Antenna (n = ',num2str(n),')'];
title(name);
grid on

