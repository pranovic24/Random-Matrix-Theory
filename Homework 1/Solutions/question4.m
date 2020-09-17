%III. Numerical Solutions
%Question 4.
close all;
clear all;

nt = 4;
nr = 4;
K = 4;
Pi = 10^(5/10);
P_dB = 0:40;
P = 10.^(P_dB/10);

ergodic_capacity = zeros(1,length(P));
mean_capacity = zeros(1,length(P));

for i = 1:length(P)
    for j = 1:20000
        H = 1/sqrt(2)*(randn(nr,nt) + 1i*randn(nr,nt));
        A = H*(H');
        Hi = 1/sqrt(2)*(randn(nr,K*nt) + 1i*randn(nr,K*nt));
        B = Hi*(Hi');
        C = log(det(eye(nr)+(P(i)/(Pi/K))*A*(inv(B))));
        ergodic_capacity(j) = C;
    end
    mean_capacity(i) = mean(ergodic_capacity);
end

N = nr;
m1 = nt;
m2 = K*nt;
alpha = m1-N;
beta = m2-N;
steps = 0.001;
lamda = (0:steps:1) + rand*1e-5;

polynomial = compute_jacobi(N-1,alpha,beta,lamda);
polynomialSum = sum(polynomial.^2);
p_lamda = (1/N).*(polynomialSum).*(lamda.^alpha).*((1-lamda).^beta);

theoretical_capacity = zeros(1,length(P));
for i = 1:length(P)
    capacity = N*sum(log(1+(P(i)/(Pi/K))*(lamda./(1-lamda))).*p_lamda.*steps);
    theoretical_capacity(i) = capacity; 
end



set(gcf,'color','w');
plot(P_dB,mean_capacity,'b-','LineWidth',2);
hold on; 
plot(P_dB,theoretical_capacity,'or','LineWidth',2);
hold off;
legend('Simulated','Theoretical','Location','NorthWest');
name = ['Ergodic Capacity of the System'];
xlabel('Transmit Power P(dB)');
ylabel('Average Ergodic Capacity');
title(name);
grid on;




