%III. Numerical Solutions
%Question 5 Part 1.
close all;
clear all;

nt = 4;
nr = 4;
%n = 4;
K = 5:5:35;
Pi = 10^(5/10);
P_dB = 0:40;
P = 10.^(P_dB/10);
P_eqv = P/Pi;

ergodic_capacity = zeros(1,length(P));
mean_capacity = zeros(1,length(P));
ergodic_capacity_su = zeros(1,length(P_eqv));
mean_capacity_su = zeros(1,length(P_eqv));

for n = 1:length(K)
    for i = 1:length(P)
        for j = 1:20000
            H = 1/sqrt(2)*(randn(nr,nt) + 1i*randn(nr,nt));
            A = H*(H');
            Hi = 1/sqrt(2)*(randn(nr,K(n)*nt) + 1i*randn(nr,K(n)*nt));
            B = Hi*(Hi');
            C = log(det(eye(nr)+(P(i)/(Pi/K(n)))*A*(inv(B))));
            ergodic_capacity(j) = C;
        end
        mean_capacity(i) = mean(ergodic_capacity);
    end
    plot(P_dB,mean_capacity,'LineWidth',1);
    hold on; 
end

for l = 1:length(P_eqv)
    for m = 1:20000
        H = 1/sqrt(2)*(randn(nr,nt) + 1i*randn(nr,nt));
        A = H*(H');
        C_su = log(det(eye(nr)+(P_eqv(l)/nr)*A));
        ergodic_capacity_su(m) = C_su;
    end
    mean_capacity_su(l) = mean(ergodic_capacity_su);
end

plot(P_dB,mean_capacity_su,'bo-','LineWidth',1);
hold off;

set(gcf,'color','w');
legend('K=5','K=10','K=15','K=20','K=25','K=30','K=35','Ergodic Capacity E[C_s_u]','Location','Northwest');
name = ['Ergodic Capacity of the System'];
xlabel('Transmit Power P(dB)');
ylabel('Average Ergodic Capacity');
title(name);
grid on;