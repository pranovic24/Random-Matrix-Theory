%III. Numerical Solutions
%Question 5 Part 2.
close all;
clear all;

nr = 4;
nt = 4;
n = nr;
Pi = 10^(5/10);
P_dB = 0:40;
P = 10.^(P_dB/10);
P_eqv = P/Pi;

ergodic_capacity_su = zeros(1,length(P));
mean_capacity_su = zeros(1,length(P));





set(gcf,'color','w');
plot(P_dB,mean_capacity_su,'bo-','LineWidth',1);
legend('Simulated','Location','NorthWest');
name = ['Ergodic Capacity E[C_s_u] of the System'];
xlabel('Transmit Power P(dB)');
ylabel('Average Ergodic Capacity');
title(name);
grid on;