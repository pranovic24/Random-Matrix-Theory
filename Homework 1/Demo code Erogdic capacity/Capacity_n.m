%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
% Demo Code for ELEC 6910h, Random Matrix Theory (Spring 2016)      %
%                                                                   %
% Description: This function will numerically plot the theoretical  %
% average ergodic capacity of a channel (with a number of receiving % 
% antennas equal to the number of transmitting antennas) and the    %
% average ergodic capacity obtained using the MP-law, vs. the SNR   %
%                                                                   %
% Authors: Nicolas Auguin and Matthew McKay                         %
%                                                                   %
% Date: 9th March, 2016                                             %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script calls for the function Laguerre.m, which computes the
% coefficients of the n-th order Laguerre polynomial

close all 
clear all


P = [0:2.5:35];     % SNR
P1 = 10.^(P/10);
asymp_capacity = zeros(1,length(P));

% Average asymptotic ergodic capacity (computed using the MP-law)
wd = 0;
for p=P1
    wd = wd +1;
    f = @(x) 1/pi*sqrt(1./x-1/4).*log(1+p*x);   % MP-law for n_r=n_t
    asymp_capacity(wd) = quadgk(f,0.01,4);       % Evaluates the integral for computing the asymptotic ergodic capacity
end


% Average true ergodic capacity
step = 0.001;
x = 0:step:4.5;

N = 5;             % Number of receiving/transmitting antennas

figure;


for n = 1:N

Laguerre_sum = exp(-n*x);       % Laguerre(0) = 1

for i=1:n-1
    
    xx = zeros(length(x),i+1);
    for j=1:length(x)
        xx(j,:) = sqrt(exp(-n*x(j))).*((n*x(j)).^[i:-1:0]);
    end
    
    Laguerre_sum = Laguerre_sum + (Laguerre(i)*xx').^2;     % Laguerre(i) calculates the coefficients of the i-th Laguerre polynomial
    
end

p_lambda = Laguerre_sum;

finite_capacity = zeros(1,length(P));

wd = 0;
for p=P1
    wd = wd + 1;
    finite_capacity(wd) = n*sum(step*p_lambda.*log(1+p*x));
    
end

if n==1
    plot(P,finite_capacity,'b-','LineWidth',2);
end
if n==2
    plot(P,finite_capacity,'g-','LineWidth',2);
end
if n==3
    plot(P,finite_capacity,'m-','LineWidth',2);
end
if n==4
    plot(P,finite_capacity,'y-','LineWidth',2);
end
if n==5
    plot(P,finite_capacity,'k-','LineWidth',2);
end
hold on;
grid on;

end


plot(P,N*asymp_capacity,'ro','LineWidth',3);
xlabel('SNR (dB)','fontsize',12);
ylabel('Average ergodic capacity','fontsize',12);
legend('True Capacity (n=1)','True Capacity (n=2)','True Capacity (n=3)','True Capacity (n=4)','True Capacity (n=5)','Asymptotic capacity -- MP-law (n=5)','Location','Northwest');
axis([0 35 0 40])


%%%%%%% Average capacity per antenna %%%%%%%%

P = [0:2.5:35];     % SNR
P1 = 10.^(P/10);
asymp_capacity = zeros(1,length(P));

% Average asymptotic ergodic capacity (computed using the MP-law)
wd = 0;
for p=P1
    wd = wd +1;
    f = @(x) 1/pi*sqrt(1./x-1/4).*log(1+p*x);   % MP-law for n_r=n_t
    asymp_capacity(wd) = quadgk(f,0,4);       % Evaluates the integral for computing the asymptotic ergodic capacity
end


% Average true ergodic capacity
step = 0.001;
x = 0:step:4.5;

N = 5;             % Number of receiving/transmitting antennas

figure;

for n = 1:N

Laguerre_sum = exp(-n*x);       % Laguerre(0) = 1

for i=1:n-1
    
    xx = zeros(length(x),i+1);
    for j=1:length(x)
        xx(j,:) = sqrt(exp(-n*x(j))).*((n*x(j)).^[i:-1:0]);
    end
    
    Laguerre_sum = Laguerre_sum + (Laguerre(i)*xx').^2;     % Laguerre(i) calculates the coefficients of the i-th Laguerre polynomial
    
end

p_lambda = Laguerre_sum;

finite_capacity = zeros(1,length(P));

wd = 0;
for p=P1
    wd = wd + 1;
    finite_capacity(wd) = sum(step*p_lambda.*log(1+p*x));   % Evaluates integral for computing the ergodic capacity
    
end

if n==1
    plot(P,finite_capacity,'b-','LineWidth',2);
end
if n==2
    plot(P,finite_capacity,'g-','LineWidth',2);
end
if n==3
    plot(P,finite_capacity,'m-','LineWidth',2);
end
if n==4
    plot(P,finite_capacity,'y-','LineWidth',2);
end
if n==5
    plot(P,finite_capacity,'k-','LineWidth',2);
end
hold on;
grid on;

end


plot(P,asymp_capacity,'ro','LineWidth',3);
xlabel('SNR (dB)','fontsize',12);
ylabel('Average ergodic capacity per antenna','fontsize',12);
legend('True Capacity (n=1)','True Capacity (n=2)','True Capacity (n=3)','True Capacity (n=4)','True Capacity (n=5)','Asymptotic capacity (MP-law)','Location','Northwest');
axis([0 35 0 10])
% axis([20 35 2 3.5])

