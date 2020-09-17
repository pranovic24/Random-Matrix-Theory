close all;
clear all;
clc;
load('SP100_2011_2013.mat');
m = size(Y,1);
n = length(Y);
Q = n/m;

Y = Y';
sampleMean = mean(Y);
for i = 1:n
    Xs(i,:) = Y(i,:) - sampleMean;
end
Xs = Xs';

%Section 2
%% Question 1.
% 200 Days
Xs_200 = Xs(:,1:200);
sigma_s_200 = (1/200)*(Xs_200)*(Xs_200');
portfolio_200 = ((sigma_s_200^-1)*ones(98,1))/(ones(1,98)*(sigma_s_200^-1)*ones(98,1));
Xs_201st = Xs(:,201);
valueReturn_201st = (portfolio_200')*Xs_201st;

%% Question 2.
[totalReturn_200,averageReturn_200,varianceReturn_200] = portfolioReturn(200); %Used in Question 5

%% Question 3. 
w = [100,120,150,180,200,300,400,600];
totalReturn = zeros(length(w),1);
averageReturn = zeros(length(w),1);
varianceReturn = zeros(length(w),1);
for i = 1:length(w)
    [totalR, averageR, varianceR] = portfolioReturn(w(i));
    totalReturn(i) = totalR;
    averageReturn(i) = averageR;
    varianceReturn(i) = varianceR;
end

figure(1);
loglog(w,varianceReturn,'-ro');
set(gcf,'color','w');
title('Variance of Daily Return')
xlabel('Window Size (w)')
ylabel('Variance Return (Log Scale)')
xlim([100 600]);
ylim([10^-5 10^-2]);
grid on;

figure(2);
plot(w,averageReturn,'-ro');
title('Average of Total Return');
xlabel('Window Size (w)')
ylabel('Average Return')
set(gcf,'color','w');
grid on;

figure(3);
plot(w,totalReturn,'-ro');
title('Total Return');
xlabel('Window Size (w)')
ylabel('Total')
set(gcf,'color','w');
grid on;


%% Question 4b.
K = [1:40,60,98];
totalReturn_clipped = zeros(length(K),1);
averageReturn_clipped = zeros(length(K),1);
varianceReturn_clipped = zeros(length(K),1);
for i = 1:length(K)
    [totalR_clipped,averageR_clipped,varianceR_clipped] = portfolioReturn_Clipped(K(i));
    totalReturn_clipped(i) = totalR_clipped;
    averageReturn_clipped(i) = averageR_clipped;
    varianceReturn_clipped(i) = varianceR_clipped;
end

figure(4);
plot(K,varianceReturn_clipped,'-ro');
title('Variance of Daily Return (Clipped)');
xlabel('K')
ylabel('Variance Return')
set(gcf,'color','w');
grid on;

% figure(5);
% plot(K,totalReturn_clipped,'-ro');
% title('Total Return (Clipped)');
% xlabel('K')
% ylabel('Total Return')
% set(gcf,'color','w');
% grid on;


%% Question 4c.
varianceReturn_clipped2 = zeros(length(K),1);
for i = 1:length(K)
    [varianceR_clipped2] = portfolioReturn_Clipped2(K(i));
    varianceReturn_clipped2(i) = varianceR_clipped2;
end
figure(6);
plot(K,varianceReturn_clipped,'-o');
hold on;
plot(K,varianceReturn_clipped2,'-o');
legend
hold off;
title('Variance of Daily Return (Clipped)');
legend('Clipped','Clipped 2');
xlabel('K')
ylabel('Variance Return')
set(gcf,'color','w');
grid on;

%% Question 6.
t = [1,5:5:60];
totalReturn_clipped_period = zeros(length(t),1);
varianceReturn_clipped_period = zeros(length(t),1);
for i = 1:length(t)
    [totalR_clipped_period,varianceR_clipped_period] = portfolioReturn_Clipped_Period(t(i));
    totalReturn_clipped_period(i) = totalR_clipped_period;
    varianceReturn_clipped_period(i) = varianceR_clipped_period;
end
figure(7);
plot(t,varianceReturn_clipped_period,'-ro');
title('Variance of Daily Return (Clipped with K*=10)');
xlabel('t (days)')
ylabel('Variance Return')
set(gcf,'color','w');
grid on;

%% Question 7.
totalReturn_directClipped = zeros(length(K),1);
varianceReturn_directClipped = zeros(length(K),1);
for i = 1:length(K)
    [totalR_directClipped,varianceR_directClipped] = portfolioReturn_DirectClipped(K(i));
    totalReturn_directClipped(i) = totalR_directClipped;
    varianceReturn_directClipped(i) = varianceR_directClipped;
end

figure(8);
plot(K,varianceReturn_clipped,'-o');
hold on; 
plot(K,varianceReturn_directClipped,'-o');
legend('Indirect Clipping','Direct Clipping');
hold off;
title('Variance of Daily Return (Direct Clipped)');
xlabel('K')
ylabel('Variance Return')
set(gcf,'color','w');
grid on;

% figure(9);
% plot(K,totalReturn_directClipped,'-ro');
% title('Total Return (Direct Clipped)');
% xlabel('K')
% ylabel('Total Return')
% set(gcf,'color','w');
% grid on;
