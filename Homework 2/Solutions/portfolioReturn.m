function [totalReturn,averageReturn,varianceReturn] = portfolioReturn(days)

load('SP100_2011_2013.mat');
n = length(Y);
Y = Y';
sampleMean = mean(Y);
for i = 1:n
    Xs(i,:) = Y(i,:) - sampleMean;
end
Xs = Xs';

lastDay = 752;
valueReturn = zeros(lastDay-days+1,1);
day_start = 1;

for i = days:lastDay
    Xs_data_days = Xs(:,day_start:days-1+day_start);
    n_days = length(Xs_data_days);
    
    sigma_s_days = (1/n_days)*(Xs_data_days)*(Xs_data_days');
    
    portfolio_days = ((sigma_s_days^-1)*ones(98,1))/(ones(1,98)*(sigma_s_days^-1)*ones(98,1));
    Xs_nextDay = Xs(:,days-1+day_start+1);
    valueReturn_nextDay = (portfolio_days')*Xs_nextDay;
    valueReturn(i-(days-1)) = valueReturn_nextDay;
    day_start = day_start+1;
end


totalReturn = sum(valueReturn);
averageReturn = mean(valueReturn); 
varianceReturn = var(valueReturn);


