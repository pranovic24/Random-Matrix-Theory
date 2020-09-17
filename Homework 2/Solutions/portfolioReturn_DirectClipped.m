function [totalReturn,varianceReturn] = portfolioReturn_DirectClipped(K)

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
value_return = zeros(553,1);
day_start = 1;

for i = 200:752
    Xs_data = Xs(:,day_start:199+day_start);
    n_days = length(Xs_data);
    sigma_s_data = (1/n_days)*(Xs_data)*(Xs_data'); 
    [sigma_s_data_eVector,sigma_s_data_eValue] = eig(sigma_s_data);
    [sigma_s_data_eValue_sorted, sigma_s_data_eValue_index] = sort(diag(sigma_s_data_eValue),'descend');
    
    sigma_s_data_eVector_sorted = sigma_s_data_eVector(:,sigma_s_data_eValue_index);
    constant = mean(sigma_s_data_eValue_sorted(K+1:98));
    
    sigma_s_data_eValue_sorted_new = sigma_s_data_eValue_sorted;
    sigma_s_data_eValue_sorted_new(K+1:98) = constant;
    
    sigma_s_data_eValue_sorted_diag = diag(sigma_s_data_eValue_sorted_new);
    sigma_s_clipped = (sigma_s_data_eVector_sorted)*(sigma_s_data_eValue_sorted_diag)*(sigma_s_data_eVector_sorted');
    
    portfolio_clipped = ((sigma_s_clipped^-1)*ones(98,1))/(ones(1,98)*(sigma_s_clipped^-1)*ones(98,1));
    Xs_next = Xs(:,199+day_start+1);
    valueReturn_next = (portfolio_clipped')*Xs_next;
    valueReturn(i-199) = valueReturn_next;
    
    day_start = day_start+1;
end

totalReturn = sum(valueReturn);
% averageReturn = mean(valueReturn); 
varianceReturn = var(valueReturn);

end











