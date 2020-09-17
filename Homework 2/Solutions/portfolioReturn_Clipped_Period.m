function [totalReturn,varianceReturn] = portfolioReturn_Clipped_Period(t)
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
valueReturn = [];
% day_start = 0;
K = 10;

for i = 200:t:752
    Xs_data = Xs(:,i-199:i);
    n_days = length(Xs_data);
    sigma_s_data = (1/n_days)*(Xs_data)*(Xs_data'); 
    D_data = diag(diag(sigma_s_data));
    Cs_data = (D_data^(-1/2))*(sigma_s_data)*(D_data^(-1/2));

    [Cs_data_eVector,Cs_data_eValue] = eig(Cs_data);
    [Cs_data_eValue_sorted, Cs_data_eValue_index] = sort(diag(Cs_data_eValue),'descend');

    Cs_data_eVector_sorted = Cs_data_eVector(:,Cs_data_eValue_index);
    constant = mean(Cs_data_eValue_sorted(K+1:98));

    Cs_data_eValue_sorted_new = Cs_data_eValue_sorted;
    Cs_data_eValue_sorted_new(K+1:98) = constant;

    Cs_data_eValue_sorted_diag = diag(Cs_data_eValue_sorted_new);

    Cs_clipped = (Cs_data_eVector_sorted)*(Cs_data_eValue_sorted_diag)*(Cs_data_eVector_sorted');
    sigma_s_clipped = (D_data^(1/2))*Cs_clipped*(D_data^(1/2));
    
    portfolio_clipped = ((sigma_s_clipped^-1)*ones(98,1))/(ones(1,98)*(sigma_s_clipped^-1)*ones(98,1));
    Xs_next = Xs(:,i+1:min(i+t,753));
    valueReturn_next = (portfolio_clipped')*Xs_next;
    valueReturn = [valueReturn valueReturn_next];
    
%     day_start = day_start+1;
end

totalReturn = sum(valueReturn);
% averageReturn = mean(valueReturn); 
varianceReturn = var(valueReturn);
end
