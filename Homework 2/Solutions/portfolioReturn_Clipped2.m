function [varianceReturn] = portfolioReturn_Clipped2(K)
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
    D_data = diag(diag(sigma_s_data));
    Cs_data = (D_data^(-1/2))*(sigma_s_data)*(D_data^(-1/2));

    [Cs_data_eVector,Cs_data_eValue] = eig(Cs_data);
    [Cs_data_eValue_sorted, Cs_data_eValue_index] = sort(diag(Cs_data_eValue),'descend');

    Cs_data_eVector_sorted = Cs_data_eVector(:,Cs_data_eValue_index);
    
    if K == 1
        constant = mean(Cs_data_eValue_sorted(1:98));

        Cs_data_eValue_sorted_new = Cs_data_eValue_sorted;
        Cs_data_eValue_sorted_new(1:98) = constant;
        
    elseif K >= 2 && K < 98
        constant = mean([Cs_data_eValue_sorted(1:1);Cs_data_eValue_sorted(K+1:98)]);
        
        Cs_data_eValue_sorted_new = Cs_data_eValue_sorted;
        Cs_data_eValue_sorted_new(K+1:98) = constant;
        Cs_data_eValue_sorted_new(1:1) = constant;
        
    elseif K == 98
        constant = mean(Cs_data_eValue_sorted(1:1));
        Cs_data_eValue_sorted_new = Cs_data_eValue_sorted;
        Cs_data_eValue_sorted_new(1:1) = constant;
    end

    Cs_data_eValue_sorted_diag = diag(Cs_data_eValue_sorted_new);

    Cs_clipped = (Cs_data_eVector_sorted)*(Cs_data_eValue_sorted_diag)*(Cs_data_eVector_sorted');
    sigma_s_clipped = (D_data^(1/2))*Cs_clipped*(D_data^(1/2));
    
    portfolio_clipped = ((sigma_s_clipped^-1)*ones(98,1))/(ones(1,98)*(sigma_s_clipped^-1)*ones(98,1));
    Xs_next = Xs(:,199+day_start+1);
    valueReturn_next = (portfolio_clipped')*Xs_next;
    valueReturn(i-199) = valueReturn_next;
    
    day_start = day_start+1;
end

varianceReturn = var(valueReturn);

end
