%III. Numerical Solutions
%Question 4.
close all;
clear all;

nt = 4;
nr = 4; 
K = 4; 
Pi = 10^(5/10);
P = 10.^([0:2:40]/10);
P_db = 0:2:40;


ergodic_capacity = zeros(1,length(P));
average_capacity = zeros(1,length(P));

for i = 1:length(P)
    for j = 1:10000

        H = randn(nr,nt) + 1*j*randn(nr,nt);
        A = H*(H');

        Hi = randn(nr,K*nt) + 1*j*randn(nr,K*nt);
        B = Hi*(Hi');

        C = log(det(eye(nr)+(P(i)/(Pi/K))*A*(inv(B))));
        ergodic_capacity(i) = C;
    end
    average_ergodic(i) = mean(ergodic_capacity);
end

plot(P_db,average_ergodic);