%III. Numerical Solutions
%Question 2.
close all;
clear all;

N = 4; %4,20,40
m1 = 5; %5,25,50
m2 = 6; %6,30,60

alpha = m1 - N; 
beta = m2 - N;
X = 0:0.01:1;
eigenC = [];


for i = 1:20000
    Y = 1/sqrt(2)*(randn(N,m1) + 1i*randn(N,m1));
    A = Y*(Y');

    Z = 1/sqrt(2)*(randn(N,m2) + 1i*randn(N,m2));
    B = Z*(Z');

    C = A*(inv(A+B));
    eigenvalues = eig(C);
    eigenC = [eigenC; eigenvalues];
end

polynomial = compute_jacobi(N-1,alpha,beta,X);

polynomialSum = sum(polynomial.^2);

p_lamda = (1/N).*(polynomialSum).*(X.^alpha).*((1-X).^beta);

set(gcf,'color','w');
histn(eigenC,0,0.01,1);
hold on;
plot(X,p_lamda,'LineWidth',2);
hold off;
legend('Simulated','Theoretical');
name = ['N = ',num2str(N),', m_1 = ',num2str(m1),', m_2 = ',num2str(m2)];
title(name);


