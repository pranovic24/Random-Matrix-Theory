%III. Numerical Solutions
%Question 3.
close all;
clear all;

N = 40; %4,20,40
m1 = 50; %5,25,50
m2 = 60; %6,30,60

c1 = N/m1;
c2 = N/m2;
lamda = 0:0.01:1;


b0 = (c1-c2)*lamda-c1+2;
b1 = (2*c1-2*c2)*(lamda.^2)+(2-3*c1+c2)*lamda+c1-1;
b2 = (c1-c2)*(lamda.^3)+((-1)*2*c1+c2)*(lamda.^2)+c1*lamda;


% alpha = m1 - N; 
% beta = m2 - N;
X = 0:0.01:1;
eigenC = [];


for i = 1:20000
    Y = 1/sqrt(2)*(randn(N,m1) + 1i*randn(N,m1));
    A = Y*(Y');

    Z = 1/sqrt(2)*(randn(N,m2) + 1i*randn(N,m2));
    B = Z*(Z');

    C = (1/m1).*A*(inv((1/m1).*A+(1/m2).*B));
    eigenvalues = eig(C);
    eigenC = [eigenC; eigenvalues];
end

%polynomial = compute_jacobi(N-1,alpha,beta,X);

%polynomialSum = sum(polynomial.^2);

%p_lamda = (1/N).*(polynomialSum).*(X.^alpha).*((1-X).^beta);

p_jac = (sqrt(4.*b0.*b2-(b1.^2)))./(2.*pi.*b2);
set(gcf,'color','w');
histn(eigenC,0,0.01,1);
hold on;
plot(X,p_jac,'LineWidth',2);
hold off;
legend('Simulated','Theoretical');
name = ['N = ',num2str(N),', m_1 = ',num2str(m1),', m_2 = ',num2str(m2)];
title(name);


