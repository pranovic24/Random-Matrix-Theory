close all;
clear all;
load('SP100_2011_2013.mat');
m = size(Y,1);
n = length(Y);
Q = n/m;

%Section 1
%% Question 1.
Y = Y';
sampleMean = mean(Y);
for i = 1:n
    Xs(i,:) = Y(i,:) - sampleMean;
end
Xs = Xs';

%% Question 2.
sigma_s = (1/n)*(Xs)*(Xs');
eigSigma_s = eig(sigma_s);

figure(1);
plot(eigSigma_s,'-r');
title('Scree-plot of Eigenvalues');
xlabel('Component Number');
ylabel('Eigenvalues');
set(gcf,'color','w');
legend('Eigenvalues');
grid on;

%% Question 3a.
D = diag(diag(sigma_s));
Cs = (D^(-1/2))*(sigma_s)*(D^(-1/2));

[Cs_eigenVector,Cs_eigenValue] = eig(Cs);
[lambdaSorted,index_lambda]=sort(diag(Cs_eigenValue),'descend');
%Cs_eigenValue_Sorted = Cs_eigenValue(index_lambda,index_lambda);
Cs_eigenVector_Sorted = Cs_eigenVector(:,index_lambda);

figure(2);
histn(lambdaSorted,0,0.1,50);
title('Eigenvalues of C_{s}');
xlabel('Eigenvalue (\lambda)');
ylabel('Probability Density P(\lambda)');
set(gcf,'color','w');
legend('Eigenvalues');
grid on;

%% Question 3b.
d_max = 1+(1/Q)+2*(sqrt(1/Q));
d_min = 1+(1/Q)-2*(sqrt(1/Q));
d = linspace(d_min,d_max);
P_rcm = (Q/(2*pi)).*(sqrt((d_max-d).*(d-d_min))./d);

figure(3);
histn(lambdaSorted,0,0.1,50);
hold on;
plot(d,P_rcm,'-r','LineWidth',2);
hold off;
title('Eigenvalues of C_{s}');
xlabel('Eigenvalue (\lambda)');
ylabel('Probability Density P(\lambda)');
legend('Eigenvalues','P_{rcm}');
set(gcf,'color','w');
grid on;

%% Question 3c.
lambdaSorted_5 = lambdaSorted(5:end);
lambdaSorted_5Norm = lambdaSorted_5*(sum(lambdaSorted)/sum(lambdaSorted_5));
figure(4);
histn(lambdaSorted_5Norm,0,0.1,max(lambdaSorted_5Norm));
hold on;
plot(d,P_rcm,'-r','LineWidth',2);
hold off;
title('Eigenvalues of C_{s} (Removal of largest 4 Eigenvalues)');
xlabel('Eigenvalue (\lambda)');
ylabel('Probability Density P(\lambda)');
legend('Eigenvalues','P_{rcm}');
set(gcf,'color','w');
grid on;

%% Question 3d.
H = randn(m,n);
W = (1/n)*(H)*(H');

[W_eigenVector, W_eigenValue] = eig(W);
[W_eigenValue_Sorted, W_eigenValue_Index] = sort(diag(W_eigenValue),'descend');
%W_eigenValue_Sorted_matrix = W_eigenValue( W_eigenValue_Index, W_eigenValue_Index);
W_eigenVector_Sorted = W_eigenVector(:,W_eigenValue_Index);

figure(5);
histn(real(W_eigenValue_Sorted),0,0.1,max(real(W_eigenValue_Sorted)));
hold on;
plot(d,P_rcm,'-r','LineWidth',2);
hold off;
title('Eigenvalues of W');
xlabel('Eigenvalue (\lambda)');
ylabel('Probability Density P(\lambda)');
legend('Eigenvalues','P_{rcm}');
set(gcf,'color','w');
grid on;

%% Question 3e.
Xs_shuffle = zeros(m,n);
for i = 1:m
    Xs_shuffle(i,:) = Xs(i,randperm(n));
end
sigma_s_shuffle = (1/n)*(Xs_shuffle)*(Xs_shuffle');
eigSigma_s_shuffle = eig(sigma_s_shuffle);
D_shuffle = diag(diag(sigma_s_shuffle));
Cs_shuffle = (D_shuffle^(-1/2))*(sigma_s_shuffle)*(D_shuffle^(-1/2));

[Cs_shuffle_eigenVector,Cs_shuffle_eigenValue] = eig(Cs_shuffle);
[Cs_shuffle_eigenValue_Sorted,Cs_shuffle_eigenValue_index]=sort(diag(Cs_shuffle_eigenValue),'descend');
%Cs_shuffle_eigenValue_Sorted_matrix = Cs_shuffle_eigenValue(Cs_shuffle_eigenValue_index,Cs_shuffle_eigenValue_index);
Cs_shuffle_eigenVector_Sorted = Cs_shuffle_eigenVector(:,Cs_shuffle_eigenValue_index);

figure(6);
histn(Cs_shuffle_eigenValue_Sorted,0,0.1,max(Cs_shuffle_eigenValue_Sorted));
hold on;
plot(d,P_rcm,'-r','LineWidth',2);
hold off;
title('Eigenvalues of C_{shuffle}');
xlabel('Eigenvalue (\lambda)');
ylabel('Probability Density P(\lambda)');
legend('Eigenvalues','P_{rcm}');
set(gcf,'color','w');
grid on;

%% Question 4a.
v = -4:0.1:4;
p_v = (1/sqrt(2*pi)).*exp((-1.*(v.^2))./2);
k = 70;

u_k = W_eigenVector_Sorted(:,k);

figure(7);
histn(sqrt(m)*(u_k),-4,0.1,4);
hold on; 
plot(v,p_v,'-r','LineWidth',2);
hold off;
title(['Eigenvector u_{',num2str(k),'} of W']);
xlabel('Eigenvector components u');
ylabel('Probability Density \rho(u)');
legend(['Eigenvector u_{',num2str(k),'}'],'Pdf');
set(gcf,'color','w');
grid on;

%% Question 4b.
%lambda_k = Cs_eigenValue_Sorted(k,k);
v_k = Cs_eigenVector_Sorted(:,k);

figure(8);
histn(sqrt(m)*(v_k),-4,0.1,4);
hold on; 
plot(v,p_v,'-r','LineWidth',2);
hold off;
title(['Eigenvector v_{',num2str(k),'} of C_{s}']);
xlabel('Eigenvector components v');
ylabel('Probability Density \rho(v)');
legend(['Eigenvector v_{',num2str(k),'}'],'Pdf');
set(gcf,'color','w');
grid on;

%% Question 4c. 
v_1 = Cs_eigenVector_Sorted(:,1);
v_2 = Cs_eigenVector_Sorted(:,2);
v_3 = Cs_eigenVector_Sorted(:,3);

figure(9);
histn(-1*sqrt(m)*(v_1),-4,0.1,4);
hold on; 
plot(v,p_v,'-r','LineWidth',2);
hold off;
title('Eigenvector v_{1} of C_{s}');
xlabel('Eigenvector components v');
ylabel('Probability Density \rho(v)');
legend('Eigenvector v_{1}','Pdf');
set(gcf,'color','w');
grid on;

figure(10);
histn(-1*sqrt(m)*(v_2),-4,0.1,4);
hold on;
plot(v,p_v,'-r','LineWidth',2);
hold off;
title('Eigenvector v_{2} of C_{s}');
xlabel('Eigenvector components v');
ylabel('Probability Density \rho(v)');
legend('Eigenvector v_{2}','Pdf');
set(gcf,'color','w');
grid on;

figure(11);
histn(-1*sqrt(m)*(v_3),-4,0.1,4);
hold on;
plot(v,p_v,'-r','LineWidth',2);
hold off;
title('Eigenvector v_{3} of C_{s}');
xlabel('Eigenvector components v');
ylabel('Probability Density \rho(v)');
legend('Eigenvector v_{3}','Pdf');
set(gcf,'color','w');
grid on;

%% Question 5b.
IPR_W = sum(W_eigenVector_Sorted.^4);
figure(12);
loglog(W_eigenValue_Sorted,IPR_W,'-o');
set(gcf,'color','w');
title('IPR of W');
xlabel('Eigenvalue');
ylabel('Inverse Participation Ratio');
xlim([10^-2 10^2]);
ylim([10^-3 10^0]);
grid on;

%% Question 5b. 
IPR_Cs = sum(Cs_eigenVector_Sorted.^4);
figure(13);
loglog(W_eigenValue_Sorted,IPR_W,'-o');
hold on;
loglog(lambdaSorted,IPR_Cs,'-o');
hold off;
set(gcf,'color','w');
title('IPR of W and IPR of C_s');
xlabel('Eigenvalue');
ylabel('Inverse Participation Ratio');
legend('IPR(W)','IPR(C_s)');
xlim([10^-2 10^2]);
ylim([10^-3 10^0]);
grid on;

lambdaSorted_2 = lambdaSorted(2:end);
lambdaSorted_2Norm = lambdaSorted_2*(sum(lambdaSorted)/sum(lambdaSorted_2));
IPR_Cs_2 = sum(Cs_eigenVector_Sorted(:,2:end).^4);
figure(14);
loglog(W_eigenValue_Sorted,IPR_W,'-o');
hold on;
loglog(lambdaSorted_2Norm,IPR_Cs_2,'-o');
hold off;
set(gcf,'color','w');
title('IPR of W and IPR of C_s (Removal of 1 Eigenvalue & Eigenvector)');
xlabel('Eigenvalue');
ylabel('Inverse Participation Ratio');
legend('IPR(W)','IPR(C_s)');
xlim([10^-2 10^2]);
ylim([10^-3 10^0]);
grid on;

%% Question 6

XsCap_shuffle = zeros(m,n);

for i = 1:n
    XsCap_shuffle(:,i) = Xs(randperm(m),i);
end

sigmaCap_s_shuffle = (1/n)*(XsCap_shuffle)*(XsCap_shuffle');
eigSigmaCap_s_shuffle = eig(sigmaCap_s_shuffle);
DCap_shuffle = diag(diag(sigmaCap_s_shuffle));
CsCap_shuffle = (DCap_shuffle^(-1/2))*(sigmaCap_s_shuffle)*(DCap_shuffle^(-1/2));

[CsCap_shuffle_eigenVector,CsCap_shuffle_eigenValue] = eig(CsCap_shuffle);
[CsCap_shuffle_eigenValue_Sorted,CsCap_shuffle_eigenValue_index]=sort(diag(CsCap_shuffle_eigenValue),'descend');
%Cs_shuffle_eigenValue_Sorted_matrix = Cs_shuffle_eigenValue(Cs_shuffle_eigenValue_index,Cs_shuffle_eigenValue_index);
CsCap_shuffle_eigenVector_Sorted = CsCap_shuffle_eigenVector(:,CsCap_shuffle_eigenValue_index);

figure(15);
histn(CsCap_shuffle_eigenValue_Sorted,0,0.1,50);
hold on;
plot(d,P_rcm,'-r','LineWidth',2);
hold off;
title('Eigenvalues of C_{shuffle}');
xlabel('Eigenvalue (\lambda)');
ylabel('Probability Density P(\lambda)');
legend('Eigenvalues','P_{rcm}');
set(gcf,'color','w');
grid on;

