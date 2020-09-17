function phi = compute_jacobi(n,a,b,x)
% x is a column vector
% n is an integer. 
% a and b are alpha and beta in the assignment
% The function compute_jacobi(n,a,b,x) computes all the Jacobi polynomials of order
% 0, 1, 2, ..., n evaluated at x. It returns a matrix phi of size (N+1, length(x)) 
% where phi(i,j) stores the value of the i-th Jacobi polynomial evaluated in x(j)
% These Jacobi polynomials are orthogonal with respect to
% J_w^{(alpha,beta)}, and they are already normalized.
% This function calls the function orthopoly_evaluate.m which evaluates the
% polynomials in x

B = diag(sqrt((a+1+(0:n))./((a+b)+(2:2:(2*n+2)))).*sqrt((a+b+1+(0:n))./((a+b)+1+(0:2:(2*n)))));
B = B - diag(sqrt((b+(1:n))./((a+b)+(2:2:(2*n)))).*sqrt((1:n)./(a+b+1+(2:2:(2*n)))),1);
c = factorial(a)*factorial(b)/factorial(a+b+1);
T = B'*B;
phi = orthopoly_evaluate(T,x);

phi = phi/sqrt(c); % Normalization: \int_0^1 J_w^{(alpha, beta)} P^{(alpha,beta)}_n(x)^2 dx = 1


