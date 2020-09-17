function phi = orthopoly_evaluate(T,x)
% Evaluate the Jacobi polynomial in x

[H, Lambda] = eig(T); Lambda = diag(Lambda); Hn = H(end,:)';
for m=1:length(x)
    xi = x(m);
    u = H*(Hn./(Lambda-xi));
    phi(:,m) = u/u(1);
end
end