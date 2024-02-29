function s = SISSY(x,A,T,lambda,alpha,Niter)
[N,D] = size(A);
E = size(T,1);
rho = 1;
s = zeros(D,1);
z = zeros(E,1);
y = zeros(D,1);
u = zeros(E,1);
v = zeros(D,1);

P = sparse(rho*(T.'*T+speye(size(T,2))));
APi = A/P;
L = chol(eye(size(A,1))+APi*A.','lower');
s1 = A.'*x;

for i=1:Niter
    b = s1+rho*(T.'*(z+u/rho)+y+v/rho);
    s = P\b - APi'*(L'\(L\(APi*b)));
    z = prox_op(T*s-1/rho*u,'L1',lambda/rho);
    y = prox_op(s-1/rho*v,'L1',lambda*alpha/rho);
    u = u + rho*(z-T*s);
    v = v+rho*(y-s);
end

s;

end