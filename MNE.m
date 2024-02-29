function SOut=MNE(X,A,lambda)

N = size(A,1);

SOut = A' * pinv(A*A' + lambda * eye(N)) * X;
end