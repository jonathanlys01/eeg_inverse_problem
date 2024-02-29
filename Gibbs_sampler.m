function [SOut,LambdaOut]=Gibbs_sampler(X,A, SNR)

%number of iterations
Niter=10;

%get number of sensors and number of dipoles
[N,D]=size(A);

% constants 
nA = sum(A.^2,1); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% implement the Gibbs sampler here according to the pseudocode
%
% store vectors q and s of each iteration in matrices Q (D x niter) and S
% (D x niter)
% use variables sigma_n2 and sigma_s2 for the variances of noise and
% signals

sigma_s2 = 0.3077; % measured in the data
sigma_n2 = SNR*0.3077; % SNR = 1
alpha = 1;
beta = 1;
Q = zeros(D,Niter);
S = zeros(D,Niter);


s = zeros(D,1);


lambda = betarnd(alpha,beta);

for j=1:Niter
    e = X - A*s;
    
    for i=1:D
        Ai = A(:,i);
        ei = e + Ai*s(i);
        
        sigma2i = sigma_n2*sigma_s2/(sigma_n2 + sigma_s2*norm(Ai)^2);
        
        
        mui = sigma2i/sigma_n2*Ai'*ei;
       
        
        if mui^2/(2*sigma2i) > log(realmax("double"))
            lambdai = 0.999999;
        else   
            nui = lambda*sqrt(sigma2i/sigma_s2)*exp(mui^2/(2*sigma2i));
           

            lambdai = nui/(nui+1-lambda);
        end
        
        qi = binornd(1,lambdai);
        
        if qi
            S(i,j) = mui+sqrt(sigma2i)*randn();
        else
            S(i,j) = 0;
        end
        
        Q(i,j) = qi;
        e = ei - Ai*s(i);
    end
    
    L = sum(Q(:,j));
    lambda = betarnd(alpha+L,beta+D-L);
     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp("Make the estimation from Q and S using the MAP criterium")
%Q(:,Niter/2:Niter)

cout = zeros(Niter/2,1);

for j = 1:Niter/2
    disp(j)
    q = Q(:,Niter/2+j); 
    idx = find(q==1); 
    R = (A(:,idx)'*A(:,idx))/sigma_n2 + eye(length(idx))/sigma_s2; 
    S_opt = (R\(A(:,idx)'*X))/sigma_n2; 
    cout(j) = - S_opt'*R*S_opt/sigma_n2; 
end

[Val, OptIdx] = min(cout); 
q = Q(:,Niter/2+OptIdx);
idx = find(q==1); 
R = (A(:,idx)'*A(:,idx))/sigma_n2 + eye(length(idx))/sigma_s2; 
SOut = zeros(D,1); 
SOut(idx) = (R\(A(:,idx)'*X))/sigma_n2;
LambdaOut = length(idx)/D; 



