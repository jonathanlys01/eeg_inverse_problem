clear;
close all;
clc;

load TP_data;

%generate linear mixture of source signals
Xs=G*S;

%determine maximum of the signal of interest (here an epileptic spike) to
%apply source loclization algorithms to this time point in the following
[~,id]=max(mean(S,1));

%visualize original source distribution
figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),S(:,id));
title('original source configuration: two source regions','FontSize',18); axis off;

%generate Gaussian random noise
Noise=randn(size(Xs));

%normalize noise
Noise=Noise/norm(Noise,'fro')*norm(Xs,'fro');

%signal to noise ratio
SNR=10;

%generate noisy data according to given SNR
X=Xs+1/sqrt(SNR)*Noise;

%visualize data (for a reduced number of sensors whose indices are
%specified by idx_electrodes)
plot_eeg(X(idx_electrodes,:),max(max(X(idx_electrodes,:))),256,channel_names);
title('noisy EEG data','FontSize',18);

% ------------ TP1 -----------------------------------
%% Gibbs sampler

% Shat = Gibbs_sampler(X(:,id),A);  %

new_id = id;
% 
[Shat, Lhat] = Gibbs_sampler(X(:,new_id),G, SNR);
disp(Lhat);
figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),Shat);
title("Gibbs sampler");

% ------------ end of TP1 ----------------------------



% -----------  TP2 / 3 -------------------------------

%% manual test of lambda
lambda = logspace(0,3,6);
% lambda: 10^0 to 10^3

Shat = zeros(size(S,1),length(lambda));
for k=1:length(lambda)
    Shat(:,k)=MNE(X(:,id),G,lambda(k));  % function to be implemented
    figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),Shat(:,k));
    title(strcat('\lambda=',num2str(lambda(k))));
end

%% Robustness to SNR

Shat=MNE(X(:,id),G,461); 
figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),Shat);
%title(strcat('SNR=',num2str(SNR)));

%% L-curve criterion

lambda=logspace(-10,10,20);

Shat = zeros(size(S,1),length(lambda));
for k=1:length(lambda)
    Shat(:,k)=MNE(X(:,id),G,lambda(k)); 
end

n = sum( (X(:,id)-G*Shat).^2 , 1)
s = sum (Shat.^2 , 1)

figure; loglog(n,s,"-o")


%% discrepancy principle

noise_gt = sum(Noise(:,id).^2);


semilogy(sum (Shat.^2, 1), "-o")
hold on
semilogy(ones(length(lambda),1)*noise_gt)
hold off



%% Generalized Cross-Validation
%  code to be added in the L-curve iteration

GCV = zeros(length(lambda),1);

for k=1:length(lambda)
    q_ = (trace(eye(size(G,1)) - G*G' * pinv(G*G'+lambda(k)*eye(size(G,1)))))^2;
    GCV(k) = n(k) / q_;
end

plot(GCV)
%% Comparison of the methods

idx_L = 8;
idx_DP = 6;
idx_GCV = 12;

disp(norm(S(:,id) - MNE(X(:,id), G, lambda(idx_L))));
disp(norm(S(:,id) - MNE(X(:,id), G, lambda(idx_DP))));
disp(norm(S(:,id) - MNE(X(:,id), G, lambda(idx_GCV))));


%% SISSY
disp("SISSY");
T=variation_operator(mesh,'face');
lambda=logspace(-2,3,6); % lambda from 0.01 to 1000
Shat = zeros(size(S,1),length(lambda));
MaxIter = 60;
for k=1:length(lambda)
    Shat(:,k)=SISSY(X(:,id),G,T,lambda(k),0.1, MaxIter);  %  function to be implemented
    figure; trisurf(mesh.f,mesh.v(:,1),mesh.v(:,2),mesh.v(:,3),Shat(:,k));
    title(strcat('\lambda=',num2str(lambda(k))));
end
