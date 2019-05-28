function u_DTV=DTV(b,N,th,lambda)



A  = 1; % forward operator for denoising problem
AT = 1; % adjoint of forward operator for denoising problem

%% solve L2-DTV problem
% estimate angle from noisy image
%lambda  = 0.18;     % regularization parameter
Ni      = 50;      % maximum number of iterations
tol     = 1e-6;     % stopping criterion tolerance
a       = 0.15;     % width parameter for DTV

[u_DTV,on_DTV,sn_DTV,n_DTV] = l2_DTV(b(:),A,AT,N,lambda,th,a,Ni,tol);