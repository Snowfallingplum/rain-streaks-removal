function [u,on,sn,n] = l2_DTV(b,A,AT,N,lambda,th,a,Ni,tol)
% Solve the L2-DTV problem:
% min_u J(u) := 1/2*|Au - b|_2^2 + lambda*DTV(u),
% using the Primal-Dual-Hybrid-Gradient method by Chambolle and Pock.
%
% Input
% b:        data (e.g. noisy) stacked in column vector
% A:        system matrix/forward operator size N^2 x length(b)
% AT:       transpose of system matrix
% N:        width and height of 2D function u (assumed to be square)
% lambda:   regularization parameter. i.e. should be non-negative scalar
% th:       angle for DTV regularizer in degrees
% a:        width parameter for DTV, should be in ]0,1]
% Ni:       maximum number of iterations
% tol:      tollerance for stopping criterion (relative change of J(u)
%
% Output
% u:        2D function solution, size N x N
% on:       relative change of J(u) for each iteration
% un:       relative change of u for each iteration
% n:        number of used iterations
%
% By Rasmus Dalgas Kongskov, 18/09/2015, DTU Compute

% finite diffrence derivative approximation matrix
Df = sparse(-eye(N)+diag(ones(N-1,1),1)); Df(end,end)=0; % forward diff.
Db = sparse(eye(N)-diag(ones(N-1,1),-1)); Db(1,1)=0; % backward diff.
Dx = kron(Df,eye(N));                           % derivative in x-direction
if th>=0
    Dy = kron(eye(N),Df);                       % derivative in y-direction
else
    Dy = kron(eye(N),Db);                       % derivative in y-direction
end

% power method to estimate |K|
ux = ones(N^2,1);
for n = 1:50
    ux = AT*(A*ux);
    ux = ux/norm(ux);
end
nA = norm(A*ux);

% scale to operator-norm 1
K  = @(ui) double(1/nA*A*ui);
KT = @(qi) double(1/nA*AT*qi);
b  = 1/nA*b;

% parameters from Chambolle-Pock
L     = 1/2*(sqrt(32) + 18);
tau   = 1 * 1/sqrt(L*1.01);
sigma = 1 * 1/sqrt(L*1.01);

% initialization
u = zeros(N^2,1); ub = u;
v = zeros(N^2,2);
q = b*0;

on = zeros(Ni,1); sn = on;

% Ellipse gradient and divergence projection operators
th = th/180*pi;
c   = cos(th);
s   = sin(th);

Exf  = @(uin) c*Dx*uin - s*Dy*uin;
Eyf  = @(uin) a*(c*Dy*uin + s*Dx*uin);

Exb  = @(pin) -Dx'*(c*pin) + Dy'*(s*pin);
Eyb  = @(pin) a*(-Dy'*(c*pin) - Dx'*(s*pin));

for n = 1:Ni
    
    % step in v
    v(:,1) = v(:,1) + sigma*Exf(ub);
    v(:,2) = v(:,2) + sigma*Eyf(ub);
    ma = max(1,sqrt(v(:,1).^2 + v(:,2).^2)/lambda);
    v(:,1) = (v(:,1))./ma;
    v(:,2) = (v(:,2))./ma;
    
    % step in q
    q = (q + sigma*(K(ub)-b))/(1+sigma);
    
    % save current u
    ub = u;
    
    % step in u
    u = u + tau*(Exb(v(:,1)) + Eyb(v(:,2))) - tau*KT(q);
    
    % objective-function- and solution-change
    on(n) = nA/2*norm(K(u)-b)^2+lambda*sum(sqrt((Exf(u)).^2+(Eyf(u)).^2));
    sn(n) = norm(u-ub)/norm(ub);
    
    % step in ub
    ub = 2*u - ub;
    
    % stop if relative objective function change is below tol
    if n>1 && (abs(on(n)-on(n-1))/on(n-1))<tol;
        break
    end
    
end

u  = reshape(u,N,N);
on = on(1:n); sn = sn(1:n);
on = abs(diff([2*on(1);on]))./on;

end