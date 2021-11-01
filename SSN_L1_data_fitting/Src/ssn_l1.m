function [x,pp,ithist] = ssn_l1(K,y,alpha,xe,opts)
% [x,pp,ithist] = SSN_L1(K,y,alpha,xe,opts)
% semi-smooth Newton method for L^1-fitting-L^2-regularization using a
%   path-following method for the regularization of dual problem
% Input:  K      -- operator to be inverted
%         y      -- given data
%         xe     -- exact solution (to compute error)
%         alpha  -- regularization parameter for primal problem
%         opts   -- options for semismooth Newton method (see below)
% Output: x      -- solution of primal problem (L^1 reconstruction)
%         pp     -- solution of dual problem
%         ithist -- iteration history in path-following method (see below)
% JIN Bangti(kimbtsing@yahoo.com.cn)
% Christian Clason(christian.clason@uni-graz.at)
% April 14, 2009

%% Set parameters
N = size(K,2);

% Default parameters
if nargin <= 4 || isempty(opts)
    beta    = 1;              % initial value for smoothing parameter
    q       = 1/5;            % reduction factor
    c       = 1e9;            % regularization parameter of box constraints
    maxit   = 10;             % maximum number of iterations
    dispopt = 0;              % plot iterates in SSN iteration
else
    beta    = opts.beta0;      % initial value for smoothing parameter
    q       = opts.q;          % reduction factor
    c       = opts.c;          % regularization parameter of box constraints
    maxit   = opts.maxit;      % maximum number of iterations
    dispopt = opts.dispopt;    % plot iterates in SSN iteration

end

% H^1-smoothing: Laplace with homogeneous Dirichlet conditions
e=ones(N,1); h=1/N;
DD=1/h^2*spdiags([-e 2*e -e],[-1 0 1],N,N);

KtK = K*K';         % precompute

pp   = zeros(N,1);  % initial guess
Ip_t = zeros(N,1);  % initialize active sets
In_t = zeros(N,1);

ithist = [];

%% path-following method for beta
while (beta>1e-16)
    update=1e99;
    pold=pp;
    M0 = 1/alpha*KtK + beta*DD; % precompute

    %% semismooth Newton iteration
    for iter=1:maxit
        % find active sets
        Ip = (pp>1); In= (pp<-1); % (upper, lower bound)

        % compute dual variable
        M  = M0+c*spdiags(Ip+In,0,N,N);
        pp = M\(y+c*(Ip-In));

        if dispopt
            % plot dual, primal iterate
            x  = (K'*pp)*1/alpha;
            figure(3), plot(1:N,pp,'r')
            figure(4), plot(1:N,xe,'k',1:N,x,'r');xlabel('ssn');drawnow
        end

        % terminate inner loop if active sets haven't changed
        update = nnz(Ip-Ip_t)+nnz(In-In_t);
        if update == 0
            break
        end

        Ip_t = Ip;
        In_t = In;
    end
    %% end ssn

    if nargout == 3
        % calculate primal iterate
        x      = (K'*pp)*1/alpha;
        % save iteration history for path-following method
        err    = norm(x-xe)/sqrt(N);  % L^2 norm of error
        res    = norm(K*x-y,1)/N;     % L^2 norm of residual
        dfval  = 1/(2*alpha*N)*norm(K'*pp)^2-pp'*y*1/N; % dual funct.val.
        pfval  = res+alpha/(2*N)*norm(x)^2;             % primal funct.val.
        dnorm  = 1/2*pp'*(DD*pp);     % H^1 seminorm of dual variable
        ithist = [ithist;beta iter update err res dfval pfval dnorm]; %#ok<AGROW>
    end

    % If SSN iteration terminated without a (nearly) feasible solution
    if (iter == maxit) && max(abs(pp)>2)
        % return iterate for last beta which gave feasible pp
        pp=pold;
        break
    end

    % reduce beta
    beta = beta*q;
end

%% compute primal solution
x=(K'*pp)*1/alpha;

save ssn_l1
