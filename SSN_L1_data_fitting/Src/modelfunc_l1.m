function [x,fphist] = modelfunc_l1(K,y,xe,model_opt,ssn_opt)
% [x,fphist] = MODELFUNC_L1(K,y,alpha0,maxit,dispopt,xe)
% model function approach for optimal regularization parameter choice in
%   L^1 fitting with L^2 regularization
% Input:  K         -- operator to be inverted
%         y         -- given data
%         xe        -- exact solution (to compute error)
%         model_opt -- options for fixed point iteration (alpha), see below
%         ssn_opt   -- options for path-following iteration (beta)
% Output: x         -- solution of L^1 fitting problem
%         fphist    -- history of iterates in fixed point iteration
% JIN Bangti(kimbtsing@yahoo.com.cn)
% Christian Clason(christian.clason@uni-graz.at)
% April 14, 2009

%% Set parameters
N = size(K,2);

% Default parameters
if nargin <= 3
    alpha   = 1e-2;                % initial guess for parameter alpha
    sigma   = 1.05;                % relaxation
    maxit   = 10;                  % maximum number of iterations
    dispopt = 0;                   % plot iterates in fixed point iteration
else
    alpha   = model_opt.alpha0;    % initial guess for parameter alpha
    sigma   = model_opt.sigma;     % relaxation
    maxit   = model_opt.maxit;     % maximum number of iterations
    dispopt = model_opt.dispopt;   % plot iterates in fixed point iteration
    if nargin <= 4
        ssn_opt = [];              % use default parameters in SSN
    end
end

fphist = [];

%% Fixed point iteration for optimal alpha

b = norm(y,1)/N;

for i = 1:maxit
    if nargout >= 2
        [x,pp,ithist] = ssn_l1(K,y,alpha,xe,ssn_opt);
    else
        x = ssn_l1(K,y,alpha,xe,ssn_opt);
    end
    
    if dispopt
        figure(2),plot(1:N,x,'r',1:N,xe,'k');drawnow
    end
    
    F  = norm(K*x-y,1)/N + alpha/2*norm(x)^2*1/N; % functional value
    Fp = 1/2*norm(x)^2*1/N;          % derivative of functional
    t  =  (b-F)/Fp - alpha;          % parameter for model function
    c  = -(F-b)^2/Fp;                % parameter for model function
    deltaest = 1*(b + c/t);          % estimate of noise level: m(0)
    
    if nargout >= 2
        fphist = [fphist;deltaest alpha ithist(end,4)]; %#ok<AGROW>
    end
    
    intcpt = F-alpha*Fp;             % y-intercept
    alpha  = c/(intcpt*sigma-b)-t;   % estimate of opt. reg. param.
    
    if dispopt
    display(sprintf('iter = %d, deltaest = %g, alpha = %g',...
        i,deltaest,alpha))
    end
end

save modelfunc_l1
