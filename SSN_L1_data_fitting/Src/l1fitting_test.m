% L1FITTING_TEST
% test script for L^1 data fitting using automatic parameter choice and
% semismooth Newton method.
%
% JIN Bangti(kimbtsing@yahoo.com.cn)
% Christian Clason(christian.clason@uni-graz.at)
% April 14, 2009

%% Parameters
n     = 100;       % problem size
d_mag = 1.0;       % magnitude of noise
d_per = 0.3;       % percentage of corrupted data points

% ... for model function method
model_opt.alpha0  = 1e-2;     % initial guess for parameter alpha
model_opt.sigma   = 1.05;     % relaxation
model_opt.maxit   = 10;       % maximum number of fixed point iterations
model_opt.dispopt = 0;        % output iterates in fixed point iteration

% ... for semismooth Newton method
ssn_opt.beta0     = 1;        % initial guess for parameter beta
ssn_opt.q         = 1/5;      % reduction factor for beta
ssn_opt.c         = 1e9;      % regularization parameter for box constraint
ssn_opt.maxit     = 10;       % maximum number of semismooth Newton steps
ssn_opt.dispopt   = 0;        % plot iterates in SSN iteration

%% Data

% operator K, exact data ye, reference solution xe
% deriv2 from Regularization Tools: http://www2.imm.dtu.dk/~pch/Regutools/
[K,ye,xe] = deriv2(n,3);

% add impulse noise to data
drnd   = rand(n,1);
ind    = find(drnd<d_per);
imprnd = zeros(n,1);
imprnd(ind) = randn(length(ind),1);
noise  = d_mag*max(abs(ye))*imprnd;

y = ye+noise;                   % noisy data

% show noisy data
figure(1), plot(1:n,ye,'k',1:n,y,'r');legend('exact','noisy');

% compute exact noise level
delta = norm(y-ye,1)/n;

%% Calculate minimizer
[x,fphist] = modelfunc_l1(K,y,xe,model_opt,ssn_opt);

err      = norm(x-xe)/sqrt(n); % L^2 norm of reconstruction error
deltaest = fphist(end,1);      % estimate of noise level
alpha    = fphist(end,2);      % estimate of optimal parameter

display(sprintf('delta = %e, deltaest = %e, alpha = %e, err = %e', ...
     delta, deltaest, alpha, err));

figure(2),plot(1:n,x,'r',1:n,xe,'k');h=legend('L1 reconstruction','exact');
