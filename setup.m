% Set global parameters
Ne = 72; % Ensemble size
Nt = 1024; % Number of assimilation cycles
locRad = 16; % Localization radius

% Set up localization
x = [0:128 -127:-1];

% Initialize ensemble
load ~/MATLAB/Wind_DA/Reference/qp_ref_01.mat
rng(0); % For reproducibility
tmp = randperm(size(qp_ref,3));
Q = qp_ref(:,:,tmp(1:Ne));
clear tmp

% Get obs
%load obs_grid.mat % or obs_rand.mat

% setup
k = [0:128 -127:-1];
[KX,KY] = meshgrid(k,k);

% Prep likelihood
theta_grid = linspace(-pi,pi,257)';
l = zeros(256,1);
tmp = linspace(-1,2,25);
tmp = tmp(2:end);
tmp = ppval(bspline([-1 0 1 2]),tmp);
l(1:13) = tmp(12:end);
l(end-10:end) = tmp(1:11);
l = circshift(l,[127 0])/max(l);
clear tmp
