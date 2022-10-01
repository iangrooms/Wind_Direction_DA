function [p,R] = get_prior_pdf(theta)
% Uses kernel density estimation with a von Mises kernel to estimate the
% pdf associated with the sample contained in eta. The von Mises kernel is 
% truncated to avoid extremely small probability densities.

% Get an ad hoc estimate of the kernel bandwidth
R = abs(mean(exp(1i*theta)));
kappa = -.445*(length(theta)^0.4)/log(R);

% Now construct KDE estimator
t = linspace(-pi,pi,257)';
t = t(2:end);
theta = reshape(theta,1,[]);
p = sum(exp(kappa*cos(t-theta)),2);

% Normalize; quadrature using trapezoid rule
p = p / ( (t(2)-t(1)) * sum(p) );

% Truncate and normalize again
p = max(1E-8,p);
p = p / ( (t(2)-t(1)) * sum(p) );