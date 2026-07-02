function [a,b] = fit_GK(x,y)
% Define the power function model: y = a * x^b
powerLawModel = @(params, x) params(1) * x .^ params(2);

% Define the error function to minimize sum of squared differences
error_fn = @(params) sum((y - powerLawModel(params, x)).^2);

% Initial guess for parameters [a, b]
params0 = [1e5, -0.5];

% Optimize parameters using fminsearch
params_opt = fminsearch(error_fn, params0);

% Extract fitted parameters
a = params_opt(1);
b = params_opt(2);