% nonlinear_exp.m
% This function fits a tree ring time series to the non-linear
% equation used by Ed Cook in his ARSTAN program.

% Function written Jan 12, 2004.
% Function last revised Jan 29, 2004.

function qq = nonlinear_exp(beta,x)

% Assign parameters from beta vector.
a=beta(1);
b=beta(2);
d=beta(3);
qq = a*exp(-b*x)+d;