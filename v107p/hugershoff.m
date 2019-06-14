% hugershoff.m
% This function fits a tree ring time series to the growth trend equation
% developed by Warren (1980) TRR.

% Function written Nov 12, 2013.
% Function last revised Nov 12, 2013.

function qq = hugershoff(beta,x)

% Assign parameters from beta vector.
a=beta(1);
b=beta(2);
c=beta(3);
k=beta(4);
qq=a*((x).^b).*exp(-c*(x))+k;
% qq=a*((x+1).^b).*exp(-c*(x+1))+k;
%q=log(a)+b*log(x)%-c*x;
%qq=exp(