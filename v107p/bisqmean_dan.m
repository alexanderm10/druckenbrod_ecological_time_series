function [ymn,varyh,df,w,ybar,se]=bisqmean_dan(y)
%
% Biweight mean for a vector of numbers.
% Last revised 2011-7-09
% Revised by Dan Druckenbrod 2012-1-11
%
% Source:  Mosteller and Tukey (1977, p. 205, p 351-352)
%		Cook and  Kairiukstis (1990, p. 125-126)
%
%
%****************  INPUT *************************
%
% y (? x 1)r  vector of data -- say, indices for ? cores in a year
% 
%
%********************  OUTPUT ************************
%
% ymn (1 x 1)r  biweight mean
% varyh (1 x 1)r   asymptotic standard dev of biweight mean - p. 208, 
%		third eqn from top of page
% df (1x1)r degrees of freedom
% w (? x 1)r  final weights on values in y
% ybar (1 x 1)r arithmetic mean corresponding to ymn
% se (1 x 1)r  standard error of ybar
%
%********************  NOTES *********************
%
% ybar and se just included in debugging  to double check
% on closeness of ybar to ymn, se to sqrt(varyh)
%
%*******************************************************************


sens = 0.001;  % hard coded theshold of sensitivity for stopping iterat
nits = 100;  % max number of allowed iterations

[n,ny]=size(y);
if ny > 1;
	error('y should be a vector')
end

if any(isnan(y));
    error('y not permitted to have NaNs');
end;

if n<6; % if fewer than 6 sample size, use median
    ymn = median(y);
    w=[];
    ybar=mean(y);
    se= sqrt(var(y)/n); % standard error of mean
    df=[];
    varyh=NaN;
    return;
end;



ww = 1/n; % weight for even average
ybar = mean(y); % arith mean
%ybar=median(y);
se= sqrt(var(y)/n); % standard error of mean

nz=0;
ymn = ybar; % initial biweight mean as arith mean

for i = 1: nits;  % iterate max of nits times
	ymnold = ymn;  % store old value of mean
	e = y-ymn; % deviations from mean
	S = median(abs(e));  % median abs deviation
	u = e / (6*S);  % scaled deviations

	w = (1 - u.^2).^2;  % compute weights
	L1 = abs(u)>=1; % flag huge errors
	L1s = sum(L1);
	if L1s>0
		nz=0;
		nz= nz(ones(L1s,1),:);
		w(L1)=nz;  % set weights on those obs to zero
	end
	w = w / sum(w); % adjust weights to sum to 1.0
	
	ymn = sum(w .* y); % compute biweight mean


	% Variance of estimate of biweight mean
	ui= e / (9*S);
	L2 = ui>1;
	ui(L2)=[];
	z =y(~L2);
	nz = length(z);
	nom1 = (z - ymn) .^2;
	nom2 = (1-ui .^2) .^4;
	nom = sum(nom1 .* nom2);

	den1 = sum((1-ui .^2) .* (1-5*ui .^2));
	
    varyh_hoaglin=(n^0.5)*(nom^0.5)/den1; % Dan: p. 417 3rd equation
    den2 = -1 + sum ((1-ui .^2) .* (1-5*ui .^2));
	% varyh = nom / (den1*den2);  % variance of biweight mean 
	% last eqn, p. 208
    
	varyh = n^.5*nom^.5 / ((den1*den2)^.5); % Dan: p.417 Kafadar approach
    
    df = 0.7 * (nz -1); % degrees of freedom 


   % if little change in mean, exit loop
	if abs (ymn - ymnold) < sens
		return
	end
end

