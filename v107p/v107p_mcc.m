% v107p_mcc.m
% This function extracts a single tree-ring time series from 
% ringwidth_import.m and places it in a vector for time series analysis.
% The function power transforms and removes the mean to create transformed
% residuals.  The function then detrends with an iterative neg. exponential
% fit, or if that does not fit or fails to find a solution, then a linear 
% regression with either a positive or negative slope is fit to the data.  
% Using the maximum entropy model solution otherwis known as the Burg 
% method, the autoregressive model that is the best fit for the series is 
% determined.  Using the best fit model, the function searches for 
% autoregressive outliers iteratively.  These outliers may either be pulse 
% events (1 yr) or CSTs (> minimum no. of yrs).  After the first pass,
% the outliers are removed and the series is reconstituted.  The best ar 
% order is then redetermined and the function searches for additional 
% outliers.  The # of iterations is set by the user (8 should be enough).  
% This version uses a power transformation to minimize
% the heteroscedastic nature of my time series. 'fig' is a flag that 
% specifies whether you want a figure (=1) or not (=0).  Missing years are
% set to the average of neighboring rings.  The central limit theorem 
% is used to search the residuals for trend outliers.  This version also
% uses Dave Meko's biweight mean code and currently runs with a window of 
% 9 to 30 yrs.  Estimated values fpr missing rinngs are removed in the 
% output series.  This version uses a modified Hugershoff curve with a 
% potentially nonzero asymptote to detrend + and - disturbance events. 
% It also returns the transformed standardized series.

% Function written Sep 10, 2002.
% Function last revised Sept 25, 2015.

function [YEARS,transformed,detrended,St,Str,Dtr,Atr,age,outs]=...
    v107p_mcc(core,fig,iter,filename,scale,kind,years,col_header,rings)
global PARAM; PARAM=0; % vector of parameters for best order AR model.
global ORDER; ORDER=0; % best order of AR model determined by AIC.
global YEARS; YEARS=0; % calendar years of tree growth from datafile

% Import tree-ring data (returns vars *col_header* and *rings*)
% [col_header,rings]=dc_ringwidth_import(0,filename);

% Find pointer to start and end of series
sos=find(rings(:,(core))>0, 1);
eos=find(rings(:,(core))>0, 1, 'last');

% Assign years and raw widths to respective vectors.
YEARS=years(sos:eos,1); %#ok<*NASGU>
raw=rings(sos:eos,core);

disp(['Core: ' char(col_header(core))])
nyrs=length(YEARS);
disp(['Total no. of measured years: ' int2str(nyrs)])
disp(['First year is ' num2str(YEARS(1))])
disp(['Last year is ' num2str(YEARS(nyrs))])

% Estimate missing ring widths using mean of neighboring rings
mss=NaN(length(raw),1);

if find(raw==0)
   m1=find(raw==0);
   disp(['Missing rings at years ' num2str([YEARS(m1)'])])
   for nm=1:length(m1)
       prior=mean(raw(find(raw(1:m1(nm)),1,'last')));
       subs=mean(raw(find(raw(m1(nm):length(raw)),1,'first')+m1(nm)-1));
       mss(m1(nm))=mean([prior subs]);
   end
   raw=nansum([raw mss],2);
end

% Power transformation.
fdiff=0;
for x=1:(length(YEARS)-1) % Calculate 1st differences
    fdiff(x,1)=raw(x+1);
    fdiff(x,2)=abs(raw(x+1)-raw(x));
end

s=1;
for q=1:(length(YEARS)-1)
    if (fdiff(q,1)~=0) && (fdiff(q,2)~=0)
        nz_fdiff(s,:)=fdiff(q,1:2);% non-zero ring widths
        s=s+1;
    end
end
log_fdiff=[log(nz_fdiff(:,1)) log(nz_fdiff(:,2))];
 
X=[ones(length(log_fdiff(:,1)),1) log_fdiff(:,1)];
bb = regress(log_fdiff(:,2), X);
optimal_line = bb(2)*log_fdiff(:,1)+bb(1);

optimal_pwr = 1-bb(2);
disp(['Optimal Power = ' num2str(optimal_pwr)])
if optimal_pwr <= 0.05
    transformed=log10(raw);
    tzero=log10(0.001);
    disp('Series was log10 transformed')
elseif optimal_pwr>1
    optimal_pwr=1;
    transformed=(raw.^(optimal_pwr));
    tzero=0.001.^(optimal_pwr);
    disp('Series was power transformed with power =1')
else
    transformed=(raw.^(optimal_pwr));
    disp(['Series was power transformed with power = ' ...
        num2str(optimal_pwr)])
    tzero=0.001.^(optimal_pwr);
end
transm=mean(transformed);
% disp(['tzero = ' num2str(tzero)])

% tresids=transformed-transm; % transformed residuals cannot be used with
% iterative age detrending since half of residuals are negative.


% warning('OFF',stats:nlinfit:IterationLimitExceeded);
% warning('OFF',stats:nlinfit:IllConditionedJacobian);
% Nonlinear detrending option.
% Function nlinfit employs nonlinear least squares data fitting by the
% Gauss-Newton Method.
crashed=zeros(nyrs,1);
wlngth=zeros(nyrs,1);
trendtype=0; % Neg exp = 1, neg linear reg = 2, or pos linear reg = 3 
minyr=30; % minimum # of yrs to fit to nlinfit
if minyr>nyrs
    disp('Insufficient # of years to fit minimum nonlinear age trend.')
end
b=zeros(nyrs,3);
mse=NaN(nyrs,1);
warning off
for i=minyr:nyrs
    try
        lastwarn('') 
        beta = [.5 .1 1];
        xyrs = 1:i; % set years from 1 to length of series
        [b(i,1:3),~,~,~,mse(i)]=nlinfit(... 
            xyrs(1:i),transformed(1:i)','nonlinear_exp',beta);
        crashed(i)=1;
        msgstr = lastwarn;
        % [msgstr, msgid] = lastwarn;
        wlngth(i)=length(msgstr);
        % disp(['Variance of the error term: ' num2str(mse(i))])
    catch % Stops code from crashing because of problems fitting exp curve
        crashed(i)=2;
        % disp('Exponential fit function failed, reverting to linear regression.')
    end
end

warning on
% warning('ON',stats:nlinfit:IterationLimitExceeded);
% warning('ON',stats:nlinfit:IllConditionedJacobian);
i_c=0;

% Dissallow curve to be concave up and make sure nlinfit
% converges by making b(2) sufficiently large.
% constant b(3) must be >=0 in original mm
i_c=find(crashed==1 & b(:,1)>=0 & b(:,2)>0.001 & b(:,3)>=tzero & wlngth==0); % & b(:,2)<0.5);
[mmse,imse]=min(mse(i_c));

if fig==1 % if you want a figure as output
    figure('Position', [10 150 600 600])
    subplot(3,1,1)
    plot(YEARS,raw,'k','LineWidth',2)
    ylabel('\bf Ring width (mm)')
    % title([{'\bf CST Intervention Detection on Core '} col_header(core+1)])
    fig1atext = {['Optimal power = ', num2str(optimal_pwr,4)]};
    text(range(YEARS)/3+YEARS(1), max(raw)/1.2,fig1atext)
end

if i_c(imse)>0
    disp(['Lowest error from fit = ' num2str(mmse)])
    disp(['Best age trend fit from years ' num2str(YEARS(1)) ' to ' ...
        num2str(YEARS(i_c(imse)))])
    disp(['Best fit extends for ' num2str(i_c(imse)) ' years'])
    best=b(i_c(imse),:);
    % best=[.82207569 .01318938 .42187578]
    trendtype=1;
    y_exp=nonlinear_exp(best,xyrs);
    % y_exp=b(1)*exp(-b(2)*xyrs)+b(3); % Equivalent code to preceeding line.
    detrended=transformed-y_exp';
    disp('Initial Age Detrending')
    disp(['Y = ', num2str(best(1),4), '*exp(-', num2str(best(2),4),...
                    '*x)+', num2str(best(3),4)]);
    if fig==1 % if you want a figure as output
        subplot(3,1,2) 
        [h312a, h312h1, h312h2] = plotyy(YEARS,[transformed y_exp'],YEARS(i_c),mse(i_c));
        set(h312h1(1),'LineWidth',2)
        set(h312h1(1),'Color',[0 0 0])
        set(h312h1(2),'Color',[.2 .2 1])
        set(h312h2(1),'LineStyle','none','Marker','.','MarkerFaceColor',[1 .2 .2])
        fig1btext = {['Y = ', num2str(best(1),4), '*exp(-', num2str(best(2),4),...
                    '*x)+', num2str(best(3),4)]};
        text(range(YEARS)/3+YEARS(1), max(transformed)/1.2,fig1btext)
        line([YEARS(i_c(imse)) YEARS(i_c(imse))],... 
            [y_exp(i_c(imse))+.2 y_exp(i_c(imse))+.2],'Color','k','Marker',...
            'v','MarkerEdgeColor',[1 .2 .2],'MarkerFaceColor',[1 .2 .2])
        set(get(h312a(1),'Ylabel'),'String','\bf Transformed width') 
        set(get(h312a(2),'Ylabel'),'String','\bf Error Term Variance')
        % ylabel('\bf Transformed width')
        subplot(3,1,3) 
        plot(YEARS,detrended,'k','LineWidth',2)%,YEARS,zeros(1,length(YEARS)),'r')
        ylabel('\bf Transformed width')
        xlabel('\bf Year')
    end
else
    trendtype=2;
    xyrs=(1:nyrs)';
    % Linear detrending option used if neg. exponential curve dissallwd.
    [b,~,~,~,stats] = regress(transformed,...
        [ones(length(YEARS),1) xyrs]); % [ones(length(YEARS),1) YEARS]);
    if b(2)>=0; trendtype=3; end % Find positive age trends
    % if b(2)<0;
        y_lin=b(2)*xyrs +b(1);
        detrended=transformed-y_lin;
        disp('Initial Age Detrending')
        disp(['Y = ', num2str(b(2)), ' * X + ', num2str(b(1))]);
        if fig==1 % if you want a figure as output
            subplot(3,1,2)
            h312b=plot(YEARS,transformed,'k',YEARS,y_lin,'k--');
            set(h312b(1),'LineWidth',2)
            fig1btext = {['Y = ', num2str(b(2)), ' * X + ', num2str(b(1))]};
            text(range(YEARS)/3+YEARS(1), max(transformed)/1.2,fig1btext)
            ylabel('\bf Transformed width')
            subplot(3,1,3)
            plot(YEARS,detrended,'k','LineWidth',2)%,YEARS,zeros(1,length(YEARS)),'r')
            ylabel('\bf Transformed width')
            xlabel('\bf Year')
        end
    % else
    %    y_lin=zeros(nyrs,1);
    %    detrended=transformed-y_lin;
    % end
end
% Output age detrending info
age={char(col_header(core)); trendtype; YEARS(i_c(imse))};

% % Plot a histogram of the data to investigate its skew.
% figure('Position', [50 25 400 300]);
% hist(detrended)
% title('\bf Histogram of Detrended Ring Widths')

% Initialize arrays.
next_iter=1; %Switch to determine whether next iteration is needed
St=detrended; % St will be the iterated series (standardized)
Atr=NaN(length(raw),1); % Age trend re-expressed in raw units
rline=NaN(length(raw),1); % Just the slope of the intervention
tline=NaN(length(raw),iter); % Slope and constant of the intervention
outs=zeros(iter,5);

for q=1:iter % Iterate AR model 'iter' times to remove all outliers
    if next_iter==1
        bckcasted=0;ar_estimate=0;residuals=0;area_t=0;
        iter_i=St; % Initial values of series for ith iteration.
        disp(' ')
        disp(['Statistics for AR model iteration ' int2str(q) ':'])
        
        % Calculate best AR model order and return in the following order:
        % residuals (white noise) and ar model estimates
        [ar_white, ar_model]=ar_order(St);
        
        % Use new coefficients to prewhiten ORIGINAL series without
        % downweighted originals.
        
        % Backcast for pth order years of AR model.
        % bckcasted=backcast(detrended);
        bckcasted=backcast(St); %I think this is correct here.
        
        for g=ORDER:length(bckcasted) % g = observation year
            ar=0; % ar model estimate for order i, year g
            for k=1:ORDER % kth parameter of order ORDER
                if (g-ORDER)>0 % ensure obs yr > model order
                    ar=PARAM(k)*(bckcasted(g-k))+ar;
                end
            end
            if g-ORDER>0 % calculate model estimate and residuals
                if detrended(g-ORDER)==0 % Set missing rings to ar estimate value
                    disp(['Missing ring at year: ' int2str(YEARS(g-ORDER))])
                    bckcasted(g)=ar;
                end
                ar_estimate(g-ORDER)=ar;
                residuals(g-ORDER)=(bckcasted(g)-ar);
            end
        end
        ar_estimate=ar_estimate';
        residuals=residuals';
        if fig==1 % if you want a figure as output
            figure('Position', [600 150 600 600])
            subplot(4,1,1)
            h411=plot(YEARS,St,'k',YEARS,ar_estimate,'k-.');
            set(h411(1),'LineWidth',2)
            title(['\bf' char(col_header(core)) ' iteration ' int2str(q)])
            ylabel('\bf Trans. width')
            xlabel('\bf Year')
            axis([min(YEARS) max(YEARS) min(St)*.9 max(St)*1.1])
        end
   
        % Find release outliers
        [downres, mres, otype]=outlier_clt(residuals,scale,kind,fig);
        f=find(downres~=0);
        
        if q>1 && length(f)>1 && outs(q-1,1)==YEARS(min(f));
            disp('Disturbance detected in same year as previous iteration.')
            disp(['Last iteration: ' num2str(outs(q-1,1))...
                ' and this iteration: ' num2str(YEARS(min(f)))]);
            disp('Terminating iteration since series is not detrending well.')
            f=[];
        end
        if otype==0 && ~isempty(f) % Pulse Outlier Detected
            St(f)=ar_estimate(f);
        elseif otype>0 && length(f)>1 % Trend Outlier Detected
            w=[ones(length(f),1) (1:length(f))'];
            slope=regress(St(f),w);
            % disp(['Constant and slope = ' num2str([slope(1) slope(2)])])
            
            % Fit Hugershoff curve to remainder of series
            lngthw=min(f):length(St);
            lngthwf=(max(f)+1):length(St);
            lngthn=1:length(lngthw);
            lngthn=lngthn(:);
            opts = statset('nlinfit');
            % opts.RobustWgtFun = 'bisquare';
            opts.FunValCheck = 'off';
            opts.MaxIter = 400;
            bw=nlinfit(lngthn,St(lngthw),@hugershoff,[.1 .5 .1 .1],opts);
            disp(['Hugershoff Parameters: ' num2str(bw)])
            ar_est=ar_estimate(f(1)); 
            rline(lngthw)=-bw(1)*(lngthn.^bw(2)).*exp(-bw(3)*lngthn)-bw(4);
            
            % If nlinfit returns NaN, then try again with diffent initial parameters.
            if find(isnan(rline(lngthw)))>0; rline(lngthw)=0;
                disp('Default initial parameters for Hugershoff curve failed')
                disp('Fitting alternate, robust initial parameters [.1 .5 .1 .1]')
                opts.RobustWgtFun = 'bisquare'
                bw=nlinfit(lngthn,St(lngthw),'hugershoff',[.1 .5 .1 .1],opts);
                disp(['Hugershoff Parameters: ' num2str(bw)])
                ar_est=ar_estimate(f(1));                
                rline(lngthw)=-bw(1)*(lngthn.^bw(2)).*exp(-bw(3)*lngthn)-bw(4);
            end
            
            % If nlinfit returns NaN, then end outlier iterations and quit.
            if find(isnan(rline(lngthw)))>0; rline(lngthw)=0;
               disp('Unable to fit Hugershoff curve')
               ar_est=0;
               next_iter=0; 
               outs(q,1:5)=[0 0 0 0 0];
            end            
            
            if f(1)>1
                St(lngthw)=rline(lngthw)+St(lngthw)+ar_est;
                % St(v)=rline(v)+St(v)-slope(1)+St(min(f)-1);
                tline(lngthw,q)=-rline(lngthw);
            elseif f(1)==1 % If trend occurs in 1st yr of series
                St(lngthw)=rline(lngthw)+St(lngthw);
                tline(lngthw,q)=-rline(lngthw);
            end            
            
            outs(q,1:5)=[YEARS(min(f)) YEARS(max(f)) slope(1) slope(2) otype];
        end
        
        if isempty(f) % Determine whether any outliers...
            next_iter=0; % were detected on this iteration
        end
        
        if q==iter && ~isempty(f)
            disp('Need to run additional iterations to resolve series!')
        end
        
        if fig==1 % if you want a figure as output
            subplot(4,1,4)
            hold on
            if ~isempty(f) && min(f)>1 % Draw detrended regression line
                line([YEARS(min(f)) YEARS(max(f))],[ar_est ar_est],...
                    'Color',[.6 .6 .6],'LineStyle','--','LineWidth',2)
                % line([YEARS(min(f)) YEARS(max(f))],[St(min(f)-1) St(min(f)-1)],...
                %     'Color',[.6 .6 .6],'LineStyle','--','LineWidth',2)
            else % Draw same line, but set to first year of series
                line([YEARS(1) YEARS(max(f))],[0 0],...
                    'Color',[.6 .6 .6],'LineStyle','--','LineWidth',2)
            end
            h414=plot(YEARS,iter_i,'k',YEARS,St,'k',YEARS,tline(:,q),'k');
            set(h414(1),'LineWidth',2)
            set(h414(3),'LineWidth',2,'Color',[.6 .6 .6])
            ylabel('\bf Trans. width')
            xlabel('\bf Year')
            ymin=min([min(iter_i) min(St)])*.9;
            ymax=max([max(iter_i) max(St)])*1.1;
            axis([min(YEARS) max(YEARS) ymin ymax])
            box on
            hold off
            
            subplot(4,1,2)
            h412=plot(YEARS,residuals,'k-.',YEARS,zeros(1,length(YEARS)),'k',...
                YEARS,mres);
            set(h412(3),'Color',[.6 .6 .6])
            set(h412(3),'LineWidth',2)
            % title(['\bf AR Residuals and running residual mean for '...
            %     int2str(length(f)) ' years'])
            axis([min(YEARS) max(YEARS) min(residuals)*1.1 max(residuals)*1.1])
            ylabel('\bf Residuals')
            xlabel('\bf Year')
        end
    end
end


if fig==1 % if you want a figure as output
    % Shows final iterated series in transformed units
    figure('Position', [1200 150 600 400])
    subplot(2,1,1)
    transDt=detrended-St; % transformed outlier series
    % disp(['transDtarea = ' num2str(sum(transDt(f)))])
    % hpentult=plot(YEARS,detrended,'k',YEARS,transDt,'LineWidth',2);
    hpentult=plot(YEARS,detrended,'k',YEARS,St,'k--');
    set(hpentult(1),'LineWidth',2)
    set(hpentult(2),'LineWidth',2)
    % set(hpentult(2),'Color',[.6 .6 .6])
    % title('\bf Outlier Series')
    legend('Age-detrended series','Standardized series',...
        'Location','NorthWest')
    legend('boxoff')
    ylabel('\bf Transformed width')
%     subplot(2,1,2)
%     plot(YEARS,St,'k','LineWidth',2)
%     title('\bf Standardized Series')
%     ylabel('\bf transformed mm')
end

% Shows final iterated series in original units (mm presumably)

% St=St+St_pos;

if trendtype==1 % negative exponential trend
    Stt=y_exp'+St; % Size trend & first detrending
    % Stt=Stt+transm; % Add transformed mean back to series

    
    if optimal_pwr <= 0.05
        Str=10.^(Stt);% Size trend in original (raw) units
        Atr=10.^(y_exp');% Age trend in original (raw) units
    else
        Stt(Stt<=0)=0; % Set neg values to zero
        Str=(Stt).^(1/optimal_pwr);
        Atr=(y_exp').^(1/optimal_pwr);
    end
elseif trendtype==2 || trendtype==3 % linear regression trend
    Stt=y_lin+St; % Size trend & first detrending
     % Stt=Stt+transm; % Add transformed mean back to series
    
    if optimal_pwr <= 0.05
        Str=10.^(Stt);% Size trend in original (raw) units
        Atr=10.^(y_lin);% Age trend in original (raw) units
    else
        Stt(Stt<=0)=0; % Set neg values to zero
        Str=(Stt).^(1/optimal_pwr);
        Atr=(y_lin).^(1/optimal_pwr);
    end
else
    disp('Error in trend type designation')
end
% outs(:,6)=(nanmean(prline,1))'; % mean prior to intervention removal
% outs(:,7)=(nanmean(arline,1))'; % mean after intervention removal

raw(mss>0)=NaN; % Remove estimated values of missing rings
Str(mss>0)=NaN; % Remove estimated values of missing rings
Dtr=raw-Str; % Remove estimated values of missing rings

if fig==1 % if you want a figure as output
    subplot(2,1,2)
    h_end=plot(YEARS,Dtr,YEARS,Str,'k--',YEARS,raw,'k');%,YEARS,Atr,'b');
    set(h_end(1),'Color',[.6 .6 .6])
    set(h_end(1),'LineWidth',2)
    set(h_end(2),'LineWidth',2)
    set(h_end(3),'LineWidth',2)
    legend('Disturbance index','Standardized series','Original series',...
        'Location','NorthWest')
    legend('boxoff')
    ylabel('\bf Ring width (mm)')
    xlabel('\bf Year')
end

% ar_order.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This subfunction is based on series_ar.m and determines the
% autoregressive parameters for the best model order as calculated using
% AIC criteria. The function returns the residuals, and AR model
% estimate of the best order found with AIC criteria.
function [out_res, out_est]=ar_order(series)
global PARAM
global ORDER
global YEARS

% Initializes variables for autoregressive modeling.
ar_param=0; residuals=0;

% Calculate Autoregressive parameters for orders 1 through 10.
for ar_order=1:10
    ar_param(ar_order,1:ar_order+1)=-arburg(series,ar_order);
end

%Remove first column of minus ones from ar_param.
ar_param(:,1)=[];

% Calculate residuals for particular AR order model.
for i=1:10 % i = ar model order
    for g=1:length(YEARS) % g = observation year
        ar=0; % ar model estimate for order i, year g
        for k=1:length(ar_param(1,:)) % kth parameter of order i
            if (g-k)>0 % ensure obs yr > model order
                ar=ar_param(i,k)*(series(g-k))+ar;
            end
        end
        if g-i>0 % calculate residuals
            residuals(g,i)=(series(g)-ar);
            ar_estimate(g,i)=ar;
        end
    end
    % Calculate the total variance of the residuals by model order
    resid_var(i)=var(residuals(:,i));
end

% Calculate variance of the residuals of a particular AR order model.
% Reference Box & Jenkins & Reinsel 1994 pp. 200-201.
% Using Akaike Information Criteria (A-key-key)
% Equation now uses natural log and simply 'n' in the denominator.
% t+1 (or # of params+1) is a penalty factor for estimating the mean.
aic=0;
for t=1:length(resid_var)
    aic(t)=log(resid_var(t))+(2*(t+1))/length(YEARS);
end
% Find the first minimum AIC order (ie first saw-tooth-shaped dip).
% If AIC values monotonically decrease, set best_order=9.
best_order=0;
for s=2:length(aic)
    if((aic(s)>=aic(s-1)) && (best_order==0))
        best_order=s-1;
    end
end

if best_order==0;
    best_order=9;
end

% disp('AIC Values:')
% disp(num2str(aic))

ORDER=best_order;
PARAM=ar_param(best_order,1:best_order);

% disp(['Best order = ' int2str(ORDER)])
disp(['AR Model Parameters: ' num2str(PARAM)])

out_res=residuals(:,best_order);
out_est=ar_estimate(:,best_order);

% outlier_clt.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                               %
% This subfunction determines the auto regressive outliers in the
% residuals that are greater than a given number of std devs using the
% central limit theorem.
% 99% of the observations lie within 2.58(std_res)
% 97.5% of the observations lie within 2.24(std_res)
% Ed recommends 3(std_res)

function [dres, rmr, type]=outlier_clt(in,zsc,kind2,fig2)
global YEARS

% initialize variables
type=0; % Type of outlier detected (1=positive, 2=negative)
lngth=length(YEARS);
a=9; b=30;
% a=9; b=30; %b=lngth/3; %b=lngth-40; % min and max of trend window
if b>lngth/4
    b=floor(lngth/4);
    disp(['Maximum outlier detection length reduced to ' num2str(b) ' due to low ring #'])
end
lt=a; % Length of trend
window=0;
% zsc=3.29; % zsc is the scale value for detection (ie zscore) 
% zsc=1.96; % 95 pct CI
% zsc=2.58; % 99 pct CI8
% zsc=3.29; % 99.9 pct CI
dres=zeros(lngth,1); % downweighted residuals
mr=zeros(lngth,1); % residuals mean in window
rmr=zeros(lngth,1);
rmu=0;
rshat=0;

% initialize masked to ones.
marker=zeros(length(YEARS),1);
masked=zeros(length(YEARS),1);

std_res = (var(in))^0.5; % Calculate std dev of residuals

% for u=1:length(YEARS) % Detect pulse outliers
%     rres=in(u)/(zsc*std_res); % calculates relative residuals
%     if  rres>=1.0
%         dres(u)=1;
%         type=0;
%         disp(['Positive Pulse in ' int2str(YEARS(u))])
%         disp(['Outlier value = ' num2str(rres)])
%     elseif rres<=-1.0
%         dres(u)=1;
%         type=0;
%         disp(['Negative Pulse in ' int2str(YEARS(u))])
%         disp(['Outlier value = ' num2str(rres)])
%         % elseif abs(rres)<1.0
%         % psi=rres;
%         % The code below simply produces psi = rres.  Why did Ed code it
%         % this way in his robar function?
%         % psi=rres*exp(-exp(3.0*(abs(rres)-3.0)));
%     else
%         disp('No pulse outliers detected')
%     end
% end

if a<=b
    for v=a:b
        z=v-a+1;
        for u=1:(lngth-v)%+1) % Changed 4-13-13
            window=in(u:(u+v-1));
            mr(u,z)=mean(window);
        end
        % [muhat(z),sigmahat(z)] = normfit(mr(:,z)) % Arithmetic mean
        % Uses Tukey's bi-weight mean instead (Hoaglin 1983, Meko's code)
        [ymn(z),varyh(z),~,~,~,se]=bisqmean_dan(mr(:,z));
        [mam(z),imax(z)]=max(mr(:,z)); % Find max. means & their locations
        [mim(z),imin(z)]=min(mr(:,z)); % Find min. means & their locations
    end
    
    poso=(mam-ymn)./varyh; % Determines # deviations from mean of means
    nego=(ymn-mim)./varyh;
    [relmam, rimax]=max(poso); % Find max. of positive dev.s & their locations
    [relmim rimin]=max(nego); % Find max of negative dev.s & their locations
    disp(['Max departure = ' num2str(relmam)])
    disp(['Min departure = ' num2str(relmim)])
    if kind2==1 % Detect positive outliers only
        if (poso(rimax)>=zsc)
            type=1;
            lt=rimax+a-1; % length of trend
            dres(imax(rimax):(imax(rimax)+lt-1))=ymn(rimax);
            disp(['Release detected in ' int2str(YEARS(imax(rimax)))...
                ' for ' int2str(lt) ' years'])
            rmu=ymn(rimax);
            rshat=varyh(rimax);
            rmr=mr(:,rimax);
        else
            rmu=ymn(1);
            rshat=varyh(1);
            rmr=mr(:,1);
            disp(['No trend outliers detected up to ' int2str(b) ' yrs'])
        end
    elseif kind2==2 % Detect negative outliers only
        if (nego(rimin)>=zsc) % Comment elseif for only + outliers
            type=2;
            lt=rimin+a-1;
            dres(imin(rimin):(imin(rimin)+lt-1))=ymn(rimin);
            disp(['Suppression detected in ' int2str(YEARS(imin(rimin)))...
                ' for ' int2str(lt) ' years'])
            rmu=ymn(rimin);
            rshat=varyh(rimin);
            rmr=mr(:,rimin);
        else
            rmu=ymn(1);
            rshat=varyh(1);
            rmr=mr(:,1);
            disp(['No trend outliers detected up to ' int2str(b) ' yrs'])
        end
    elseif kind2==3 % Detect both positive and negative outliers
        if (poso(rimax)>=zsc) && (relmam>=relmim)
            type=1;
            lt=rimax+a-1; % length of trend
            dres(imax(rimax):(imax(rimax)+lt-1))=ymn(rimax);
            disp(['Release detected in ' int2str(YEARS(imax(rimax)))...
                ' for ' int2str(lt) ' years'])
            rmu=ymn(rimax);
            rshat=varyh(rimax);
            rmr=mr(:,rimax);
        elseif (nego(rimin)>=zsc)
            type=2;
            lt=rimin+a-1;
            dres(imin(rimin):(imin(rimin)+lt-1))=ymn(rimin);
            disp(['Suppression detected in ' int2str(YEARS(imin(rimin)))...
                ' for ' int2str(lt) ' years'])
            rmu=ymn(rimin);
            rshat=varyh(rimin);
            rmr=mr(:,rimin);
        else
            rmu=ymn(1);
            rshat=varyh(1);
            rmr=mr(:,1);
            disp(['No trend outliers detected up to ' int2str(b) ' yrs'])
        end
    end
end

if fig2==1 % if you want a figure as output
    subplot(4,1,3)
    hold on
    hist(rmr)
    h413 = findobj(gca,'Type','patch');
    set(h413,'FaceColor','k')
    box on
    % title('\bf Histogram of Running AR Residual Means')
    line([rmu-zsc*rshat rmu+zsc*rshat],[10 10],'Color',[.6 .6 .6])
    plot(rmu,10,'o','Color',[.6 .6 .6])
    ylabel('\bf Frequency')
    xlabel(['\bf' int2str(lt) '\bf Yr Residual Means'])
    hold off
end

% backcast.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                              
% This subfunction estimates the first elements of a series for which the   
% residuals could not be calculated owing to the use of ar modeling.       
function bckcasted=backcast(seriesb)
global PARAM
global ORDER
global YEARS

% Invert time series for backcasting.
flipped=flipud(seriesb); 

% Add in backcasted AR estimates as new values at end of inverted series.
for g=(length(YEARS)+1):(length(YEARS)+ORDER) % g = backcasted years
    ar=0; % ar model estimate for order i, year g
    for k=1:ORDER % kth parameter of order ORDER
        ar=PARAM(k)*(flipped(g-k))+ar;
    end
    flipped(g)=ar;
end

% Re-invert series and return as output.
bckcasted=flipud(flipped);
% disp('Backcasted Values:')
% for h=ORDER:-1:1
%     disp(['Year -' int2str(h) ': ' num2str(bckcasted(h))])
% end