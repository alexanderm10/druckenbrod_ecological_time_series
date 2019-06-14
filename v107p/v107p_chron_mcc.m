% v107p_chron_mcc.m
% This function calculates each autoregressive outlier for each
% core in a dataset and lumps those results by tree.  This version also
% returns the transformed standardized series for each core in St.

% Function written Mar 10, 2004.
% Function last revised Sep 28, 2015.

function v107p_chron_mcc

disp('Disturbance Growth Chronology ver107pn')
disp('Compiled Sept 30 2015')
disp('This version detrends positive and negative events using a Hugershoff curve')
disp('Contact ddruckenbrod@rider.edu for distribution or support')
disp('Select ringwidth file with values ending in either -9999 or 999.')

% Import tree-ring data (returns vars *col_header* and *rings*)
hdr=0;
hdr=input('How many header rows preceed ringwidth data? Enter 0 for none: ');
[col_header,rings,filename]=ringwidth_import_v107(hdr);

iter=8; % Set maximum number of iterations to 20.
years=rings(:,1); % years for entire chronology
rings(:,1)=[]; % remove year column
col_header(1)=[]; % remove year label from array
ncores=size(rings,2); % # cores in group
nyrs=size(rings,1); % total # of years in chronology
expval=NaN(1,ncores); % value of last year that neg exp curve fits
yrs=NaN(nyrs,ncores); % years for each cores
tres=NaN(nyrs,ncores); % Transformed residuals for each core
det=NaN(nyrs,ncores); % Detrended series for each core
St=NaN(nyrs,ncores); % Undisturbed series for each core in transformed units
Straw=NaN(nyrs,ncores); % Undisturbed series for each core
Atraw=NaN(nyrs,ncores); % Age series for each core
Dtraw=NaN(nyrs,ncores); % Disturbance series for each core
% agestats=NaN(2,ncores); % power and trend type for transformed core
agestats=cell(3,ncores); % power and trend type for transformed core
% agestats(3,:)={'0000'};
out=NaN(iter,7,ncores); % Outlier statistics

diary(strcat(filename,'.txt')) % Record screen output
disp(['Processing file ' filename])
disp(['Chronology contains ' num2str(ncores) ' cores.']);
disp(['Chronology extends from ' num2str(min(years)) ' to ' num2str(max(years))]);


ts=3.29; tk=1; tyr=2; figs=2;% Default Values for Menu
menu_options=1;
while menu_options==1
    disp(' ')
    disp('Set Disturbance Growth Chronology Parameters.')
    disp('Current parameters are listed after arrow. Select using numbers')
    disp('0: Keep current parameters and start run')
    disp(['1: Scale for Detection (ie z-score) -> ' num2str(ts)]);
    disp(['2: Disturbance events: pos=1, neg=2, both=3 -> ' num2str(tk)]);
    disp(['3: Mark expected disturbance year: yes=1, no=2 -> ' num2str(tyr)]);
    disp(['4: Show figures of individual series: yes=1, no=2 -> ' num2str(figs)]);
    t0=0; t0=input('Select option: ');
    disp(' ')
    
    if t0==1
        ts=input('Select new scale value: ');
        disp(['Scale set to: ' num2str(ts)]);
    elseif t0==2
        disp('1: Positive events only')
        disp('2: Negative events only')
        disp('3: Both positive and negative events')
        tk=input('Select disturbance category: ');
        if tk==1|2|3
        else
            tk=1; disp('Category chosen not permissible. Reverting to default.')
        end
        disp(['Disturbance category selected: ' num2str(tk)]);
    elseif t0==3
        disp('Enter expected disturbance year in XXXX format with year > 0).')
        tyr=input('Enter 2 if you do not want to mark a expected disturbance year: ');
        disp(['Disturbance marker year set to: ' num2str(tyr)]);
    elseif t0==4
        disp('Show figures of individual series?')
        disp('1: Yes')
        disp('2: No')
        figs=input('Enter selection: ');
        if figs==1|2
        else
            figs=1; disp('Category chosen not permissible. Reverting to default.')
            enddisp(['Disturbance marker year set to: ' num2str(figs)]);
        end
    elseif t0==0
        menu_options=0;
    end
end

for i=1:ncores
    disp(' ')
    disp(['Series #' num2str(i) '---------------------------------------'])
    
    % Find pointer to start and end of series
    s=find(rings(:,i)>0,1);
    e=find(rings(:,i)>0,1,'last');
    [yrs(s:e,i),tres(s:e,i),det(s:e,i),St(s:e,i),Straw(s:e,i),Dtraw(s:e,i),...
        Atraw(s:e,i),agestats(:,i),out(1:iter,1:5,i)]=...
        v107p_mcc(i,figs,iter,filename,ts,tk,years,col_header,rings);
    
    b=find(out(:,5,i)==2);% find all releases for a core
    if b>0
        for c=1:length(b) % iterate through each release
            startyr=find(years==out(b(c),1,i));
            dbh_rel(startyr,i)=sum(rings(s:startyr,i))/1000;
            age_rel(startyr,i)=length(rings(s:startyr,i));
        end
    end
    
    if figs==1
        disp('Continue to display figures for individual series? ')
        disp('1: Yes')
        disp('2: No')
        figs=input('Enter selection: ');
    end
    close all % Close all figures
    clear b
    clear startyr
end

% % Figure shows DBH and age at time of disturbance relative to other samples
% figure('Position', [10 5 700 800]) 
% subplot(2,1,1)
% nanrings=cumsum(rings,1); % cumulative dbh
% coreage=rings;
% coreage(find(coreage>0))=1;
% coreage=cumsum(coreage,1); % count age of each core
% coreage(coreage==0)=NaN;
% av_ca=nanmean(coreage,2);
% v_ca=nanstd(coreage,1,2);
% nanrings(nanrings==0)=NaN;
% av_dbh=nanmean(nanrings,2)/1000;
% v_dbh=nanstd(nanrings,1,2)/1000;
% hold on
% fill([years; flipud(years)],[v_dbh+av_dbh; flipud(av_dbh-v_dbh)],...
% [.7 .7 .7],'EdgeColor','none')
% plot(years,av_dbh,'k',years,dbh_rel,'k--o')
% ylabel('\bf Av. Inside Diameter (m)')
% xlabel('\bf Year')
% hold off
% 
% subplot(2,1,2)
% hold on
% fill([years; flipud(years)],[v_ca+av_ca; flipud(av_ca-v_ca)],...
% [.7 .7 .7],'EdgeColor','none')
% plot(years,av_ca,'k',years,age_rel,'k--o')
% ylabel('\bf Av. Age')
% xlabel('\bf Year')

rel=[0 0 0];
sup=[0 0 0];
d=1; % release counter
f=1; % suppression counter

% Summarize outlier descriptive statistics
for a=1:size(Dtraw,2)% # of cores
    b=find(out(:,5,a)==1); % find all releases for a core
    if b>0
        for c=1:length(b) % iterate through each release
            startyr=out(b(c),1,a);
            endyr=out(b(c),2,a);
            rel(d,1:3)=[a startyr endyr]; 
            d=d+1;
        end
    end
    clear b
    
    g=find(out(:,5,a)==2); % find all suppressions
    if g>0
        for h=1:length(g)
            startyr2=out(g(h),1,a);
            endyr2=out(g(h),2,a);
            sup(f,1:3)=[a startyr2 endyr2];
            f=f+1;
        end
    end
    clear g
end
rings(find(~rings))=NaN; % Convert rings matrix zeros to NaNs

% mraw=nanmean(rings,2);
% mStraw=nanmean(Straw,2);

% Find and average all cores with pos or neg outliers
subset=unique([rel(:,1); sup(:,1)]);
subset=subset(find(subset));% Find & remove nonzeros if no pos or neg outliers found
if subset % Only graph if interventions found.
    sigDtraw=Dtraw(:,subset);
    sigDtm=nanmean(Dtraw(:,subset),2);
    
    depth=size(Dtraw,2)-sum(isnan(Dtraw),2); % Total sample depth
    % samples with outliers
    subdepth=size(Dtraw(:,subset),2)-sum(isnan(Dtraw(:,subset)),2);
    
    figure('Position', [10 5 700 800])
    subplot(4,1,1)
    h_end=plot(years,nanmean(rings-Atraw,2),'k',years,nanmean(Straw-Atraw,2),'k--');
    set(h_end(1),'LineWidth',2)
    set(h_end(2),'LineWidth',2)
    legend('Mean Ct + Dt','Mean Ct','Location','NorthWest')
    legend('boxoff')
    ylabel('\bf Residuals (mm)')
    xlabel('\bf Year')
    
    subplot(4,1,2)
    [AX,H1,H2] = plotyy(years,sigDtm,years,depth);
    set(H1,'LineWidth',2)
    set(H1,'Color','k')
    set(H2,'Color','k')
    set(AX(1),'ycolor','k')
    set(AX(2),'ycolor','k')
    set(get(AX(1),'Ylabel'),'String',{'\bf Mean Dt', '(mm)'})
    set(get(AX(2),'Ylabel'),'String','\bf Sample Size')
    set(AX(2),'ylim',[0 ceil((ncores+1)/10)*10])
    set(AX(1),'XTickLabel',[])
    bounds=xlim;
    box off
%     h_end=plot(years,sigDtm,years,mStraw,'k--',years,mraw,'k');
%     set(h_end(1),'Color',[.6 .6 .6])
%     set(h_end(1),'LineWidth',2)
%     set(h_end(2),'LineWidth',2)
%     set(h_end(3),'LineWidth',2)
%     legend('Disturbance index','Standardized series','Original series',...
%         'Location','NorthWest')
%     legend('boxoff')
%     ylabel('\bf Ring width (mm)')
%     xlabel('\bf Year')

    subplot(4,1,3)
    hold on;
    edges=[bounds(1):bounds(2)]; % Bin outliers annually
    if tk==1
        pcounts=0;pcounts=hist(rel(:,2),bounds(1):bounds(2));
        ncounts=0;
        bar(edges,pcounts,'k')
    elseif tk==2
        pcounts=0;
        ncounts=0;ncounts=-hist(sup(:,2),bounds(1):bounds(2));
        bar(edges,ncounts,'k')
    elseif tk==3
        pcounts=0;pcounts=hist(rel(:,2),bounds(1):bounds(2));
        ncounts=0;ncounts=-hist(sup(:,2),bounds(1):bounds(2));
        bar(edges,pcounts,'k')
        bar(edges,ncounts,'k')
    end
    
    if tyr~=2;
    line([tyr tyr],[min(ncounts) max(pcounts)],'LineWidth',1,'Color',[.8 .2 .2]);
    end
    xlim([bounds(1) bounds(2)]);
    ylabel('\bf Dt Initiation Yrs')
    hold off
    
    subplot(4,1,4) % show cores that are open grown initially
    agenum=cell2mat(agestats(2,:));
    %       disp('Trees that likely established in open conditions')
    %    disp(col_header(find(agenum==1)));% | agenum==2))));
    %     if ~isempty(opngrwn>1)
    %         hold on
    %         plot(years,Atraw,'r',years,opngrwn,'b')
    %         if size(agestats,1)==3
    %             for m=1:ncores
    %                 if cell2mat(agestats(3,m))
    %                     yrval(m)=cell2mat(agestats(3,m));
    %                     expval(m)=Atraw(find(years==yrval(m)),m);
    %                 end
    %             end
    %             yrval=yrval(yrval>0);
    %             expval=expval(expval>0);
    %             scatter(yrval,expval,'ob')
    %         end
    %     ylabel('\bf Gt (mm)')
    %     xlabel('\bf Year')
    %     hold off
    
    % edges=mindecade:10:maxdecade; % Bin outliers by year
    for i=1:size(yrs,2); firstyr(i)=yrs(find(rings(:,i)>0,1,'first'),i);end
    opngrwn=firstyr(agenum==1);
    % disp('Trees that likely established in open conditions')
    % disp(col_header(agenum==1));
    clsdcan=firstyr(agenum==2|agenum==3);
    % Establishment dates binned by year
    o=histc(opngrwn,edges); o=o(:);% open grown establishment
    u=histc(clsdcan,edges); u=u(:);% understory (shaded) establishment
    Y=[u o];
    bar(edges,Y, 'stacked')
    xlim([bounds(1) bounds(2)])
    ylim([0 max(o+u)])
    % axis([mindecade maxdecade 0 max(o+u+5)])
    ylabel('\bf Tree Recruitment')
    xlabel('\bf Year')
    ColorOrder2=...
        [0 0 0; .6 .6 .6];
    colormap(ColorOrder2)

%     figure(3); hold on;
%     counts=hist(rel(:,2),1910:2010);hist(rel(:,2),1910:2010)
%     axis([min(years) max(years) 0 max(counts)])
%     line([1990 1990],[0 max(counts)],'LineWidth',2,'Color',[.8 .2 .2]); hold off;
else
    sigDtraw=0;
    sigDtm=0;
end

p_out=NaN(iter,ncores);
n_out=NaN(iter,ncores);

rr=0; % Organize release years for output file
for i=1:ncores
    rr=find(rel(:,1)==i);
    if rr
    p_out(1:length(rr),i)=rel(rr,2);
    end
    rr=0;
end

rs=0; % Organize suppression years for output file
for i=1:ncores
    rs=find(sup(:,1)==i);
    if rs
    n_out(1:length(rs),i)=rel(rs,2);
    end
end

warning('off', 'MATLAB:xlswrite:AddSheet');
dist_total=size(find(rel(:,1)),1)+size(find(sup(:,1)),1);
disp(['A total of ' num2str(dist_total) ' disturbance events were detected in '...
    num2str(size(subset,1)) ' series out of ' num2str(ncores)...
    ' total series in chronology.']);
disp('Saving output to Excel File')
disp('Sheet 1 = Series Names and Growth Curve Categories')
disp('Sheet 2 = At series')
disp('Sheet 3 = Dt series')
disp('Sheet 4 = Disturbance standardized series')
disp('Sheet 5 = Initial years of release by series')
disp('Sheet 6 = Initial years of suppression by series')
disp('Sheet 7 = Average Chronology Values')

% Output data to an excel file
% State shape of age/size trend for each series
ttype=cellstr(['NegExp';'NegLin';'PosLin']);
txls=ttype(cell2mat(agestats(2,:)))';
xlswrite(strcat(filename,'.xlsx'),...
    ['Series Name' col_header; 'Growth Curve' txls], 1)
xlswrite(strcat(filename,'.xlsx'),['At' col_header;... 
    num2cell([years(:) round(Atraw,3)])], 2)
% xlswrite(strcat(filename,'.xlsx'),[years(:) round(Atraw,3)], 2)
xlswrite(strcat(filename,'.xlsx'),['Dt' col_header;...
    num2cell([years(:) round(Dtraw,3)])], 3)
xlswrite(strcat(filename,'.xlsx'),['Disturbance Standardized Series' col_header;...
    num2cell([years(:) round(Straw,3)])], 4)
pxls=[col_header; num2cell(p_out)];
nxls=[col_header; num2cell(n_out)];
xlswrite(strcat(filename,'.xlsx'),pxls, 5)
xlswrite(strcat(filename,'.xlsx'),nxls, 6)
% avtxt=(cellstr(['Years            '; 'Av Raw Ringwidths'; 'Av Age/Size Trend';...
% 'Av Dist Chron    '; 'Av Dist Std Chron'; 'Av Cores w Intvns']))';
avtxt=(cellstr(['Years              '; 'Avg Raw            ';...
    'Avg At             '; 'Avg Dt             '; 'Avg Dist Stndrdized';...
    'Avg Ct             '; 'Avg Series w Events']))';
avs=num2cell(round([years(:) nanmean(rings,2) nanmean(Atraw,2)...
    nanmean(rings-Straw,2) nanmean(Straw,2) nanmean(Straw-Atraw,2)...
    (sigDtm(:))]*1000)/1000);
xlswrite(strcat(filename,'.xlsx'),[avtxt; avs], 7)
diary off % stop recording screen output