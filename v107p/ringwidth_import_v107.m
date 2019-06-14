% ringwidth_import_v107.m
% This function imports decadal format tree ring data for manipulation 
% as a matrix in Matlab.  The number of header lines must be specified
% as an input by the user.  The end of each series must be flagged 
% with 999.  Measurements are stored as one thousandth of a 
% millimeter.  The filename can either be specified as an input or 
% found using a gui. The LAST LINE of the input text file must also
% be blank! This version can handle core labels up to 8 characters and 
% file types other than txt.

% Function written Mar 4, 2004.
% Function last revised Sept 25, 2015.

function [col_header,rings,filename]=ringwidth_import_v107(header,varagin)

if nargin==1
    [filename,path]=uigetfile('*.*',  'All Files (*.*)');
elseif nargin==2
    filename=varagin;
else
    disp('Too many parameters entered.')
end

% Read in header lines 
headers=textread(filename, '%q',10)';
% disp([headers]) % Display 1st 10 words as screen output.
[label yr y0 y1 y2 y3 y4 y5 y6 y7 y8 y9]=...
    textread(filename,...
    '%8s %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d',...
    'headerlines',header);
% Place decadal format widths into one matrix
widths=[y0 y1 y2 y3 y4 y5 y6 y7 y8 y9];

% Determine whether ring widths use 999 or -9999 as end of series 
precisn=100;markr=999; % Default values for ring wdths measured to .01 precision
if find(widths==-9999)
    markr=-9999; precisn=1000;
end

% Extract unique labels of each core.
importedrows=length(label);
a=1;b=1;
while(a<=importedrows)
    core(b)=label(a);
    corestr(:,b)=strcmp(label(a),label);
    a=max(find(corestr(:,b)==1))+1;
    b=b+1;
end

% Find range of years over all cores and set as col 1 in rings.
% As it is difficult to know how many years are in the last row
% of measurments for a core, assume that the last decade has 10
% measurements.
rings=(min(yr):(max(yr)+10))';

% Transfer widths into vectors by core
for i=1:length(core)
    core_rows=find(corestr(:,i));
    core_yr=yr(core_rows);
    core_widths=widths(core_rows,:);
    % Find # of measurements in a row and assign to vector series.
    k=1;series=0;
    for j=1:length(core_rows)
        % Look for end of series flag
        flag=find(core_widths(j,:)==markr);
        if flag>0
            msmts=flag-1;
        elseif (ceil(core_yr(j)/10)*10)-core_yr(j)==0
            msmts=10;
        else
            msmts=(ceil(core_yr(j)/10)*10)-core_yr(j);
        end
        series(k:(k+msmts-1))=core_widths(j,(1:msmts));
        k=msmts+k;
    end
    % Determine start and end of series
    sos=find(rings(:,1)==min(core_yr));
    length(series)+sos-1;
    eos=length(series)+sos-1;
    % Assign series to output matrix and convert to 1/1000th of a mm
    rings(sos:eos,i+1)=(series./precisn)';
end
% Construct column headers
col_header(2:(length(core)+1))=core;
col_header(1)={'Year'};
        