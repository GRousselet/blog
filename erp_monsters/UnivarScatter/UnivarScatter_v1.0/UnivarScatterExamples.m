% Generate 3 normal COLUMN VECTORS of different lengths
smallset=normrnd(0,3,[10,1]);
averageset=normrnd(0,3,[50,1]);
bigset=normrnd(0,3,[150,1]);

%Concatenate them, and fill the empty ones with nans, keep in mind that you
%can do this much faster with the function padcat(available in mathworks, Copyright (c) 2009, Jos van der Geest)
DataArray=nan(size(bigset,1),3);
DataArray(1:length(smallset),1)=smallset;
DataArray(1:length(averageset),2)=averageset;
DataArray(1:length(bigset),3)=bigset;

answer=questdlg('Do you want to choose the colors from the palette?','','Yes','Use Default','Use Default');
if strcmp('Yes',answer)
    Colors=ColorCoder(3);
else
    Colors=[0.9 0 0;0 0.9 0; 0 0.9 0.9];
end
%
figure
UnivarScatter(DataArray,'Label',{'Small','Average','Big'},'MarkerFaceColor',Colors);
title('example 1')
%
figure
[~,~,~,RangeCut]=UnivarScatter(DataArray,'Label',{'Small','Average','Big'},'MarkerFaceColor',Colors,'SEMColor',Colors/1.5,'StdColor',Colors/2);
title('example 2')
%
figure
UnivarScatter(DataArray,'Label',{'Small','Average','Big'},'MarkerFaceColor',Colors,'SEMColor',Colors/1.5,'StdColor',Colors/2);
title('example 2 changing pbaspect')
pbaspect([6,4,1])
%
figure
UnivarScatter(DataArray,'Label',{'Small','Average','Big'},'MarkerFaceColor',Colors,'Whiskers','lines');
title('example 2 with lines')
pbaspect([4,4,1])
%
figure
UnivarScatter(DataArray,'Label',{'Small','Average','Big'},'MarkerFaceColor',Colors,'SEMColor',Colors/1.5,'StdColor',Colors/2,'RangeCut',RangeCut+3);
title('example 2 with changed RangeCut')


% In case we are working with tables
a=[1 2 3 2 5 6 7 7 7 8 12 10 14 17 22 10 10 10 20 21 22];
b={'a' 'a' 'a' 'a' 'a' 'a' 'a' 'a' 'a' 'b' 'b' 'b' 'b' 'b' 'b' 'b' 'b' 'b' 'b' 'b' 'b'};
c=table(a',b');
figure
[xPoints,yPoints,Label]=UnivarScatter(c);

% This last one is just to give an idea of possible use of the extracted
% data
figure
hold on
title({'Points coloured in the same colour in groups of consecutive 3 points,'; 'this could be useful in a table in which every 3 values were measurements for the same individual in a group'})
UnivarScatter(c);
for i=1:3:length(xPoints(:))
    scatter(xPoints([i i+1 i+2]),yPoints([i i+1 i+2]),'MarkerEdgeColor','k','MarkerFaceColor','flat')
end
set(gca,'xtick',1:size(yPoints,2))
set(gca,'xticklabel',Label)
