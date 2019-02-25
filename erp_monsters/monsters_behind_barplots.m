% Code to reproduce the figures in my blog post.
%
% Execute each cell in turn by pressing cmd+enter (mac) or ctrl+enter (pc)
% You can navigate between cells by pressing cmd/ctrl + up or down arrow.
%
% Copyright (C) 2016 Guillaume Rousselet - University of Glasgow

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.

%% dependencies ===========================================================
%% structured jitter for scatterplots
% http://www.mathworks.com/matlabcentral/fileexchange/54243-univarscatter
% The code is also available in the blog folder, for convenience, and
% because I've created an alternative version which only outputs the
% jittered x-axis values: function `UnivarScatter_nofig.m`
%% =========================================================================

%% get data
load monsters.mat

%% Figure 1 - ERP bar graph equivalent
figure;set(gcf,'Color','w');hold on
plot(xf,mean(erp1,2),'Color',[1 0.5 0.2],'LineWidth',2)
plot(xf,mean(erp2,2),'Color',[0 .5 0],'LineWidth',2)
plot([-300 600],[0 0],'k--')
xlabel('Time in ms','FontSize',16)
ylabel('ERPs in \mu-volts','FontSize',16)
set(gca,'XLim',[-300 600],'XTick',-300:100:600,'YLim',[-15 10],'FontSize',14,'Layer','Top')
box on
legend('Condition 1','Condition 2')
plot(156,-6,'r*','MarkerSize',12,'LineWidth',2)

%% Figure 2 - ERP bar graph equivalent + difference

figure;set(gcf,'Color','w');hold on
plot(xf,mean(erp1,2),'Color',[1 0.5 0.2],'LineWidth',2)
plot(xf,mean(erp2,2),'Color',[0 .5 0],'LineWidth',2)
plot(xf,diff,'Color',[.1 .1 .1],'LineWidth',3)
plot([-300 600],[0 0],'k--')
xlabel('Time in ms','FontSize',16)
ylabel('ERPs in \mu-volts','FontSize',16)
set(gca,'XLim',[-300 600],'XTick',-300:100:600,'YLim',[-15 10],'FontSize',14,'Layer','Top')
box on
legend('Condition 1','Condition 2','difference')
plot(156,-6,'r*','MarkerSize',12,'LineWidth',2)

%% Figure 3 - difference + confidence interval

figure;set(gcf,'Color','w');hold on
plot(xf,mean(erp1,2),'Color',[1 0.5 0.2],'LineWidth',2)
plot(xf,mean(erp2,2),'Color',[0 .5 0],'LineWidth',2)
plot(xf,diff,'Color',[.1 .1 .1],'LineWidth',3)
plot(xf,CImean,'Color',[.1 .1 .1],'LineWidth',1)
plot([-300 600],[0 0],'k--')
xlabel('Time in ms','FontSize',16)
ylabel('ERPs in \mu-volts','FontSize',16)
set(gca,'XLim',[-300 600],'XTick',-300:100:600,'YLim',[-15 10],'FontSize',14,'Layer','Top')
box on

%% Figure 4 - we open the hood

figure;set(gcf,'Color','w');hold on
plot(xf,mean(erp1,2),'Color',[1 0.5 0.2],'LineWidth',2)
plot(xf,mean(erp2,2),'Color',[0 .5 0],'LineWidth',2)
% plot results from every participant
plot(xf,erp1-erp2,'Color',[.7 .7 .7],'LineWidth',1)
plot(xf,diff,'Color',[.1 .1 .1],'LineWidth',3)
plot(xf,CImean,'Color',[.1 .1 .1],'LineWidth',1)
plot([-300 600],[0 0],'k--')
xlabel('Time in ms','FontSize',16)
ylabel('ERPs in \mu-volts','FontSize',16)
set(gca,'XLim',[-300 600],'XTick',-300:100:600,'YLim',[-15 10],'FontSize',14,'Layer','Top')
box on

%% Figure 5 - scatterplot at latency of max difference

% find latency of max absolute significant difference after stimulus onset
tmp = diff(pval<0.05 & xf>0);
F = find(abs(diff)==max(abs(tmp)));

dt1 = erp1(F,:)-erp2(F,:); mdt1 = mean(dt1);

figure('Color','w','NumberTitle','off')
hold on
plot([0 2],[0 0],'k:')
xpts = UnivarScatter_nofig(dt1');
scatter(xpts,dt1,70,'k','filled')
plot([0.6 1.4],[mdt1 mdt1],'Color',[0 0 0],'LineWidth',2)
set(gca,'XLim',[0 2],'XTick',[],'FontSize',14)
box on
title('Individual differences','FontSize',20)
ylabel('Differences in \mu-volts','FontSize',16)

if sign(diff) == 1
    p1 = round(mean(dt1>0)*100)/100;
    text(0.1,.5,['P(diff>0)=',num2str(p1)],'FontSize',16)
else
    p1 = round(mean(dt1<0)*100)/100;
    text(0.1,.5,['P(diff<0)=',num2str(p1)],'FontSize',16)
end

%% Figure 1b - bar graph
% I'm making the bar graphs by hand, line by line, because this got to
% hurt. If to make bar graphs, people had to carve them in big pieces of
% oak wood, we would see very few of them.

tmp1 = erp1(F,:);tmp2 = erp2(F,:);
m1 = mean(tmp1);m2 = mean(tmp2); 
sem1 = std(tmp1) / sqrt(numel(tmp1));
sem2 = std(tmp2) / sqrt(numel(tmp2));

figure('Color','w','NumberTitle','off')
hold on

plot([0 3],[0 0],'k','LineWidth',1)

plot([0.7 1.3],[m1 m1],'k','LineWidth',2)
plot([0.7 0.7],[0 m1],'k','LineWidth',2)
plot([1.3 1.3],[0 m1],'k','LineWidth',2)
plot([1 1],[m1 m1-sem1],'k','LineWidth',2)
plot([1 1],[m1 m1+sem1],'k','LineWidth',2)
plot([0.9 1.1],[m1+sem1 m1+sem1],'k','LineWidth',2)
plot([0.9 1.1],[m1-sem1 m1-sem1],'k','LineWidth',2)

plot([1.7 2.3],[m2 m2],'k','LineWidth',2)
plot([1.7 1.7],[0 m2],'k','LineWidth',2)
plot([2.3 2.3],[0 m2],'k','LineWidth',2)
plot([2 2],[m2 m2-sem2],'k','LineWidth',2)
plot([2 2],[m2 m2+sem2],'k','LineWidth',2)
plot([1.9 2.1],[m2+sem2 m2+sem2],'k','LineWidth',2)
plot([1.9 2.1],[m2-sem2 m2-sem2],'k','LineWidth',2)

set(gca,'XLim',[0 3],'YLim',[-6 0],'FontSize',14,'YTick',-6:0)
box on

labels = {'Condition 1';'Condition 2'};
set(gca,'XTick',1:2,'XTickLabel',labels,'XTickLabelRotation',0)

plot(1.5,-5,'r*','MarkerSize',20,'LineWidth',2)

