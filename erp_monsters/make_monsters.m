% Code to generate data similar to those used in my blog post.
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
%% pink noise generator
% http://www.mathworks.com/matlabcentral/fileexchange/42919-pink--red--blue-and-violet-noise-generation-with-matlab-implementation/content/pinknoise.m
%% LIMO EEG toolbox:
% www.hindawi.com/journals/cin/2011/831409/
% https://gforge.dcn.ed.ac.uk/gf/project/limo_eeg/
% https://github.com/LIMO-EEG-Toolbox/limo_eeg
%% =========================================================================

%% get data
load erps

%% reformat
pool = zeros(40,451,192);
count = 0;
for P = 1:20
   for C = 1:2
       count = count+1;
       pool(count,:,:) = erps{P,C};
   end
end

pool = reshape(permute(pool, [2 3 1]),[451 192*40]);

% figure;hold on
% plot(xf,mean(pool,2),'g','LineWidth',3)
% plot(xf,mean(pool,2) + pinknoise(451)')

%% make monsters

Np = 30; % 30 participants
Nt = 100; % 100 trials
Nf = 451; % time points
total = size(pool,2);
pval = ones(Nf,1);

erp1 = zeros(Nf,Np);
erp2 = zeros(Nf,Np);

while sum(pval(xf>0) < 0.05) < 5
    
    for P = 1:Np
       erp1(:,P) = mean(pool(:,randi(total,Nt,1)) + reshape(pinknoise(Nf*Nt),[Nf,Nt]),2); 
       erp2(:,P) = mean(pool(:,randi(total,Nt,1)) + reshape(pinknoise(Nf*Nt),[Nf,Nt]),2); 
    end
    
    todo1 = zeros(1,Nf,Np); % reformat to electrode x time x participant
    todo1(1,:,:) = erp1;
    todo2 = zeros(1,Nf,Np); % reformat to electrode x time x participant
    todo2(1,:,:) = erp2;
    [~,diff,~,CImean,pval,~,~] = limo_yuend_ttest(todo1,todo2,0,0.05);
    
end

CImean = squeeze(CImean);
% save monsters pool erp1 erp2 diff CImean pval xf

%%
figure;set(gcf,'Color','w');hold on
plot(xf,erp1-erp2,'Color',[.7 .7 .7],'LineWidth',1)
plot(xf,mean(erp1,2),'Color',[1 0.5 0.2],'LineWidth',2)
plot(xf,mean(erp2,2),'Color',[0 .5 0],'LineWidth',2)
plot(xf,diff,'Color',[.1 .1 .1],'LineWidth',4)
plot(xf,CImean,'Color',[.1 .1 .1],'LineWidth',2)
plot([-300 600],[0 0],'k--')
xlabel('Time in ms','FontSize',16)
ylabel('ERPs in \mu-volts','FontSize',16)
set(gca,'XLim',[-300 600],'XTick',-300:100:600,'YLim',[-15 10],'FontSize',14,'Layer','Top')
box on
title('monsters in the closet','FontSize',20)

tmp = diff(pval<0.05 & xf>0);
F = find(abs(diff)==max(abs(tmp)));
plot([xf(F) xf(F)],[-15 10],'k:')

% scatterplot
dt1 = erp1(F,:)-erp2(F,:); mdt1 = mean(dt1);

figure('Color','w','NumberTitle','off')
hold on

plot([0 2],[0 0],'k:')

xpts = UnivarScatter_nofig(dt1');
scatter(xpts,dt1,70,'k','filled')
plot([0.6 1.4],[mdt1 mdt1],'Color',[0 0 0],'LineWidth',2)

% xpts = UnivarScatter_nofig(dt2');
% scatter(xpts+1,dt2,70,'k','filled')
% plot([1.6 2.4],[mdt2 mdt2],'Color',[0 0 0],'LineWidth',2)

set(gca,'XLim',[0 2],'XTick',[],'FontSize',14)
box on

title('Differences','FontSize',20)
% labels = {'Time 1';'Time 2'};
% set(gca,'XTickLabel',labels,'XTickLabelRotation',0)

if sign(diff) == 1
p1 = round(mean(dt1>0)*100)/100;
text(0.1,.5,['P(diff>0)=',num2str(p1)],'FontSize',16)
else
    p1 = round(mean(dt1<0)*100)/100;
text(0.1,.5,['P(diff<0)=',num2str(p1)],'FontSize',16)
end
