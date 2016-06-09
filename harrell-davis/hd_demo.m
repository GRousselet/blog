%% hd_demo

% Code to reproduce the figures in my blog post on
% the Harrell-Davis quantile estimator.
%
% Execute each cell in turn by pressing cmd+enter (mac) or ctrl+enter (pc)
% You can navigate between cells by pressing cmd/ctrl + up or down arrow.
%
% Copyright (C) 2016 Guillaume Rousselet - University of Glasgow

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.

%% illustrate age distribution

age = load('age.txt');
sage = sort(age);

xpts = UnivarScatter_nofig(sage');
c = linspace(1,10,numel(age));

figure;set(gcf,'Color','w');
subplot(2,1,1);hold on
hist(age,20)
h = findobj(gca,'Type','patch');
h.FaceColor = [0 .5 .5];
h.EdgeColor = 'w';
box on
set(gca,'LineWidth',2,'XLim',[0 100],'XTick',20:10:80,'YLim',[0 20],'YTick',0:5:20)
set(gca,'FontSize',14,'Layer','Top','Color',[.9 .9 .9])
ylabel('Frequency','FontSize',14)

subplot(2,1,2);hold on
scatter(sage,xpts,70,c,'filled')
box on
set(gca,'LineWidth',2,'XLim',[0 100],'XTick',20:10:80,'YLim',[0.4 1.6],'YTick',[])
set(gca,'FontSize',14,'Layer','Top','Color',[.9 .9 .9])
xlabel('Participants'' age','FontSize',14)

%% illustrate weights
cc = parula(9);

figure;set(gcf,'Color','w');hold on

for R = 1:9
    q=R/10; % quantile
    n=length(age);
    m1=(n+1).*q;
    m2=(n+1).*(1-q);
    vec=1:length(age);
    w=betacdf(vec./n,m1,m2)-betacdf((vec-1)./n,m1,m2);
    y=sort(age);
    ageq(R)=sum(w.*y);

    title('Weights of the Harrell-Davis estimator','FontSize',16)
    plot(w,'Color',cc(R,:),'LineWidth',4)
    box on;set(gca,'LineWidth',2)
    xlabel('Participants 1 to 120','FontSize',14)
    ylabel('weights','FontSize',14)
    set(gca,'FontSize',14,'Layer','Top','Color',[.9 .9 .9])
end

hl = legend('10th quantile','20th quantile','30th quantile','40th quantile',...
    '50th quantile','60th quantile','70th quantile','80th quantile','90th quantile',...
    'Location','North');
set(hl,'Box','off')

%% age results

thetaq = deciles(age); % same results as those stored in ageq from previous cell
fprintf('age deciles (hd) = %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f\n',thetaq)

% compare to Matlab's percentile function:
fprintf('age deciles (Matlab) = %.0f, %.0f, %.0f, %.0f, %.0f, %.0f, %.0f, %.0f, %.0f\n',prctile(age,10:10:90))

% age deciles (hd) = 21.1, 23.3, 29.7, 37.0, 45.3, 56.1, 63.3, 66.6, 70.4
% age deciles (Matlab) = 21, 23, 30, 36, 45, 57, 64, 66, 70

%% age results illustrated - scatterplot

figure;set(gcf,'Color','w'); hold on
scatter(sage,xpts,70,c,'filled')
box on
set(gca,'LineWidth',2,'XLim',[0 100],'XTick',20:10:80,'YLim',[0.4 1.6],'YTick',[])
set(gca,'FontSize',14,'Layer','Top','Color',[.9 .9 .9])
xlabel('Participants'' age','FontSize',14)

for d = 1:9
    h = plot([thetaq(d) thetaq(d)],[0.4 1.6],'k-');
    if d == 5 % highlight the median
        set(h,'LineWidth',2,'LineStyle','-')
    end
end

%% median age + CI

rng(1) % set seed to reproduce analysis
q = 0.5; % median

% percentile bootstrap estimation of standard error method:
nboot = 200; % default is 100
[xhd CI] = hdci(age,q,nboot); % 45.31 [35.89, 54.73]
fprintf('----------------------------\n')
fprintf('Median = %2.2f [%2.2f, %2.2f]\n',xhd,CI)

% percentile bootstrap of hd approach:
nboot = 2000; % that's the default
[xhd CIpb] = hdpbci(age,q,nboot); % 45.31 [38.49, 54.40]
fprintf('Median = %2.2f [%2.2f, %2.2f]\n',xhd,CIpb)

%% =======================================================================
%% SIMULATION ============================================================
%% median CI probability coverage

rng(1) % set seed to reproduce simulation
Npop = 1000000;
pop = chi2rnd(10,Npop,1); % skewed population

figure('Color','w','NumberTitle','off')
[N,X] = hist(pop,50);
% N = number of elements in each of nbin containers
% X = position of the bin centers
N = N./Npop; % convert histogram frequencies to proportions
h = bar(X,N);
% h = findobj(gca,'Type','patch');
h.FaceColor = [0 .5 .5];
h.EdgeColor = 'w';
box on
set(gca,'LineWidth',2,'YLim',[0 0.11],'YTick',0:0.01:0.1)
set(gca,'FontSize',14,'Layer','Top','Color',[.9 .9 .9])
xlabel('Arbitrary units','FontSize',14)
ylabel('Proportion','FontSize',14)

q = 0.5; % median
pop_med = median(pop); % population median = ground truth = 9.3415
pop_mean = mean(pop); % population median = ground truth = 10.0031
% Nsim = 10000; % number of simulations - used in post
Nsim = 1000; % start with a much smaller number
Nt = 50; % number of trials per simulation

% WARNING: CAN TAKE A VERY LONG TIME!
% LOAD RESULTS:
load simres

% OR UNCOMMENT NEXT SECTION TO RUN SIMULATION:

% xhd = zeros(Nsim,1);
% xmean = zeros(Nsim,1);
% xmedian = zeros(Nsim,1);
% CIpbse_hd = zeros(Nsim,2);
% CIpbci_hd = zeros(Nsim,2);
% CIpbci_mean = zeros(Nsim,2);
% CIpbci_median = zeros(Nsim,2);
% simres_pbse_hd = zeros(Nsim,1);
% simres_pbci_hd = zeros(Nsim,1);
% simres_pbci_median = zeros(Nsim,1);
% simres_pbci_mean = zeros(Nsim,1);
%
% for sim = 1:Nsim
%
%     samp = pop(randi(Npop,Nt,1));
%
%     [xhd(sim), CIpbse_hd(sim,:)] = hdci(samp,q,200);
%     simres_pbse_hd(sim) = pop_med > CIpbse_hd(sim,1) && pop_med < CIpbse_hd(sim,2);
%
%     [~, CIpbci_hd(sim,:)] = hdpbci(samp,q,2000);
%     simres_pbci_hd(sim) = pop_med > CIpbci_hd(sim,1) && pop_med < CIpbci_hd(sim,2);
%
%     [~,xmedian(sim),CIpbci_median(sim,:),~] = pbci(samp,1000,0.05,'median');
%     simres_pbci_median(sim) = pop_med > CIpbci_median(sim,1) && pop_med < CIpbci_median(sim,2);
%
%     [~,xmean(sim),CIpbci_mean(sim,:),~] = pbci(samp,1000,0.05,'mean');
%     simres_pbci_mean(sim) = pop_mean > CIpbci_mean(sim,1) && pop_mean < CIpbci_mean(sim,2);
%
% end

fprintf('----------------------------\n')
fprintf('probability coverage:\n')
fprintf('hd: pb(se(hd)) = %2.4f\n',mean(simres_pbse_hd))
fprintf('hd: pbci(hd) = %2.4f\n',mean(simres_pbci_hd))
fprintf('median: pbci(median) = %2.4f\n',mean(simres_pbci_median))
fprintf('mean: pbci(mean) = %2.4f\n',mean(simres_pbci_mean))

% save simres xhd xmean xmedian CIpbse_hd CIpbci_hd CIpbci_mean CIpbci_median ...
%          simres_pbse_hd simres_pbci_hd simres_pbci_median simres_pbci_mean

% probability coverage:
% hd: pb(se(hd)) = 0.9530
% hd: pbci(hd) = 0.9473
% median: pbci(median) = 0.9452
% mean: pbci(mean) = 0.9394

%% histogram of sample medians & sample means
% this part is brilliant: we get to look at the results of 10,000
% experiments! These bootstrap disitributions provide an estimation of the sampling
% distribution of hd, the median, and the mean.

labels = {'xhd';'xmedian';'xmean'};

figure('Color','w','NumberTitle','off')

for sub = 1:3

    subplot(3,1,sub); hold on
    eval(['todo = ',labels{sub},';'])
    [N,X] = hist(todo,50);
    N = N./numel(todo); % convert histogram frequencies to proportions
    h = bar(X,N);
    h.FaceColor = [0 .5 .5];
    h.EdgeColor = 'w';
    box on
    set(gca,'XLim',[6 13],'YLim',[0 0.1])
    set(gca,'FontSize',14,'Layer','Top','Color',[.9 .9 .9])
    ylabel('Proportion','FontSize',14)
    title(labels{sub}(2:end),'FontSize',18)

end

%% histogram of confidence intervals' width
% how do CIs vary in size across our 10,000 experiments?
% The size of the CIs is the difference between their upper & lower bounds.

labels = {'CIpbse_hd';'CIpbci_hd';'CIpbci_median';'CIpbci_mean'};
figure('Color','w','NumberTitle','off')

for sub = 1:4

    subplot(4,1,sub)
    eval(['todo = ',labels{sub},';'])
    diff = todo(:,2) - todo(:,1);
    [N,X] = hist(diff,50);
    N = N./numel(diff); % convert histogram frequencies to proportions
    h = bar(X,N);
    h.FaceColor = [0 .5 .5];
    h.EdgeColor = 'w';
    box on
    set(gca,'LineWidth',2,'YLim',[0 0.08],'YTick',0:0.04:0.08,'XLim',[0 6],'XTick',1:5)
    set(gca,'FontSize',14,'Layer','Top','Color',[.9 .9 .9])
    ylabel('Proportion','FontSize',14)
    title([labels{sub}(3:6),' ',labels{sub}(8:end)],'FontSize',18)

end

% pbse_hd is right shifted compared to pbci_hd, with increasing differences
% from left to right
% pbci_median starts slightly to left of pbci_hd, but then progressively
% shifts toward the right
%% =======================================================================

%% =======================================================================
%% onsets in ms - data from Bieniek et al. EJN 2015

onset = load('onset.txt');
onset = onset(:,3);

figure;set(gcf,'Color','w');
subplot(4,1,1:3);hold on
hist(onset,20)
h = findobj(gca,'Type','patch');
h.FaceColor = [0 .5 .5];
h.EdgeColor = 'w';
box on
set(gca,'LineWidth',2,'XLim',[0 300],'XTick',0:50:300,'YLim',[0 30],'YTick',0:5:25)
set(gca,'FontSize',14,'Layer','Top','Color',[.9 .9 .9])
ylabel('Frequency','FontSize',14)

subplot(4,1,4);hold on
for P = 1:numel(onset)
    plot([onset(P) onset(P)],[0 1],'Color','k')
    axis([0 250 0 1]);
end
box on
set(gca,'LineWidth',2,'XLim',[0 300],'XTick',0:50:300,'YLim',[0 1],'YTick',[])
set(gca,'FontSize',14,'Layer','Top','Color',[.9 .9 .9])
xlabel('Onsets in ms','FontSize',14)

%% onset deciles + confidence intervals - data from Bieniek et al. EJN 2015

nboot = 200;
plotCI = 0; % no plot thank you, we're going to make a different one
[xhd CI] = decilesci(onset,nboot,plotCI);

figure;set(gcf,'Color','w');
for d = 1:9
    line([CI(d,1) CI(d,2)],[d d],'Color','k','Linewidth',1);hold on
end
h=plot(xhd,1:9,'ko','Linewidth',1);hold on
set(h,'MarkerFaceColor','k')
line([xhd(5) xhd(5)],[0 10],'LineStyle','--','Color','k','LineWidth',1);hold on
title('Deciles + confidence intervals','Fontsize',20)
box on
ylabel('Onset deciles','FontSize',16);
xlabel('Onsets in ms','FontSize',16)
axis([50 140 0 10]);
set(gca,'Fontsize',14,'YTick',1:9)%,'YTickLabel',{'1' '' '' '' '5' '' '' '' '9'});
text(xhd(5)+1,4.5,sprintf('median = %2.1f [%2.1f, %2.1f]',xhd(5),CI(5,:)),'FontSize',14)

fprintf('----------------------------\n')
for d = 1:9
    fprintf('decile %i = %2.2f [%2.2f, %2.2f]\n',d,xhd(d),CI(d,:))
end

%% ---------------------------------------------
%% shape of the bootstrap decile distributions?

% get bootstrap estimates
rng(7)
nboot = 2000;
n = numel(onset);

% CI boundaries
alpha = 0.05;
lo = round(nboot*(alpha/2));
hi = nboot - lo;
CI = zeros(9,2);
HDI = zeros(9,2);

list = randi(n,nboot,n); % use same bootstrap samples for all CIs
xboot = zeros(nboot,9);
for d = 1:9
    q = d/10;
    % percentile bootstrap estimate of the standard error of hd
    for B = 1:nboot
        xboot(B,d) = hd(onset(list(B,:)),q);
    end
    HDI(d,:) = hdi(xboot(:,d),1-alpha);
end

xboot = sort(xboot,1);
% confidence intervals
CI(:,1) = xboot(lo+1,:)';
CI(:,2) = xboot(hi,:)';

%% make figure

figure;set(gcf,'Color','w');

for d = 1:9
    subplot(3,3,d); hold on
    hist(xboot(:,d),50)
    h = findobj(gca,'Type','patch');
    h.FaceColor = [0 .5 .5];
    h.EdgeColor = 'w';
    box on
    % set(gca,'LineWidth',2,'XLim',[0 300],'XTick',0:50:300,'YLim',[0 30],'YTick',0:5:25)
    set(gca,'FontSize',14,'Layer','Top','Color',[.9 .9 .9])
    ylabel('Frequency','FontSize',14)
    skewness(xboot(:,d))
    plot([CI(d,1) CI(d,2)],[100 100],'k','LineWidth',3)
    plot([HDI(d,1) HDI(d,2)],[103 103],'g','LineWidth',3)
end
