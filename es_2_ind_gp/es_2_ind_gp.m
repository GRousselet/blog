%% es_2_ind_gp
% This script demonstrates a few examples of robust measures of effect size
% for two independent groups of observations.
% Copyright (C) 2016 Guillaume Rousselet - University of Glasgow

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.

%% Cohen's d - 3 examples with samples from a normal distribution

rng(1) % set seed for reproducible results

Nt = 20; % number of trials in each group
g1 = randn(Nt,1); % group 1 = sample from a normal distribution
g2 = randn(Nt,1); % group 2 = sample from a normal distribution

cst = [0 1 2]; % 3 conditions: constant added to group 2
Nc = numel(cst);

figure('Color','w','NumberTitle','off')

for C = 1:Nc

    subplot(Nc,1,C)
    hold on

    gg2 = g2+cst(C);
    cod = abs(cohend(g1,gg2));

    title(sprintf('Cohen''s d = %.2f',cod),'FontSize',20)
    scatter(g1,ones(Nt,1),[],'filled','MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
    scatter(gg2,ones(Nt,1)*2,[],'filled','MarkerEdgeColor',[1 0.1 0.1],'MarkerFaceColor',[1 0.5 0.2],'LineWidth',1.5)

    set(gca,'YColor','k','FontSize',14,'YLim',[0 3],'YTick',[1 2],'YTickLabel',{'g1';'g2'},'XLim',[-2.5,4])

end

% fileName = 'cohend_3ex.jpg';
% set(gcf,'PaperPositionMode','auto')
% print(gcf, '-djpeg','-r150', fileName);
% close(gcf)

%% Cohen's d - systematic mapping

rng(1) % set seed for reproducible results

Nt = 20; % number of trials in each group
simval = 0:0.01:6; % repeat simulations for various shifts between the two groups

Nst = numel(simval); % number of simulations
hcod = zeros(Nst,1); % vector holding Cohen's d values

for st = 1:Nst % for each simulation

    g1 = randn(Nt,1); % group 1 = sample from a normal distribution
    g2 = randn(Nt,1); % group 2 = sample from a normal distribution
    g1 = g1 - mean(g1); g2 = g2 - mean(g2); % centre each group
    g2 = g2 + simval(st); % shift group 2 by a constant
    hcod(st) = abs(cohend(g1,g2));

end

figure('Color','w','NumberTitle','off')
hold on
plot([0 7],[0 7],'k:','Linewidth',1)
cm = linspace(1,10,Nst); % colour map for scatterplot
scatter(simval,hcod,[],cm,'filled')
box on
set(gca,'FontSize',14,'Layer','Top','XLim',[0 6],'YLim',[0 7])
xlabel('Difference in means','FontSize',16)
ylabel('Cohen''s d','FontSize',16)

% fileName = 'cohend_sysmap.jpg';
% set(gcf,'PaperPositionMode','auto')
% print(gcf, '-djpeg','-r150', fileName);
% close(gcf)

%% Cohen's d - examples with outliers

rng(1) % set seed for reproducible results

Nt = 20; % number of trials in each group
g1 = randn(Nt,1); % group 1 = sample from a normal distribution
g2 = randn(Nt,1); % group 2 = sample from a normal distribution

cst = 2; % constant added to group 2
outliers = [0 1 2 5 9]; % constant added to most extreme value from group 2
Nc = numel(outliers);

figure('Color','w','NumberTitle','off')

for C = 1:Nc

    subplot(Nc,1,C)
    hold on

    gg2 = sort(g2 + cst);
    gg2(Nt) = gg2(Nt) + outliers(C);
    cod = abs(cohend(g1,gg2));

    title(sprintf('Cohen''s d = %.2f',cod),'FontSize',20)
    scatter(g1,ones(Nt,1),[],'filled','MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
    scatter(gg2,ones(Nt,1)*2,[],'filled','MarkerEdgeColor',[1 0.1 0.1],'MarkerFaceColor',[1 0.5 0.2],'LineWidth',1.5)

    set(gca,'YColor','k','FontSize',14,'YLim',[0 3],'YTick',[1 2],'YTickLabel',{'g1';'g2'},'XLim',[-2.5,13])

end

% fileName = 'cohend_outliers.jpg';
% set(gcf,'PaperPositionMode','auto')
% print(gcf, '-djpeg','-r150', fileName);
% close(gcf)

%% Cohen'd - systematic relationship to outliers' size

rng(1) % set seed for reproducible results

Nt = 20; % number of trials in each group
outliers = 0:0.1:10; % constant added to most extreme value from group 2
Nc = numel(outliers); % number of simulations
hcod = zeros(Nc,1); % vector holding Cohen's d values
meandiff = zeros(Nc,1); % vector holding differences between group means
cst = 2; % constant added to group 2

g1 = randn(Nt,1); % group 1 = sample from a normal distribution
g2 = randn(Nt,1); % group 2 = sample from a normal distribution
g1 = g1 - mean(g1); g2 = g2 - mean(g2); % centre each group -> mean = 0
g1 = g1 * sqrt( 1 / var(g1)); % set SD to 1
g2 = g2 * sqrt( 1 / var(g2)); % set SD to 1

g2 = sort(g2 + cst); % shift group 2 by constant

for C = 1:Nc % for outliers of increasing size

    gg2 = g2;
    gg2(Nt) = gg2(Nt) + outliers(C);
    hcod(C) = abs(cohend(g1,gg2));
    meandiff(C) = abs(mean(g1) - mean(gg2));

end

figure('Color','w','NumberTitle','off')
cm = linspace(1,10,Nc); % colour map for scatterplot
scatter(meandiff,hcod,[],cm,'filled')
box on
set(gca,'FontSize',14,'Layer','Top','XTick',2:.1:2.5,'XLim',[1.9, 2.6])
xlabel('Difference between group means','FontSize',16)
ylabel('Cohen''s d','FontSize',16)

% fileName = 'cohend_sysout.jpg';
% set(gcf,'PaperPositionMode','auto')
% print(gcf, '-djpeg','-r150', fileName);
% close(gcf)

%% Cohen's d - normal distribution vs. skewed distribution

rng(20) % set seed for reproducible results

Nt = 100;
cst = 2; % constant added to group 2
g1 = randn(Nt,1);
g2 = randn(Nt,1);

% centre each group, then add constant to group 2
g1 = g1 - mean(g1); % mean = 0
g2 = g2 - mean(g2) + cst; % mean = cst

% mixture of samples from 2 normal distributions that differ
% in variance & skewness. We replace the top 10% of observations in g2:
mg2 = sort(g2);
tmp = sort(abs(randn(Nt,1)*4));
mg2(Nt*.9+1:end) = tmp(Nt*.9+1:end);

% kernel density estimates
kg1 = akerd(g1);
kg2 = akerd(g2);
kmg2 = akerd(mg2);

figure('Color','w','NumberTitle','off','Units','Normalized','Position',[0 1 0.4 0.2])

subplot(1,2,1);hold on % 2 normal distributions

cod = abs(cohend(g1,g2));

title(sprintf('Cohen''s d = %.2f',cod),'FontSize',20)
plot(sort(g1),kg1,'LineWidth',2)
plot(sort(g2),kg2,'LineWidth',2)
set(gca,'FontSize',14)
set(gca,'FontSize',14,'XLim',[-8 8],'XTick',-8:2:8)
xlabel('x','FontSize',16)
ylabel('Density function','FontSize',16)

subplot(1,2,2);hold on % 1 normal distribution, 1 mixed normal distribution

cod = abs(cohend(g1,mg2));

title(sprintf('Cohen''s d = %.2f',cod),'FontSize',20)
plot(sort(g1),kg1,'LineWidth',2)
plot(sort(mg2),kmg2,'LineWidth',2)
set(gca,'FontSize',14,'XLim',[-8 8],'XTick',-8:2:8)
xlabel('x','FontSize',16)
ylabel('Density function','FontSize',16)

% fileName = 'cohend_mixed.jpg';
% set(gcf,'PaperPositionMode','auto')
% print(gcf,'-djpeg','-r150', fileName);
% close(gcf)

%% revisit outlier example using Cliff's delta

rng(1) % set seed for reproducible results

Nt = 20; % number of trials in each group
g1 = randn(Nt,1); % group 1 = sample from a normal distribution
g2 = randn(Nt,1); % group 2 = sample from a normal distribution

cst = 2; % constant added to group 2
outliers = [0 1 2 5 9]; % constant added to most extreme value from group 2
Nc = numel(outliers);

figure('Color','w','NumberTitle','off')

for C = 1:Nc

    subplot(Nc,1,C)
    hold on

    gg2 = sort(g2 + cst);
    gg2(Nt) = gg2(Nt) + outliers(C);
    cliffd = abs(cliffdelta(g1,gg2));

    title(sprintf('Cliff''s delta = %.2f',cliffd),'FontSize',20)
    scatter(g1,ones(Nt,1),[],'filled','MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
    scatter(gg2,ones(Nt,1)*2,[],'filled','MarkerEdgeColor',[1 0.1 0.1],'MarkerFaceColor',[1 0.5 0.2],'LineWidth',1.5)

    set(gca,'YColor','k','FontSize',14,'YLim',[0 3],'YTick',[1 2],'YTickLabel',{'g1';'g2'},'XLim',[-2.5,13])

end

%% revisit outlier example using mutual information

rng(1) % set seed for reproducible results

Nt = 20; % number of trials in each group
g1 = randn(Nt,1); % group 1 = sample from a normal distribution
g2 = randn(Nt,1); % group 2 = sample from a normal distribution

cst = 2; % constant added to group 2
outliers = [0 1 2 5 9]; % constant added to most extreme value from group 2
Nc = numel(outliers);

categ = [zeros(Nt,1);ones(Nt,1)]; % declare labels for the two groups

figure('Color','w','NumberTitle','off')

for C = 1:Nc

    subplot(Nc,1,C)
    hold on

    gg2 = sort(g2 + cst);
    gg2(Nt) = gg2(Nt) + outliers(C);
    mi = gcmi_cd([g1;gg2], categ, 2);

    title(sprintf('Mutual information = %.2f',mi),'FontSize',20)
    scatter(g1,ones(Nt,1),[],'filled','MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
    scatter(gg2,ones(Nt,1)*2,[],'filled','MarkerEdgeColor',[1 0.1 0.1],'MarkerFaceColor',[1 0.5 0.2],'LineWidth',1.5)

    set(gca,'YColor','k','FontSize',14,'YLim',[0 3],'YTick',[1 2],'YTickLabel',{'g1';'g2'},'XLim',[-2.5,13])

end

%% when units matter: all pairwise differences

rng(1) % set seed for reproducible results

Nt = 20; % number of trials in each group
g1 = randn(Nt,1); % group 1 = sample from a normal distribution
g2 = lognrnd(0,1,Nt,1); % group 2 = sample from a normal distribution
% g1 = g1 - mean(g1); g2 = g2 - mean(g2); % centre groups
g1 = sort(g1);

% all pairwise differences - as implemented in wmwloc
yy = repmat(g2,[1 length(g1)])';
xx = repmat(g1,[1 length(g2)]);
alldiff = xx-yy;

% Harrell-Davis estimate of the median of all pairwise differences
md_alldiff = hd(alldiff(:));
q1_alldiff = hd(alldiff(:),0.25); % 1st quartile
q3_alldiff = hd(alldiff(:),0.75); % 3rd quartile

% kernel density estimation of the distribution of differences
kalldiff = akerd(alldiff(:));

% add jitter for illustration
jt = 0.1; % amount of jitter for scatter plot
jg1 = ones(Nt,1)+randn(Nt,1)*jt;
jg2 = ones(Nt,1)*2+randn(Nt,1)*jt;

figure('Color','w','NumberTitle','off')

subplot(1,3,1); hold on
for sub=1:Nt
    plot([1 jg2(sub)],[max(g1) g2(sub)],'r')
end
scatter(jg1(1:Nt-1),g1(1:Nt-1),[],'filled')
scatter(1,max(g1),[],'filled','MarkerFaceColor','r','LineWidth',1.5)
scatter(jg2,g2,[],'filled')
set(gca,'FontSize',14,'XTick',1:2,'XLim',[0 3])
box on
xlabel('Groups','FontSize',16)
ylabel('Observations in arbitrary units','FontSize',16)

subplot(1,3,[2 3]); hold on
plot(sort(alldiff(:)),kalldiff,'k','LineWidth',3)
set(gca,'FontSize',14)
box on
xlabel('Pairwise differences','FontSize',16)
ylabel('Density','FontSize',16)
v=axis;
plot([md_alldiff md_alldiff],[v(3) v(4)],'k') % plot median
plot([q1_alldiff q1_alldiff],[v(3) v(4)],'k--') % plot 1st quartile
plot([q3_alldiff q3_alldiff],[v(3) v(4)],'k--') % plot 3rd quartile

% fileName = 'all_pairwise_differences.jpg';
% set(gcf,'PaperPositionMode','auto')
% print(gcf, '-djpeg','-r150', fileName);
% close(gcf)

%% example with two samples differing in variance

rng(20) % set seed for reproducible results

Nt = 100;
g1 = randn(Nt,1);
g2 = randn(Nt,1)*5;

% centre each group
g1 = g1 - mean(g1); % mean = 0
g2 = g2 - mean(g2); % mean = 0

% kernel density estimates
kg1 = akerd(g1);
kg2 = akerd(g2);

% Cohen's d
cod = abs(cohend(g1,g2)); % obviously Cohen's d = 0 here!

% Cliff's delta
cliffd = abs(cliffdelta(g1,g2)); % delta is very near 0

% mutual information
categ = [zeros(Nt,1);ones(Nt,1)]; % declare labels for the two groups
mi = gcmi_cd([g1;g2], categ, 2);

% Kolmogorov-Smirnov statistic
ks = ksstat(g1,g2);

% Q
q = qhat(g1,g2);

figure('Color','w','NumberTitle','off')
hold on

title(sprintf('Cd=%.0f, delta=%.0f\n MI=%.2f, KS=%.2f, Q=%.2f',cod,cliffd,mi,ks,q),'FontSize',20)
plot(sort(g1),kg1,'LineWidth',2)
plot(sort(g2),kg2,'LineWidth',2)
set(gca,'FontSize',14)
xlabel('x','FontSize',16)
ylabel('Density function','FontSize',16)
% set(gca,'FontSize',14,'XLim',[-8 8],'XTick',-8:2:8)

% fileName = 'vardiffexample.jpg';
% set(gcf,'PaperPositionMode','auto')
% print(gcf, '-djpeg','-r150', fileName);
% close(gcf)

%% effect sizes as a function of variance differences
% WARNING: Q's bias correction is very slow to compute
% you can save time by loading the pre-computed results

load_res = 1; % load results

rng(20) % set seed for reproducible results

Nt = 100;
g1 = randn(Nt,1);
g2 = randn(Nt,1);

% centre each group
g1 = g1 - mean(g1); % mean = 0
g2 = g2 - mean(g2); % mean = 0

varm = 1:.1:10; % variance modifier
Nc = numel(varm);

if load_res == 0
    
    hvar = zeros(Nc,1);
    mi = zeros(Nc,1);
    ks = zeros(Nc,1);
    q = zeros(Nc,1);
    
    categ = [zeros(Nt,1);ones(Nt,1)]; % declare labels for the two groups
    
    for C = 1:Nc
        
        gg2 = g2 * varm(C); % increase variance of group 2
        hvar(C) = abs(var(g1) - var(gg2));
        
        mi(C) = gcmi_cd([g1;gg2], categ, 2); % mutual information
        ks(C) = ksstat(g1,gg2); % Kolmogorov-Smirnov statistic
        q(C) = qhat(g1,gg2); % Q
        
    end
    % save varres hvar mi ks q
else
    load varres
end

cm = linspace(1,10,Nc); % colour map for scatterplot

figure('Color','w','NumberTitle','off','Units','Normalized','Position',[0 1 0.5 0.3])

subplot(1,3,1);hold on % KS
scatter(hvar,ks,[],cm,'filled')
box on
set(gca,'FontSize',14,'Layer','Top','XTick',0:10:100,'XLim',[-1, 100])
xlabel('Absolute variance difference','FontSize',16)
ylabel('Kolmogorov-Smirnov statistics','FontSize',16)

subplot(1,3,2);hold on % Q
scatter(hvar,q,[],cm,'filled')
box on
set(gca,'FontSize',14,'Layer','Top','XTick',0:10:100,'XLim',[-1, 100])
xlabel('Absolute variance difference','FontSize',16)
ylabel('Q','FontSize',16)
% apparent probability of correct classification

subplot(1,3,3);hold on % MI
scatter(hvar,mi,[],cm,'filled')
box on
set(gca,'FontSize',14,'Layer','Top','XTick',0:10:100,'XLim',[-1, 100])
xlabel('Absolute variance difference','FontSize',16)
ylabel('Mutual information','FontSize',16)

% fileName = 'vardiff_map.jpg';
% set(gcf,'PaperPositionMode','auto')
% print(gcf, '-djpeg','-r150', fileName);
% close(gcf)

%% Comparisons of effect sizes: examples with kernel density estimates

rng(20) % set seed for reproducible results

Nt = 100;
g1 = randn(Nt,1);
g2 = randn(Nt,1);

% centre each group
g1 = g1 - mean(g1); % mean = 0
g2 = g2 - mean(g2); % mean = 0

cst = [0 1 2 3 5]; % constant added to group 2
Nc = numel(cst);
categ = [zeros(Nt,1);ones(Nt,1)]; % declare labels for the two groups

% kernel density estimates
kg1 = akerd(g1);

figure('Color','w','NumberTitle','off','Units','Normalized','Position',[0 1 0.5 0.2])

for C = 1:Nc

    gg2 = g2 + cst(C);

    % kernel density estimates
    kg2 = akerd(gg2);

    cliffd = abs(cliffdelta(g1,gg2)); % Cliff's delta
    mi = gcmi_cd([g1;gg2], categ, 2); % mutual information
    ks = ksstat(g1,gg2); % Kolmogorov-Smirnov statistic
    q = qhat(g1,gg2); % Q

    subplot(1,Nc,C);hold on

    title(sprintf('Mean diff = %.0f\n delta=%.2f MI=%.2f\n KS=%.2f Q=%.2f',cst(C),cliffd,mi,ks,q),'FontSize',20)

    plot(sort(g1),kg1,'LineWidth',2)
    plot(sort(gg2),kg2,'LineWidth',2)
    set(gca,'FontSize',14)
    set(gca,'FontSize',14,'XLim',[-4 8],'XTick',-4:2:8)
    xlabel('x','FontSize',16)
    ylabel('Density function','FontSize',16)

end

% fileName = 'escomp_kde.jpg';
% set(gcf,'PaperPositionMode','auto')
% print(gcf, '-djpeg','-r150', fileName);
% close(gcf)

%% Systematic comparisons of effect sizes
%% PART 1 - generate data - jump to next cell to generate figure
% WARNING: very slow because of Q

rng(20) % set seed for reproducible results
Nt = 100;
categ = [zeros(Nt,1);ones(Nt,1)];
simval = 0:0.01:6;
Nc = numel(simval);
hmi = zeros(Nc,1);
hks = zeros(Nc,1);
hcliffd = zeros(Nc,1);
hq = zeros(Nc,1);

for C = 1:Nc
    g1 = randn(Nt,1);
    g2 = randn(Nt,1);
    g1 = g1 - mean(g1); g2 = g2 - mean(g2); % centre each group
    g2 = g2 + simval(C); % shift group 2 by a constant

    hcliffd(C) = abs(cliffdelta(g1,g2)); % Cliff's delta
    hmi(C) = gcmi_cd([g1;g2], categ, 2); % mutual information
    hks(C) = ksstat(g1,g2); % Kolmogorov-Smirnov statistic
    hq(C) = qhat(g1,g2); % Q
end

save sysres hcliffd hmi hks hq Nc

%% PART 2 - make figure: estimator as a function of difference in means

load sysres
cm = linspace(1,10,Nc); % colour map for scatterplot

figure('Color','w','NumberTitle','off','Units','Normalized','Position',[0 1 0.5 0.2])

subplot(1,4,1); hold on % Cliff's delta
scatter(simval,hcliffd,[],cm,'filled')
ylabel('Cliff''s delta','FontSize',16)

subplot(1,4,2); hold on % mutual information
scatter(simval,hmi,[],cm,'filled')
ylabel('Mutual information','FontSize',16)

subplot(1,4,3); hold on % Kolmogorov-Smirnov statistics
scatter(simval,hks,[],cm,'filled')
ylabel('Kolmogorov-Smirnov statistics','FontSize',16)

subplot(1,4,4); hold on % Q
scatter(simval,hq,[],cm,'filled')
ylabel('Q','FontSize',16)

for sub = 1:4
    subplot(1,4,sub)
    box on
    set(gca,'FontSize',14,'Layer','Top','XLim',[0 6],'YLim',[0 1])
    xlabel('Difference in means','FontSize',16)
    if sub == 4
        set(gca,'YLim',[.45 1])
    end
end

% fileName = 'escomp_diffmean.jpg';
% set(gcf,'PaperPositionMode','auto')
% print(gcf, '-djpeg','-r150', fileName);
% close(gcf)


%% PART 3 - make figure: estimator comparisons

nfun = {'hmi';'hcliffd';'hks';'hq'};
nnam = {'MI';'Delta';'KS';'Q'};
nc1 = numel(nfun);
nc2 = nc1-1;

load sysres
cm = linspace(1,10,Nc); % colour map for scatterplot

figure('Color','w','NumberTitle','off','Units','Normalized','Position',[0 1 0.25 0.5])

for cond1 = 1:nc1

    nsub = setdiff(1:nc1,cond1); % 3 remaining conditions

    for cond2 = 1:nc2

        subplot(nc1,nc2,cond1*nc2-nc2+cond2); hold on % MI
        eval(['datx=',nfun{nsub(cond2)},';'])
        eval(['daty=',nfun{cond1},';'])
        scatter(datx,daty,[],cm,'filled')
        xlabel(nnam{nsub(cond2)},'FontSize',16)
        box on
        set(gca,'FontSize',14,'Layer','Top','XLim',[0 1],'YLim',[0 1])
        set(gca,'XTick',0:.2:1,'YTick',0:.2:1)

        if cond2 == 1
            ylabel(nnam{cond1},'FontSize',16,'FontWeight','bold')
        end

        if cond2 == 3 && cond1 < 4
            set(gca,'XLim',[.5 1])
            set(gca,'XTick',.5:.1:1)
        end

        if cond1 == 4
            set(gca,'YLim',[.5 1])
            set(gca,'YTick',.5:.1:1)
        end

        if cond1 == 1
            set(gca,'YLim',[-0.01 1])
            set(gca,'YTick',0:.2:1)
        end

        if cond1 > 1 && cond2 == 1
            set(gca,'XLim',[-0.01 1])
            set(gca,'XTick',0:.2:1)
        end

        ax = gca; ax.Color = [.9 .9 .9];

    end

end

% fileName = 'escomp_sys.jpg';
% set(gcf,'PaperPositionMode','auto')
% print(gcf, '-djpeg','-r150', fileName);
% close(gcf)
