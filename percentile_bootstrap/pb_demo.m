% pb_demo

% Code to illustrate the percentile bootstrap procedure and to
% reproduce the figures in my blog post.
%
% Execute each cell in turn by pressing cmd+enter (mac) or ctrl+enter (pc)
% You can navigate between cells by pressing cmd/ctrl + up or down arrow.
%
% Copyright (C) 2016 Guillaume Rousselet - University of Glasgow

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.

%% dependencies
% certain sections of this script require functions available here:
% https://github.com/GRousselet/matlab_stats

% to add structured jitter to the scatterplots, I've used this toolbox:
% http://www.mathworks.com/matlabcentral/fileexchange/54243-univarscatter
% The code is also available in the blog folder, for convenience, and
% because I've created an alternative version which only outputs the
% jittered x-axis values: function `UnivarScatter_nofig.m`

%% data from Harvey Motulsky - 2014, figure 5
% http://www.ncbi.nlm.nih.gov/pubmed/25204545

dat = [7.43899808325858
    14.1382674659797
    10.315154950167
    5.03032058297202
    11.2291896180685
    6.1092033148039
    18.82393464998
    21.011591225343
    17.459750710849
    12.7951869578105
    8.37435586003159
    18.0576113571032
    8.51925309940913
    9.828832278548
    5.72571965977189
    11.4624386473698
    9.9857126117517
    7.45283469145073
    9.2041853529859
    9.13869067011597];

nobs = numel(dat); % number of observations

%% percentile bootstrap demo
% this code could be vectorised to speed calculations,
% but having the loop helps understand the procedure I think.
% Similar code is implemented in the pbci & bootse functions.

Nb = 1000; % number of samples with replacement
boot_est = zeros(Nb,1); % vector that will collect bootstrap samples
rng(4); % set seed to reproduce figure from blog

% loop that implements the core of the bootstrap mechanism
for B = 1:Nb % bootstrap sampling

    % sample nobs with replacement among 1:nobs alternatives

    % in three steps -----------------------------------------------
    %     bootindex = randi(nobs,nobs,1); % create bootstrap indices
    %     bootsample = dat(bootindex); % get bootstrap samples
    %     boot_est(B) = median(bootsample); % compute bootstrap estimate
    % --------------------------------------------------------------

    % in one step
    boot_est(B) = median(dat(randi(nobs,nobs,1)));
%     boot_est(B) = mean(dat(randi(nobs,nobs,1)));
%     boot_est(B) = hd(dat(randi(nobs,nobs,1))); Harrell-Davis estimate of the 0.5 quantile

end

% We estimate the standard error of the test statistic by the
% standard deviation of the bootstrap replications
% Efron & Tibshinari 1993, chapter 6
% Wilcox 2005, p.44-45
boot_se = std(boot_est,0); % normalize by (n-1)

% Estimate of the sampling distribution of the median
dhat = akerd(boot_est); % kernel density estimation
figure('Color','w','NumberTitle','off')
plot(sort(boot_est),dhat,'LineWidth',2,'Color',[1 0.5 0.2])
set(gca,'FontSize',14)
xlabel('Bootstrap median estimates','FontSize',16)
ylabel('Density','FontSize',16)

% the distribution is particularly irregular, because of the
% non-linearities introduced by the median. If you repeat the
% calculations using the mean, the distribution is much smoother.

% a plot of the sorted bootstrap estimates helps understand
% what happens with the median:
figure('Color','w','NumberTitle','off')
plot(sort(boot_est),'LineWidth',2,'Color',[1 0.5 0.2])
set(gca,'FontSize',14)
xlabel('Sorted bootstrap samples','FontSize',16)
ylabel('Bootstrap median estimates','FontSize',16)

% The Harrell-Davis estimate of the 0.5 quantile also gives a much smoother
% distribution - but that's a topic for another post.

% get 95% confidence interval:
alpha = .05;
lo = round(Nb.*alpha./2);
hi = Nb - lo;
lo = lo+1;
bootsort = sort(boot_est); % sort in ascending order
ci1 = bootsort(lo);
ci2 = bootsort(hi);

% illustrate 95% confidence interval
figure('Color','w','NumberTitle','off'); hold on
tmp = median(dat);
plot([tmp tmp],[0 akerd(boot_est,tmp)],'k','LineWidth',3) % plot median
text(tmp,akerd(boot_est,tmp)+0.05,'sample median','FontSize',14)
plot([ci1 ci1],[0 akerd(boot_est,ci1)],'k','LineWidth',3) % plot CI lower bound
text(ci1-1.5,akerd(boot_est,ci1),{'CI';'lower';'bound'},'FontSize',14)
plot([ci2 ci2],[0 akerd(boot_est,ci2)],'k','LineWidth',3) % plot CI upper bound
text(ci2,akerd(boot_est,ci2)+0.05,{'CI';'upper';'bound'},'FontSize',14)
plot(sort(boot_est),dhat,'LineWidth',4,'Color',[1 0.5 0.2])
set(gca,'FontSize',14)
xlabel('Bootstrap median estimates','FontSize',16)
ylabel('Density','FontSize',16)
box on

% get bootstrap p value:
mu=0; % null hypothesis
pval = mean(boot_est > mu) + mean(boot_est == mu).*0.5;
pval = 2.*min(pval,1-pval);
% in that case, the p value is zero, because there is no overlap
% whatsoever between the null hypothesis and the bootstrap sampling
% distribution. Also, the boostrap assumes that no values other than
% those in the sample at hand can ever be observed. This is a clear
% limitation of the bootstrap - see more details here:
% http://www.sumsar.net/blog/2015/04/the-non-parametric-bootstrap-as-a-bayesian-model/

%% percentile bootstrap CI of median can be unstable

% We compute a median CI several times using the same data.
% Because of strong linearities, even with 1000 bootstrap
% samples, the CI can differ a lot among repetitions.
% Increasing the number of observations (try data = dat) or increasing the number of
% bootrstrap samples (try Nb=10000) seem to improve the results.

data = dat(1:2:end); % data = dat
Nb = 1000; % Nb = 10000;
est = 'median';
% rng(4); % set seed to reproduce figure from blog

% If you repeat the example several times (execute the cell multiple times),
% you will see that the confidence intervals can vary a lot.
% If you use the mean instead of the median
% you will see that the confidence intervals are much more stable.
% More generally, reliable CIs depend on 4 factors:
%   - technique to compute the CI
%   - estimator (mean, median etc.)
%   - shape of the distribution
%   - number of observations
% So when you read an article in which a CI is reported without details,
% it implies that a classic CI formula on means was used. It also implies
% that the CI is probably incorrect...

figure('Color','w','NumberTitle','off')
hold on

for out = 1:7

    [BEST,EST,CI,p] = pbci(data,Nb,alpha,est);
    mdata = EST;

    % scatter
    xpts = UnivarScatter_nofig(data);
    scatter(xpts+out-1,data,70,[.7 .7 .7],'filled')
%     scatter(out-1+ones(n,1)+jitter/10,data,70,[.7 .7 .7],'filled')

    % median & quartiles
    plot(out,mdata,'ko','MarkerSize',10,'MarkerFaceColor','k')
    plot([out out],[mdata CI(2)],'k','LineWidth',2) % whisker 1
    plot([out out],[CI(1) mdata],'k','LineWidth',2) % whisker 2

     % reference lines
    if out == 1
    plot([0 8],[CI(1) CI(1)],'k--')
    plot([0 8],[mdata mdata],'k:')
    plot([0 8],[CI(2) CI(2)],'k--')
    end

    set(gca,'XLim',[0 8],'YLim',[0 30],'FontSize',14)
    box on

end

xlabel('Repetitions of CI calculations','FontSize',16)

%% compute summary statistics for Motulsky's data

alpha = 0.05;

mdat = mean(dat);
y = prctile(dat,[25 50 75]);
q1dat = y(1); % 1st quartile
q2dat = y(2); % median
q3dat = y(3); % 3rd quartile
sddat = std(dat); % standard deviation
semdat = std(dat) / sqrt(nobs); % standard error of the mean

[H,P,CI,STATS] = ttest(dat,alpha); % classic 95% confidence interval
% using the formula instead:
% df = nobs-1;
% CI(1) = mean(dat) + tinv(.025,df) .* std(dat)./sqrt(nobs); % tinv = inverse Student's T cdf
% CI(2) = mean(dat) - tinv(.025,df) .* std(dat)./sqrt(nobs);

mindat = min(dat);
maxdat = max(dat);

Nb = 1000;
est = 'mean';
[BEST,EST,bootCI,p] = pbci(dat,Nb,alpha,est);

%% replicate figure 5 from Motulsky's 2014 paper
% I've added a percentile bootstrap CI + reference line, which shows
% that mean +/- SEM does not include the median.
% Common Misconceptions about Data Analysis and Statistics
% http://www.ncbi.nlm.nih.gov/pubmed/25204545

% I don't use error bars so i don't have functions for that.
% There is an errorbar function in Matlab, but the name freaks me out.
% I don't use boxplots much either, and once you put a boxplot in a figure,
% the formatting gets messy.
% So, we're going to make that figure by hand!
% Messy code but at least you'll see how it can be done.

rng(4); % set seed to reproduce figure from blog

figure('Color','w','NumberTitle','off')
hold on

% scatter + median

% simple way to add jitter to x-axis:
% scatter(ones(nobs,1)+randn(nobs,1)/10,dat,70,'k','filled')

% better way to add jitter to x-axis:
xpts = UnivarScatter_nofig(dat);
scatter(xpts,dat,70,'k','filled')
plot([0.6 1.4],[q2dat q2dat],'Color',[0 0 0],'LineWidth',2)
plot([0 8],[q2dat q2dat],'Color',[0 0 0],'LineWidth',1,'LineStyle',':')

% boxplot - by hand! (no outlier flagged)
plot([1.7 2.3],[q1dat q1dat],'k','LineWidth',2) % 1st quartile
plot([1.7 2.3],[q2dat q2dat],'k','LineWidth',2) % 2nd quartile
plot([1.7 2.3],[q3dat q3dat],'k','LineWidth',2) % 3rd quartile
plot([1.7 1.7],[q1dat q3dat],'k','LineWidth',2) % box
plot([2.3 2.3],[q1dat q3dat],'k','LineWidth',2)
plot([2 2],[q3dat maxdat],'k','LineWidth',2) % whisker 1
plot([1.9 2.1],[maxdat maxdat],'k','LineWidth',2)
plot([2 2],[mindat q1dat],'k','LineWidth',2) % whisker 2
plot([1.9 2.1],[mindat mindat],'k','LineWidth',2)

% quartiles
plot(3,q2dat,'ko','MarkerSize',10,'MarkerFaceColor','k')
plot([3 3],[q2dat q3dat],'k','LineWidth',2) % whisker 1
plot([2.9 3.1],[q3dat q3dat],'k','LineWidth',2)
plot([3 3],[q1dat q2dat],'k','LineWidth',2) % whisker 2
plot([2.9 3.1],[q1dat q1dat],'k','LineWidth',2)

% mean +/- SD
plot(4,mdat,'ko','MarkerSize',10,'MarkerFaceColor','k')
plot([4 4],[mdat mdat+sddat],'k','LineWidth',2) % whisker 1
plot([3.9 4.1],[mdat+sddat mdat+sddat],'k','LineWidth',2)
plot([4 4],[mdat mdat-sddat],'k','LineWidth',2) % whisker 2
plot([3.9 4.1],[mdat-sddat mdat-sddat],'k','LineWidth',2)

% mean with CI
plot(5,mdat,'ko','MarkerSize',10,'MarkerFaceColor','k')
plot([5 5],[mdat CI(2)],'k','LineWidth',2) % whisker 1
plot([4.9 5.1],[CI(2) CI(2)],'k','LineWidth',2)
plot([5 5],[mdat CI(1)],'k','LineWidth',2) % whisker 2
plot([4.9 5.1],[CI(1) CI(1)],'k','LineWidth',2)

% mean with bootCI
plot(6,mdat,'ko','MarkerSize',10,'MarkerFaceColor','k')
plot([6 6],[mdat bootCI(2)],'k','LineWidth',2) % whisker 1
plot([5.9 6.1],[bootCI(2) bootCI(2)],'k','LineWidth',2)
plot([6 6],[mdat bootCI(1)],'k','LineWidth',2) % whisker 2
plot([5.9 6.1],[bootCI(1) bootCI(1)],'k','LineWidth',2)

% mean +/- SEM
plot(7,mdat,'ko','MarkerSize',10,'MarkerFaceColor','k')
plot([7 7],[mdat mdat+semdat],'k','LineWidth',2) % whisker 1
plot([6.9 7.1],[mdat+semdat mdat+semdat],'k','LineWidth',2)
plot([7 7],[mdat mdat-semdat],'k','LineWidth',2) % whisker 2
plot([6.9 7.1],[mdat-semdat mdat-semdat],'k','LineWidth',2)

set(gca,'XLim',[0 8],'YLim',[0 25],'FontSize',14)
box on

labels = {'Scatter';'Box & whiskers';'Quartiles';'Mean +/- SD';...
    'Mean with CI';'Mean with bootCI';'Mean +/- SEM'};
set(gca,'XTick',1:7,'XTickLabel',labels,'XTickLabelRotation',45)

%% Example with outliers: quartiles

% Now, let's look at the behaviour of different estimators in the presence
% of outliers. To emphasis the effect of one outlier, we're going to
% consider only a subsample of the data.
% We start by considering the quartiles, which, as expected, are robust to
% outliers...

data = dat(1:2:end); % take every other observation
data = sort(data);
n = numel(data);

% rng(21); jitter = randn(n,1);

figure('Color','w','NumberTitle','off')
hold on

for out = 1:7

    data(n) = data(n)+out-1; % add outliers
    n = numel(data);
    qs = prctile(data,[25 50 75]);

    % scatter
    xpts = UnivarScatter_nofig(data);
    scatter(xpts+out-1,data,70,[.7 .7 .7],'filled')
%     scatter(out-1+ones(n,1)+jitter/10,data,70,[.7 .7 .7],'filled')

    % quartiles
    plot(out,qs(2),'ko','MarkerSize',10,'MarkerFaceColor','k')
    plot([out out],[qs(2) qs(3)],'k','LineWidth',2) % whisker 1
    plot([out out],[qs(1) qs(2)],'k','LineWidth',2) % whisker 2

    % reference lines
    if out == 1
    plot([0 8],[qs(1) qs(1)],'k--')
    plot([0 8],[qs(2) qs(2)],'k:')
    plot([0 8],[qs(3) qs(3)],'k--')
    end

    set(gca,'XLim',[0 8],'YLim',[0 30],'FontSize',14)
    box on

end

%% Example with outliers: classic confidence interval

% Contrary to the quartiles, the classic CI of the mean is not robust, so it
% provides very inaccurate results.

data = dat(1:2:end);
data = sort(data);
n = numel(data);

figure('Color','w','NumberTitle','off')
hold on

for out = 1:7

    data(n) = data(n)+out-1; % add outliers
    [H,P,CI,STATS] = ttest(data,alpha);
    mdata = mean(data);

    % scatter
    xpts = UnivarScatter_nofig(data);
    scatter(xpts+out-1,data,70,[.7 .7 .7],'filled')
%     scatter(out-1+ones(n,1)+jitter/10,data,70,[.7 .7 .7],'filled')

    % median & quartiles
    plot(out,mdata,'ko','MarkerSize',10,'MarkerFaceColor','k')
    plot([out out],[mdata CI(2)],'k','LineWidth',2) % whisker 1
    plot([out out],[CI(1) mdata],'k','LineWidth',2) % whisker 2

     % reference lines
    if out == 1
    plot([0 8],[CI(1) CI(1)],'k--')
    plot([0 8],[mdata mdata],'k:')
    plot([0 8],[CI(2) CI(2)],'k--')
    end

    set(gca,'XLim',[0 8],'YLim',[0 30],'FontSize',14)
    box on

end

%% Example with outliers: percentile bootstrap confidence interval of the mean

% The percentile bootstrap CI of the mean also gives misleading results, but
% it performs better than the classic CI: notice how, contrary to the
% classic CI, the lower bound is not affected by outliers.

rng(4); % set seed to reproduce figure from blog

data = dat(1:2:end);
data = sort(data);
n = numel(data);

Nb = 1000;
est = 'mean';

% rng(21); jitter = randn(n,1);

figure('Color','w','NumberTitle','off')
hold on

for out = 1:7

    data(n) = data(n)+out-1; % add outliers
    if out == 1
        s = rng;
    else
        rng(s) % use same seed for all bootstrap CI
    end
    [BEST,EST,CI,p] = pbci(data,Nb,alpha,est);
    mdata = EST;

    % scatter
    xpts = UnivarScatter_nofig(data);
    scatter(xpts+out-1,data,70,[.7 .7 .7],'filled')
%     scatter(out-1+ones(n,1)+jitter/10,data,70,[.7 .7 .7],'filled')

    % median & quartiles
    plot(out,mdata,'ko','MarkerSize',10,'MarkerFaceColor','k')
    plot([out out],[mdata CI(2)],'k','LineWidth',2) % whisker 1
    plot([out out],[CI(1) mdata],'k','LineWidth',2) % whisker 2

     % reference lines
    if out == 1
    plot([0 8],[CI(1) CI(1)],'k--')
    plot([0 8],[mdata mdata],'k:')
    plot([0 8],[CI(2) CI(2)],'k--')
    end

    set(gca,'XLim',[0 8],'YLim',[0 30],'FontSize',14)
    box on

end

%% Example with outliers: percentile bootstrap confidence interval of the median

% The percentile bootstrap CI of the median is robust to outliers.
% But it can be unstable, as seen in a previous cell.

rng(4); % set seed to reproduce figure from blog

data = dat(1:2:end);
data = sort(data);
n = numel(data);

Nb = 1000;
est = 'median';

% rng(21); jitter = randn(n,1);

figure('Color','w','NumberTitle','off')
hold on

for out = 1:7

    data(n) = data(n)+out-1; % add outliers
    if out == 1
        s = rng;
    else
        rng(s) % use same seed for all bootstrap CI
    end
    [BEST,EST,CI,p] = pbci(data,Nb,alpha,est);
    mdata = EST;

    % scatter
    xpts = UnivarScatter_nofig(data);
    scatter(xpts+out-1,data,70,[.7 .7 .7],'filled')
%     scatter(out-1+ones(n,1)+jitter/10,data,70,[.7 .7 .7],'filled')

    % median & quartiles
    plot(out,mdata,'ko','MarkerSize',10,'MarkerFaceColor','k')
    plot([out out],[mdata CI(2)],'k','LineWidth',2) % whisker 1
    plot([out out],[CI(1) mdata],'k','LineWidth',2) % whisker 2

     % reference lines
    if out == 1
    plot([0 8],[CI(1) CI(1)],'k--')
    plot([0 8],[mdata mdata],'k:')
    plot([0 8],[CI(2) CI(2)],'k--')
    end

    set(gca,'XLim',[0 8],'YLim',[0 30],'FontSize',14)
    box on

end
