%% dependencies

% Spearman's correlation is computed using the Robust Correlation Toolbox
% https://sourceforge.net/projects/robustcorrtool/
% That toolbox also contains functions to compute the winsorized
% correlation, the skipped correlation and a few other goodies.

%% CASE 1 ================================================================
%% create & illustrate data

rng(2)

Np = 100; % sample size / number of participants
a = randn(1,Np);
b = a * .1 + rand(1,Np);
c = a * .4 + rand(1,Np);

figure('Color','w')
subplot(1,2,1)
scatter(a,b,'MarkerEdgeColor','k','MarkerFaceColor',[.7 .7 .7],'MarkerEdgeAlpha',.5,'MarkerFaceAlpha',.5)
xlabel('a','FontSize',20,'FontWeight','bold')
ylabel('b','FontSize',20,'FontWeight','bold')

subplot(1,2,2)
scatter(a,c,'MarkerEdgeColor','k','MarkerFaceColor',[.7 .7 .7],'MarkerEdgeAlpha',.5,'MarkerFaceAlpha',.5)
xlabel('a','FontSize',20,'FontWeight','bold')
ylabel('c','FontSize',20,'FontWeight','bold')

for sub = 1:2
    subplot(1,2,sub)
    axis square
    set(gca,'FontSize',14,'XLim',[-5 5],'YLim',[-1 1.5])
    lsline
end

corr1 = Spearman(a,b);
corr2 = Spearman(a,c);

%% bootstrap

Nb = 500;
bootcorr1 = zeros(Nb,1);
bootcorr2 = zeros(Nb,1);

for B = 1:Nb
    
    bootsample = randi(Np,1,Np);
    bootcorr1(B) = Spearman(a(bootsample),b(bootsample),0);
    bootcorr2(B) = Spearman(a(bootsample),c(bootsample),0);
    
end

%% bootstrap without loop
% Cyril Pernet rightly pointed out on Twitter that the loop is unnecessary:
bootsamples = randi(Np,Np,Nb);
bc1 = Spearman(a(bootsamples),b(bootsamples),0);
bc2 = Spearman(a(bootsamples),c(bootsamples),0);

% The bootstrap loop does make the bootstrap procedure more intuitive
% though, especially for new users. 

%% confidence intervals

alpha = 0.05; % probability coverage - 0.05 for 95% CI
hi = floor((1-alpha/2)*Nb+.5);
lo = floor((alpha/2)*Nb+.5);

% for each correlation
boot1sort = sort(bootcorr1);
boot2sort = sort(bootcorr2);
boot1ci = [boot1sort(lo) boot1sort(hi)]; 
boot2ci = [boot2sort(lo) boot2sort(hi)]; 

% for the difference between correlations
bootdiff = bootcorr1 - bootcorr2;
bootdiffsort = sort(bootdiff);
diffci = [bootdiffsort(lo) bootdiffsort(hi)]; 
pvalue = mean(bootdiffsort<0);
pvalue = 2*min(pvalue,1-pvalue);

fprintf('==========================\n')
fprintf('corr(a,b) = %.2f [%.2f %.2f]\n',corr1,boot1ci(1),boot1ci(2))
fprintf('corr(a,c) = %.2f [%.2f %.2f]\n',corr2,boot2ci(1),boot2ci(2))
fprintf('difference = %.2f [%.2f %.2f]\n',corr1-corr2,diffci(1),diffci(2))

%% illustrate bootstrap distribution + confidence interval

cdiff = corr1-corr2;

figure('Color','w'); hold on
hist(bootdiff,50)
h = findobj(gca,'Type','patch');
h.FaceColor = [0 .5 .5];
h.EdgeColor = 'w';
v = axis;
plot([cdiff cdiff],[v(3) v(4)],'k','LineWidth',5)
plot([diffci(1) diffci(1)],[v(3) v(4)],'k')
plot([diffci(2) diffci(2)],[v(3) v(4)],'k')
set(gca,'FontSize',14,'Layer','Top')

%% CASE 2 ================================================================
%% create & illustrate data

rng(2)

Np = 100; % sample size / number of participants
a1 = randn(1,Np);
b1 = a * .1 + rand(1,Np);
a2 = a1 + rand(1,Np);
b2 = a2 * .2 + rand(1,Np);

figure('Color','w')
subplot(1,2,1)
scatter(a1,b1,'MarkerEdgeColor','k','MarkerFaceColor',[.7 .7 .7],'MarkerEdgeAlpha',.5,'MarkerFaceAlpha',.5)
xlabel('a1','FontSize',20,'FontWeight','bold')
ylabel('b1','FontSize',20,'FontWeight','bold')

subplot(1,2,2)
scatter(a2,b2,'MarkerEdgeColor','k','MarkerFaceColor',[.7 .7 .7],'MarkerEdgeAlpha',.5,'MarkerFaceAlpha',.5)
xlabel('a2','FontSize',20,'FontWeight','bold')
ylabel('b2','FontSize',20,'FontWeight','bold')

for sub = 1:2
    subplot(1,2,sub)
    axis square
    set(gca,'FontSize',14,'XLim',[-5 5],'YLim',[-1 1.5])
    lsline
end

corr1 = Spearman(a1,b1);
corr2 = Spearman(a2,b2);

%% bootstrap

Nb = 500;
bootcorr1 = zeros(Nb,1);
bootcorr2 = zeros(Nb,1);

for B = 1:Nb
    
    bootsample = randi(Np,1,Np);
    bootcorr1(B) = Spearman(a1(bootsample),b1(bootsample),0);
    bootcorr2(B) = Spearman(a2(bootsample),b2(bootsample),0);
    
end

%% confidence intervals

alpha = 0.05; % probability coverage - 0.05 for 95% CI
hi = floor((1-alpha/2)*Nb+.5);
lo = floor((alpha/2)*Nb+.5);

% for each correlation
boot1sort = sort(bootcorr1);
boot2sort = sort(bootcorr2);
boot1ci = [boot1sort(lo) boot1sort(hi)]; 
boot2ci = [boot2sort(lo) boot2sort(hi)]; 

% for the difference between correlations
bootdiff = bootcorr1 - bootcorr2;
bootdiffsort = sort(bootdiff);
diffci = [bootdiffsort(lo) bootdiffsort(hi)]; 
pvalue = mean(bootdiffsort<0);
pvalue = 2*min(pvalue,1-pvalue);

fprintf('==========================\n')
fprintf('corr(a1,b1) = %.2f [%.2f %.2f]\n',corr1,boot1ci(1),boot1ci(2))
fprintf('corr(a2,b2) = %.2f [%.2f %.2f]\n',corr2,boot2ci(1),boot2ci(2))
fprintf('difference = %.2f [%.2f %.2f]\n',corr1-corr2,diffci(1),diffci(2))

%% illustrate bootstrap distribution + confidence interval

cdiff = corr1-corr2;

figure('Color','w'); hold on
hist(bootdiff,50)
h = findobj(gca,'Type','patch');
h.FaceColor = [0 .5 .5];
h.EdgeColor = 'w';
v = axis;
plot([cdiff cdiff],[v(3) v(4)],'k','LineWidth',5)
plot([diffci(1) diffci(1)],[v(3) v(4)],'k')
plot([diffci(2) diffci(2)],[v(3) v(4)],'k')
plot([0 0],[v(3) v(4)],'k--','LineWidth',3)
set(gca,'FontSize',14,'Layer','Top')

%% highest-density interval
% A HDI can be calculated suing the hdi() function
% the hdi() function is available on github:
% <https://github.com/GRousselet/matlab_stats>

credmass = 0.8; % mass within interval = scalar between 0 and 1
bootdiff_hdi = hdi(bootdiff, credmass);
