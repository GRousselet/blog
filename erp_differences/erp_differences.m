%% erp_differences
% This script demonstrates how to explore the time-course of ERP effects
% between two repeated-measure conditions.
% Copyright (C) 2016 Guillaume Rousselet - University of Glasgow

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.

%% BEFORE YOU START

% Each cell (the beige area) can be evaluated by making a right click >
% Evaluate Current Cell, or by using these shortcuts:

% mac shortcut = command + enter
% pc shortcut = ctrl + enter

% You can navigate from cell to cell by clicking command+up arrow, or down arrow.

%% ------------------------------------------------------------------------
%% ADD LIMO_EEG TO YOUR PATH

% You will need LIMO EEG to run parts of this script.

% The stable version is here:
% https://gforge.dcn.ed.ac.uk/gf/project/limo_eeg/frs/

% Go here for the version in development:
% https://github.com/LIMO-EEG-Toolbox/limo_eeg

% add the folders to your path:
addpath(genpath('/yourpath/limo_eeg/'));

% Alternatively go to >Home >Set Path >Add with subfolders and choose the
% limo_eeg folder

% Similarly, add the erp_differences folder to your path.

%% ------------------------------------------------------------------------
%% LOAD DATA

% To load the data, you can:
%
% drag and drop the file erpm.mat into the command window;
%
% use the Open menu;
%
% or directly through the command window:
%
% >> load erpm % if the file is in the current directory
% >> cd /yourpath/ % to change the current directory
% >> load /yourpath/erpm % if you need to specify the path

%% COMPUTE AVERAGES

C1=mean(erpm(:,:,1),2); % average across participants for condition 1
C2=mean(erpm(:,:,2),2); % average across participants for condition 2

Xf = -300:2:600; % time vector
Nf = numel(Xf);

%% MAKE STANDARD FIGURE

figure('Color','w','NumberTitle','off') % make blank figure
hold on % we're going to plot several elements

plot(Xf,C1,'Color',[0 0 0],'LineStyle','-','LineWidth',3)
plot(Xf,C2,'Color',[.7 .7 .7],'LineStyle','-','LineWidth',3)

set(gca,'LineWidth',1,'FontSize',14,'XLim',[-300 600])
xlabel('Time in ms','FontSize',16)
ylabel('Amplitude in \muV','FontSize',16)
box on

% you could change many plot properties -- see guidelines here:
% http://www.mathworks.co.uk/help/matlab/creating_plots/using-high-level-plotting-functions.html#f6-35125

% To change text properties, see:
% http://www.mathworks.co.uk/help/matlab/ref/text_props.html

%% SAVE FIGURE

% in the figure panel > File > Save as

% in a script, we can generate figures automatically:

% JPEG FORMAT
fileName = 'standard_figure.jpg'; 
set(gcf,'PaperPositionMode','auto')
print(gcf,'-djpeg','-r150', fileName); 

% TIFF FORMAT
% fileName = 'standard_figure.tif'; 
% set(gcf,'PaperPositionMode','auto')
% print(gcf,'-dtiff','-r150', fileName); 

% EPS FORMAT (for publications)
% fileName = 'standard_figure.eps';
% set(gcf,'PaperPositionMode','auto')
% print(gcf, '-dpsc','-r150', fileName);


%% STATS

% We do a t-test at each time point & plot the results on the figure.
% This is done in one step using the function limo_yuend_ttest
% from the LIMO_EEG toolbox. It performs a t-test on trimmed means for
% dependent groups.
% >> help limo_yuend_ttest % to learn more about this function

c1=zeros(1,451,20);c2=c1; % reformat data so that dimension 1 = 1 electrode
% if we had say 100 electrodes, c1 would have dimensions 100,451,20

c1(1,:,:) = erpm(:,:,1); c2(1,:,:) = erpm(:,:,2);
percent = 0; % with 0% trimming, this is a regular t-test on means
alpha = 0.05;
[Ty,diff,se,CI,p,tcrit,df] = limo_yuend_ttest(c1,c2,percent,alpha);

% We can also do a ttest on means using the function limo_ttest:
% [m,dfe,ci,sd,n,t,p] = limo_ttest(1,c1,c2,alpha);

% We could also use limo_rep_anova with one 2-level factor:
% data=zeros(size(c1,2),size(c1,3),2);
% data(:,:,1)=c1;data(:,:,2)=c2;
% result=limo_rep_anova(data,[],[2]);
% You can check that `result.F` is the same as either `t` or `Ty` squared.

%% ADD STATS TO STANDARD FIGURE

figure('Color','w','NumberTitle','off') % make blank figure
hold on % we're going to plot several elements
plot(Xf,C1,'Color',[0 0 0],'LineStyle','-','LineWidth',3)
plot(Xf,C2,'Color',[.7 .7 .7],'LineStyle','-','LineWidth',3)
set(gca,'LineWidth',1,'FontSize',14,'XLim',[-300 600])
xlabel('Time in ms','FontSize',16)
ylabel('Amplitude in \muV','FontSize',16)
box on

% We add a horizontal line showing time points at which a significant mean
% difference was observed, according to a paired t-test.
v=axis; % get current axis limits
alpha=0.05; % set alpha level
plot(Xf,(p<=alpha)*100-100+v(3)*.95,'r.','MarkerSize',15)
axis(v) % restore axis limits to hide non-significant points

%% SAVE FIGURE

fileName = 'standard_figure_with_stats.jpg'; 
set(gcf,'PaperPositionMode','auto')
print(gcf,'-djpeg','-r150', fileName); 

%% ADD CONFIDENCE INTERVALS

figure('Color','w','NumberTitle','off') % make blank figure
hold on % we're going to plot several elements
plot(Xf,C1,'Color',[0 0 0],'LineStyle','-','LineWidth',3)
plot(Xf,C2,'Color',[.7 .7 .7],'LineStyle','-','LineWidth',3)
set(gca,'LineWidth',1,'FontSize',14,'XLim',[-300 600])
xlabel('Time in ms','FontSize',16)
ylabel('Amplitude in \muV','FontSize',16)
box on

% get confidence intervals for each condition
percent=0;
alpha=0.05; % to get a 95% confidence interval
nullvalue=0;
[t,tmdata,trimci_c1,p,tcrit,df] = limo_trimci(c1, percent, alpha, nullvalue);
[t,tmdata,trimci_c2,p,tcrit,df] = limo_trimci(c2, percent, alpha, nullvalue);

% plot confidence intervals
plot(Xf,squeeze(trimci_c1),'Color',[0 0 0],'LineStyle','-','LineWidth',1)
plot(Xf,squeeze(trimci_c2),'Color',[.7 .7 .7],'LineStyle','-','LineWidth',1)

% plot results of t-test
[Ty,diff,se,CI,p,tcrit,df] = limo_yuend_ttest(c1,c2,percent,alpha);
v=axis; % get current axis limits
alpha=0.05; % set alpha level 
plot(Xf,(p<=alpha)*100-100+v(3)*.95,'r.','MarkerSize',15)
axis(v) % restore axis limits to hide non-significant points

set(gca,'Layer','Top') % so the box surrounding the plot stays on top

%% OTHER CONFIDENCE INTERVAL OPTIONS

% In addition to `limo_trimci`, you can compute confidence intervals using:
% `limo_pbci` - percentile bootstrap
% `limo_bootttest1` - percentile-t bootstrap

% Other functions return confidence intervals for the difference, e.g.:
% `limo_yuend_ttest` & `limo_yuen_ttest`

% On github you will also find:
% `limo_central_estimator` - returns the highest density interval of
%                            Bayesian boostrap samples

%% SAVE FIGURE

fileName = 'standard_figure_with_ci.jpg'; 
set(gcf,'PaperPositionMode','auto')
print(gcf,'-djpeg','-r150', fileName); 

%% USE SHADED AREAS INSTEAD

figure('Color','w','NumberTitle','off') % make blank figure
hold on % we're going to plot several elements
plot(Xf,C1,'Color',[0 0 0],'LineStyle','-','LineWidth',3)
plot(Xf,C2,'Color',[.7 .7 .7],'LineStyle','-','LineWidth',3)
set(gca,'LineWidth',1,'FontSize',14,'XLim',[-300 600])
xlabel('Time in ms','FontSize',16)
ylabel('Amplitude in \muV','FontSize',16)
box on

% plot confidence intervals
x = Xf; % time vector
y = squeeze(trimci_c1(:,:,1)); % CI lower bound
z = squeeze(trimci_c1(:,:,2)); % CI upper bound
c = [1 0.5 0.2]; % set colour
t = .1; % set transparency [0, 1]
x = x(:)';y=y(:)';z=z(:)';
X = [x,fliplr(x)]; 
Y = [y,fliplr(z)]; 
hf = fill(X,Y,c); % plot filled area
set(hf,'FaceAlpha',t,'EdgeColor',[0 0 0]);

y = squeeze(trimci_c2(:,:,1));
z = squeeze(trimci_c2(:,:,2));
c = [0 .7 .1];
x=x(:)';y=y(:)';z=z(:)';
X = [x,fliplr(x)];
Y = [y,fliplr(z)];
hf = fill(X,Y,c); % plot filled area
set(hf,'FaceAlpha',t,'EdgeColor',[.7 .7 .7]);

set(gca,'Layer','Top') % so the box surrounding the plot stays on top

%% SAVE FIGURE

fileName = 'standard_figure_with_ci2.jpg'; 
set(gcf,'PaperPositionMode','auto')
print(gcf,'-djpeg','-r150', fileName); 

%% ADD DIFFERENCE AND ITS CONFIDENCE INTERVAL

scrsz = get(groot,'Screensize'); % get the size of the screen
figure('Color','w','NumberTitle','off','Position',[1 scrsz(4) scrsz(3)/2.5 scrsz(4)])

subplot(2,1,1);hold on  % plot conditions 1 and 2 **********************

plot(Xf,C1,'Color',[0 0 0],'LineStyle','-','LineWidth',3)
plot(Xf,C2,'Color',[.7 .7 .7],'LineStyle','-','LineWidth',3)

% get confidence intervals for each condition
percent=0;
alpha=0.05; % to get a 95% confidence interval
nullvalue=0;
[t,tmdata,trimci_c1,p,tcrit,df]=limo_trimci(c1, percent, alpha, nullvalue);
[t,tmdata,trimci_c2,p,tcrit,df]=limo_trimci(c2, percent, alpha, nullvalue);

% plot confidence intervals
plot(Xf,squeeze(trimci_c1),'Color',[0 0 0],'LineStyle','-','LineWidth',1)
plot(Xf,squeeze(trimci_c2),'Color',[.7 .7 .7],'LineStyle','-','LineWidth',1)

% add lines marking the 2 difference peaks
F1 = Xf(diff(Xf<200)==max(diff(Xf<200)));
F2 = Xf(diff==max(diff));
v = axis;
plot([F1 F1],[v(3) v(4)],'r')
plot([F2 F2],[v(3) v(4)],'r')

subplot(2,1,2);hold on % plot difference ********************************
% we use the outputs diff and CI from `limo_yuend_ttest`
plot(Xf,diff,'Color',[1 0 0],'LineStyle','-','LineWidth',3)
plot(Xf,squeeze(CI),'Color',[1 0 0],'LineStyle','-','LineWidth',1)

v=axis; % get current axis limits

% add lines marking the 2 difference peaks
plot([F1 F1],[v(3) v(4)],'r')
plot([F2 F2],[v(3) v(4)],'r')

[Ty,diff,se,CI,p,tcrit,df]=limo_yuend_ttest(c1,c2,percent,alpha);
alpha=0.05; % set alpha level
plot(Xf,(p<=alpha)*100-100+v(3)*.95,'r.','MarkerSize',10)
axis(v) % restore axis limits to hide non-significant points

plot([v(1) v(2)],[0 0],'k--') % add horizontal line at zero

for sub=1:2 % format 2 subplots
    subplot(2,1,sub)
    set(gca,'LineWidth',1,'FontSize',14,'XLim',[-300 600])
    xlabel('Time in ms','FontSize',16)
    ylabel('Amplitude in \muV','FontSize',16)
    box on
end

%% SAVE FIGURE

fileName = 'figure_with_difference.jpg'; 
set(gcf,'PaperPositionMode','auto')
print(gcf,'-djpeg','-r150', fileName); 

%% HOW MANY PARTICIPANTS SHOW AN EFFECT?

scrsz = get(groot,'ScreenSize'); % get the size of the screen
figure('Color','w','NumberTitle','off','Position',[1 scrsz(4) scrsz(3)/2.5 scrsz(4)]) 

% --------------------------------------------------------------------------
subplot(5,1,[1 2]);hold on  % plot conditions 1 and 2 **********************

plot(Xf,C1,'Color',[0 0 0],'LineStyle','-','LineWidth',3)
plot(Xf,C2,'Color',[.7 .7 .7],'LineStyle','-','LineWidth',3)

% get confidence intervals for each condition
percent=0;
alpha=0.05; % to get a 95% confidence interval
nullvalue=0;
[t,tmdata,trimci_c1,p,tcrit,df]=limo_trimci(c1, percent, alpha, nullvalue);
[t,tmdata,trimci_c2,p,tcrit,df]=limo_trimci(c2, percent, alpha, nullvalue);

% plot confidence intervals
plot(Xf,squeeze(trimci_c1),'Color',[0 0 0],'LineStyle','-','LineWidth',1)
plot(Xf,squeeze(trimci_c2),'Color',[.7 .7 .7],'LineStyle','-','LineWidth',1)

[Ty,diff,se,CI,p,tcrit,df]=limo_yuend_ttest(c1,c2,percent,alpha);

set(gca,'LineWidth',1,'FontSize',14,'XLim',[-300 600])
%     xlabel('Time in ms','FontSize',16)
    ylabel('Amplitude in \muV','FontSize',16)
    box on

% --------------------------------------------------------------------------
subplot(5,1,[3 4]);hold on % plot difference ********************************

plot(Xf,squeeze(c1-c2),'Color',[.7 .7 .7]) % PLOT DIFFERENCE FOR EVERY PARTICIPANT <<<<

plot(Xf,diff,'Color',[1 0 0],'LineStyle','-','LineWidth',3)
plot(Xf,squeeze(CI),'Color',[1 0 0],'LineStyle','-','LineWidth',1)

v=axis; % get current axis limits
alpha=0.05; % set alpha level 
plot(Xf,(p<=alpha)*100-100+v(3)*.95,'r.','MarkerSize',10)
axis(v) % restore axis limits to hide non-significant points

plot([v(1) v(2)],[0 0],'k--') % add horizontal line at zero

set(gca,'LineWidth',1,'FontSize',14,'XLim',[-300 600])
%     xlabel('Time in ms','FontSize',16)
    ylabel('Amplitude in \muV','FontSize',16)
    box on

% --------------------------------------------------------------------------
subplot(5,1,5);hold on % plot proportion of participants with difference > 0

prop = sum(squeeze(c1-c2)>0,2)./size(c1,3);

plot(Xf,prop,'Color',[1 0 0],'LineStyle','-','LineWidth',3)

set(gca,'LineWidth',1,'FontSize',14,'XLim',[-300 600])
    xlabel('Time in ms','FontSize',16)
    ylabel({'Proportion';'of participants'},'FontSize',16)
    box on
plot([v(1) v(2)],[.5 .5],'k--') % add horizontal line at 0.5

%% SAVE FIGURE

fileName = 'full_figure.jpg'; 
set(gcf,'PaperPositionMode','auto')
print(gcf,'-djpeg','-r150', fileName); 

%% ---------------------------------------
%% EXPLORE RESULTS FROM EVERY PARTICIPANT:

% first we need to load erps.mat
load erps
Np = length(erps);

%% T-TESTS ON MEANS

alpha=0.05;
percent=0;

% we declare variables that will hold the results
diff_all = zeros(Np,Nf);
CI_all = zeros(Np,Nf,2);
p_all = zeros(Np,Nf);

for P = 1:Np % for each participant
    
   [Ty,diff_all(P,:),se,CI_all(P,:,:),p_all(P,:),tcrit,df] = limo_yuen_ttest(erps{P,1},erps{P,2},percent,alpha);
   % We use `limo_yuen_ttest`, which performs t-tests on independent groups.
   % For each participant, distributions of single-trials are independent.
   % Previously, we used `limo_yuend_ttest` for group level analyses
   % between 2 dependent groups.
    
end

%% SINGLE-PARTICIPANT FIGURES

scrsz = get(groot,'ScreenSize'); % get the size of the screen
figure('Color','w','NumberTitle','off','Position',[1 1 scrsz(3) scrsz(4)]) 

for P = 1:Np 
    
    subplot(4,5,P);hold on
    title(['P',num2str(P)],'FontSize',20)
    plot(Xf,diff_all(P,:),'Color','r','LineStyle','-','LineWidth',3)
    plot(Xf,squeeze(CI_all(P,:,:)),'Color','r','LineStyle','-','LineWidth',1)
   
    axis tight
    v=axis;
    plot([v(1) v(2)],[0 0],'k--')
    plot([0 0],[v(3) v(4)],'k--')
    box on
    set(gca,'LineWidth',1,'FontSize',14,'XLim',[-300 600],'YLim',[v(3)-1 v(4)+1])
    
end

%% SAVE FIGURE

fileName = 'every_participant.jpg'; 
set(gcf,'PaperPositionMode','auto')
print(gcf,'-djpeg','-r150', fileName); 
