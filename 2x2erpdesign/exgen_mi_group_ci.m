%% load data
load gpmi
% gpmi contains variables times, eo, ey, go, gy
% letters stand for:
%  - o = old
%  - y = young
%  - e = expressive task
%  - g = gender task

% exclude last 2 columns
eo = eo(:,1:18);
go = go(:,1:18);
ey = ey(:,1:24);
gy = gy(:,1:24);

Ntimes = numel(times);

%% compute group medians using the Harrell-Davis estimate of the 50th quantile
% the hd() function is available on github
% <https://github.com/GRousselet/matlab_stats>

mdey = zeros(1,Ntimes);
mdeo = zeros(1,Ntimes);
mdgy = zeros(1,Ntimes);
mdgo = zeros(1,Ntimes);

% task differences
mdty = zeros(1,Ntimes);
mdto = zeros(1,Ntimes); 

for T = 1:Ntimes
    mdey(T) = hd(ey(T,:));
    mdeo(T) = hd(eo(T,:));
    mdgy(T) = hd(gy(T,:));
    mdgo(T) = hd(go(T,:));
    
    % task differences: expressive - gender
    mdty(T) = hd(ey(T,:) - gy(T,:));
    mdto(T) = hd(eo(T,:) - go(T,:));
end

% group differences for each task
diff_md_g = mdgy - mdgo;
diff_md_e = mdey - mdeo;

% interaction
md_interaction = mdty - mdto;

%% compute group means

mey = mean(ey,2);
meo = mean(eo,2);
mgy = mean(gy,2);
mgo = mean(go,2);

% plot individual MI time-courses for each group and condition
% superimpose means and medians

figure('Color','w','NumberTitle','off','Position', [1000 1000 800 600])

subplot(2,2,1); hold on
plot(times,ey)
plot(times,mey,'Color',[.5 .5 .5],'LineWidth',3)
plot(times,mdey,'Color',[.2 .2 .2],'LineWidth',3)
title('Young participants','FontSize',20)
ylabel('Expressive task','FontSize',16)

subplot(2,2,3); hold on
plot(times,gy)
plot(times,mgy,'Color',[.5 .5 .5],'LineWidth',3)
plot(times,mdgy,'Color',[.2 .2 .2],'LineWidth',3)
xlabel('Time in ms','FontSize',16)
ylabel('Gender task','FontSize',16)

subplot(2,2,2); hold on
plot(times,eo)
plot(times,meo,'Color',[.5 .5 .5],'LineWidth',3)
plot(times,mdeo,'Color',[.2 .2 .2],'LineWidth',3)
title('Older participants','FontSize',20)

subplot(2,2,4); hold on
plot(times,go)
plot(times,mgo,'Color',[.5 .5 .5],'LineWidth',3)
plot(times,mdgo,'Color',[.2 .2 .2],'LineWidth',3)
xlabel('Time in ms','FontSize',16)

for sub = 1:4
    subplot(2,2,sub)
    set(gca,'XLim',[-300 600],'YLim',[0 0.12],'XTick',-300:100:600)
    set(gca,'FontSize',14)
end

%% bootstrap

Nb = 1000; % number of bootstrap samples
Ny = size(ey,2);
No = size(eo,2);

% declare bootstrap variables:
% original conditions
boot_mdey = zeros(Ntimes,Nb);
boot_mdeo = zeros(Ntimes,Nb);
boot_mdgy = zeros(Ntimes,Nb);
boot_mdgo = zeros(Ntimes,Nb);
% pairwise task differences (expressive - gender)
boot_mdty = zeros(Ntimes,Nb);
boot_mdto = zeros(Ntimes,Nb);

for B = 1:Nb
    
    disp(['bootstrap = ',num2str(B),'/',num2str(Nb)])
    
   % sample participants with replacements, independently in each group
   booty = randi(Ny,1,Ny);
   booto = randi(No,1,No);
   
   % get bootstrap samples
   boot_ey = ey(:,booty);
   boot_gy = gy(:,booty);
   boot_eo = eo(:,booto);
   boot_go = go(:,booto);
   
   % compute medians of bootstrap samples --------------
   for T = 1:Ntimes
       
       % original conditions
       boot_mdey(T,B) = hd(boot_ey(T,:));
       boot_mdeo(T,B) = hd(boot_eo(T,:));
       boot_mdgy(T,B) = hd(boot_gy(T,:));
       boot_mdgo(T,B) = hd(boot_go(T,:));
       
       % medians of pairwise task differences
       boot_mdty(T,B) = hd(boot_ey(T,:) - boot_gy(T,:));
       boot_mdto(T,B) = hd(boot_eo(T,:) - boot_go(T,:));
   end
   
end

save bootres boot_mdey boot_mdeo boot_mdgy boot_mdgo boot_mdty boot_mdto

%% get confidence intervals

alpha = 0.05;
lo = round(Nb*alpha/2);
hi = Nb - lo;
lo = lo + 1;

% for each task and group
sort_boot_mdey = sort(boot_mdey,2);
sort_boot_mdeo = sort(boot_mdeo,2);
sort_boot_mdgy = sort(boot_mdgy,2);
sort_boot_mdgo = sort(boot_mdgo,2);

mdey_ci = sort_boot_mdey(:,[lo hi]);
mdeo_ci = sort_boot_mdeo(:,[lo hi]);
mdgy_ci = sort_boot_mdgy(:,[lo hi]);
mdgo_ci = sort_boot_mdgo(:,[lo hi]);
       
% for group differences for each task
sort_boot_diff_md_g = sort(boot_mdgy - boot_mdgo,2);
sort_boot_diff_md_e = sort(boot_mdey - boot_mdeo,2);

diff_md_g_ci = sort_boot_diff_md_g(:,[lo hi]);
diff_md_e_ci = sort_boot_diff_md_e(:,[lo hi]);

% for task differences for each group
sort_boot_mdty = sort(boot_mdty,2);
sort_boot_mdto = sort(boot_mdto,2);

mdty_ci = sort_boot_mdty(:,[lo hi]);
mdto_ci = sort_boot_mdto(:,[lo hi]);

% for the interaction between groups and tasks (difference of difference)
sort_boot_md_interaction = sort(boot_mdty - boot_mdto,2);
md_interaction_ci = sort_boot_md_interaction(:,[lo hi]);

save gpci mdey_ci mdeo_ci mdgy_ci mdgo_ci diff_md_g_ci diff_md_e_ci mdty_ci mdto_ci md_interaction_ci

% If you want to report highest density intervals instead of confidence
% intervals, see the script exgen_mi_group_hdi.m, which uses a Matlab
% adaptation of R code from John Kruschke in the function hdi().

%% ===========================================================
%% plot results
%% ===========================================================
%% plot individual MI time-courses for each group and condition
% superimpose medians and confidence intervals

% load gpci

figure('Color','w','NumberTitle','off','Position', [1000 1000 800 600])

subplot(2,2,1); hold on
for P = 1:Ny
h1 = plot(times,ey(:,P)); h1.Color(4) = 0.5;
end
plot(times,mdey_ci,'Color',[.2 .2 .2],'LineWidth',1)
plot(times,mdey,'Color',[.2 .2 .2],'LineWidth',3)
title('Young participants','FontSize',20)
ylabel('Expressive task','FontSize',16)

subplot(2,2,3); hold on
for P = 1:Ny
h1 = plot(times,gy(:,P)); h1.Color(4) = 0.5;
end
plot(times,mdgy_ci,'Color',[.2 .2 .2],'LineWidth',1)
plot(times,mdgy,'Color',[.2 .2 .2],'LineWidth',3)
xlabel('Time in ms','FontSize',16)
ylabel('Gender task','FontSize',16)

subplot(2,2,2); hold on
for P = 1:No
h1 = plot(times,eo(:,P)); h1.Color(4) = 0.5;
end
plot(times,mdeo_ci,'Color',[.2 .2 .2],'LineWidth',1)
plot(times,mdeo,'Color',[.2 .2 .2],'LineWidth',3)
title('Older participants','FontSize',20)

subplot(2,2,4); hold on
for P = 1:No
h1 = plot(times,go(:,P)); h1.Color(4) = 0.5;
end
plot(times,mdgo_ci,'Color',[.2 .2 .2],'LineWidth',1)
plot(times,mdgo,'Color',[.2 .2 .2],'LineWidth',3)
xlabel('Time in ms','FontSize',16)

for sub = 1:4
    subplot(2,2,sub)
    set(gca,'XLim',[-300 600],'YLim',[0 0.12],'XTick',-300:100:600)
    set(gca,'FontSize',14)
end

%% define colours

% gender
gcol = [0 .5 .5]; % soft green blue
% task effect
tcol = [.2 .2 .2]; % dark grey
% group difference
dcol = [.2 .2 .2]; % dark grey
% expressive
ecol = [1 0.5 0.2]; % orange
% young
ycol = [1 0.1 0.1]; % red
% older
ocol = [0 .5 0]; % dark green

%% task effects for each group - expressive minus gender

figure('Color','w','NumberTitle','off','Position', [1000 1000 800 600])

subplot(2,2,1); hold on
plot(times,mdey,'Color',ecol,'LineWidth',3)
plot(times,mdgy,'Color',gcol,'LineWidth',3)
legend('Expressive','Gender')
plot(times,mdey_ci,'Color',ecol,'LineWidth',1)
plot(times,mdgy_ci,'Color',gcol,'LineWidth',1)
title('Young participants','FontSize',20)
ylabel('Mutual information','FontSize',16)

subplot(2,2,3); hold on % young: task differences
title('Expressive - gender','FontSize',16)
for P = 1:Ny
h1 = plot(times,ey(:,P)-gy(:,P)); h1.Color(4) = 0.5;
end
plot([-300 600],[0 0],'k')
plot(times,mdty,'Color',tcol,'LineWidth',3)
plot(times,mdty_ci,'Color',tcol,'LineWidth',1)
xlabel('Time in ms','FontSize',16)
ylabel('Mutual information','FontSize',16)

subplot(2,2,2); hold on
plot(times,mdeo,'Color',ecol,'LineWidth',3)
plot(times,mdgo,'Color',gcol,'LineWidth',3)
legend('Expressive','Gender')
plot(times,mdeo_ci,'Color',ecol,'LineWidth',1)
plot(times,mdgo_ci,'Color',gcol,'LineWidth',1)
title('Older participants','FontSize',20)

subplot(2,2,4); hold on % older: task differences
title('Expressive - gender','FontSize',16)
for P = 1:No
h1 = plot(times,eo(:,P)-go(:,P)); h1.Color(4) = 0.5;
end
plot([-300 600],[0 0],'k')
plot(times,mdto,'Color',tcol,'LineWidth',3)
plot(times,mdto_ci,'Color',tcol,'LineWidth',1)
xlabel('Time in ms','FontSize',16)

for sub = 1:4
    subplot(2,2,sub)
    set(gca,'XLim',[-300 600],'XTick',-300:100:600)
    set(gca,'FontSize',14)
    if sub < 3
        set(gca,'YLim',[0 0.06])
    else
        set(gca,'YLim',[-0.06 0.09])
    end
end

%% group differences for each task - young minus older
% with bootstrap samples

figure('Color','w','NumberTitle','off','Position', [1000 1000 800 600])

subplot(2,2,1); hold on
plot(times,mdey,'Color',ycol,'LineWidth',3)
plot(times,mdeo,'Color',ocol,'LineWidth',3)
legend('Young','Older')
plot(times,mdey_ci,'Color',ycol,'LineWidth',1)
plot(times,mdeo_ci,'Color',ocol,'LineWidth',1)
title('Expressive','FontSize',20)
ylabel('Mutual information','FontSize',16)

subplot(2,2,3); hold on 
title('Young - older','FontSize',16)
plot(times,boot_mdey - boot_mdeo,'Color',[.7 .7 .7])
plot([-300 600],[0 0],'k')
plot(times,diff_md_e,'Color',dcol,'LineWidth',3)
plot(times,diff_md_e_ci,'Color',dcol,'LineWidth',1)
xlabel('Time in ms','FontSize',16)
ylabel('Mutual information','FontSize',16)

subplot(2,2,2); hold on
plot(times,mdgy,'Color',ycol,'LineWidth',3)
plot(times,mdgo,'Color',ocol,'LineWidth',3)
legend('Young','Older')
plot(times,mdgy_ci,'Color',ycol,'LineWidth',1)
plot(times,mdgo_ci,'Color',ocol,'LineWidth',1)
title('Gender','FontSize',20)

subplot(2,2,4); hold on 
title('Young - older','FontSize',16)
plot(times,boot_mdgy - boot_mdgo,'Color',[.7 .7 .7])
plot([-300 600],[0 0],'k')
plot(times,diff_md_g,'Color',dcol,'LineWidth',3)
plot(times,diff_md_g_ci,'Color',dcol,'LineWidth',1)
xlabel('Time in ms','FontSize',16)

for sub = 1:4
    subplot(2,2,sub)
    set(gca,'XLim',[-300 600],'XTick',-300:100:600)
    set(gca,'FontSize',14)
    if sub < 3
        set(gca,'YLim',[0 0.06])
    else
%         set(gca,'YLim',[-0.02 0.03]) % zoom
        set(gca,'YLim',[-0.03 0.04])
    end
end

%% group x task interaction
% with bootstrap samples
    
figure('Color','w','NumberTitle','off','Position', [1000 1000 400 600])

subplot(2,1,1); hold on
plot([-300 600],[0 0],'k')
h1 = plot(times,mdty,'Color',ycol,'LineWidth',3);
h2 = plot(times,mdto,'Color',ocol,'LineWidth',3);
legend([h1 h2],'Young','Older')
plot(times,mdty_ci,'Color',ycol,'LineWidth',1)
plot(times,mdto_ci,'Color',ocol,'LineWidth',1)
title('Expressive - gender','FontSize',20)
ylabel('Mutual information','FontSize',16)

subplot(2,1,2); hold on % group differences
title({'Young - Older';'(differences of differences)'},'FontSize',16)
plot(times,boot_mdty - boot_mdto,'Color',[.7 .7 .7])
plot([-300 600],[0 0],'k')
plot(times,md_interaction,'Color',gcol,'LineWidth',3)
plot(times,md_interaction_ci,'Color',gcol,'LineWidth',1)
xlabel('Time in ms','FontSize',16)
ylabel('Mutual information','FontSize',16)

for sub = 1:2
    subplot(2,1,sub)
    set(gca,'XLim',[-300 600],'XTick',-300:100:600)
    set(gca,'FontSize',14)
    if sub < 2
        set(gca,'YLim',[-0.01 0.02])
    else
%         set(gca,'YLim',[-0.015 0.015]) % zoom
    set(gca,'YLim',[-0.02 0.025])
    end
end
 
