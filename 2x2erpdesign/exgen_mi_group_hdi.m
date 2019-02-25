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

%% load bootstrap results

load bootres % computed in the script exgen_mi_group_ci.m

%% ==========================================================
%% get highest-density intervals - hdi() function
% the hdi() function is available on github
% <https://github.com/GRousselet/matlab_stats>

credmass = 0.8; % mass within interval = scalar between 0 and 1

mdey_hdi = zeros(Ntimes,2);
mdeo_hdi = zeros(Ntimes,2);
mdgy_hdi = zeros(Ntimes,2);
mdgo_hdi = zeros(Ntimes,2);
diff_md_g_hdi = zeros(Ntimes,2);
diff_md_e_hdi = zeros(Ntimes,2);
mdty_hdi = zeros(Ntimes,2);
mdto_hdi = zeros(Ntimes,2);
md_interaction_hdi = zeros(Ntimes,2);

for T = 1:Ntimes
    
    % for each task and group
    mdey_hdi(T,:) = hdi(boot_mdey(T,:), credmass);
    mdeo_hdi(T,:) = hdi(boot_mdeo(T,:), credmass);
    mdgy_hdi(T,:) = hdi(boot_mdgy(T,:), credmass);
    mdgo_hdi(T,:) = hdi(boot_mdgo(T,:), credmass);
    
    % for group differences for each task
    diff_md_g_hdi(T,:) = hdi(boot_mdgy(T,:) - boot_mdgo(T,:), credmass);
    diff_md_e_hdi(T,:) = hdi(boot_mdey(T,:) - boot_mdeo(T,:), credmass);
    
    % for task differences for each group
    mdty_hdi(T,:) = hdi(boot_mdty(T,:), credmass);
    mdto_hdi(T,:) = hdi(boot_mdto(T,:), credmass);
    
    % for the interaction between groups and tasks (difference of difference)
    md_interaction_hdi(T,:) = hdi(boot_mdty(T,:) - boot_mdto(T,:), credmass);
    
end

save gphdi mdey_hdi mdeo_hdi mdgy_hdi mdgo_hdi diff_md_g_hdi diff_md_e_hdi mdty_hdi mdto_hdi md_interaction_hdi

%% ===========================================================
%% plot individual MI time-courses for each group and condition
% superimpose medians and HDIs

figure('Color','w','NumberTitle','off','Position', [1000 1000 800 600])

subplot(2,2,1); hold on
plot(times,ey)
plot(times,mdey_hdi,'Color',[.2 .2 .2],'LineWidth',1)
plot(times,mdey,'Color',[.2 .2 .2],'LineWidth',2)
title('Young participants','FontSize',20)
ylabel('Expressive task','FontSize',16)

subplot(2,2,3); hold on
plot(times,gy)
plot(times,mdgy_hdi,'Color',[.2 .2 .2],'LineWidth',1)
plot(times,mdgy,'Color',[.2 .2 .2],'LineWidth',2)
xlabel('Time in ms','FontSize',16)
ylabel('Gender task','FontSize',16)

subplot(2,2,2); hold on
plot(times,eo)
plot(times,mdeo_hdi,'Color',[.2 .2 .2],'LineWidth',1)
plot(times,mdeo,'Color',[.2 .2 .2],'LineWidth',2)
title('Older participants','FontSize',20)

subplot(2,2,4); hold on
plot(times,go)
plot(times,mdgo_hdi,'Color',[.2 .2 .2],'LineWidth',1)
plot(times,mdgo,'Color',[.2 .2 .2],'LineWidth',2)
xlabel('Time in ms','FontSize',16)

for sub = 1:4
    subplot(2,2,sub)
    set(gca,'XLim',[-300 600],'YLim',[0 0.12],'XTick',-300:100:600)
    set(gca,'FontSize',14)
end

%% group differences for each task - young minus older
% superimpose medians and HDIs
% with bootstrap samples

figure('Color','w','NumberTitle','off','Position', [1000 1000 800 600])

subplot(2,2,1); hold on
plot(times,mdey,'Color',ycol,'LineWidth',3)
plot(times,mdeo,'Color',ocol,'LineWidth',3)
legend('Young','Older')
plot(times,mdey_hdi,'Color',ycol,'LineWidth',1)
plot(times,mdeo_hdi,'Color',ocol,'LineWidth',1)
title('Expressive','FontSize',20)
ylabel('Mutual information','FontSize',16)

subplot(2,2,3); hold on 
title('Young - older','FontSize',16)
plot(times,boot_mdey - boot_mdeo,'Color',[.7 .7 .7])
plot([-300 600],[0 0],'k')
plot(times,diff_md_e,'Color',dcol,'LineWidth',3)
plot(times,diff_md_e_hdi,'Color',dcol,'LineWidth',1)
xlabel('Time in ms','FontSize',16)
ylabel('Mutual information','FontSize',16)

subplot(2,2,2); hold on
plot(times,mdgy,'Color',ycol,'LineWidth',3)
plot(times,mdgo,'Color',ocol,'LineWidth',3)
legend('Young','Older')
plot(times,mdgy_hdi,'Color',ycol,'LineWidth',1)
plot(times,mdgo_hdi,'Color',ocol,'LineWidth',1)
title('Gender','FontSize',20)

subplot(2,2,4); hold on 
title('Young - older','FontSize',16)
plot(times,boot_mdgy - boot_mdgo,'Color',[.7 .7 .7])
plot([-300 600],[0 0],'k')
plot(times,diff_md_g,'Color',dcol,'LineWidth',3)
plot(times,diff_md_g_hdi,'Color',dcol,'LineWidth',1)
xlabel('Time in ms','FontSize',16)

for sub = 1:4
    subplot(2,2,sub)
    set(gca,'XLim',[-300 600],'XTick',-300:100:600)
    set(gca,'FontSize',14)
    if sub < 3
        set(gca,'YLim',[0 0.06])
    else
        set(gca,'YLim',[-0.02 0.03])
    end
end
