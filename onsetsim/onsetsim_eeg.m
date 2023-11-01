% onsetsim

corepath = pwd;
addpath([corepath,'/functions'])
ath = 0.05; % arbitrary threshold value for ttest
pth = 1 - ath; % permutation threshold
Nperm = 2000;

%% Make template

% true onset = 160 ms, F=17, max at F=26
true_onset = 160;
Xf = 0:10:500;
Nf = length(Xf);
temp1 = zeros(1,Nf);
%temp2 = [zeros(1,15), linspace(0,1,21), ones(1,15)];
erp = normpdf(linspace(-1.5,1.5,21), 0, 1);
erp = erp - min(erp);
erp = erp / max(erp);
temp2 = [zeros(1,15), erp, zeros(1,15)];
% erp = minmax_scale(norm.pdf(np.linspace(-1.5, 1.5, 21)))

figure('Color', 'w', 'NumberTitle', 'off')
hold on
plot(Xf, temp1, 'k', 'LineWidth', 3)
plot(Xf, temp2, 'r', 'LineWidth', 2)

%% Illustrate noise

% Xf = 0:1:500; % 1000 Hz
Xf = 0:10:500; % 100 Hz
Nf = length(Xf);
outvar = 1;

figure('Color', 'w', 'NumberTitle', 'off')
hold on
% plot(Xf, pinknoise(Nf), 'r', 'LineWidth',2)
plot(Xf, noise(Nf, 1, 100, outvar), 'Color', [0 0.6 0], 'LineWidth',1)
plot(Xf, noise(Nf, 1, 100, outvar), 'Color', [0.6 0 0], 'LineWidth',1)
plot(Xf, noise(Nf, 1, 100, outvar), 'Color', [0 0 0.6], 'LineWidth',1)
plot(Xf, noise(Nf, 1, 100, outvar), 'Color', [1 0 0], 'LineWidth',1)
plot(Xf, noise(Nf, 1, 100, outvar), 'Color', [0 0 0], 'LineWidth',1)


%% Template + noise -- noise.m version

Nt = 50; % number of trials per condition
cond1 = zeros(Nf, Nt);
cond2 = zeros(Nf, Nt);

outvar = 1;

for T = 1:Nt
    cond1(:,T) = temp1 + noise(Nf, 1, 100, outvar);
    cond2(:,T) = temp2 + noise(Nf, 1, 100, outvar);
end

es = (mean(cond2, 2) - mean(cond1, 2)) ./ sqrt((var(cond1,0,2) + var(cond2,0,2)) / 2);
disp(['Mean effect size = ', num2str(round(mean(es(25:27)), 2))])
%figure;plot(Xf, abs(es))

% Illustrate single-trials + means

figure('Color', 'w', 'NumberTitle', 'off')
hold on
plot(Xf, cond1, 'Color', [.7 .7 .7], 'LineWidth', 0.5)
plot(Xf, cond2, 'Color', [.6 0 0], 'LineWidth', 0.5)
plot(Xf, mean(cond1, 2), 'k', 'LineWidth', 3)
plot(Xf, mean(cond2, 2), 'r', 'LineWidth', 3)
plot(Xf, zeros(1, Nf), 'k')

%% Illustrate variance

figure;plot(Xf, var(cond2, 0, 2))

%% Illustrate means only

figure('Color', 'w', 'NumberTitle', 'off')
hold on
% plot(Xf, cond1, 'k', 'LineWidth', 1)
% plot(Xf, cond2, 'r', 'LineWidth', 1)
plot(Xf, mean(cond1, 2), 'k', 'LineWidth', 2)
plot(Xf, mean(cond2, 2), 'r', 'LineWidth', 2)
plot(Xf, zeros(1, Nf), 'k')

%% Generate data + permutation distribution + cluster test

outvar = 1;

cond1 = zeros(Nf, Nt);
cond2 = zeros(Nf, Nt);

for T = 1:Nt
    cond1(:,T) = temp1 + noise(Nf, 1, 100, outvar);
    cond2(:,T) = temp2 + noise(Nf, 1, 100, outvar);
end

% [~, ~, ~, ~, ~, tval, ~] = limo_ttest(2, cond1, cond2, ath);
tval = limo_ttest_light(2, cond1, cond2);
t2 = tval.^2;

% get permutation estimates
t2_perm = zeros(Nf,Nperm);
erpall = cat(2,cond1,cond2);
pval_perm = zeros(Nf,1);

for perm_iter = 1:Nperm
    
    perm_trials = randperm(Nt + Nt);
    perm_cond1 = erpall(:,perm_trials(1:Nt));
    perm_cond2 = erpall(:,perm_trials(Nt+1:Nt+Nt));
    %     [~, ~, ~, ~, ~, tval, ~] = limo_ttest(2,perm_cond1,perm_cond2,ath);
    tval = limo_ttest_light(2,perm_cond1,perm_cond2);
    t2_perm(:,perm_iter) = tval.^2;
    
end
maxt2_perm = max(t2_perm, [], 1);
for F = 1:Nf
    pval_perm(F) = (sum(t2_perm(F,:) > t2(F)) + 1) / (Nperm + 1);
end

% get threshold
perm_th = prctile(t2_perm, pth*100, 2); % univariate thresholds
% because there is no p value for max t^2,
% we use the univariate permutation distributions to threshold themselves
% fake p values: 0 if above threshold, 1 otherwise
pval = t2_perm < repmat(perm_th, 1, Nperm);
% threshold permutation distribution
tmp = t2_perm;
th = limo_ecluster_make(tmp, pval, ath);
% threshold T2 results
sigcluster = limo_ecluster_test(t2, t2 < perm_th, th, ath);
% find onset
onset = find_onset(sigcluster.elec_mask, Xf, 1);

% graphical check
figname = 'onsetsim: T2 CLUSTER onset estimate';
figure('Name', figname,'NumberTitle','off','Color','white');
plot(Xf, t2_perm, 'Color', [.7 .7 .7])
hold on
plot(Xf, zeros(Nf,1), 'k', 'LineWidth', 1) % zero reference line
plot(Xf, t2, 'r', 'LineWidth', 2)
plot(Xf, perm_th, 'k', 'LineWidth', 2)
v = axis;
plot(Xf, 2000.*(sigcluster.elec_mask - 1), 'go', 'MarkerSize', 4)
plot([onset onset], [v(3) v(4)], 'k:', 'LineWidth', 2)
plot([true_onset true_onset], [v(3) v(4)], 'k', 'LineWidth', 1)
set(gca, 'Layer', 'top')
ylim([v(3) v(4)])
xlim([0 500])
title(['onset = ',num2str(round(onset)),' ms'])
xlabel('Time in ms')
ylabel('T^2')
ax = gca;
ax.FontSize = 14;

% save('orit2.txt', 't2', '-ascii')
% save('permt2.txt', 't2_perm', '-ascii')

%% MAX onset estimate

max_perm_th = prctile(maxt2_perm, pth*100); % MAX stat threshold
onset = find_onset(t2 > max_perm_th, Xf, 1);

% graphical check
figname = 'onsetsim: T2 MAX onset estimate';
figure('Name', figname,'NumberTitle','off','Color','white');
plot(Xf, t2_perm, 'Color', [.7 .7 .7])
hold on
plot(Xf, zeros(Nf,1), 'k', 'LineWidth', 1) % zero reference line
plot(Xf, t2, 'r', 'LineWidth', 2)
plot([0 500], [max_perm_th max_perm_th], 'k', 'LineWidth', 2)
v = axis;
plot(Xf, 2000.*((t2 > max_perm_th) - 1), 'go', 'MarkerSize', 4)
plot([onset onset], [v(3) v(4)], 'k:', 'LineWidth', 2)
plot([true_onset true_onset], [v(3) v(4)], 'k', 'LineWidth', 1)
set(gca, 'Layer', 'top')
ylim([v(3) v(4)])
xlim([0 500])
title(['onset = ',num2str(round(onset)),' ms'])
xlabel('Time in ms')
ylabel('T^2')
ax = gca;
ax.FontSize = 14;

%% FDR onset estimate

pID = limo_FDR(pval_perm, ath);
onset = find_onset(pval_perm < pID, Xf, 1);

% graphical check
figname = 'onsetsim: T2 FDR onset estimate';
figure('Name', figname,'NumberTitle','off','Color','white');
plot(Xf, t2_perm, 'Color', [.7 .7 .7])
hold on
plot(Xf, zeros(Nf,1), 'k', 'LineWidth', 1) % zero reference line
plot(Xf, t2, 'r', 'LineWidth', 2)
plot(Xf, perm_th, 'k', 'LineWidth', 2)
v = axis;
plot(Xf, 2000.*((pval_perm < pID) - 1), 'go', 'MarkerSize', 4)
plot([onset onset], [v(3) v(4)], 'k:', 'LineWidth', 2)
plot([true_onset true_onset], [v(3) v(4)], 'k', 'LineWidth', 1)
set(gca, 'Layer', 'top')
ylim([v(3) v(4)])
xlim([0 500])
title(['onset = ',num2str(round(onset)),' ms'])
xlabel('Time in ms')
ylabel('T^2')
ax = gca;
ax.FontSize = 14;

%% Test change point function

%findchangepts(t2,'Statistic','rms', 'MaxNumChanges', 2)
res = findchangepts(t2,'Statistic','rms', 'MaxNumChanges', 2);
disp(['RMS: onset = ',num2str(Xf(res(1))),' ms'])

%% ---------------------------------------------------------
%% Simulation: true positives, t^2, Nt=50, var=1

rng(666) % set random seed

Nsim = 10000;
Nperm = 2000;

Nt = 50; % sample size per condition
outvar = 1; % single-trial noise variance

onset_cluster = zeros(1, Nsim);
onset_max = zeros(1, Nsim);
onset_fdr = zeros(1, Nsim);
cohend = zeros(1, Nsim);
onset_cp = zeros(1, Nsim);

for iter = 1:Nsim
    
    if rem(iter,100) == 0
        disp(['onsetsim tp t^2, Nt=50, var=1, EEG noise: ',num2str(iter),' / ',num2str(Nsim)])
    end
    
    % template + noise -- max sample size
    cond1 = zeros(Nf, Nt);
    cond2 = zeros(Nf, Nt);
    
    for T = 1:Nt
        cond1(:,T) = temp1 + noise(Nf, 1, 100, outvar);
        cond2(:,T) = temp2 + noise(Nf, 1, 100, outvar);
    end
    
    % save max Cohen'd
    allcohen = (mean(cond2, 2) - mean(cond1, 2)) ./ sqrt((var(cond1,0,2) + var(cond2,0,2)) / 2);
    cohend(iter) = mean(allcohen(25:27));
    
    % t-test
    tval = limo_ttest_light(2, cond1, cond2);
    t2 = tval.^2;
    
    % change point
    % find the two points at which the mean and SD change the most
    % onset = first of these two points.
    res = findchangepts(t2, 'Statistic', 'std', 'MaxNumChanges', 2);
    try
        onset_cp(iter) = Xf(res(1));
    catch
        onset_cp(iter) = NaN;
    end
    
    % get permutation estimates
    t2_perm = zeros(Nf,Nperm);
    erpall = cat(2,cond1,cond2);
    
    for perm_iter = 1:Nperm
        
        perm_trials = randperm(Nt + Nt);
        perm_cond1 = erpall(:,perm_trials(1:Nt));
        perm_cond2 = erpall(:,perm_trials(Nt+1:Nt+Nt));
        tval = limo_ttest_light(2,perm_cond1,perm_cond2);
        t2_perm(:,perm_iter) = tval.^2;
        
    end
    
    maxt2_perm = max(t2_perm, [], 1);
    pval_perm = zeros(Nf,1);
    for F = 1:Nf
        pval_perm(F) = (sum(t2_perm(F,:) > t2(F)) + 1) / (Nperm + 1);
    end
    
    % CLUSTER onset estimate
    
    % get threshold
    perm_th = prctile(t2_perm, pth*100, 2); % univariate thresholds
    % because there is no p value for max t^2,
    % we use the univariate permutation distributions to threshold themselves
    % fake p values: 0 if above threshold, 1 otherwise
    pval = t2_perm <= repmat(perm_th, 1, Nperm);
    % threshold permutation distribution
    tmp = t2_perm;
    th = limo_ecluster_make(tmp, pval, ath);
    % threshold T2 results
    sigcluster = limo_ecluster_test(t2, t2 < perm_th, th, ath);
    % find onset
    try
        onset_cluster(iter) = find_onset(sigcluster.elec_mask, Xf, 0);
    catch
        onset_cluster(iter) = NaN;
    end
    
    % MAX onset estimate
    max_perm_th = prctile(maxt2_perm, pth*100); % MAX stat threshold
    try
        onset_max(iter) = find_onset(t2 > max_perm_th, Xf, 0);
    catch
        onset_max(iter) = NaN;
    end
    
    % FDR onset estimate
    [pID,pN] = FDR(pval_perm, ath);
    try
        onset_fdr(iter) = find_onset(pval_perm < pID, Xf, 0);
    catch
        onset_fdr(iter) = NaN;
    end
    
end

save([corepath,'/onsetsim_t2_n50_var1_eegnoise'], ...
    'onset_cluster', 'onset_max', 'onset_fdr', 'cohend', ...
    'onset_cp', 'Nsim', 'Nperm')

%% Simulation results: true positives, t^2 -- one gamma/outvar/Nt
% Check results for one combination of variables.

Nsim = 10000;

true_onset = 160;
load([corepath,'/onsetsim_t2_n50_var1_eegnoise'])

% remove NaNs
oc = onset_cluster(isfinite(onset_cluster));
om = onset_max(isfinite(onset_max));
of = onset_fdr(isfinite(onset_fdr));
ocp = onset_cp(isfinite(onset_cp));
coh = cohend(isfinite(cohend));

% check true positives
tpc = length(oc) / Nsim;
tpm = length(om) / Nsim;
tpf = length(of) / Nsim;
tpcp = length(ocp) / Nsim;

% median onsets
mdo_c = median(oc);
mdo_m = median(om);
mdo_f = median(of);
mdo_cp = median(ocp);
mdo = [mdo_c mdo_m mdo_f mdo_cp];

disp('-----------------------------------------------------')
disp('True positives:')
disp(['cluster = ',num2str(round(tpc*100,1)),' %'])
disp(['max = ',num2str(round(tpm*100,1)),' %'])
disp(['fdr  = ',num2str(round(tpf*100,1)),' %'])
disp(['cp  = ',num2str(round(tpcp*100,1)),' %'])
disp('-----------------------------------------------------')
disp('Max Cohen d:')
disp(['min = ',num2str(round(min(coh),1))])
disp(['max = ',num2str(round(max(coh),1))])
disp(['median  = ',num2str(round(median(coh),1))])
disp('-----------------------------------------------------')
disp('Median onsets:')
disp(['cluster = ',num2str(round(mdo_c)),' ms'])
disp(['max = ',num2str(round(mdo_m)),' ms'])
disp(['fdr = ',num2str(round(mdo_f)),' ms'])
disp(['cp = ',num2str(round(mdo_cp)),' ms'])
disp('-----------------------------------------------------')
disp('Mean absolute error (MAE):')
disp(['cluster = ',num2str(round(mean(abs(oc-true_onset)))),' ms'])
disp(['max = ',num2str(round(mean(abs(om-true_onset)))),' ms'])
disp(['fdr = ',num2str(round(mean(abs(of-true_onset)))),' ms'])
disp(['cp = ',num2str(round(mean(abs(ocp-true_onset)))),' ms'])
disp('-----------------------------------------------------')
disp('Bias:')
disp(['cluster = ',num2str(round(median(oc)-true_onset)),' ms'])
disp(['max = ',num2str(round(median(om)-true_onset)),' ms'])
disp(['fdr = ',num2str(round(median(of)-true_onset)),' ms'])
disp(['cp = ',num2str(round(median(ocp)-true_onset)),' ms'])
disp('-----------------------------------------------------')
disp('Underestimations of at least 40 ms:')
disp(['cluster = ',num2str(round(100*mean((oc-true_onset)<= -40),1)),' %'])
disp(['max = ',num2str(round(100*mean((om-true_onset)<= -40),1)),' %'])
disp(['fdr = ',num2str(round(100*mean((of-true_onset)<= -40),1)),' %'])
disp(['cp = ',num2str(round(100*mean((ocp-true_onset)<= -40),1)),' %'])
disp('-----------------------------------------------------')
disp('Proportion too early:')
disp(['cluster = ',num2str(round(100*mean((oc-true_onset)< 0),1)),' %'])
disp(['max = ',num2str(round(100*mean((om-true_onset)< 0),1)),' %'])
disp(['fdr = ',num2str(round(100*mean((of-true_onset)< 0),1)),' %'])
disp(['cp = ',num2str(round(100*mean((ocp-true_onset)< 0),1)),' %'])
disp('-----------------------------------------------------')
disp('Variance:')
disp(['cluster = ',num2str(round(var(oc,1))),' ms'])
disp(['max = ',num2str(round(var(om,1))),' ms'])
disp(['fdr = ',num2str(round(var(of,1))),' ms'])
disp(['cp = ',num2str(round(var(ocp,1))),' ms'])

figure('NumberTitle','off', 'Name', ['nt=',num2str(Nt),', var=',num2str(outvar),', EEG noise'])

subplot(1,4,1)
[f,x] = ksdensity(oc);
plot(x,f,'k')
% histogram(oc, 'Normalization', 'count')
title(['CLUSTER median onset = ',num2str(round(mdo_c))])

subplot(1,4,2)
[f,x] = ksdensity(om);
plot(x,f,'k')
% histogram(om, 'Normalization', 'count')
title(['MAX median onset = ',num2str(round(mdo_m))])

subplot(1,4,3)
[f,x] = ksdensity(of);
plot(x,f,'k')
% histogram(of, 'Normalization', 'count')
title(['FDR median onset = ',num2str(round(mdo_f))])

subplot(1,4,4)
[f,x] = ksdensity(ocp);
plot(x,f,'k')
title(['CP median onset = ',num2str(round(mdo_cp))])

for sub = 1:4
    subplot(1,4,sub)
    hold on
    xlim([0 500])
    v = axis;
    plot([150 150], [v(3) v(4)], 'k', 'LineWidth', 1)
    plot([mdo(sub) mdo(sub)], [v(3) v(4)], 'k--', 'LineWidth', 1)
end




