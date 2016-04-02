%% fake_erp_differences
% This script illustrates the relationship between 2 ERP conditions and
% their difference in different situations, using fake data.
% Copyright (C) 2016 Guillaume Rousselet - University of Glasgow

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.

stdev = 10;
maxHalfSize = stdev * 3; %cuts the gaussian after 3 * stdev
[y,x] = meshgrid(-maxHalfSize:maxHalfSize, -maxHalfSize:maxHalfSize);
gauss = exp(-(x.^2/stdev^2)-(y.^2/stdev^2));
clear x y
gauss = gauss(31,:);

stdev = 17;
maxHalfSize = stdev * 3; %cuts the gaussian after 3 * stdev
[y,x] = meshgrid(-maxHalfSize:maxHalfSize, -maxHalfSize:maxHalfSize);
gauss_tmp = exp(-(x.^2/stdev^2)-(y.^2/stdev^2));
clear x y
gauss17 = gauss_tmp(round(size(gauss_tmp,1)/2),:);

%%

xf = 1:200;
htodo = zeros(length(xf),1);
ntp = length(gauss);
ntp17 = length(gauss17);

figure('Color','w','NumberTitle','off','Units','Normalized','Position',[0 0 1 .4])

% amplitude difference <<<<<<<<<<<<<<<<<<<<<<<<
subplot(2,4,1) 
hold on
title('Amplitude difference','FontSize',18)
gauss1 = htodo; ofs=30; gauss1(ofs:ntp+ofs-1) = gauss;
plot(xf,gauss1,'k','LineWidth',2)
gauss2 = htodo; ofs=30; gauss2(ofs:ntp+ofs-1) = gauss .* 0.9;
plot(xf,gauss2,'Color',[.7 .7 .7],'LineWidth',2)

subplot(2,4,1+4) 
hold on
gauss1 = htodo; ofs=30; gauss1(ofs:ntp+ofs-1) = gauss;
plot(xf,gauss1,'k','LineWidth',2)
gauss2 = htodo; ofs=30; gauss2(ofs:ntp+ofs-1) = gauss .* 0.9;
plot(xf,gauss2,'Color',[.7 .7 .7],'LineWidth',2)
diff = gauss1-gauss2;
plot(xf,diff,'r','LineWidth',2)
F = find(diff == max(diff));
plot([F F],[-2 2],'r--')

% latency difference <<<<<<<<<<<<<<<<<<<<<<<<
subplot(2,4,2) 
hold on
title('Latency difference','FontSize',18)
gauss1 = htodo; ofs=30; gauss1(ofs:ntp+ofs-1) = gauss;
plot(xf,gauss1,'k','LineWidth',2)
gauss2 = htodo; ofs=35; gauss2(ofs:ntp+ofs-1) = gauss;
plot(xf,gauss2,'Color',[.7 .7 .7],'LineWidth',2)

subplot(2,4,2+4) 
hold on
gauss1 = htodo; ofs=30; gauss1(ofs:ntp+ofs-1) = gauss;
plot(xf,gauss1,'k','LineWidth',2)
gauss2 = htodo; ofs=35; gauss2(ofs:ntp+ofs-1) = gauss;
plot(xf,gauss2,'Color',[.7 .7 .7],'LineWidth',2)
diff = gauss1-gauss2;
plot(xf,diff,'r','LineWidth',2)
F = find(diff == max(diff));
plot([F F],[-2 2],'r--')

% latency & amplitude difference (I) <<<<<<<<<<<<<<<<<<<<<<<<
subplot(2,4,3) 
hold on
title('Latency & amplitude difference (I)','FontSize',18)
gauss1 = htodo; ofs=30; gauss1(ofs:ntp+ofs-1) = gauss;
plot(xf,gauss1,'k','LineWidth',2)
gauss2 = htodo; ofs=35; gauss2(ofs:ntp+ofs-1) = gauss .* 0.8;
plot(xf,gauss2,'Color',[.7 .7 .7],'LineWidth',2)

subplot(2,4,3+4) 
hold on
gauss1 = htodo; ofs=30; gauss1(ofs:ntp+ofs-1) = gauss;
plot(xf,gauss1,'k','LineWidth',2)
gauss2 = htodo; ofs=35; gauss2(ofs:ntp+ofs-1) = gauss .* 0.8;
plot(xf,gauss2,'Color',[.7 .7 .7],'LineWidth',2)
diff = gauss1-gauss2;
plot(xf,diff,'r','LineWidth',2)
F = find(diff == max(diff));
plot([F F],[-2 2],'r--')

% latency & amplitude difference (II) <<<<<<<<<<<<<<<<<<<<<<<<
subplot(2,4,4) 
hold on
title('Latency & amplitude difference (II)','FontSize',18)
gauss1 = htodo; ofs=30; gauss1(ofs:ntp+ofs-1) = gauss .* 0.8;
plot(xf,gauss1,'k','LineWidth',2)
gauss2 = htodo; ofs=31; gauss2(ofs:ntp+ofs-1) = gauss;
plot(xf,gauss2,'Color',[.7 .7 .7],'LineWidth',2)

subplot(2,4,4+4) 
hold on
gauss1 = htodo; ofs=30; gauss1(ofs:ntp+ofs-1) = gauss.* 0.8;
plot(xf,gauss1,'k','LineWidth',2)
gauss2 = htodo; ofs=31; gauss2(ofs:ntp+ofs-1) = gauss;
plot(xf,gauss2,'Color',[.7 .7 .7],'LineWidth',2)
diff = gauss1-gauss2;
plot(xf,diff,'r','LineWidth',2)
F = find(diff == min(diff));
plot([F F],[-2 2],'r--')

for sub = 1:4
    subplot(2,4,sub)
   set(gca,'FontSize',14,'Layer','Top','YLim',[0 1.1],'XLim',[0 200])
   box on
   
   if sub == 1
      xlabel('Time','FontSize',16)
      ylabel('Amplitude','FontSize',16)
   end
   
    subplot(2,4,sub+4)
   set(gca,'FontSize',14,'Layer','Top','YLim',[-0.5 1.1],'XLim',[0 200])
   box on
end

