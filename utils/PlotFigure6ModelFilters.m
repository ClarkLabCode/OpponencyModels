function PlotFigure6ModelFilters(filters)
% Visualize filters for Figure 6 models - JZV, 20180824

% Load colormap
load('blueRedColorMap.mat', 'cmpBlueRed');

% Plot temporal filters
figure('Position',[200,500,500,700],'WindowStyle','docked');
plot(filters.t, filters.lp, 'linewidth', 2);
hold on;
plot([0; filters.t], [0; filters.hp], 'linewidth', 2);
xlabel('time in past (s)');
ylabel('filter magnitude');
legend({'low-pass','high-pass'});
localConfAxis(16);

% Plot spatial filters
figure('Position',[200,500,500,700],'WindowStyle','docked');
plot(filters.x, [filters.w1; filters.w2; filters.w3], 'linewidth', 2);
xlabel('spatial position (\circ)');
ylabel('filter magnitude');
legend({'h_1','h_2', 'h_3'});
localConfAxis(16);

%% Three-input model effective filters

% Plot three-input effective filter
ca = max(max(max(abs(filters.num), abs(filters.den))));
figure('Position',[200,500,500,700],'WindowStyle','docked');
imagesc(filters.x,filters.t,filters.num);
axis('square','tight');
colormap(cmpBlueRed);
cbar = colorbar;
cbar.Ticks = [-ca, 0, ca];
caxis([-ca ca]);
ylim([0 1]);
cbar.TickLabels = {'-','0','+'};
title('numerator');
xlabel('spatial location (\circ)');
ylabel('time in past (s)');
localConfAxis(16);

% Plot three-input denominator effective filter
figure('Position',[200,500,500,700],'WindowStyle','docked');
imagesc(filters.x,filters.t,filters.den);
axis('square','tight');
colormap(cmpBlueRed);
cbar = colorbar;
cbar.Ticks = [-ca, 0, ca];
caxis([-ca ca]);
ylim([0 1]);
cbar.TickLabels = {'-','0','+'};
title('denominator');
xlabel('spatial location (\circ)');
ylabel('time in past (s)');
localConfAxis(16);

%% Adaptive nonlinearity filters

% Plot DS filter
ca = max(max(max(abs(filters.dsFilt), abs(filters.adaptFilt))));
figure('Position',[200,500,500,700],'WindowStyle','docked');
imagesc(filters.x,filters.t,filters.dsFilt);
axis('square','tight');
colormap(cmpBlueRed);
cbar = colorbar;
cbar.Ticks = [-ca, 0, ca];
caxis([-ca ca]);
ylim([0 1]);
cbar.TickLabels = {'-','0','+'};
title('Adaptive nonlinearity DS filter');
xlabel('spatial location (\circ)');
ylabel('time in past (s)');
localConfAxis(16);

% Plot adaptation filter
ca = max(max(max(abs(filters.dsFilt), abs(filters.adaptFilt))));
figure('Position',[200,500,500,700],'WindowStyle','docked');
imagesc(filters.x,filters.t,-filters.adaptFilt);
axis('square','tight');
colormap(cmpBlueRed);
cbar = colorbar;
cbar.Ticks = [-ca, 0, ca];
caxis([-ca ca]);
cbar.TickLabels = {'-','0','+'};
title('Adaptive nonlinearity adaptation filter');
xlabel('spatial location (\circ)');
ylabel('time in past (s)');
ylim([0 1]);
localConfAxis(16);

end

