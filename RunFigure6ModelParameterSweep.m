function RunFigure6ModelParameterSweep()
%RunFigure6ModelParameterSweep: Sweep the parameters of the two models in
%Figure 6 for which there are free parameters. - JZV 

%% Set parameter grids

% Define parameters to test for adaptive gain model
a = (10:10:400)';
b = (10:10:400)';
c = (10:10:400)';

% Set fixed values of b and c as we sweep a & b together with c fixed and
% then sweep a & c together with b fixed
bFix = 100;
cFix = 50;

% Define conductances to test
gInh = (0.1:0.1:4)';
gExc = (0.1:0.1:4)';

%% Set other parameters

% Set temporal frequency
tf = 1;

% Set contrast
cont = 1/2;

% Define model names
modelNames = {'Dynamic gain nonlinearity', 'Conductance nonlinearity'};

% Set options for overlaid contours
contourOpts = {-1:0.1:1,'EdgeColor', 'k', 'linewidth', 2};

% Set options for points indicating chosen parameter values
pointOpts = {'o', 'MarkerSize', 10, 'Color', 'k' 'MarkerFaceColor', 'k'};

% Load colormap
load('utils/blueRedColorMap.mat', 'cmpBlueRed');

%% Compute conductance sweeps

% Add utils folder to MATLAB path
addpath('utils');

% Set parameters for models
[ params ] = SetFigure6ModelParameters();

% Make the filters
[ filters ] = MakeFigure6Filters(params, tf);

tic;
[ dsi, opi, odi ] = ComputeFigure6ModelParameterSweep(params, tf, a, b, c, bFix, cFix, gInh, gExc, cont, filters);
fprintf('Completed parameter sweeps in %f seconds\n', toc);

%% Plot DSIs

figure('Position',[200,500,500,700],'WindowStyle','docked');
imagesc(a, b, dsi(:,:,1)');
hold on;
contour(a, b, dsi(:,:,1)', contourOpts{:});
plot(params.alpha, params.beta, pointOpts{:});
axis('xy','square','tight');
caxis([-1 1]);
cbar = colorbar;
ylabel(cbar, 'DSI');
xlabel('\alpha');
ylabel('\beta');
localConfAxis(16);
colormap(cmpBlueRed);
title(modelNames{1})

figure('Position',[200,500,500,700],'WindowStyle','docked');
imagesc(a, c, dsi(:,:,2)');
hold on;
contour(a, c, dsi(:,:,2)', contourOpts{:});
plot(params.alpha, params.gamma, pointOpts{:});
axis('xy','square','tight');
caxis([-1 1]);
cbar = colorbar;
ylabel(cbar, 'DSI');
xlabel('\alpha');
ylabel('\gamma');
localConfAxis(16);
colormap(cmpBlueRed);
title(modelNames{1})

figure('Position',[200,500,500,700],'WindowStyle','docked');
imagesc(gInh, gExc, dsi(:,:,3)');
hold on;
contour(gInh, gExc, dsi(:,:,3)', contourOpts{:});
plot(params.g1, params.g2, pointOpts{:});
axis('xy','square','tight');
caxis([-1 1]);
cbar = colorbar;
ylabel(cbar, 'DSI');
xlabel('g_{inh} / g_{leak}');
ylabel('g_{exc} / g_{leak}');
localConfAxis(16);
colormap(cmpBlueRed);
title(modelNames{2});

%% Plot opponency indicies

figure('Position',[200,500,500,700],'WindowStyle','docked');
imagesc(a, b, opi(:,:,1)');
hold on;
contour(a, b, opi(:,:,1)', contourOpts{:});
plot(params.alpha, params.beta, pointOpts{:});
axis('xy','square','tight');
caxis([-0.5 0.5]);
cbar = colorbar;
ylabel(cbar, 'opponency index');
xlabel('\alpha');
ylabel('\beta');
localConfAxis(16);
colormap(cmpBlueRed);
title(modelNames{1});

figure('Position',[200,500,500,700],'WindowStyle','docked');
imagesc(a, b, opi(:,:,2)');
hold on;
contour(a, b, opi(:,:,2)', contourOpts{:});
plot(params.alpha, params.gamma, pointOpts{:});
axis('xy','square','tight');
caxis([-0.5 0.5]);
cbar = colorbar;
ylabel(cbar, 'opponency index');
xlabel('\alpha');
ylabel('\gamma');
localConfAxis(16);
colormap(cmpBlueRed);
title(modelNames{1});

figure('Position',[200,500,500,700],'WindowStyle','docked');
imagesc(gInh, gExc, opi(:,:,3)');
hold on;
contour(gInh, gExc, opi(:,:,3)', contourOpts{:});
plot(params.g1, params.g2, pointOpts{:});
axis('xy','square','tight');
caxis([-0.5 0.5]);
cbar = colorbar;
ylabel(cbar, 'opponency index');
xlabel('g_{inh} / g_{leak}');
ylabel('g_{exc} / g_{leak}');
localConfAxis(16);
colormap(cmpBlueRed);
title(modelNames{2});

%% Plot OD indicies

figure('Position',[200,500,500,700],'WindowStyle','docked');
imagesc(a, b, odi(:,:,1)');
hold on;
contour(a, b, odi(:,:,1)', contourOpts{:});
plot(params.alpha, params.beta, pointOpts{:});
axis('xy','square','tight');
caxis([-0.5 0.5]);
cbar = colorbar;
ylabel(cbar, 'orthogonal direction index');
xlabel('\alpha');
ylabel('\beta');
localConfAxis(16);
colormap(cmpBlueRed);
title(modelNames{1})

figure('Position',[200,500,500,700],'WindowStyle','docked');
imagesc(a, b, odi(:,:,2)');
hold on;
contour(a, b, odi(:,:,2)', contourOpts{:});
plot(params.alpha, params.gamma, pointOpts{:});
axis('xy','square','tight');
caxis([-0.5 0.5]);
cbar = colorbar;
ylabel(cbar, 'orthogonal direction index');
xlabel('\alpha');
ylabel('\gamma');
localConfAxis(16);
colormap(cmpBlueRed);
title(modelNames{1});

figure('Position',[200,500,500,700],'WindowStyle','docked');
imagesc(gInh, gExc, odi(:,:,3)');
hold on;
contour(gInh, gExc, odi(:,:,3)', contourOpts{:});
plot(params.g1, params.g2, pointOpts{:});
axis('xy','square','tight');
caxis([-0.5 0.5]);
cbar = colorbar;
ylabel(cbar, 'orthogonal direction index');
xlabel('g_{inh} / g_{leak}');
ylabel('g_{exc} / g_{leak}');
localConfAxis(16);
colormap(cmpBlueRed);
title(modelNames{2});

end

