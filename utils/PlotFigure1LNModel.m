function PlotFigure1LNModel
%Make cartoon LN model for Figure 1 - JZV, 20180824

%% Set parameters

dx = 0.1;
dt = 1;
tauLp = 75;
tauHp = 75;
lambda = 45;
eyeDist = 5;
spatialStd = 5 / (2*sqrt(2*log(2)));
numShift = 8;

tf = 1;
numPeriods = 2;

% Set the contrasts to be tested
cont = 1/2;

inputNames = {'PD','ND','PD+ND'};

%% Compute needed spacetime vectors

% Compute spatial position vector
x = -15:dx:15;

% Compute the total simulation time
tEnd = numPeriods*ceil(tf*8*max(tauLp, tauHp)/1000)/tf;

% Compute temporal position vector
t = (0:dt/1000:tEnd-dt/1000)';

% Define phase offsets for PD + ND and PD + OD
phi = permute((0:1:numShift-1)*2*pi/numShift, [1,3,2]);

%% Define temporal filters

[Blp, Alp] = butter(1, 1/(pi * tauLp / dt), 'low');
[Bhp, Ahp] = butter(1, 1/(pi * tauHp / dt), 'high');

filtT = zeros(length(t),1);
filtT(1) = 1;

% Expand to the time duration needed
hp = filter(Bhp, Ahp, filtT);
lp =  filter(Blp, Alp, filtT);

% Smooth responses by composing with second low-pass filter
hp = filter(Blp, Alp, hp);
lp = filter(Blp, Alp, lp);

% Precompute FFTs of temporal filters
HP = fft(hp, [], 1);
LP = fft(lp, [], 1);

%% Define spatial filters

normfactor = sqrt(2*pi*spatialStd^2) / dx;
w1 = exp(-(x+eyeDist).^2/(2*spatialStd.^2)) / normfactor;
w2 = exp(-(x).^2/(2*spatialStd.^2)) / normfactor;
w3 = exp(-(x-eyeDist).^2/(2*spatialStd.^2)) / normfactor;

%% Compute input stimuli

numStim = 3;
xtPlot = cell(numStim,1);

% PD
xtPlot{1} = sin(2*pi*tf*t - (2 * pi * x / lambda));

% ND
xtPlot{2} = sin(2*pi*tf*t + (2 * pi * x / lambda + phi));

% PD + ND (CIS)
xtPlot{3} = xtPlot{1} + xtPlot{2};

% Adjust contrasts
for ind = 1:length(xtPlot), xtPlot{ind} = cont*xtPlot{ind}; end

%% Filter the inputs (in real space)

% Allocate containers
resp1 = cell(numStim, 1);
resp2 = cell(numStim, 1);
resp3 = cell(numStim, 1);

% Iterate over stimuli
for stimInd = 1:numStim
    % Note that we compute the convolutions for each component of each
    % separable filter independently to improve efficiency. Note also that
    % the spatial filter is not reversed in correctly computing the
    % convolution, as the temporal filter is effectively twice-inverted,
    % since its temporal axis is time in past, whereas the temporal axis of
    % the stimulus is time in future. 

    % Compute spatial responses
    r1 = sum(w1 .* xtPlot{stimInd}, 2);
    r2 = sum(w2 .* xtPlot{stimInd}, 2);
    r3 = sum(w3 .* xtPlot{stimInd}, 2);
    
    % Compute convolutions in time
    resp1{stimInd} = squeeze(ifft(LP .* fft(r1, [], 1), [], 1));
    resp2{stimInd} = squeeze(ifft(HP .* fft(r2, [], 1), [], 1));
    resp3{stimInd} = squeeze(ifft(LP .* fft(r3, [], 1), [], 1));
end

clearvars r1 r2 r3;


%% Compute filtered sinewaves

filtInputs = cellfun(@(r1,r2,r3) r1 + r2 - r3, resp1, resp2, resp3, 'uni', false);

% Compute half-squared input to cartoon LN model response
meanResp = cellfun(@(f) mean((f(:).*(f(:)>0)).^2), filtInputs);

%% Visualize STRF

% Compute separable filters for each input
s1 = lp*w1;
s2 = hp*w2;
s3 = lp*w3;

% Compute STRF
strf = s2 - s3;
% strf = s1 + s2 - s3;

% Plot STRF
ca = max(abs(strf(:)));
figure('Position',[200,500,500,700],'WindowStyle','docked');
imagesc(x,t,strf);
axis('square','tight');
colormap(b2r(-ca, ca));
cbar = colorbar;
cbar.Ticks = [-ca, 0, ca];
cbar.TickLabels = {'-','0','+'};
title('STRF');
xlabel('spatial location (\circ)');
ylabel('time in past (s)');
localConfAxis(16);
ylim([0 0.5]);
xlim([-10 10]);

%% Plot sinusoidal traces

% Set color order
corder = [21 114 186; 216 85 39; 129 47 140; 60, 181, 74]/256;

% Select phase of ND component of CIS
ndPhase = 3;

% Plot time traces
figure('Position',[200,500,500,700],'WindowStyle','docked');
hold on;
set(gca, 'colororder', corder);
plot(t, [filtInputs{1}(:,1), filtInputs{2}(:,ndPhase), filtInputs{3}(:, ndPhase)], 'linewidth', 2);
legend(inputNames, 'location','northwest', 'fontsize', 16);
legend('boxoff');
axis('square', 'off');

%% Plot nonlinearities

r = -1:0.01:1;

figure('Position',[200,500,500,700],'WindowStyle','docked');

% Full quadratic
subplot(2,2,1);
hold on;
plot([0 0], [0 1], 'k', 'linewidth', 1);
plot([-1 1], [0 0], 'k', 'linewidth', 1);
plot(r, r.^2, 'k','linewidth', 2);
hold off;
axis('equal','square','off');
xlim([-1 1]);
ylim([0 1]);

% Exponential
subplot(2,2,2);
hold on;
plot([0 0], [0 1], 'k', 'linewidth', 1);
plot([-1 1], [0 0], 'k', 'linewidth', 1);
plot(r, exp(2*(r-1)), 'k','linewidth', 2);
hold off;
axis('equal','square','off');
xlim([-1 1]);
ylim([0 1]);

% Half quadratic
subplot(2,2,3);
hold on;
plot([0 0], [0 1], 'k', 'linewidth', 1);
plot([-1 1], [0 0], 'k', 'linewidth', 1);
plot(r, (r.*(r>0)).^2, 'k','linewidth', 2);
hold off;
axis('equal','square','off');
xlim([-1 1]);
ylim([0 1])

% Soft rectifier
subplot(2,2,4);
hold on;
plot([0 0], [0 1], 'k', 'linewidth', 1);
plot([-1 1], [0 0], 'k', 'linewidth', 1);
plot(r, log( (1+exp(10*r)).^2 ) / 20,'k','linewidth', 2);
hold off;
axis('equal','square','off');
xlim([-1 1]);
ylim([0 1])

%% Cartoon bar plot of mean responses

% Set color order
corder = [21 114 186; 216 85 39; 129 47 140; 60, 181, 74]/256;

% Make bar plot
figure('Position',[200,500,500,700],'WindowStyle','docked');
bar(meanResp / meanResp(1), 'FaceColor', 'flat', 'CData', corder(1:3,:));
xticks(1:3);
xticklabels(inputNames)
ylabel('normalized response (arb. units)');
axis('square');
yticks(0:0.5:2);
localConfAxis(16);

end

