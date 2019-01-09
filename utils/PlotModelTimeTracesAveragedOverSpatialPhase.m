function PlotModelTimeTracesAveragedOverSpatialPhase(params, modelResp, meanResp, tf, cont, modelName, inputNames)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

% Set bar color order
corder = [21 114 186; 216 85 39; 129 47 140; 60, 181, 74]/256;

selTF = 1;
selStim = 1:3; % Exclude PD+OD
% selStim = 1:4; % Include PD+OD

% Plot bar graph at selected TF and contrast 1/2
selectedTf = (tf==selTF);
selectedCont = (cont==1/2);

if any(selectedTf) && any(selectedCont)
    pdResp = meanResp{selectedTf, selectedCont}(1);
    
    % Filter in time to mimic calcium indicator dynamics
    [Blp, Alp] = butter(1, 1/(pi * params.GC6fTimeConstant / params.dt), 'low');
    plotData = modelResp{selectedTf, selectedCont}(:,:,selStim);
    plotData = filter(Blp, Alp, plotData , [], 1);
    
    plotData = plotData / pdResp;
    
    meanTimeSeries = squeeze(mean(plotData,2));
    sdTimeSeries = squeeze(std(plotData,0,2)) / sqrt(params.numShift);
    
    % Time in seconds
    t = (0:params.dt/1000:params.tEnd-params.dt/1000)';
    t = t-params.tOff;
    
    figure('Position',[200,500,500,700],'WindowStyle','docked');
    area([0;t(end)-params.tOff], [2;2], 'FaceColor', [1 1 1]/2, 'FaceAlpha', 1/2, 'EdgeColor','None');
    hold on;
    PlotErrorPatch(t, meanTimeSeries, sdTimeSeries, corder);
    hold off;
    ylim([0 2]);
    yticks(0:0.5:2);
    
    title(sprintf('%s, 1 Hz, contrast 1/2', modelName));
    xlabel('time (s)');
    ylabel('response averaged over spatial phase (arb. units)');
    legend(inputNames(selStim));
    axis('square');
    localConfAxis(16);
    
end

end

