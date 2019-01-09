function RunFigure6Models()
%Run models from Figure 6 and S6 - JZV

%% Set parameters

% Add utils folder to MATLAB path
addpath('utils');

% Set parameters for models
[ params ] = SetFigure6ModelParameters();

% For TF sweep:
% tf = 2.^(-3:0.5:7)';

% For timetraces:
tf = 1;

% Set the contrasts to be tested
cont = [1/4, 1/2];

% Define stimulus names
numStim = 6;
inputNames = {'PD','ND','PD+ND','PD+OD', 'CISSumPD', 'CISSumND'};

% Define model names
modelNames = {'Rectified multiplier', 'Dynamic gain nonlinearity', 'Conductance nonlinearity', 'Sigmoidal LN'};

%% Run the models

% Get the number of TFs, the number of contrasts, and the number of models
numTf = length(tf);
numCont = length(cont);
numModel = length(modelNames);

% Allocate containers
meanResp = cell(numTf,numCont);
voltResp = cell(numTf,numCont);
modelResp = cell(numTf,numCont);
numResp = nan(numTf, numCont, numStim);
denResp = nan(numTf, numCont, numStim);

% Iterate over temporal frequencies
for tfInd = 1:numTf
    
    % Iterate over contrasts
    for contInd = 1:numCont
        tic;
        
        % Make filters
        [ filters ] = MakeFigure6Filters(params, tf(tfInd));
        
        % Run models
        [ meanResp{tfInd, contInd}, voltResp{tfInd, contInd}, modelResp{tfInd, contInd}, numResp(tfInd, contInd, :), denResp(tfInd, contInd, :) ] = ...
            ComputeFigure6ModelResponses(params, tf(tfInd), cont(contInd), filters);
        
        % Print a status update to the terminal
        fprintf('tf %d of %d, contrast %d of %d, %f\n', tfInd, numTf, contInd, numCont, toc);
    end
end

%% Visualize filters

% Make filters
params.numPeriods = 1;
[ filters ] = MakeFigure6Filters(params, 1);

% Visualize filters
PlotFigure6ModelFilters(filters);

%% Make bar plots

% Iterate over models
for ind = 1:numModel
    
    % Get the response for the selected model
    selResp = cellfun(@(x) x(:,ind), meanResp, 'uni', false);
    
    % Make bar plots
    if contains(modelNames{ind}, {'three', 'Three'})
        MakeModelBarPlots(selResp, tf, cont, modelNames{ind}, inputNames, numResp, denResp);
    else
        MakeModelBarPlots(selResp, tf, cont, modelNames{ind}, inputNames);
    end
end

%% Plot TF sweeps

% Only plot if multiple TFs were tested
if length(tf) > 1
    
    % Iterate over models
    for ind = 1:numModel
        
        % Get the response for the selected model
        selResp = cellfun(@(x) x(:,ind), meanResp, 'uni', false);
        
        % Plot TF sweeps
        if contains(modelNames{ind}, {'three', 'Three'})
            PlotModelTFSweep(selResp, tf, cont, modelNames{ind}, inputNames, numResp, denResp);
        else
            PlotModelTFSweep(selResp, tf, cont, modelNames{ind}, inputNames);
        end
    end
end

%% Plot time traces (averaged over spatial phase)

% Iterate over models
for ind = 1:numModel
    
    % Get the response for the selected model
    selResp = cellfun(@(x) x(:,:,:,ind), modelResp, 'uni', false);
    selMeanResp = cellfun(@(x) x(:,ind), meanResp, 'uni', false);
    
    % Make bar plots
    PlotModelTimeTracesAveragedOverSpatialPhase(params, selResp, selMeanResp, tf, cont, modelNames{ind}, inputNames);
end


%% Plot linearity analysis from Weinecke et al. 2018
ind = 3;

% Get the response for the selected model
selResp = cellfun(@(x) x(:,:,:,ind), voltResp, 'uni', false);

% Plot linearity analysis
PlotLinearityAnalysis(selResp, tf, cont, modelNames{ind});

end
