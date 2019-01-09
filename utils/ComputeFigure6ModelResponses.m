function [ meanResp, modelVoltageResp, modelCalciumResp, meanNumResp, meanDenResp ] = ComputeFigure6ModelResponses(params, tf, cont, filters)
% Compute responses for Figure 6 models - JZV

%% Compute needed spacetime vectors

% Compute spatial position vector
x = (-15:params.dx:15)';
numX = length(x);

% Compute temporal position vector
t = (0:params.dt/1000:params.tEnd-params.dt/1000)';
numT = length(t);

% Convert time and space vectors into units of radians and reformat to
% achieve the desired [time x phase x numStim x space] index ordering
tVec = 2 * pi * tf * t;
xVec = 2 * pi * permute(x, [4,3,2,1]) / params.lambda;

% Define stimulus presentation mask
stimMask = (t > params.tOff) & (t < (params.tEnd - params.tOff));

% Define phase offsets for PD + ND and PD + OD
numShift = params.numShift;
if params.useRandomShifts
    phi = 2 * pi * [0, rand(1, numShift-1)];
    phi2 = 2 * pi * [0, rand(1, numShift-1)];
else
    %     phi = (0:1:params.numShift-1)*2*pi/params.numShift;
    %     phi2 = phi;
    
    % Sample numShift points from the grid of phase shifts
    nn = floor(sqrt(numShift));
    numShift = nn^2;
    [phi, phi2] = meshgrid((0:1:nn-1)*2*pi/nn);
    phi = phi(:)';
    phi2 = phi2(:)';
end

%% Compute input stimuli
% The format of the xtPlot array is [time x phase x numStim x space]

% Set the number of stimuli
numStim = 6;

% Allocate a container
xtPlot = zeros(numT, numShift, numStim, numX);

% PD
xtPlot(:,:,1,:) = sin(tVec - xVec - phi);

% ND
xtPlot(:,:,2,:) = sin(tVec + xVec + phi);

% PD + ND (CIS)
xtPlot(:,:,3,:) = sin(tVec - xVec - phi) + sin(tVec + xVec + phi2);

% PD + OD
xtPlot(:,:,4,:) = sin(tVec - xVec - phi) + sin(tVec + phi2);

% CIS for linearity analysis
phiLin = (0:1:7)/8*pi; % Phase offsets used in Wienecke et al. 2018
xtPlot(:,1:8,5,:) = sin((tVec + phiLin) - (xVec + phiLin)) + sin((tVec + phiLin) + (xVec + phiLin));
xtPlot(:,1:8,6,:) = sin((tVec - phiLin) - (xVec + phiLin)) + sin((tVec - phiLin) + (xVec + phiLin));

% Adjust contrasts and apply mask
xtPlot = cont .* stimMask .* xtPlot;

%% Filter the inputs

% Note that we compute the convolutions for each component of each
% separable filter independently to improve efficiency. Note also that
% the spatial filter is not reversed in correctly computing the
% convolution, as the temporal filter is effectively twice-inverted,
% since its temporal axis is time in past, whereas the temporal axis of
% the stimulus is time in future.

% Note also that the format of the resp arrays is [time x phase x numStim]

% Compute spatial responses
r1 = sum(permute(filters.w1, [4,3,1,2]) .* xtPlot, 4);
r2 = sum(permute(filters.w2, [4,3,1,2]) .* xtPlot, 4);
r3 = sum(permute(filters.w3, [4,3,1,2]) .* xtPlot, 4);

% Compute convolutions in time
resp1 = squeeze(ifft(filters.LP .* fft(r1, [], 1), [], 1));
resp2 = squeeze(ifft(filters.HP .* fft(r2, [], 1), [], 1));
resp3 = squeeze(ifft(filters.LP .* fft(r3, [], 1), [], 1));

% Clear temporary containers
clearvars r1 r2 r3;

%% Allocate containers to store model responses

numModel = 3;
% The format of the response arrays is [time x phase x numStim x numModel]

modelVoltageResp = nan(numT, numShift, numStim, numModel);
modelCalciumResp = nan(numT, numShift, numStim, numModel);

% modelResp = cell(3, 1);
% meanResp = cell(3,1);

%% Compute two-input multiplicative model response

% Hard rectifier
relu = @(f) (f .* (f>0));

% Compute raw response
modelVoltageResp(:,:,:,1) = resp1 .* resp2;

% Rectify and compute mean response
modelCalciumResp(:,:,:,1) = relu(modelVoltageResp(:,:,:,1));

%% Compute adaptive nonlinearity model response

% Hard rectifier
relu = @(f) (f .* (f>0));

% Compute raw linear filter response
modelVoltageResp(:,:,:,2) = params.alpha * resp2 - params.beta * resp3;

% Compute adaptation linear filter response
adaptFiltResp = -params.gamma * resp1;

% Compute the response for the full adaptive nonlinearity model
modelCalciumResp(:,:,:,2) = (relu(modelVoltageResp(:,:,:,2)).^2) ./ (1 + relu(adaptFiltResp).^2);

%% Compute three-input biophysical model response

% Define the input rectifiers
if isinf(params.inputRectBeta)
    relu = @(f) (f .* (f>0));
else
    relu = @(f) f.*(erf(params.inputRectBeta * f)+1)/2;
end

% Define the output half-quadratic
if isinf(params.outputRectBeta)
    halfsquare = @(f) (f .* (f>0)).^2;
else
    halfsquare = @(f) (f.*(erf(params.outputRectBeta*f)+1)/2).^2;
end

% Compute the numerator and denominator of the three-input model
numResp = params.V1*params.g1*relu(-resp1) + params.V2*params.g2*relu(resp2) + params.V3*params.g3*relu(resp3);
denResp = params.g1*relu(-resp1) + params.g2*relu(resp2) + params.g3*relu(resp3);

% Compute the voltage response of the full model
modelVoltageResp(:,:,:,3) = numResp ./ (params.gleak + denResp);

% Model the transformation from membrane voltage to calcium concentration
% as a half-quadratic
modelCalciumResp(:,:,:,3) = halfsquare(modelVoltageResp(:,:,:,3));

% Compute averaged numerator and denominator LN responses
numResp = halfsquare(-params.V1*params.g1*resp1 + params.V2*params.g2*resp2 + params.V3*params.g3*resp3);
denResp = halfsquare(-params.g1*resp1 + params.g2*resp2 + params.g3*resp3);
meanNumResp = squeeze(nanmean(nanmean(numResp(stimMask, :, :),1),2));
meanDenResp = squeeze(nanmean(nanmean(denResp(stimMask, :, :),1),2));

%% Compute sigmoidal nonlinearity LN model response

% Compute the linear response
modelVoltageResp(:,:,:,4) = resp1 + resp2 - resp3;

% Apply the sigmoidal nonlinearity
modelCalciumResp(:,:,:,4) = 1 ./ (1 + exp(-params.sigmoidK1*(modelVoltageResp(:,:,:,4) - params.sigmoidK2)));

%% Average model responses over time and phase

meanResp = squeeze(nanmean(nanmean(modelCalciumResp(stimMask, :, :, :),1),2));

end
