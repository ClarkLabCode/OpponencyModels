function [ dsi, opi, odi ] = ComputeFigure6ModelParameterSweep(params, tf, a, b, c, bFix, cFix, g1, g2, cont, filters)
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
    % Sample numShift points from the grid of phase shifts
    nn = floor(sqrt(numShift));
    numShift = nn^2;
    [phi, phi2] = meshgrid((0:1:nn-1)*2*pi/nn);
    phi = phi(:)';
    phi2 = phi2(:)';
end

%% Set up parallel pool

% % Check if a pool is open
% poolObj = gcp('nocreate');
%
% % If no pool is open, create a pool
% if isempty(poolObj)
%     poolObj = parpool('local');
% end

%% Compute input stimuli
% The format of the xtPlot array is [time x phase x numStim x space]

% Set the number of stimuli
numStim = 4;

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
% phiLin = (0:1:7)/8*pi; % Phase offsets used in Wienecke et al. 2018
% xtPlot(:,1:8,5,:) = sin((tVec + phiLin) - (xVec + phiLin)) + sin((tVec + phiLin) + (xVec + phiLin));
% xtPlot(:,1:8,6,:) = sin((tVec - phiLin) - (xVec + phiLin)) + sin((tVec - phiLin) + (xVec + phiLin));

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

n1 = length(g1);
n2 = length(g2);
dsi = nan(n1,n2, numModel);
opi = nan(n1,n2, numModel);
odi = nan(n1,n2, numModel);

%% Compute adaptive nonlinearity model response

% Hard rectifier
relu = @(f) (f .* (f>0));

% Iterate over alpha and beta with gamma fixed
tic;
for ind1 = 1:n1
    for ind2 = 1:n2
        % Compute raw linear filter response
        modelResp = (relu(1+ a(ind1) * resp2 - b(ind2)*resp3) ./ (1 + relu(-cFix*resp1))).^2;
        
        % Average responses over time and spatial phase
        meanResp = squeeze(mean(mean(modelResp,1),2));
        
        % Compute indices
        dsi(ind1,ind2,1) = (meanResp(1) - meanResp(2)) / (meanResp(1) + meanResp(2));
        opi(ind1,ind2,1) = (meanResp(3) - meanResp(1)) / (meanResp(3) + meanResp(1));
        odi(ind1,ind2,1) = (meanResp(4) - meanResp(1)) / (meanResp(4) + meanResp(1));
    end
end
fprintf('Completed adaptive nonlinearity model alpha-beta parameter sweep in %f seconds\n', toc);

% Iterate over alpha and gamma with beta fixed
tic;
for ind1 = 1:n1
    for ind2 = 1:n2
        % Compute raw linear filter response
        modelResp = (relu(a(ind1) * resp2 - bFix*resp3).^2) ./ (1 + relu(-c(ind2)*resp1).^2);
        
        % Average responses over time and spatial phase
        meanResp = squeeze(mean(mean(modelResp,1),2));
        
        % Compute indices
        dsi(ind1,ind2,2) = (meanResp(1) - meanResp(2)) / (meanResp(1) + meanResp(2));
        opi(ind1,ind2,2) = (meanResp(3) - meanResp(1)) / (meanResp(3) + meanResp(1));
        odi(ind1,ind2,2) = (meanResp(4) - meanResp(1)) / (meanResp(4) + meanResp(1));
    end
end
fprintf('Completed adaptive nonlinearity model alpha-gamma parameter sweep in %f seconds\n', toc);

%% Compute three-input biophysical nonlinearity model response
tic;

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

% Compute the needed rectified response
rect1 = relu(-resp1);
rect2 = relu(resp2);
rect3 = relu(resp3);

% Iterate over conductances
for ind1 = 1:n1
    for ind2 = 1:n2
        % Compute the numerator and denominator of the three-input model
        numResp = params.V1*g1(ind1)*rect1 + params.V2*g2(ind2)*rect2 + params.V3*g1(ind1)*rect3;
        denResp = g1(ind1)*rect1 + g2(ind2)*rect2 + g1(ind1)*rect3;
        
        % Compute the response of the full model
        modelResp = halfsquare(numResp ./ (params.gleak + denResp));
        
        % Average responses over time and spatial phase
        meanResp = squeeze(mean(mean(modelResp,1),2));
        
        % Compute indices
        dsi(ind1,ind2,3) = (meanResp(1) - meanResp(2)) / (meanResp(1) + meanResp(2));
        opi(ind1,ind2,3) = (meanResp(3) - meanResp(1)) / (meanResp(3) + meanResp(1));
        odi(ind1,ind2,3) = (meanResp(4) - meanResp(1)) / (meanResp(4) + meanResp(1));
    end
end
fprintf('Completed conductance model parameter sweep in %f seconds\n', toc);

end
