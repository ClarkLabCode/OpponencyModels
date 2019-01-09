function [ params ] = SetFigure6ModelParameters()
% Set parameters for Figure 6 models - JZV, 20190109

%% Set timing parameters

% Pre-stimulus interleave 
params.tOff = 1; % s

% Stimulus duration
params.tOn = 4; % s

% Total simulation time
params.tEnd = 6;

%% Set filter parameters

% Spatial resolution (degrees)
params.dx = 0.1;

% Temporal resolution (ms)
params.dt = 1;

% Filter time constants (ms)
params.tauLp = 75;
params.tauHp = 75;

% Spatial wavelength of input sinusoids (degrees)
params.lambda = 45;

% Spatial separation of inputs (degrees)
params.eyeDist = 5;

% Standard deviation of Gaussian spatial filter (degrees)
params.spatialStd = 5 / (2*sqrt(2*log(2)));

% Time constant of GC6f lowpass filter (for timeseries only, in ms)
params.GC6fTimeConstant = 200;

%% Set spatial phase shifts
% Note that these phase shifts are not used in the linearity analysis,
% where the phase shifts are fixed to be (0:1:7)/8*pi as in Wienecke et al.

params.numShift = 256;
params.useRandomShifts = false;

%% Set hardness of rectifiers for three-input model
% For soft rectification:
% params.inputRectBeta = 2;
% params.outputRectBeta = 32;

% For hard rectification:
params.inputRectBeta = Inf;
params.outputRectBeta = Inf;

%% Set reversal potentials & conductances for three-input model

% Input 1 ~ Mi9
% Input 2 ~ Mi1
% Input 3 ~ Mi4
params.V1 = - 30;
params.V2 = + 60;
params.V3 = - 30;
params.gleak = 1;

% Used with hard rectification
params.g1 = 3;
params.g2 = 2;
params.g3 = 3;

% Used with soft rectification
% params.g1 = 1.2;
% params.g2 = 1;
% params.g3 = 1.2;

%% Set parameters for adaptive nonlinearity model

params.alpha = 300;
params.beta = 100;
params.gamma = 50;

%% Set parameters for sigmoidal nonlinearity LN model

params.sigmoidK1 = 20;
params.sigmoidK2 = 0.4;


end

