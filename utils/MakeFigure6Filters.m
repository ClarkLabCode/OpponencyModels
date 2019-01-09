function [ filters ] = MakeFigure6Filters(params, tf)
% Make filters for Figure 6 models - JZV, 20180824

%% Compute needed spacetime vectors

% Compute and store spatial position vector
x = -15:params.dx:15;
filters.x = x;

% Compute and store temporal position vector
t = (0:params.dt/1000:params.tEnd-params.dt/1000)';
filters.t = t;

%% Define temporal filters

[Blp, Alp] = butter(1, 1/(pi * params.tauLp / params.dt), 'low');
[Bhp, Ahp] = butter(1, 1/(pi * params.tauHp / params.dt), 'high');

filters.Blp = Blp;
filters.Alp = Alp;
filters.Bhp = Bhp;
filters.Ahp = Ahp;

filtT = zeros(length(t),1);
filtT(1) = 1;

% Expand to the time duration needed
hp = filter(Bhp, Ahp, filtT);
lp =  filter(Blp, Alp, filtT);

% Smooth responses by composing with second low-pass filter
filters.hp = filter(Blp, Alp, hp);
filters.lp = filter(Blp, Alp, lp);

% Precompute FFTs of temporal filters
filters.HP = fft(filters.hp, [], 1);
filters.LP = fft(filters.lp, [], 1);

%% Define spatial filters

normfactor = sqrt(2*pi*params.spatialStd^2) / params.dx;
filters.w1 = exp(-(x+params.eyeDist).^2/(2*params.spatialStd.^2)) / normfactor;
filters.w2 = exp(-(x).^2/(2*params.spatialStd.^2)) / normfactor;
filters.w3 = exp(-(x-params.eyeDist).^2/(2*params.spatialStd.^2)) / normfactor;

%% Compute equivalent STRFs

% Compute separable filters for each input
s1 = filters.lp*filters.w1;
s2 = filters.hp*filters.w2;
s3 = filters.lp*filters.w3;

% Grab the parameters
V1 = params.V1;
V2 = params.V2;
V3 = params.V3;
g1 = params.g1;
g2 = params.g2;
g3 = params.g3;

% Compute the full filters
num = -V1 * g1 * s1 + V2 * g2 * s2 + V3 * g3 * s3;
den = -g1 * s1 + g2 * s2 + g3 * s3;

% Normalize the filters
filters.num = num / vecnorm(num(:));
filters.den = den / vecnorm(den(:));

%% Define filters for adaptive nonlinearity model

% Define the response function for the linear filter
dsFilt = V2*g2*s2 + V3*g3*s3;

% Define the response function for the linear filter
adaptFilt = V1*g1*s1;

% Normalize the filters
filters.dsFilt = dsFilt / vecnorm(dsFilt(:));
filters.adaptFilt = adaptFilt / vecnorm(adaptFilt(:));

end