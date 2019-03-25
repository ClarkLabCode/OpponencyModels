function SpatialCorrFigure(localPath)
addpath('nsutils');

if nargin < 1
    localPath = pwd;
end

%% load images
allImages = load(fullfile(localPath, '\savedData\combinedFiltered2D.mat'));
allImages = allImages.scenes;

imageInd = 1:421;

% organizing the images as (x,y) instead of (y,x) might improve
% performance
imageIn = bsxfun(@rdivide,allImages(:,:,imageInd),max(max(allImages(:,:,imageInd))));
imageSize = size(imageIn);

vel = (50:50:200)';
numVel = length(vel);
sampleFreq = 1000;
phi = 5;

xRes = 360/imageSize(2);

% x = (xRes:2*xRes:360)';
x = xRes;
y = (xRes:1:imageSize(1)*xRes)';

[xMat,yMat] = meshgrid(x,y);

xMatCol = Columnize(xMat);
yMatCol = Columnize(yMat);
xCorrMat = cell(numVel,1);
xCorrMatStd = cell(numVel,1);
dt = cell(numVel,1);

for vv = 1:numVel
    velCol = repmat(vel(vv),[length(xMatCol) 1]);
    tMin = 360/vel(vv);
    xtPlot = MakeXtPlotSmall(imageIn,xMatCol,yMatCol,velCol,sampleFreq,tMin,phi);
    xtPlot = bsxfun(@minus,xtPlot,mean(xtPlot,1));
    xtPlot = bsxfun(@rdivide,xtPlot,std(xtPlot,[],1));
    fftLeft = fft(xtPlot(:,1,:,:),[],1)/sqrt(size(xtPlot,1));
    fftRight = fft(xtPlot(:,2,:,:),[],1)/sqrt(size(xtPlot,1));
    xCorrMatAll = fftshift(ifft(conj(fftLeft).*fftRight,[],1));
    xCorrMat{vv} = mean(mean(xCorrMatAll,4),3);
    xCorrMatStd{vv} = std(mean(xCorrMatAll,4),[],3);
    dt{vv} = (((1:sampleFreq*tMin)-floor(sampleFreq*tMin/2)-1)/sampleFreq)';
end

%% plot
plotColor = lines(numVel);

MakeFigure;
hold on;
for vv = 1:numVel
    %         PlotXvsY(dt{vv}*1000,xCorrMat{vv},'error',xCorrMatStd{vv},'plotColor',plotColor(vv,:));
    plot(dt{vv}*1000,xCorrMat{vv},'MarkerEdgeColor',plotColor(vv,:));
    %     PlotConstLine(phi/vel(vv),2);
end
xlim([-0.3 0.3]*1000);
ylim([0 1]);
PlotConstLine(0,2);
hold off;
velLegend = cellfun(@num2str,(num2cell(vel)),'UniformOutput',0);
ConfAxis('labelX','dt (ms)','labelY','correlation','figLeg',velLegend);
end