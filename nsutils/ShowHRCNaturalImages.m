function [modelOutput,vel] = ShowHRCNaturalImages(modelName,velStd,numSamples,localPath,saveVar)

    forceFilter = 0;
    
    if nargin<1 || isempty(modelName)
        modelName = 'FullReichardt';
    end

    % rectify the model output?
    if nargin<2
        velStd = 100;
    end
    
    modelFunc = str2func(modelName);

    if forceFilter
        NaturalScenesToContrast;
    end

    imageIn = load(fullfile(localPath,'/savedData/combinedFiltered2D.mat'));
    filteredScenes = imageIn.finalContrast;
%         filteredScenes = filteredScenes(:,:,1:20);
    numScenes = size(filteredScenes,3);
    imageSize = size(filteredScenes);
    xRes = 360/imageSize(2);
    imageSizeDeg = imageSize(1:2)'*xRes;
    
    %% load in the model we'll use for opponency
    sampleFreq = 100; % hz
    [modelStructure,tMin,phi] = modelFunc(sampleFreq);

    %% decide which x and y points to use and what velocity

    % generate xy points and double it so we get +- all the same
    % velociites
    xPoints = repmat(imageSizeDeg(2)*rand(numSamples/2,1),[2 1]);
    yPoints = repmat(imageSizeDeg(1)*rand(numSamples/2,1),[2 1]);

    velHalf = randn(numSamples/2,1)*velStd;
    vel = [velHalf; -velHalf];

    %% generate xt plots from the image samples
    xtPlots = MakeXtPlotSmall(filteredScenes,xPoints,yPoints,vel,sampleFreq,tMin,phi);

%     xtPlots = bsxfun(@rdivide,xtPlots,sqrt(sum(xtPlots.^2,1)));
    
    %% get model output
    % model ouput in the form of scenes by samples
    modelOutput = zeros(numScenes,numSamples);

    for ss = 1:numSamples
        modelTimeStart = tic;

        modelOutput(:,ss) = modelStructure(xtPlots(:,:,:,ss));

        modelTimeEnd = toc(modelTimeStart);

        if mod(ss,round(numSamples/10))==0
            disp([num2str(ss) '/' num2str(numSamples) ' models computed and took ' num2str(modelTimeEnd) ' seconds per model']);
            disp([num2str(modelTimeEnd*(numSamples-ss)) ' seconds remaining']);
        end
    end
    
%     save(['reichResp_' modelName '.mat'], 'vel', 'modelOutput', 'velStd')
  if saveVar
    save(['savedData/reichResp_' modelName '.mat'],'vel','modelOutput','velStd');
  end
end