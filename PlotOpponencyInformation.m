function PlotOpponencyInformation(localPath, modelNames, forceReichCalc, saveV)
addpath('nsutils');

if nargin < 1
    localPath = pwd;
end

if nargin<2 || isempty(modelNames)
    modelNames = {'HalfReichardt' ,'FullReichardt', 'RectHalfReichardt','RectFullReichardt' };
end

if ~iscell(modelNames)
    modelNames = {modelNames};
end

if nargin<3
    forceReichCalc = true;
    
end
saveVar = false;
if nargin > 4
    saveVar = saveV;
end

numModels = length(modelNames);
modelOutput = cell(numModels,1);
vel = cell(numModels,1);
velStd = cell(numModels,1);
rect = false(numModels,1);

for mm = 1:numModels
    if isequal(modelNames{mm}(1:4),'Rect')
        rect(mm) = true;
        modelNames{mm} = modelNames{mm}(5:end);
    end
end

for mm = 1:numModels
    if forceReichCalc
        velStd{mm} = 100;
        numSamples = 1000;
        [modelOutput{mm},vel{mm}] = ShowHRCNaturalImages(modelNames{mm},velStd{mm},numSamples,localPath,saveVar);
    else
        savedData = load(fullfile(localPath,['/savedData/reichResp_' modelNames{mm} '.mat']));
        vel{mm} = savedData.vel;
        modelOutput{mm} = savedData.modelOutput;
        velStd{mm} = savedData.velStd;
    end
    
    %% rectify the output
    if rect(mm)
        modelOutput{mm}(modelOutput{mm}<0) = 0;
    end
end

%% general values
numScenes = size(modelOutput{1},1);
numBinsVel = 11;
numBinsResp = 11;
numPlotPoints = 11;

%% calculate generalized correlation across natural images
genCorrRtoV = zeros(numScenes,numModels);

pRV = cell(numModels,1);
pR = cell(numModels,1);
pV = cell(numModels,1);
pVGivenR = cell(numModels,1);
pRGivenV = cell(numModels,1);
eVGivenR = cell(numModels,1);
eRGivenV = cell(numModels,1);
semRGivenV = cell(numModels,1);
stdVGivenR = cell(numModels,1);
stdRGivenV = cell(numModels,1);

rAll = cell(numModels,1);
vAll = cell(numModels,1);

tickV = cell(numModels,1);
tickR = cell(numModels,1);

genCorrCi = zeros(numModels,2);
nBoot = 1000;

responseMult = zeros(numModels,1);

for mm = 1:numModels
    %% make histogram for each individual scene
    %         maxVel = std(vel{mm})*2;
    maxVel = velStd{mm}*2;
    percentToKeep = 95;
    maxAbsResp = prctile(abs(Columnize(modelOutput{mm})),percentToKeep);
    
    rDiff = 2*maxAbsResp/(numBinsResp-1);
    vDiff = 2*maxVel/(numBinsVel-1);
    
    rEdges = linspace(-maxAbsResp-rDiff/2,maxAbsResp+rDiff/2,numBinsResp+1)';
    vEdges = linspace(-maxVel-vDiff/2,maxVel+vDiff/2,numBinsVel+1)';
    
    r = rEdges(1:end-1)+diff(rEdges(1:2))/2;
    v = vEdges(1:end-1)+diff(vEdges(1:2))/2;
    
    for sc = 1:numScenes
        velocityVectThisScene = vel{mm};
        modelOutVectThisScene = modelOutput{mm}(sc,:)';
        countsThisScene = histcounts2(modelOutVectThisScene,velocityVectThisScene,rEdges,vEdges);
        pRVThisScene = countsThisScene/sum(Columnize(countsThisScene));
        genCorrRtoV(sc,mm) = CalcGeneralCorr(pRVThisScene',v);
    end
    
    bootFun = @(x)BootstrapModelsOnNaturalImages(x,vel{mm},velStd{mm});
    %
    tic;
    %         genCorrCi(mm,:) = bootci(nBoot,{bootFun,modelOutput{mm}},'alpha',0.01);
    toc;
    
    %% make histogram for all scenes together
    velocityVect = Columnize(repmat(vel{mm},[1 numScenes]));
    modelOutVect = Columnize(modelOutput{mm}');
    
    % define v and r
    maxVel = velStd{mm}*2;
    percentToKeep = 95;
    maxAbsResp = prctile(abs(modelOutVect),percentToKeep);
    
    rEdges = linspace(-maxAbsResp-rDiff/2,maxAbsResp+rDiff/2,numBinsResp+1)';
    vEdges = linspace(-maxVel-vDiff/2,maxVel+vDiff/2,numBinsVel+1)';
    
    rAll{mm} = rEdges(1:end-1)+diff(rEdges(1:2))/2;
    vAll{mm} = vEdges(1:end-1)+diff(vEdges(1:2))/2;
    
    tickV{mm} = vAll{mm}(round(linspace(1,numBinsVel,numPlotPoints)));
    tickR{mm} = rAll{mm}(round(linspace(1,numBinsResp,numPlotPoints)));
    
    % get join probability distribution
    counts = histcounts2(modelOutVect,velocityVect,rEdges,vEdges);
    
    %% get information metrics
    pRV{mm} = counts/sum(Columnize(counts));
    pR{mm} = sum(pRV{mm},2);
    pV{mm} = sum(pRV{mm},1)';
    pVGivenR{mm} = bsxfun(@rdivide,pRV{mm},pR{mm});
    pRGivenV{mm} = bsxfun(@rdivide,pRV{mm},pV{mm}');
    eVGivenR{mm} = sum(bsxfun(@times,pVGivenR{mm},vAll{mm}'),2);
    stdVGivenR{mm} = sqrt(sum(bsxfun(@times,pVGivenR{mm},(vAll{mm}.^2)'),2));
    
    eRGivenV{mm} = sum(bsxfun(@times,pRGivenV{mm},rAll{mm}),1)';
    %         eRGivenVSplit{mm} = [eRGivenV{mm}(ceil(end/2):end) flipud(eRGivenV{mm}(1:ceil(end/2)))];
    stdRGivenV{mm} = sqrt(sum(bsxfun(@times,pRGivenV{mm},rAll{mm}.^2),1))';
    semRGivenV{mm} = stdRGivenV{mm}/sqrt(numScenes);
    %         semRGivenVSplit{mm} = [semRGivenV{mm}(ceil(end/2):end) flipud(semRGivenV{mm}(1:ceil(end/2)))];
    %         vAllSplit{mm} = vAll{mm}(ceil(end/2):end);
    
    responseMult(mm) = 1/max(eRGivenV{mm});
end

%% plotting

meanGenCorr = mean(genCorrRtoV,1);
%     genCorrCi(:,1) = prctile(genCorrRtoV,95,1);
%     genCorrCi(:,2) = prctile(genCorrRtoV,5,1);

pRVVect = log(Columnize(cell2mat(pRV)))/log(10);
minLogPRV = min(pRVVect(~isinf(pRVVect)));
maxLogPRV = max(pRVVect(~isinf(pRVVect)));
%     maxConditional = max(Columnize([cell2mat(pVGivenR) cell2mat(pRGivenV)]));
maxConditional = max(Columnize([cell2mat(pVGivenR)]));

pRVect = log(Columnize(cell2mat(pR)))/log(10);
minLogPR = min(pRVect(~isinf(pRVect)));
maxLogPR = max(pRVect(~isinf(pRVect)));

pVVect = log(Columnize(cell2mat(pV)))/log(10);
minLogPV = min(pVVect(~isinf(pVVect)));
maxLogPV = max(pVVect(~isinf(pVVect)));

figure;
% plot general correlation
bar((1:length(meanGenCorr))',meanGenCorr,'b');
hold on;
PlotErrBars((1:length(meanGenCorr))',meanGenCorr,[],{genCorrCi(:,1) genCorrCi(:,2)},'k.');
hold off;
ConfAxis();

MakeFigure;
for ii = 1:numModels
    % plot p(v|r)
    subplot(numModels,2,2*ii-1);
    %         imagesc(vAll{ii},rAll{ii},pRGivenV{ii});
    %         set(gca,'YDir','normal');
    %         colorbar;
    %         caxis([0 ceil(maxConditional*20)/20]);
    title('p(r|v)');
    %         PlotXvsY(vAll{ii},[eRGivenV{ii} eRGivenV{ii}]*responseMult(ii),'error',[semRGivenV{ii} stdRGivenV{ii}]*responseMult(ii));
    PlotXvsY(vAll{ii},eRGivenV{ii}*responseMult(ii),'error',semRGivenV{ii}*responseMult(ii));
    hold on;
    PlotConstLine(0,1);
    PlotConstLine(0,2);
    hold off;
    ConfAxis('labelX',['velocity (' char(186) '/s)'],'labelY','response');
    %         ConfAxis('tickX',tickV{ii},'tickLabelX',round(tickV{ii}),'labelX','velocity','tickY',tickR{ii},'tickLabelY',round(tickR{ii}*responseRes),'labelY','response');
    
    % plot p(v|r)
    subplot(numModels,2,2*ii);
    imagesc(rAll{ii},vAll{ii},pVGivenR{ii}');
    set(gca,'YDir','normal');
    colorbar;
    caxis([0 ceil(maxConditional*20)/20]);
    title('p(v|r)');
    hold on;
    plot(rAll{ii},eVGivenR{ii},'k');
    %         plot(rAll{ii},eVGivenR{ii}+stdVGivenR{ii},'color',[0.25 0.25 0.25]);
    %         plot(rAll{ii},eVGivenR{ii}-stdVGivenR{ii},'color',[0.25 0.25 0.25]);
    PlotConstLine(0,1);
    PlotConstLine(0,2);
    hold off;
    ConfAxis('tickX',tickR{ii},'tickLabelX',round(tickR{ii}*responseMult(ii)),'labelX','response','tickY',tickV{ii},'tickLabelY',round(tickV{ii}),'labelY','velocity');
end
colormap(b2r(0, ceil(maxConditional*20)/20));

for ii = 1:numModels
    MakeFigure;
    % plot p(r,v)
    subplot(4,4,2:4);
    plot(vAll{ii},log(pV{ii})/log(10));
    xlim([vAll{ii}(1) vAll{ii}(end)]);
    ConfAxis();
    ylim([floor(minLogPV) ceil(maxLogPV)]);
    subplot(4,4,5:4:16);
    plot(log(pR{ii})/log(10),rAll{ii}*responseMult(ii));
    xlim([floor(minLogPR) ceil(maxLogPR)]);
    set(gca,'XDir','reverse');
    ylim([rAll{ii}(1) rAll{ii}(end)]*abs(responseMult(ii)));
    ConfAxis();
    subplot(4,4,[6:8 10:12 14:16]);
    imagesc(vAll{ii},rAll{ii},log(pRV{ii})/log(10));
    caxis([minLogPRV maxLogPRV]);
    set(gca,'YDir','normal');
    colorbar;
    title('p(v,r)');
    hold on;
    PlotConstLine(0,1);
    PlotConstLine(0,2);
    hold off;
    ConfAxis('tickX',tickV{ii},'tickLabelX',round(tickV{ii}),'labelX','velocity','tickY',tickR{ii},'tickLabelY',round(tickR{ii}*responseMult(ii)),'labelY','response');
    caxis([floor(minLogPRV) ceil(maxLogPRV)]);
    colormap(b2r(minLogPRV, maxLogPRV));
end

end

