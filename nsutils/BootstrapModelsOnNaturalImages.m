function genCorr = BootstrapModelsOnNaturalImages(modelOutput,vel,velStd)
    maxVel = velStd*2;
    percentToKeep = 95;
    numBinsVel = 21;
    numBinsResp = 21;
    numScenes = size(modelOutput,1);
    
    maxAbsResp = prctile(abs(Columnize(modelOutput)),percentToKeep);

    rDiff = 2*maxAbsResp/(numBinsResp-1);
    vDiff = 2*maxVel/(numBinsVel-1);

    rEdges = linspace(-maxAbsResp-rDiff/2,maxAbsResp+rDiff/2,numBinsResp+1)';
    vEdges = linspace(-maxVel-vDiff/2,maxVel+vDiff/2,numBinsVel+1)';

%     r = rEdges(1:end-1)+diff(rEdges(1:2))/2;
    v = vEdges(1:end-1)+diff(vEdges(1:2))/2;
    genCorrRtoV = zeros(numScenes,1);

    for sc = 1:numScenes
        velocityVectThisScene = vel;
        modelOutVectThisScene = modelOutput(sc,:)';
        countsThisScene = histcounts2(modelOutVectThisScene,velocityVectThisScene,rEdges,vEdges);
        pRVThisScene = countsThisScene/sum(Columnize(countsThisScene));
        genCorrRtoV(sc) = CalcGeneralCorr(pRVThisScene',v);
    end
    
    genCorr = mean(genCorrRtoV);
end