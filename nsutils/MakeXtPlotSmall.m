function xtPlot = MakeXtPlotSmall(xyPlot,xPoints,yPoints,vel,sampleFreq,tMin,phi)
    % shift the image over time at velocity vel to make an xt plot
    numSamples = size(xPoints,1);
    numScenes = size(xyPlot,3);

    tRes = 1/sampleFreq;

    xRes = 360/size(xyPlot,2);
    yPointsInd = ceil(yPoints/xRes);

    % add one extra value so you can interpolate circularly
    xyPlot = [xyPlot xyPlot(:,1,:,:)];
    x = (0:xRes:360)';
    
    numT = length(0:tRes:tMin-tRes);
    xtPlot = zeros(numT,2,numScenes,numSamples);
    
    for ss = 1:numSamples
        samTimeStart = tic;

        for sc = 1:numScenes
        
            baseSample = (0:tRes:tMin-tRes)'*-vel(ss);
            
            phis = -1:1:1;
            
            for pp = 1:length(phis)
                xSam = round(mod(baseSample+xPoints(ss)+phis(pp),360));
                xtPlot(:,pp,sc,ss) = interp1(x,xyPlot(yPointsInd(ss),:,sc),xSam);
            end
            
        end
        
        samTimeEnd = toc(samTimeStart);
        
        if mod(ss,round(numSamples/10))==0
            disp([num2str(ss) '/' num2str(numSamples) ' xt plots made and took ' num2str(samTimeEnd) ' seconds per plot']);
            disp([num2str(samTimeEnd*(numSamples-ss)) ' seconds remaining']);
        end
    end
end