function [modelStructure,tMin,phi]=FullReichardt(sampleFreq)
    % arm1 is a delta function
    % arm2 is a exponential
    
    tau = 0.15; % in seconds
    phi = 5; % in degrees
    tMin = 8*tau;
    modelStructure = @(x)HalfCorrStructure(x,sampleFreq,tau);
end

function response = HalfCorrStructure(xtPlot,sampleFreq,tau)    
    tStep = 1/sampleFreq;
    lengthT = size(xtPlot,1)*tStep;
    t = (0:tStep:lengthT-tStep)';
    
    tFiltD = heaviside(t).*exp(-t/tau);
    tFiltD = tFiltD/sum(tFiltD);
    
    arm1ND = xtPlot(end,1,:);
    arm1D = sum(bsxfun(@times,xtPlot(:,1,:),flipud(tFiltD)),1);
    arm2ND = xtPlot(end,2,:);
    arm2D = sum(bsxfun(@times,xtPlot(:,2,:),flipud(tFiltD)),1);

    response = arm1D.*arm2ND-arm1ND.*arm2D;
    response = permute(response,[3 1 2]);
end