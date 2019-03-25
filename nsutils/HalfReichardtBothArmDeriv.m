function [modelStructure,tMin,phi]=HalfReichardtDeriv(sampleFreq)
    % arm1 is a delta function
    % arm2 is a exponential
    
    tau = 0.15; % in seconds
    phi = 5; % in degrees
    tMin = 8*tau;
    modelStructure = @(x)HalfCorrDerivStructure(x,sampleFreq,tau);
end

function response = HalfCorrDerivStructure(xtPlot,sampleFreq,tau)    
    tStep = 1/sampleFreq;
    lengthT = size(xtPlot,1)*tStep;
    t = (0:tStep:lengthT-tStep)';
    
    tFiltD = heaviside(t).*t.*exp(-t/tau/2);
    tFiltD = tFiltD/sum(tFiltD);
    tFiltD = [diff(tFiltD); 0]/(t(2)-t(1));
    
    tFiltNd = heaviside(t).*t.*exp(-t/tau);
    tFiltNd = tFiltNd/sum(tFiltNd);
    tFiltNd = [diff(tFiltNd); 0]/(t(2)-t(1));
    
    arm1D = sum(bsxfun(@times,xtPlot(:,1,:),flipud(tFiltD)),1);
    arm2ND = sum(bsxfun(@times,xtPlot(:,2,:),flipud(tFiltNd)),1);

    response = arm1D.*arm2ND;
    response = permute(response,[3 1 2]);
end
