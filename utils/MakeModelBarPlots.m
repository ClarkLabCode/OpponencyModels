function MakeModelBarPlots(meanResp, tf, cont, modelName, inputNames, numResp, denResp)
% Make bar plots for Figure 6 models - JZV, 20180815

% Set bar color order
corder = [21 114 186; 216 85 39; 129 47 140; 60, 181, 74]/256;

selTF = 1;

% Plot bar graph at selected TF and contrast 1/2
selectedTf = (tf==selTF);
selectedCont = (cont==1/2);
if any(selectedTf) && any(selectedCont)
    figure('Position',[200,500,500,700],'WindowStyle','docked');
    plotData = meanResp{selectedTf, selectedCont}(1:4);
    bar(plotData / plotData(1), 'FaceColor', 'flat', 'CData', corder);
    xticks(1:4);
    xticklabels(inputNames)
    ylabel('normalized response (arb. units)');
    axis('square');
    yticks(0:0.5:2);
    localConfAxis(16);
    title(sprintf('%s, 1 Hz, contrast 1/2', modelName));
end

% Plot bar graph at selected TF and contrast 1/4
selectedTf = (tf==selTF);
selectedCont = (cont==1/4);
if any(selectedTf) && any(selectedCont)
    figure('Position',[200,500,500,700],'WindowStyle','docked');
    plotData = meanResp{selectedTf, selectedCont}(1:4);
    bar(plotData / plotData(1), 'FaceColor', 'flat', 'CData', corder);
    xticks(1:4);
    xticklabels(inputNames)
    ylabel('normalized response (arb. units)');
    axis('square');
    yticks(0:0.5:2);
    localConfAxis(16);
    title(sprintf('%s, 1 Hz, contrast 1/4', modelName));
end

if nargin > 5
    % Plot bar graph of numerator response at selected TF and contrast 1/2
    selectedTf = (tf==selTF);
    selectedCont = (cont==1/2);
    if any(selectedTf) && any(selectedCont)
        figure('Position',[200,500,500,700],'WindowStyle','docked');
        plotData = squeeze(numResp(selectedTf, selectedCont, 1:4));
        bar(plotData, 'FaceColor', 'flat', 'CData', corder);
        xticks(1:4);
        xticklabels(inputNames)
        ylabel('numerator response (arb. units)');
        axis('square');
        % yticks(0:0.5:2);
        localConfAxis(16);
        title(sprintf('%s, 1 Hz, contrast 1/2', modelName));
    end
    
    % Plot bar graph of numerator response at selected TF and contrast 1/4
    selectedTf = (tf==selTF);
    selectedCont = (cont==1/4);
    if any(selectedTf) && any(selectedCont)
        figure('Position',[200,500,500,700],'WindowStyle','docked');
        plotData = squeeze(numResp(selectedTf, selectedCont, 1:4));
        bar(plotData, 'FaceColor', 'flat', 'CData', corder);
        xticks(1:4);
        xticklabels(inputNames)
        ylabel('numerator response (arb. units)');
        axis('square');
        % yticks(0:0.5:2);
        localConfAxis(16);
        title(sprintf('%s, 1 Hz, contrast 1/4', modelName));
    end
    
    % Plot bar graph of denominator response at selected TF and contrast 1/2
    selectedTf = (tf==selTF);
    selectedCont = (cont==1/2);
    if any(selectedTf) && any(selectedCont)
        figure('Position',[200,500,500,700],'WindowStyle','docked');
        plotData = squeeze(denResp(selectedTf, selectedCont, 1:4));
        bar(plotData, 'FaceColor', 'flat', 'CData', corder);
        xticks(1:4);
        xticklabels(inputNames)
        ylabel('denominator response (arb. units)');
        axis('square');
        % yticks(0:0.5:4);
        localConfAxis(16);
        title(sprintf('%s, 1 Hz, contrast 1/2', modelName));
    end
    
    % Plot bar graph of denominator response at selected TF and contrast 1/4
    selectedTf = (tf==selTF);
    selectedCont = (cont==1/4);
    if any(selectedTf) && any(selectedCont)
        figure('Position',[200,500,500,700],'WindowStyle','docked');
        plotData = squeeze(denResp(selectedTf, selectedCont, 1:4));
        bar(plotData, 'FaceColor', 'flat', 'CData', corder);
        xticks(1:4);
        xticklabels(inputNames)
        ylabel('denominator response (arb. units)');
        axis('square');
        % yticks(0:0.5:10);
        localConfAxis(16);
        title(sprintf('%s, 1 Hz, contrast 1/4', modelName));
    end
end

end

