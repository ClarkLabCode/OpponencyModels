function PlotModelTFSweep(meanResp, tf, cont, modelName, inputNames, numResp, denResp)
% Plot TF sweeps for Figure 6 models - JZV, 20180815

corder = [21 114 186; 216 85 39; 129 47 140; 60, 181, 74]/256;

% TF sweep of full model at contrast 1/2
selectedCont = (cont==1/2);
if any(selectedCont)
    figure('Position',[200,500,500,700],'WindowStyle','docked');
    hold on;
    set(gca, 'colororder', corder);
    
    plotData = cat(2, meanResp{:, selectedCont});
    plotData = plotData(1:4,:);
    plot(log2(tf), plotData, 'linewidth', 2);
    legend(inputNames(1:4));
    xticks(-3:1:7);
    xticklabels(cellstr(num2str((-3:1:7)', '2^{%d}')));
    xlabel('tf (Hz)');
    ylabel('model response (arb. units)');
    title(sprintf('%s, contrast 1/2', modelName));
    localConfAxis(16);
end

% TF sweep of full model at contrast 1/4
selectedCont = (cont==1/4);
if any(selectedCont)
    figure('Position',[200,500,500,700],'WindowStyle','docked');
    hold on;
    set(gca, 'colororder', corder);
    plotData = cat(2, meanResp{:, selectedCont});
    plotData = plotData(1:4,:)';
    plot(log2(tf), plotData, 'linewidth', 2);
    legend(inputNames(1:4));
    xticks(-3:1:7);
    xticklabels(cellstr(num2str((-3:1:7)', '2^{%d}')));
    xlabel('tf (Hz)');
    ylabel('model response (arb. units)');
    title(sprintf('%s, contrast 1/4', modelName));
    localConfAxis(16);
end

if nargin > 5
    % TF sweep of numerator at contrast 1/2
    selectedCont = (cont==1/2);
    if any(selectedCont)
        figure('Position',[200,500,500,700],'WindowStyle','docked');
        hold on;
        set(gca, 'colororder', corder);
        plotData = squeeze(numResp(:, selectedCont, 1:4));
        plot(log2(tf), plotData, 'linewidth', 2);
        legend(inputNames(1:4));
        xticks(-3:1:7);
        xticklabels(cellstr(num2str((-3:1:7)', '2^{%d}')));
        xlabel('tf (Hz)');
        ylabel('numerator response (arb. units)');
        title(sprintf('%s, contrast 1/2', modelName));
        localConfAxis(16);
    end
    
    % TF sweep of numerator at contrast 1/4
    selectedCont = (cont==1/4);
    if any(selectedCont)
        figure('Position',[200,500,500,700],'WindowStyle','docked');
        hold on;
        set(gca, 'colororder', corder);
        plotData = squeeze(numResp(:, selectedCont, 1:4));
        plot(log2(tf), plotData, 'linewidth', 2);
        legend(inputNames(1:4));
        xticks(-3:1:7);
        xticklabels(cellstr(num2str((-3:1:7)', '2^{%d}')));
        xlabel('tf (Hz)');
        ylabel('numerator response (arb. units)');
        title(sprintf('%s, contrast 1/4', modelName));
        localConfAxis(16);
    end
    
    % TF sweep of denominator at contrast 1/2
    selectedCont = (cont==1/2);
    if any(selectedCont)
        figure('Position',[200,500,500,700],'WindowStyle','docked');
        hold on;
        set(gca, 'colororder', corder);
        plotData = squeeze(denResp(:, selectedCont, 1:4));
        plot(log2(tf), plotData, 'linewidth', 2);
        legend(inputNames(1:4));
        xticks(-3:1:7);
        xticklabels(cellstr(num2str((-3:1:7)', '2^{%d}')));
        xlabel('tf (Hz)');
        ylabel('denominator response (arb. units)');
        title(sprintf('%s, contrast 1/2', modelName));
        localConfAxis(16);
    end
    
    % TF sweep of denominator at contrast 1/4
    selectedCont = (cont==1/4);
    if any(selectedCont)
        figure('Position',[200,500,500,700],'WindowStyle','docked');
        hold on;
        set(gca, 'colororder', corder);
        plotData = squeeze(denResp(:, selectedCont, 1:4));
        plot(log2(tf), plotData, 'linewidth', 2);
        legend(inputNames(1:4));
        xticks(-3:1:7);
        xticklabels(cellstr(num2str((-3:1:7)', '2^{%d}')));
        xlabel('tf (Hz)');
        ylabel('denominator response (arb. units)');
        title(sprintf('%s, contrast 1/4', modelName));
        localConfAxis(16);
    end
end

end

