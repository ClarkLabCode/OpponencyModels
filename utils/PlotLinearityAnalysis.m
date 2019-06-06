function PlotLinearityAnalysis(modelResp, tf, cont, modelName)
% Plot the linearity analysis used in Wienecke et al. 2018 - JZV, 20180815

% Set line colors
corder2 = [0.3467    0.5360    0.6907; 0.9153    0.2816    0.2878];

%% Contrast 1/2
selectedTf = (tf==1);
selectedCont = (cont==1/2);
if any(selectedTf) && any(selectedCont)
    
    % Extract the PD and ND responses
    pd = nanmean(modelResp{selectedTf, selectedCont}(:, 1, 1),2);
    nd = nanmean(modelResp{selectedTf, selectedCont}(:, 1, 2),2);
    
    % Predict the PD and ND response from the CIS responses
    cp = sum(modelResp{selectedTf, selectedCont}(:, 1:8, 5),2)/4;
    cn = sum(modelResp{selectedTf, selectedCont}(:, 1:8, 6),2)/4;
    
    figure('Position',[200,500,500,700],'WindowStyle','docked');
    sp(1) = subplot(2,1,1);
    hold on;
    set(gca, 'colororder', corder2);
    plot(pd, '-', 'linewidth', 2);
    plot(cp, '-', 'linewidth', 2);
    
    axis('off');
    plot([0, 0, 500], [-5, -15, -15], '-k', 'linewidth', 2);
    h = text(-20, -8, '10 mV', 'horiz','right','vert','middle', 'fontsize', 16);
    set(h, 'rotation', 90);
    text(50, -15, '500 ms', 'horiz','center','vert','top', 'fontsize', 16);
    ylim([-30 30]);
    legend({'Model PD response','Linear prediction of PD response'},'location','northwest');
    legend boxoff;
    localConfAxis(16);
    title(sprintf('%s, 1 Hz, contrast 1/2', modelName));
    
    sp(2) = subplot(2,1,2);
    hold on;
    set(gca, 'colororder', corder2);
    plot(nd, '-', 'linewidth', 2);
    plot(cn, '-', 'linewidth', 2);
    
    axis('off');
    plot([0, 0, 500], [-5, -15, -15], '-k', 'linewidth', 2);
    h = text(-20, -8, '10 mV', 'horiz','right','vert','middle', 'fontsize', 16);
    set(h, 'rotation', 90);
    text(50, -15, '500 ms', 'horiz','center','vert','top', 'fontsize', 16);
    ylim([-30 30]);
    legend({'Model ND response','Linear prediction of ND response'},'location','northwest');
    legend boxoff;
    
    localConfAxis(16);
    linkaxes(sp);
end

%% Contrast 1/4
selectedTf = (tf==1);
selectedCont = (cont==1/4);

if any(selectedTf) && any(selectedCont)
    % Extract the PD and ND responses
    pd = nanmean(modelResp{selectedTf, selectedCont}(:, 1, 1),2);
    nd = nanmean(modelResp{selectedTf, selectedCont}(:, 1, 2),2);
    
    % Predict the PD and ND response from the CIS responses
    cp = nanmean(modelResp{selectedTf, selectedCont}(:, 1:8, 5),2);
    cn = nanmean(modelResp{selectedTf, selectedCont}(:, 1:8, 6),2);
    
    figure('Position',[200,500,500,700],'WindowStyle','docked');
    sp(1) = subplot(2,1,1);
    hold on;
    set(gca, 'colororder', corder2);
    plot(pd, '-', 'linewidth', 2);
    plot(cp, '-', 'linewidth', 2);
    
    axis('off');
    plot([0, 0, 500], [-5, -15, -15], '-k', 'linewidth', 2);
    h = text(-20, -8, '10 mV', 'horiz','right','vert','middle', 'fontsize', 16);
    set(h, 'rotation', 90);
    text(50, -15, '500 ms', 'horiz','center','vert','top', 'fontsize', 16);
    ylim([-30 30]);
    legend({'Model PD response','Linear prediction of PD response'},'location','northwest');
    legend boxoff;
    localConfAxis(16);
    title(sprintf('%s, 1 Hz, contrast 1/4', modelName));
    
    sp(2) = subplot(2,1,2);
    hold on;
    set(gca, 'colororder', corder2);
    plot(nd, '-', 'linewidth', 2);
    plot(cn, '-', 'linewidth', 2);
    
    axis('off');
    plot([0, 0, 500], [-5, -15, -15], '-k', 'linewidth', 2);
    h = text(-20, -8, '10 mV', 'horiz','right','vert','middle', 'fontsize', 16);
    set(h, 'rotation', 90);
    text(50, -15, '500 ms', 'horiz','center','vert','top', 'fontsize', 16);
    ylim([-30 30]);
    legend({'Model ND response','Linear prediction of ND response'},'location','northwest');
    legend boxoff;
    
    localConfAxis(16);
    linkaxes(sp);
end

end

