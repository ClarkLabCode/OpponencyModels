function g = PlotErrorPatch(x,y,err,corder)
% Plot lines with error patches. Individual lines are assumed to be row
% vectors in the arrays x, y, and err

if ~isequal(size(y), size(err)) || (size(x,1)~=size(y,1))
    error('Input arrays must be of equal size!');
end

if size(x,2) ~= size(y,2)
    x = repmat(x, 1, size(y,2));
end

if nargin < 4
    corder = lines(size(x,2));
end

yErrorTop = y + err;
yErrorBottom = y - err;

%Make the bottom run the opposite direction to plot around the eventual
%shape of the error patch clockwise
yErrorBottom = yErrorBottom(end:-1:1,:);
ye=[yErrorTop; yErrorBottom];

%Similarily run the x back
xe = [x; x(end:-1:1,:)];
xe = repmat(xe,[1 size(ye,2)/size(xe,2)]);
x = repmat(x,[1 size(y,2)/size(x,2)]);

corder = repmat(corder, [ceil(size(x,2)/size(corder,1)) 1]);
corder = corder(1:size(x,2),:);

% Get the current hold status
hStat = ishold;

set(gca, 'ColorOrder', corder, 'NextPlot', 'replacechildren');
if size(x,1) < 10
%     g=errorbar(x,y,err,'marker','o','LineWidth',2,'MarkerSize',8);
    g=errorbar(x,y,err,'marker','o','LineWidth',2,'MarkerSize',8);
else
%     g=plot(x,y,'LineWidth',2);
    g=plot(x,y,'LineWidth',1);
end

if any(err(:))
    
    hold on;
    
        colormap(corder);
        h = fill(xe,ye,repmat(0:size(xe,2)-1,[size(xe,1) 1]),'linestyle','none','FaceAlpha',0.25, 'FaceColor', 'flat');
%     h = patch(xe, ye, permute(corder, [1 3 2]), 'FaceAlpha', 0.25, 'LineStyle','None', 'FaceColor', 'flat');
    
    hAnnotation = get(h,'Annotation');
    
    if ~iscell(hAnnotation)
        hAnnotation = {hAnnotation};
    end
    
    for ii = 1:length(h)
        hLegendEntry = get(hAnnotation{ii},'LegendInformation');
        set(hLegendEntry,'IconDisplayStyle','off');
    end
end

if hStat
    hold on;
end

end