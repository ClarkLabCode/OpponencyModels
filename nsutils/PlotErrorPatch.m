function g=PlotErrorPatch(x,y,e,color)
    %x=x(:)';
    %y=y(:)';
    %e=e(:)';
    if nargin<4
        color = [];
    end
    
    yErrorTop=y+e;
    yErrorBottom=y-e;
    
    %Make the bottom run the opposite direction to plot around the eventual
    %shape of the error patch clockwise
    yErrorBottom=yErrorBottom(end:-1:1,:);
    ye=[yErrorTop;yErrorBottom];
    
    %Similarily run the x back
    xe=[x;x(end:-1:1,:)];
    xe = repmat(xe,[1 size(ye,2)/size(xe,2)]);
    x = repmat(x,[1 size(y,2)/size(x,2)]);

    color = repmat(color,[ceil(size(x,2)/size(color,1)) 1]);
    color = color(1:size(x,2),:);
    
    hStat = ishold;
    
    set(gca, 'ColorOrder', color, 'NextPlot', 'replacechildren');
    if size(x,1) < 20 % Previously <50
        g=errorbar(x,y,e,'marker','o','LineWidth',2,'MarkerSize',8);
    else
        g=plot(x,y,'LineWidth',2);
    end
    
    if all(e==0)
        return;
    end
    
    hold on;
    
    colormap(color);
    h=fill(xe,ye,repmat(0:size(xe,2)-1,[size(xe,1) 1]),'linestyle','none','FaceAlpha',0.25);
    
    hAnnotation = get(h,'Annotation');
    
    if ~iscell(hAnnotation)
        hAnnotation = {hAnnotation};
    end
    
    for ii = 1:length(h)
        hLegendEntry = get(hAnnotation{ii},'LegendInformation');
        set(hLegendEntry,'IconDisplayStyle','off');
    end
    
    if hStat, hold on; end
end