function [fig,han1,han2,han3,c]=plot_img(h,X,Y,Z,map,axetop,axerev,cb,fs,xtitle,ytitle,ztitle,...
    xlimit,ylimit,zlimit,xticks,yticks,zticks,xline,yline,isoline,sizefig,sizeax,vertex,blocky)

%%% S. Pasquet - V17.05.15
%
% plot_img(h,X,Y,Z,map,axetop,axerev,cb,fs,xtitle,ytitle,ztitle,xlimit,ylimit,zlimit,...
%    xticks,yticks,zticks,xline,yline,isoline,sizefig,sizeax,vertex,blocky)
%
% h -> Figure handle
% X, Y -> vector of length N & M of meshgrid matrix of size N-by-M
% Z -> matrix N-by-M
% map -> colormap
% axetop -> =1 draw X axis on top (default), =0 bottom
% axerev -> =1 draw Y axis pointing downward (default), =0 upward
% cb -> =0 no colorbar, =1 colorbar to the right, =2 colorbar at the bottom
% fs -> font size
% xtitle, ytitle, ztitle -> axis title
% xlimit, ylimit, zlimit -> axis min and max values
% xticks, yticks, zticks -> axis ticks
% xline, yline -> 2 values vector to draw straight horizontal or vertical line
% isoline -> Vector to draw contour lines
% sizefig -> 4 values vector [left, bottom, width, height] for figure
% sizeax -> 4 values vector [left, bottom, width, height] for axis
% vertex -> Vertical exageration: =0 for square axis, =1 for equal axis, >1
% increase Y axis size, <1 decrease Y axis size
% blocky -> =0 for blocky image, =1 for smooth interp, =2 for smooth contour

% figHandles = findall(0,'Type','figure')

% Figure handle
if isempty(h)==0 && h~=0
    fig=figure(h);
elseif h==0
    fig=figure('visible','off');
else
    fig=figure;
end
set(fig,'Units','centimeters');

if exist('fs','var')==1 && isempty(fs)~=1
    set(fig,'DefaultTextFontSize',fs,'DefaultAxesFontSize',fs);
    set(gca,'FontSize',fs);
end

% Colormap
if exist('map','var')==1 && isempty(map)~=1
    colormap(map);
else
    map=colormap;
end

if exist('blocky','var')==1 && isempty(blocky)~=1 && blocky==2
    if exist('zlimit','var')==1 && isempty(zlimit)~=1
        isolevels=linspace(zlimit(1),zlimit(2),length(map)+1);
        Z(Z<zlimit(1))=zlimit(1);
    else
        isolevels=linspace(min(min(Z)),max(max(Z)),length(map)+1);
    end
end

% Plot image
if exist('blocky','var')==1 && isempty(blocky)~=1 && blocky==1
    han1=pcolor(X,Y,Z); shading interp;
elseif exist('blocky','var')==1 && isempty(blocky)~=1 && blocky==2
    han1=contourf(X,Y,Z,isolevels,'edgecolor','none');
else
    if isempty(Z(isnan(Z)))==1
        han1=imagesc(X,Y,Z);
    else
        han1=imagescnan(X,Y,Z);
    end
end
hold on
grid off; view(0,90);


% Vertical exageration
if exist('vertex','var')==1 && isempty(vertex)~=1
    if vertex==1
        axis equal
    elseif vertex==0
        axis square
    else
        daspect([vertex 1 1])
    end
end

% Plot X axis on top or bottom
if exist('axetop','var')==1 && isempty(axetop)~=1 && axetop==0
    set(gca,'XAxisLocation','bottom');
else
    set(gca,'XAxisLocation','top');
end

% Plot Y axis pointing up or down
if exist('axerev','var')==1 && isempty(axerev)~=1 && axerev==1
    set(gca,'YDir','reverse');
else
    set(gca,'YDir','normal');
end

% Write titles
if exist('xtitle','var')==1 && isempty(xtitle)~=1
    xlabel(xtitle);
end
if exist('ytitle','var')==1 && isempty(ytitle)~=1
    ylabel(ytitle);
end

% Change axis limits
if exist('xlimit','var')==1 && isempty(xlimit)~=1
    xlim(xlimit)
else
    xlim([min(X(:))-abs(median(unique(diff(X(:))))) max(X(:))+abs(median(unique(diff(X(:)))))]);
end
if exist('ylimit','var')==1 && isempty(ylimit)~=1
    ylim(ylimit);
else
    ylim([min(Y(:))-abs(median(unique(diff(Y(:))))) max(Y(:))+abs(median(unique(diff(Y(:)))))]);
end

% Change ticks
if exist('xticks','var')==1 && isempty(xticks)~=1
    set(gca,'XTick',xticks);
end
if exist('yticks','var')==1 && isempty(yticks)~=1
    set(gca,'YTick',yticks);
end

if exist('zlimit','var')==1 && isempty(zlimit)~=1
    caxis(zlimit);
else
    caxis([min(Z(:)) max(Z(:))]);
end

% Plot colorbar
if exist('cb','var')==1 && isempty(cb)~=1 && cb~=0
    c=colorbar; % Colorbar
    if exist('ztitle','var')==1 && isempty(ztitle)~=1
        if cb==2
            if axetop==1
                set(cbhandle,'location','southoutside');
            else
                set(cbhandle,'location','northoutside');
            end
            cblabel(ztitle,'Rotation', 0);
        elseif cb==1
            cblabel(ztitle,'Rotation', 270,'VerticalAlignment','Bottom');
        end
    end
    if exist('zticks','var')==1 && isempty(zticks)~=1
        if cb==2
            ticklength=get(cbhandle,'TickLength');
            set(cbhandle,'XTick',zticks,'TickLength',[ticklength(1)/3 ticklength(2)]);
        else
            set(cbhandle,'YTick',zticks);
        end
    end
    ticklength=get(cbhandle,'TickLength');
    set(cbhandle,'linewidth',1.5,'box','on','TickLength',ticklength*2);    
else
    c=[];
end

% Plot vertical or horizontal line
if exist('xline','var')==1 && isempty(xline)~=1
    xL=get(gca,'XLim');
    han2=dashline([xline xline],xL,3,3,3,3,'color',[1 0 0],'linewidth',5);
else
    han2=[];
end
if exist('yline','var')==1 && isempty(yline)~=1
    yL=get(gca,'YLim');
    han3=dashline([yline yline],yL,3,3,3,3,'color',[1 0 0],'linewidth',5);
%     han3=line([yline yline],yL,'linestyle','--','color',[0.5 0.5 0.5],'linewidth',2);
else
    han3=[];
end

if exist('isoline','var')==1 && isempty(isoline)~=1
    if length(isoline)==1
        isoline=[isoline isoline];
    end
    [cs,hc]=contour(X,Y,Z,isoline,'k','linewidth',1);
    clabel(cs, hc, 'Color', 'k', 'Rotation', 0);
end

set(gca,'TickDir','out','linewidth',1.5,'XMinorTick','on','YMinorTick','on');

% Change figure size
if exist('sizeax','var')==1 && isempty(sizeax)~=1
    set(gca,'Position',[sizeax(1),sizeax(2),sizeax(3),sizeax(4)]);
end
if exist('sizefig','var')==1 && isempty(sizefig)~=1
    set(gcf,'Position',[sizefig(1),sizefig(2),sizefig(3),sizefig(4)])
else
    sizefig=get(gcf,'Position');
end
set(gcf,'PaperSize',[sizefig(3) sizefig(4)]);
hold off
end