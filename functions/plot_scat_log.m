function [fig,han1,han2,han3,c]=plot_scat_log(h,X,Y,Z,marker,markersize,filled,map,axetop,axerev,cb,fs,xtitle,ytitle,ztitle,...
    xlimit,ylimit,zlimit,xticks,yticks,zticks,xline,yline,sizefig,sizeax,vertex)

%%% S. Pasquet - V18.08.14
%
% plot_scat_log(h,X,Y,Z,marker,markersize,filled,map,axetop,axerev,cb,fs,...
% xtitle,ytitle,ztitle,xlimit,ylimit,zlimit,xticks,yticks,zticks,xline,yline,sizefig,sizeax,vertex)
%
% f -> Figure handle
% X, Y -> vector of length N & M of meshgrid matrix of size N-by-M
% Z -> matrix N-by-M
% marker -> marker type
% markersize -> maker size
% filled -> marker filled (=1) or not (=0)
% map -> colormap
% axetop -> =1 draw X axis on top (default), =0 bottom
% axerev -> =1 draw Y axis pointing downward (default), =0 upward
% cb -> =0 no colorbar, =1 colorbar to the right, =2 colorbar at the bottom
% fs -> font size
% xtitle, ytitle, ztitle -> axis title
% xlimit, ylimit, zlimit -> axis min and max values
% xticks, yticks, zticks -> axis ticks
% xline, yline -> 2 values vector to draw straight horizontal or vertical line
% sizefig -> 4 values vector [left, bottom, width, height] for figure
% sizeax -> 4 values vector [left, bottom, width, height] for axis
% vertex -> Vertical exageration: =0 for square axis, =1 for equal axis, >1
% increase Y axis size, <1 decrease Y axis size

% figHandles = findall(0,'Type','figure')

matrelease=version('-release'); % Get matlab release

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

if exist('marker','var')==0 || isempty(marker)==1
    marker='o';
end

if exist('markersize','var')==0 || isempty(markersize)==1
    markersize=20;
end

% Convert data for log representation
Z(Z==0)=NaN;
Z(isinf(Z))=NaN;
Zlog=log10(Z);
mn = min(Zlog(:));
rng = max(Zlog(:))-mn;
Zlog = 1+(length(map)-1)*(Zlog-mn)/rng; % Self scale data

if exist('zlimit','var')==1 && isempty(zlimit)~=1
    zloglim=1+(length(map)-1)*(log10(zlimit)-mn)/rng;
    Zlog(Zlog<zloglim(1))=zloglim(1);
end

% Plot image
if exist('filled','var')==0 || isempty(filled)==1 || filled==1
    han1=scatter(X,Y,markersize,Zlog,marker,'filled');
else
    han1=scatter(X,Y,markersize,Zlog,marker);
end
hold on
grid on; box on;

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
    xlim([min(X(:))-min(abs(diff(X(~isnan(X))))) max(X(:))+min(abs(diff(X(~isnan(X)))))]);
end
if exist('ylimit','var')==1 && isempty(ylimit)~=1
    ylim(ylimit);
else
    ylim([min(Y(:))-min(abs(diff(Y(~isnan(Y))))) max(Y(:))+min(abs(diff(Y(~isnan(Y)))))]);
end

% Change ticks
if exist('xticks','var')==1 && isempty(xticks)~=1
    set(gca,'XTick',xticks);
end
if exist('yticks','var')==1 && isempty(yticks)~=1
    set(gca,'YTick',yticks);
end

if exist('zlimit','var')==1 && isempty(zlimit)~=1
    if zlimit(1)<=0
        zlimit(1)=0.0001;
    end
    zlimit=1+(length(map)-1)*(log10(zlimit)-mn)/rng;
    if isnan(zlimit(1))==0 && isnan(zlimit(2))==0
        caxis(zlimit);
    end
else
    zlimit=[min(min(Z)) max(max(Z))];
    if zlimit(1)<=0
        zlimit(1)=0.0001;
    end
    zlimit=1+(length(map)-1)*(log10(zlimit)-mn)/rng;
    if isnan(zlimit(1))==0 && isnan(zlimit(2))==0
        caxis(zlimit);
    end
end

% Plot colorbar
if exist('cb','var')==1 && isempty(cb)~=1 && cb==1
    c=colorbar; % Colorbar
    if str2double(matrelease(1:4))<=2014
        c=cbhandle();
    end
    if exist('ztitle','var')==1 && isempty(ztitle)~=1
        if cb==2
            if axetop==1
                set(c,'location','southoutside');
            else
                set(c,'location','northoutside');
            end
            if str2double(matrelease(1:4))>2014
                c.Label.String = ztitle;
                c.Label.Rotation = 0;
            else
                cblabel(ztitle,'Rotation', 0);
            end
        elseif cb==1
            if str2double(matrelease(1:4))>2014
                c.Label.String = ztitle;
                c.Label.Rotation = 270;
                c.Label.VerticalAlignment = 'Bottom';
            else
                cblabel(ztitle,'Rotation', 270,'VerticalAlignment','Bottom');
            end
        end
    end
    if exist('zticks','var')==1 && isempty(zticks)~=1
        ztickslog=1+(length(map)-1)*(log10(zticks)-mn)/rng; % Tick mark positions
    else
        if cb==2
            ztickslog=get(c,'Xtick');
        else
            ztickslog=get(c,'Ytick');
        end
        zticks=10.^(mn+rng*(ztickslog-1)/(length(map)-1));
        if abs(min(zticks))<1
            prec = -log10(abs(min(zticks)));
        else
            prec = 0;
        end
        zticks = round(10^prec*zticks)/10^prec;
    end
    if str2double(matrelease(1:4))<=2014
                ticklength=get(c,'TickLength');
                set(c,'XTick',zticks,'TickLength',[ticklength(1)/3 ticklength(2)]);
            else
                set(c,'XTick',zticks);
            end
    if cb==2
        if str2double(matrelease(1:4))<=2014
            ticklength=get(c,'TickLength');
            set(c,'Xtick',ztickslog,'XTicklabel',zticks,'TickLength',[ticklength(1)/3 ticklength(2)]);
        else
            set(c,'XTick',ztickslog,'YTicklabel',zticks);
        end
    else
        ztickslog(1) = ztickslog(1) + 1e-12;
        ztickslog(end) = ztickslog(end) - 1e-12;
        set(c,'Ytick',ztickslog,'YTicklabel',zticks);
    end
    ticklength=get(c,'TickLength');
    if str2double(matrelease(1:4))<=2014
        set(c,'LineWidth',1.5,'box','on','TickLength',ticklength*2);
    else
        set(c,'LineWidth',1.5,'box','on');
    end
else
    c=[];
end

% Plot vertical or horizontal line
if exist('xline','var')==1 && isempty(xline)~=1
    xL=get(gca,'XLim');
    han2=dashline([xline xline],xL,2,2,2,2,'color','k','linewidth',1.5);
else
    han2=[];
end
if exist('yline','var')==1 && isempty(yline)~=1
    yL=get(gca,'YLim');
    han3=line([yline yline],yL,'linestyle','--','color','k','linewidth',1.5);
else
    han3=[];
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