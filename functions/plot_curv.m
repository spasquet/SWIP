function [fig,han1,han2,han3,c]=plot_curv(h,X,Y,errbar,symb,col,lw,axetop,axerev,cb,fs,xtitle,ytitle,ztitle,...
    xlimit,ylimit,zlimit,xticks,yticks,zticks,xline,yline,sizefig,sizeax,vertex)

%%% S. Pasquet - V16.6.29
%
% [fig,han1,han2,han3,c]=plot_curv(h,X,Y,errbar,symb,col,lw,axetop,axerev,cb,fs,xtitle,ytitle,ztitle,...
%    xlimit,ylimit,zlimit,xticks,yticks,zticks,xline,yline,sizefig,sizeax,vertex)
%
% h -> Figure handle
% X, Y -> vectors of length N
% errbar -> error bar (=1) or not (=0)
% symb -> symbol (e.g. 'x')
% col -> color in rgb (e.g. [0 0 0] for black)
% lw -> linewidth
% axetop -> =1 draw axis on top (default), =0 bottom
% axerev -> =1 draw Y axis pointing downward (default), =0 upward
% cb -> =0 no colorbar, =1 colorbar to the right, =2 colorbar at the bottom
% fs -> font size
% xtitle, ytitle, ztitle -> axis title and colorbar title
% xlimit, ylimit, zlimit -> axis and colorbar min and max values (zlimit
% required to display colorbar)
% xticks, yticks, zticks -> axis and colorbar ticks
% xline, yline -> 2 values vector to draw straight horizontal or vertical line
% sizefig -> 4 values vector [left, bottom, width, height] for figure
% sizeax -> 4 values vector [left, bottom, width, height] for axis
% vertex -> Vertical exageration: =0 for square axis, =1 for equal axis, >1
% increase Y axis size, <1 decrease Y axis size

matrelease=version('-release');

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
if (size(X,1)>1 && size(Y,1)>1) && (size(X,2)>1 && size(Y,2)>1)
   multiplot=1;
else
    multiplot=0;
end

if exist('col','var')==0 || isempty(col)==1
    col=[0 0 0];
end

if exist('symb','var')==0 || isempty(symb)==1
    if multiplot==1 
        symb='none';
    else
        symb='x';
    end
end

if exist('lw','var')==0 || isempty(lw)==1
    lw=1.5;
end
if exist('errbar','var')==1 && isempty(errbar)~=1
    % Plot errorbars
    if multiplot==0
        if str2double(matrelease(1:4))>2014
            han1=plot(X,Y,symb,'Color',col,'linewidth',lw,'markersize',10);
            hold on;
            herrorbars=terrorbar(X,Y,errbar,1,'units');
            set(herrorbars,'Color',col,'linewidth',lw);
        else
            han1=errorbar(X,Y,errbar,symb,'Color',col,'linewidth',lw,'markersize',10);
            errorbar_tick(han1,1,'units');
        end
    else
        for i=1:size(X,2)
            if str2double(matrelease(1:4))>2014
                han1=line(X(i),Y(i),'color',col(i,:),'markerfacecolor',col(i,:),'markersize',5,'marker',symb,'linestyle','none');
                hold on;
                herrorbars=terrorbar(X(:,i),Y(:,i),errbar(:,i),1,'units');
                set(herrorbars,'Color',col(i,:),'linewidth',lw);
            else
                han1=errorbar(X(:,i),Y(:,i),errbar(:,i),'Color',col(i,:),'linewidth',lw,'markersize',10,'marker',symb);
                errorbar_tick(han1,1,'units');
            end
            if i==1
                hold on
            end
        end
    end
else
    % Plot curves
    if multiplot==1 && min(size(col))>1
        set(gca,'ColorOrder',col);
        if size(X,1)>1 && size(X,2)>1 && size(Y,1)>1 && size(Y,2)>1
            han1=line(X,Y,'linewidth',lw,'markersize',10,'marker',symb);
        else
            for i=1:length(X)
                han1=line(X(i),Y(i),'color',col(i,:),'markerfacecolor',col(i,:),'markersize',5,'marker',symb,'linestyle','none');
            end
        end
    else
        han1=plot(X,Y,symb,'Color',col,'linewidth',lw,'markersize',10);    
    end
end
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
if exist('axerev','var')==1 && isempty(axerev)~=1 && axerev==0
    set(gca,'YDir','normal');
else
    set(gca,'YDir','reverse');
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
end
if exist('ylimit','var')==1 && isempty(ylimit)~=1
    ylim(ylimit);
end

% Change ticks
if exist('xticks','var')==1 && isempty(xticks)~=1
    set(gca,'XTick',xticks);
end
if exist('yticks','var')==1 && isempty(yticks)~=1
    set(gca,'YTick',yticks);
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
            set(cbhandle,'XTick',zticks);
        else
            set(cbhandle,'YTick',zticks);
        end
    end
    if exist('zlimit','var')==1 && isempty(zlimit)~=1
        caxis(zlimit);
    end
    set(cbhandle,'linewidth',1.5);
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