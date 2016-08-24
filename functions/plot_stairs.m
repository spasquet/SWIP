function [fig,han1,han2,han3]=plot_stairs(h,X,Y,symb,col,lw,axetop,fs,xtitle,ytitle,...
    xlimit,ylimit,xticks,yticks,xline,yline,sizefig,sizeax,vertex)

%%% S. Pasquet - V16.4.25
%
% [fig,han]=plot_stairs(h,X,Y,errbar,symb1,symb2,col,axetop,fs,xtitle,ytitle,
%     xlimit,ylimit,xticks,yticks,xline,yline,sizefig,sizeax,axeq)
%
% h -> Figure handle
% X, Y -> vectors of length N
% symb -> symbol for Y (e.g. '-')
% col -> color in rgb (e.g. [0 0 0] for black)
% axetop -> =1 draw axis on top (default), =0 bottom
% fs -> font size
% xtitle, ytitle, -> axis title
% xlimit, ylimit, -> axis min and max values 
% xticks, yticks -> axis ticks
% xline, yline -> 2 values vector to draw straight horizontal or vertical line
% sizefig -> 4 values vector [left, bottom, width, height] for figure
% sizeax -> 4 values vector [left, bottom, width, height] for axis
% axeq -> Vertical exageration: =0 for square axis, =1 for equal axis, >1
% increase Y axis size, <1 decrease Y axis size

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

if exist('col','var')==0 || isempty(col)==1
    col=[0 0 0];
end

if exist('symb','var')==0 || isempty(symb)==1
    symb='-';
end

if exist('lw','var')==0 || isempty(lw)==1
    lw=1.5;
end

% Plot curves
han1=stairs(X,Y,symb,'Color',col,'linewidth',lw);
hold on
grid on; view(90,90);

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

set(gca,'XAxisLocation','bottom');
% Plot X axis on top or bottom
if exist('axetop','var')==1 && isempty(axetop)~=1 && axetop==0
    set(gca,'YAxisLocation','left');
else
    set(gca,'YAxisLocation','right');
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