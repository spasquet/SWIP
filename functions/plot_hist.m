function [fig,han1,han2,han3]=plot_hist(h,X,N,colbin,coledge,axetop,axerev,fs,xtitle,ytitle,...
    xlimit,ylimit,xticks,yticks,xline,yline,sizefig,sizeax,vertex,xlogscale)

%%% S. Pasquet - V17.09.18
%
% [fig,han1,han2,han3,c]=plot_hist(h,X,Y,col,axetop,axerev,fs,xtitle,ytitle,...
%    xlimit,ylimit,xticks,yticks,xline,yline,sizefig,sizeax,xlogscale)
%
% h -> Figure handle
% X -> vectors of data
% N -> Number of bins in histogram
% col -> color in rgb (e.g. [0 0 0] for black)
% axetop -> =1 draw axis on top (default), =0 bottom
% axerev -> =1 draw Y axis pointing downward (default), =0 upward
% fs -> font size
% xtitle, ytitle -> axis title
% xlimit, ylimit -> axis min and max values
% xticks, yticks -> axis ticks
% xline, yline -> 2 values vector to draw straight horizontal or vertical line
% sizefig -> 4 values vector [left, bottom, width, height] for figure
% sizeax -> 4 values vector [left, bottom, width, height] for axis
% vertex -> Vertical exageration: =0 for square axis, =1 for equal axis, >1
% increase Y axis size, <1 decrease Y axis size
% xlogscale

% Figure handle
if isempty(h)==0 && h~=0
    fig=figure(h);
elseif h==0
    fig=figure('visible','off');
else
    fig=figure;
end
set(fig,'Units','centimeters');

if exist('N','var')==0 || isempty(N)==1
    N=20;
end

if exist('fs','var')==1 && isempty(fs)~=1
    set(fig,'DefaultTextFontSize',fs,'DefaultAxesFontSize',fs);
    set(gca,'FontSize',fs);
end

if exist('colbin','var')==0 || isempty(colbin)==1
    colbin=[1 0 0];
end

if exist('coledge','var')==0 || isempty(coledge)==1
    coledge=[0 0 0];
end

if exist('xlogscale','var')==0 || isempty(xlogscale)==1
    xlogscale=0;
end

% Plot curves

hist(X,N);
if xlogscale==1
    set(gca,'xscale','log');
end
hi=findobj(gca,'Type','patch');
set(hi,'FaceColor',colbin,'EdgeColor',coledge);
han1=get(gca);

grid on; box on;
hold on;

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
    set(gca,'XAxisLocation','top');
else
    set(gca,'XAxisLocation','bottom');
end

% Plot Y axis pointing up or down
if exist('axerev','var')==1 && isempty(axerev)~=1 && axerev==0
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
    han2=dashline(xL,[xline xline],2,2,2,2,'color','k','linewidth',1.5);
else
    han2=[];
end
if exist('yline','var')==1 && isempty(yline)~=1
    yL=get(gca,'YLim');
    han3=line([yline yline],yL,'linestyle','--','color','k','linewidth',1.5);
else
    han3=[];
end

ticklength=get(gca,'ticklength');
set(gca,'TickDir','out','linewidth',2,'XMinorTick','on','YMinorTick','on','ticklength',ticklength*2);

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