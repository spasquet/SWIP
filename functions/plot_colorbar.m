function fig = plot_colorbar(h, dims, orientation, title_string, cmap, Clogscale, cticks, val, fs, sizefig, sizeax, prec)
% PLOT_COLORBAR plot a standalone colorbar for inclusion in a publication
%
% Matt Foster <ee1mpf@bath.ac.uk>
% Modified by S. Pasquet - V16.11.18
%
% fig = plot_colorbar(h, dims, orientation, title_string, cmap, Clogscale, cticks, val, fs, sizefig, sizeax, prec)

error(nargchk(1, 12, nargin, 'struct'))

% Extract the width froms dims, if there is one.
if length(dims) < 2
    width = 1;
else
    width = dims(2);
end

if nargin < 5
    cmap = jet(32);
end

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
map = [cmap(1,:);cmap;cmap(end,:)] ;

if Clogscale == 1
    cticksplot = log10(cticks);
    val = log10(val);
else
    cticksplot = cticks;
end
valplot = linspace(val(1),val(2),length(map));
dv = unique(round(diff(valplot)*1e6)/1e6);
valplot = [valplot(1)-dv valplot valplot(end)+dv];
valplot = valplot-dv/2;

% test=10;
% while ismember(0,diff(round(test*cticks)/test))==1
%     test=test*10;
% end
if abs(min(cticks))<1
    prec = ceil(-log10(abs(min(cticks))))+1;
else
    prec = 0;
end
cticks = round(10^prec*cticks)/10^prec;

switch orientation
    case {1}
        image(1:width,valplot,repmat(cat(3, map(:,1), map(:,2), map(:,3)), 1, width));
        ylabel(title_string);
        ylim([val(1)-dv val(2)]);

        % Remove ticks we dont want.
        set(gca, 'xtick', 0);
        xlabel('      ');
        
        set(gca,'position',[sizeax(1:2) sizeax(3)/20 sizeax(4)]);
        set(gca,'YTick',cticksplot,'YTickLabel',cticks,...
            'Yaxislocation','right','ydir','normal');
        set(get(gca,'YLabel'),'Rotation',270,'VerticalAlignment','bottom');
    case {2}
        image(valplot,1:width,repmat(cat(3, map(:,1)', map(:,2)', map(:,3)'),width, 1));
        xlabel(title_string);
        xlim([val(1)-dv val(2)]);

        % Remove ticks we dont want.
        set(gca, 'ytick', 0);
        
        set(gca,'position',[sizeax(1) sizeax(2)*2 sizeax(4)*0.75 sizeax(3)/18]);
        set(gca,'XTick',cticksplot,'XTickLabel',cticks,...
            'Xaxislocation','bottom');
        
    otherwise
        error('unknown colorbar orientation');
end

set(gca,'TickDir','in','linewidth',1.5);

if exist('sizefig','var')==1 && isempty(sizefig)~=1
    set(gcf,'Position',[sizefig(1),sizefig(2),sizefig(3),sizefig(4)])
else
    sizefig=get(gcf,'Position');
end
set(gcf,'PaperSize',[sizefig(3) sizefig(4)]);
hold off