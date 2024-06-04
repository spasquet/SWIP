function [fig,han1]=plot_wiggle(h,seismomat,xseis,tseis,scal,clip,perc,fs,xtitle,...
    ttitle,xlimit,tlimit,xticks,tticks,sizefig,sizeax)

%%% S. Pasquet - V16.11.18
% plot_wiggle: Plot seismic data using wiggles
% plot_wiggle(h,seismomat,xseis,tseis,scal,clip,perc,fs,xtitle,ttitle,xlimit,tlimit,xticks,tticks,sizefig,sizeax)

% Figure handle
if isempty(h)==0 && h~=0
    fig=figure(h);
elseif h==0
    fig=figure('visible','off');
else
    fig=figure;
end
set(fig,'Units','centimeters');
grid on;

if exist('fs','var')==1 && isempty(fs)~=1
    set(fig,'DefaultTextFontSize',fs,'DefaultAxesFontSize',fs);
    set(gca,'FontSize',fs);
end

[nt,nx]=size(seismomat);

trmx= max(abs(seismomat));
amx=trmx;
if exist('clip','var')==0 || isempty(clip)==1
    clip=1;
end
if exist('perc','var')==0 || isempty(perc)==1
    perc=100;
end
if (nargin <= 2)
    xseis=1:nx;
    tseis=[1:nt]*0.002;
end
if (nargin <= 4)
    scal=1;
end

if nx <= 1
    disp(' ERR:PlotWig: nx has to be more than 1');
    return
end

% take the average as dx
dx1 = abs(xseis(2:nx)-xseis(1:nx-1));
dx = median(dx1);

dz=tseis(2)-tseis(1);
% xmx=max(max(seismomat));
% xmn=min(min(seismomat));

if isempty(scal) || scal == 0
    scal=1;
end
seismomat=bsxfun(@rdivide,seismomat*dx/2,amx);
% seismomat=bsxfun(@rdivide,seismomat,amx);


npts=round(size(seismomat,1)*size(seismomat,2)*perc/100);
allpts=sort(abs(reshape(seismomat,1,size(seismomat,1)*size(seismomat,2))));
p=allpts(npts);

% seismomat=seismomat/p;
% seismomat = seismomat*dx;

% seismomat(seismomat>dx)=dx; 
% seismomat(seismomat<-dx)=-dx;

seismomat = seismomat*scal;

if clip==1
   seismomat(seismomat>dx/1.5)=dx/1.5; 
   seismomat(seismomat<-dx/1.5)=-dx/1.5;
end

% set display range
x1=min(xseis)-1.0*dx; x2=max(xseis)+1.0*dx;
z1=min(tseis)-dz; z2=max(tseis)+dz;

set(gca,'NextPlot','add','Box','on', ...
    'XLim', [x1 x2], 'YDir','reverse', ...
    'YLim',[z1 z2]);

fillcolor = [0 0 0];
linecolor = [0 0 0];
linewidth = 0.5;

tseis=tseis'; 	% input as row vector
zstart=tseis(1);
zend  =tseis(nt);

for i=1:nx
    if trmx(i) ~= 0;    % skip the zero traces
        tr=seismomat(:,i); 	% --- one scale for all section
        s = sign(tr) ;
        i1= find( s(1:nt-1) ~= s(2:nt) );	% zero crossing points
        if isempty(i1)
            if isempty(find(s>0, 1))
                i1 = 1;
            else
                i1 = nt-1;
            end
        end
        
        try
            zadd = i1 + tr(i1) ./ (tr(i1) - tr(i1+1)); %locations with 0 amplitudes
        catch
            keyboard
        end
        aadd = zeros(size(zadd));
        
        [zpos] = find(tr >0);
        [zz,iz] = sort([zpos; zadd]); 	% indices of zero point plus positives
        aa = [tr(zpos); aadd];
        aa = aa(iz);
        
        % be careful at the ends
        if tr(1)>0
            a0=0; z0=1.00;
        else
            try
            a0=0; z0=zadd(1);
            catch
                keyboard
            end
        end;
        if tr(nt)>0
            a1=0; z1=nt;
        else
            a1=0; z1=max(zadd);
        end;
        
        zz = [z0; zz; z1; z0];
        aa = [a0; aa; a1; a0];
        
        
        zzz = zstart + zz*dz -dz;
        
        patch( aa+xseis(i) , zzz,  fillcolor);
        
        line( 'Color',[1 1 1], ...
            'Xdata', xseis(i)+[0 0], 'Ydata',[zstart zend]); % remove zero line 
        
        line( 'Color',linecolor,...
            'LineWidth',linewidth, ...
            'Xdata', tr+xseis(i), 'Ydata',tseis);	% negatives line
        
    else % zeros trace
        line( 'Color',linecolor,...
            'LineWidth',linewidth, ...
            'Xdata', [xseis(i) xseis(i)], 'Ydata',[zstart zend]);
    end
end

% Plot X axis on bottom
    set(gca,'XAxisLocation','top');

% Plot Y axis pointing down
    set(gca,'YDir','reverse');

% Write titles
if exist('xtitle','var')==1 && isempty(xtitle)~=1
    xlabel(xtitle);
end
if exist('ttitle','var')==1 && isempty(ttitle)~=1
    ylabel(ttitle);
end

% Change axis limits
if exist('xlimit','var')==1 && isempty(xlimit)~=1
    xlim(xlimit)
end
if exist('tlimit','var')==1 && isempty(tlimit)~=1
    ylim(tlimit);
end

% Change ticks
if exist('xticks','var')==1 && isempty(xticks)~=1
    set(gca,'XTick',xticks);
end
if exist('tticks','var')==1 && isempty(tticks)~=1
    set(gca,'YTick',tticks);
end

set(gca,'TickDir','out','linewidth',1.5,'YMinorTick','on');
han1=gca;

% Change figure size
if exist('sizeax','var')==1 && isempty(sizeax)~=1
    set(gca,'Position',[sizeax(1),sizeax(2),sizeax(3),sizeax(4)]);
end
if exist('sizefig','var')==1 && isempty(sizefig)~=1
    set(gcf,'Position',[sizefig(1),sizefig(2),sizefig(3),sizefig(4)])
end
hold off;
end
