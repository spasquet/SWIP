
function [xx,zz,xzv,nodes,MP,input,xsrc,xrec] = genmesh(input,topo,dx,dz,maxdepth,dl,min_vel,max_vel,grad,makeplot)
%
% Generate a mesh, velocity model, and raytracing nodes
%
% Modified by S. Pasquet (03/2016)
% Modified for crosshole tomography (10/2018)
%
% Usage:
%
%   [xx,zz,xzv,nodes,MP,input,xsrc,xrec] = genmesh(input,topo,dx,dz,maxdepth,dl,min_vel,max_vel,grad,makeplot)
%
%
% Inputs:  input- tomo_input file
%          topo- topography file
%          dx- x dimension of cells
%          dz- min and max z dimension of cells (increases linearly)
%          maxdepth- maximum depth of model
%          dl- raytracing node spacing
%          min_vel- minimum velocity
%          max_vel- maximum velocity
%          grad- starting gradient
%          makeplot- 1 = plot the model 0= don't;
%
% To get a layered starting model min_vel must list the velocity for the
%       top of each layer and grad must look like this:
%       grad = [grad1 depth_to_bottom_of_layer1;
%               grad2 depth_to_bottom_of_layer2;
%               etc....]
%
%
%
% Outputs:
%          xx,zz -mesh
%          xvz- velocity model
%          nodes- raytracing nodes
%          MP- model parameter id
%          input- slightly modified input array
%          xsrc- x-coordinates of sources
%          xrec- x-coordinates of geophones
%


% if numel(dz) == 1;
%     dz = [dz dz];
% end


% generate x vector and topography of model

if length(dx) == 1;
    x = min(topo(:,1))-4*dx:dx:max(topo(:,1))+4*dx;
else
    x = [dx max(dx)+mode(diff(dx))];
end
z = interp1(topo(:,1),topo(:,2),x,'linear','extrap');

% generate z vector
if length(dz) == 1;
    dz = [dz dz + 0.001*dz];
end
Nz =ceil( maxdepth / (dz(1) + 0.5*(dz(2)-dz(1)) ));

ddz = (dz(2)-dz(1))/Nz;

zm = (1:Nz)';
zm = zm*dz(1) + 0.5*ddz*zm.^2;

ddz = diff([0; zm]);


zz = z;
xx = x;

if length(grad) == 1;
    
    vv = min_vel*ones(size(xx));
    
    for i = 1:length(zm);
        
        zz = [zz; z-zm(i)];
        xx = [xx; x];
        vv = [vv; vv(i,:) + grad*(ddz(i))];
        
    end
    
else
    
    vv = min_vel(1)*ones(size(xx));
    
    cnt = 1;
    
    for i = 1:length(ddz);
        
        if abs(sum(ddz(1:i))) > abs(grad(cnt,2));
            cnt = cnt+1;  vv(i,:) = min_vel(cnt)*ones(size(x));
        end
        
        zz = [zz; z-zm(i)];
        xx = [xx; x];
        vv = [vv; vv(i,:) + grad(cnt,1) * ddz(i)];
        
    end
    
end

vv(end,:) = []; vv(:,end) = [];
vv(vv>max_vel) = max_vel;
max_vel = max(max(vv));

% get coordinates of cell center and provide velocity

[m,n] = size(xx);


nmodel = (m-1)*(n-1);

m = m-1;

xzv(nmodel,1) = struct('xz',zeros(1,2),'v',0,'u',0,'att',0,'X',zeros(4,1),'Z',zeros(4,1));

for i = 1:nmodel;
    
    cc = ceil(i/m);
    
    xzv(i,1).xz = [(xx(i-(cc-1)*m,cc) + xx(i-(cc-1)*m,cc+1)) / 2 ...
        (zz(i-(cc-1)*m,cc) + zz(i-(cc-1)*m,cc+1)) / 2 - ddz(i-(cc-1)*m)/2];
    xzv(i,1).v = vv(i);
    xzv(i,1).u = vv(i)^-1;
    
    xzv(i,1).X = [xx(i-(cc-1)*m,cc)  xx(i-(cc-1)*m,cc+1) ...
        xx(1+i-(cc-1)*m,cc+1)  xx(1+i-(cc-1)*m,cc)];
    xzv(i,1).Z = [zz(i-(cc-1)*m,cc)  zz(i-(cc-1)*m,cc+1) ...
        zz(1+i-(cc-1)*m,cc+1)  zz(1+i-(cc-1)*m,cc)];
    
end

% generate nodes

xnodes = min(x)+dx/5:dl:max(x)-dx/5;
znodes = (max(zz(:)):-dl:min(zz(:)));
nodes = makegrid(xnodes,znodes);


% determine neighbors
upbound = interp1(xx(1,:),zz(1,:),xnodes,'linear','extrap');
bbound = interp1(xx(1,:),zz(end-1,:),xnodes,'linear','extrap');


% remove nodes that live outside of model
for i = 1:length(xnodes);
    
    I = find(nodes(:,1) == xnodes(i));
    J = (nodes(I,2) < bbound(i));
    nodes(I(J),:) = [];
    
    
    I = find(nodes(:,1) == xnodes(i));
    J = (nodes(I,2) > upbound(i));
    nodes(I(J),:) = [];
    
end


% fix receiver and shot locations to interpolated topography

xsurf = unique(input(:,2));
zsurf = interp1(x,z,xsurf,'linear','extrap');
topospline=spline(topo(:,1),topo(:,2));

for i = 1:length(input(:,1));
    diffelev=ppval(topospline,input(i,2))-input(i,3);
    I = input(i,2) == xsurf;
    if diffelev>1
        input(i,3) = zsurf(I)-diffelev;
    else
        input(i,3) = zsurf(I);
    end
end

%
% I = input(:,1) ~= 0;
% xzrec = unique(input(I,2:3),'rows');
%
%
% % put receivers in nodes
%
% nodes = [xzrec; nodes];


% determine which model cell each node lives in

MP = zeros(nmodel,1);

if length(dx) > 1; x(end) = []; end;

for i = 1:length(nodes(:,1));
    
    JJ = floor((nodes(i,1)-min(x))/dx)+1; 
    JJ = [JJ JJ+1];
    
    testx = xx(:,JJ(1):JJ(2));
    testz = zz(:,JJ(1):JJ(2));
    
    testx = testx - testx(1,1);
    
    px = nodes(i,1) - testx(1,1);
    
    sl =( testz(1,2) - testz(1,1) )/ ( testx(1,2) - testx(1,1) );
    
    pz = sl*-px + nodes(i,2); pz = round(pz*1e9)*1e-9;
    if all(testz(:,1) < pz)
        pz = max(testz(:,1));
    end
    
    testz(:,1) = round(testz(:,1)*1e9)*1e-9;
    
    II = find(testz(:,1) >= pz,1,'last');
    if II > Nz
        II = Nz;
    end
    
    try
        MP(i) = sub2ind(size(vv),II,JJ(1));
    catch
        continue
    end
    
end

% find out where the sources live
srcid = unique(input(:,1)); srcid(1) = [];

for i = 1:sum(input(:,1) == 0);
    
    KK = input(:,4) == srcid(i);
    JJ = floor((input(KK,2)-min(x))/dx)+1;
    JJ = [JJ JJ+1];
    
    testx = xx(:,JJ(1):JJ(2));
    testz = zz(:,JJ(1):JJ(2));
    
    testx = testx - testx(1,1);
    
    px = input(KK,2) - testx(1,1);
    
    sl =( testz(1,2) - testz(1,1) )/ ( testx(1,2) - testx(1,1) );
    
    pz = sl*-px + input(KK,3); pz = round(pz*1e9)*1e-9;
    
    if all(testz(:,1) < pz)
        pz = max(testz(:,1));
    end
    
    testz(:,1) = round(testz(:,1)*1e9)*1e-9;
    
    II = find(testz(:,1) >= pz,1,'last');
    
    if II > Nz
        II = Nz;
    end
    
    try
        input(KK,5) = sub2ind(size(vv),II,JJ(1));
    catch
        continue
    end
    
    
end

input = input;

if makeplot;
    
    min_vel = min(min_vel);
    
    % generate colormap using default (jet(64));
    
    % defined velocity bins
    
    cmap = colormap;
    
    ci = linspace(min_vel,max_vel,length(cmap)-1); ci = [ci(:); max_vel];
    
    figure(makeplot); close(makeplot); figure(makeplot);
    set(gcf,'Units','centimeters');
    set(gcf,'Position',[30 20 24 12]);
    axis equal; box on;
    xlim([min(min(xx))-min(abs(unique(diff(xx')))) max(max(xx))+min(abs(unique(diff(xx'))))]);
    ylim([min(min(zz))-abs(median(unique(diff(zz)))) max(max(zz))+abs(median(unique(diff(zz))))]);
    xlabel('X (m)','fontsize',14); ylabel('Z (m)','fontsize',14);
    hold on;
    
    [m,n] = size(xx); m = m-1;
    
    for i = 1:nmodel;
        CI = find(ci>=xzv(i,1).v,1);
        patch(xzv(i,1).X,xzv(i,1).Z,cmap(CI,:),'edgecolor','none');
    end
    
    hold off; drawnow
    caxis([min(ci) max(ci)]);
    colormap(cmap);
    c=colorbar;
    set(get(c,'ylabel'),'String','Velocity (m/s)','Rotation', 270,...
        'VerticalAlignment','Bottom','fontsize',14);
    set(gca,'TickDir','out','linewidth',1,'XMinorTick','on','YMinorTick','on');
    h=findall(gcf,'Type','Axes'); set(h,'FontSize',14);
    
end

end
