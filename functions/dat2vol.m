function [XX,YY,ZZ,DATA_depth]=dat2vol(filename,dx,dy,dz,utm_start,utm_end)

aa = load(filename);

X = aa(:,1);
Y = aa(:,2);
Z = aa(:,3);
DATA = aa(:,4);

% plot_scat_log3(2,X,Y,Z,DATA); drawnow;

if exist('utm_start','var')==0 || isempty(utm_start)==1
    [~, seg_len] = arclength(X,Y);
    utm_start = [X(end) Y(end)];
    utm_end = [X(find(seg_len>median(seg_len)*2,1,'last')+1) Y(find(seg_len>median(seg_len)*4,1,'last')+1)];
%     figure;
%     plot(X,Y,'x');
%     hold on;
%     plot([utm_start(1) utm_end(1)],[utm_start(2) utm_end(2)],'r');
   
%    Calculate local coordinates and new grid
    [X, Y] = utm2local(X,Y,utm_start,utm_end);
%     figure;
%     plot(X,Y,'x');
end

x = [min(unique(X)):dx/9:min(unique(X))+dx min(unique(X))+2*dx:dx:max(unique(X))-2*dx max(unique(X))-dx:dx/9:max(unique(X))];
y = [min(unique(Y)):dy/9:min(unique(Y))+dy min(unique(Y))+2*dy:dy:max(unique(Y))-2*dy max(unique(Y))-dy:dy/9:max(unique(Y))];

[uni_XY,~,J] = unique([X Y],'rows');
topo = zeros(size(uni_XY,1),1);
maxdepth = zeros(size(uni_XY,1),1);
for i = 1:length(uni_XY)
    topo(i) = max(Z(ismember(J,i)));
    maxdepth(i) = topo(i) - min(Z(ismember(J,i)));
end
maxdepth = max(maxdepth);

[Xt,Yt] = ndgrid(x,y);
F1 = scatteredInterpolant(uni_XY(:,1),uni_XY(:,2),topo,'natural','linear');
topoI = F1(Xt,Yt);

if length(dz) == 1;
%     z = [min(unique(Z)):dz:max(unique(Z))];
    z = [0:dz:maxdepth];
else
    z = gen_zvec(maxdepth,dz);
    z(end) = -maxdepth;
end
topoI = repmat(topoI,1,1,length(z));

[XX,YY,ZZ] = ndgrid(x,y,z);
ZZ = ZZ + topoI;
F2 = scatteredInterpolant(X,Y,Z,DATA,'natural','nearest');
DATA_depth=F2(XX,YY,ZZ);

if exist('utm_start','var') && ~isempty(utm_start)
    % Calculate UTM coordinates and new grid
    [X_UTM, Y_UTM] = local2utm(XX(:),YY(:),utm_start,utm_end);
    XX = reshape(X_UTM,length(x),length(y),length(z));
    YY = reshape(Y_UTM,length(x),length(y),length(z));
end


