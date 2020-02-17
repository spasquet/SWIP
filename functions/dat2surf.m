function [XX,YY,ZZ,DATA_depth,XX_loc,topo_interp]=dat2surf(filename,dx,dz,utm_start,utm_end,zshift)

if exist('zshift','var')==0 || isempty(zshift)==1
    zshift = 0;
end

if exist('utm_start','var')==0 || isempty(utm_start)==1
    data_type = 1;
    % Read file
    if isstr(filename)
        data = load(filename);
    else
       data = filename; 
    end
    X = data(:,1);
    Y = data(:,2);
    Z = data(:,3);
    DATA = data(:,4);
    utm_start = [X(1) Y(1)];
    utm_end = [X(end) Y(end)];

    % Calculate local coordinates and new grid
    X = utm2local(X,Y,utm_start,utm_end);
    if exist('dx','var')==0 || isempty(dx)==1
        x = unique(X)';
    else
        x = min(unique(X)):dx:max(unique(X));
    end
    
    % Get topography
    uniq_x = unique(X);
    topo = zeros(size(uniq_x));
    z_tmp = Z;
    
    for i=1:length(uniq_x)
        topo(i) = max(Z(X==uniq_x(i)));
        z_tmp(X==uniq_x(i)) = z_tmp(X==uniq_x(i))-topo(i);
    end
    if exist('dz','var')==0 || isempty(dx)==1
        depth = flipud(unique(z_tmp))';
    elseif length(dz)==1
        depth = 0:-dz:-(min(topo)-min(Z(:)));
    else
        depth = dz;
    end
    topo_interp = interp1(uniq_x,topo,x);
    
    [XX,ZZ] = meshgrid(x,depth);
    F2 = scatteredInterpolant(X,z_tmp,DATA,'natural','none');
    DATA_depth=F2(XX,ZZ);   
else
    if exist('utm_end','var')==0 || isempty(utm_end)==1
        data_type = 2;
    else
        data_type = 0;
    end
    % Read file
    [DATA,X,Z] = readtomo(filename,1,[],[],100);
    
    if all(isnan(DATA(1,:)))
        newdata = 0;
    else
        newdata = 1;
    end
    
    if newdata == 0
        % Get topography
        [r,c] = find(~isnan(DATA'));
        firstIndex = accumarray(r,c,[size(DATA',1),1],@max,max(c));
        topo = Z(firstIndex,1);
        z_tmp = Z - repmat(topo',size(Z,1),1);        
    else
        topo = Z(1,:);
        for ii = 1:size(DATA,2)
            test = find(~isnan(DATA(:,ii)),1);
            if ~isempty(test)
                topo(ii) = Z(test,ii);
            end
        end
        z_tmp = Z - repmat(topo,size(Z,1),1);
    end
    
    % Flatten and interpolate
    if isempty(dx)==1
        x = unique(X)';
    else
        x = min(unique(X)):dx:max(unique(X));
    end
    if isempty(dz)==1
        depth = flipud(unique(z_tmp))';
    elseif length(dz)==1
        depth = 0:-dz:min(z_tmp(:));
    else
        depth = [0 -dz' -max(dz)-mean(diff(dz)):-mean(diff(dz)):min(z_tmp(:))];
        %         depth = -dz';
    end
    topo_interp = interp1(unique(X),topo,x);

    [XX,ZZ] = meshgrid(x,depth);
    F2 = scatteredInterpolant(X(z_tmp<=0),z_tmp(z_tmp<=0),DATA(z_tmp<=0),'natural','none');
    DATA_depth=F2(XX,ZZ);
end
XX_loc = XX;

if data_type==2
    aa = load(utm_start);
    if size(aa,2)>3
        pt = curvspace([aa(:,2),aa(:,3)],length(x));
    else
        pt = curvspace([aa(:,1),aa(:,2)],length(x));
    end
    X_UTM = pt(:,1);
    Y_UTM = pt(:,2);
    if zshift == 1
        Z_UTM_interp = interp1(aa(:,1),aa(:,4),x);
        shift_z = mean(Z_UTM_interp - topo_interp);% + 3;
    end
else
    % Calculate UTM coordinates and new grid
    [X_UTM, Y_UTM] = local2utm(x,[],utm_start,utm_end);
end

[XX,~] = meshgrid(X_UTM,depth);
[YY,ZZ] = meshgrid(Y_UTM,depth);
ZZ = repmat(topo_interp,length(depth),1)+ZZ;
if zshift == 1
    ZZ = ZZ+shift_z;
end
