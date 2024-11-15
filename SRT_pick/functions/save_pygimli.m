function save_pygimli(pickfile_GIMLI,Gx,Gz,Sx,Sz,Gxpick,Gzpick,Sxpick,Szpick,Tpick,err_pc,err_val,err_val_min,err_val_max)

if ~exist('err_pc','var') || isempty(err_pc)
    err_pc = 0;
end

if ~exist('err_val','var') || isempty(err_val)
    if err_pc == 1
        err_val = 0.1;
    else 
        err_val = 1;
    end
end

if ~exist('err_val_min','var') || isempty(err_val_min)
    err_val_min = 0;
end
if ~exist('err_val_max','var') || isempty(err_val_max)
    err_val_max = 1e6;
end

% Convert arrays to table
T = table(Sx, Gx, Sz, Gz);
Tpick = table(Sxpick, Gxpick, Szpick, Gzpick, Tpick);

% Sort by Sx and then by Gx
T = sortrows(T, {'Sx', 'Gx'});

% Sort by Sxpick and then by Gxpick
Tpick = sortrows(Tpick, {'Sxpick', 'Gxpick'});

% Convert back to arrays
Sx = T.Sx; Gx = T.Gx; Sz = T.Sz; Gz = T.Gz;
Sxpick = Tpick.Sxpick; Gxpick = Tpick.Gxpick; Szpick = Tpick.Szpick; Gzpick = Tpick.Gzpick; Tpick = Tpick.Tpick;

% Save pick file for pyGIMLI
fidb=fopen(pickfile_GIMLI,'wt'); % File opening %
array = check_depth_array(Gx,Gz,Sx,Sz);

if ~isempty(array{2}) && ~isempty(array{3}) && ~isempty(array{4})
    all_points_surf = unique([array{1}.G_sing; array{1}.S_sing; array{2}.S_sing; array{4}.G_sing],'rows');
    all_points_depth = unique([array{3}.G_sing; array{3}.S_sing; array{4}.S_sing; array{2}.G_sing],'rows');
    all_points = [all_points_surf; flipud(all_points_depth)];
else
    all_points = unique([[Gx;Sx],[Gz;Sz]],'rows');
end

% err_val = 0.125; % strengbach percentage
% err_val = 1; % Error in ms

npoints = size(all_points,1);
fprintf(fidb,'%d\n',npoints);
fprintf(fidb,'%s\n','# x y z');
for ii = 1:npoints
    fprintf(fidb,'%f\t%f\t%f\n',all_points(ii,1),all_points(ii,2),0);
end
fprintf(fidb,'%d\n',length(Gxpick));
fprintf(fidb,'%s\n','# g s err t');
for jj = 1:length(Gxpick)
    [~,ind_s] = ismember([Sxpick(jj),Szpick(jj)],all_points,'rows');
    [~,ind_g] = ismember([Gxpick(jj),Gzpick(jj)],all_points,'rows');
    
    if err_pc == 1
        err_ok = err_val*(Tpick(jj))/1000;
        if err_ok > err_val_max/1000
            err_ok = err_val_max/1000;
        elseif err_ok < err_val_min/1000
            err_ok = err_val_min/1000;
        end
    else
        err_ok = err_val/1000;
    end
    fprintf(fidb,'%d\t%d\t%f\t%f\n',ind_g,ind_s,err_ok,(Tpick(jj))/1000);
end
fclose(fidb);
