function array  = check_depth_array(Gx,Gz,Sx,Sz)


G_sing = unique([Gx Gz],'rows');
S_sing = unique([Sx Sz],'rows');

try
    nG = hist(G_sing(:,1),unique(G_sing(:,1)));
    nS = hist(S_sing(:,1),unique(S_sing(:,1)));
catch
   keyboard
end

G_sing_X = unique(G_sing(:,1));
S_sing_X = unique(S_sing(:,1));

if any(nG>1)
    ind_nG = find(nG>1);
    for i = 1:length(ind_nG)
        ind_GX = (G_sing(:,1) == G_sing_X(ind_nG(i)));
        ind_GZ_depth = (G_sing(:,2) == min(G_sing(ind_GX,2)));
        ind_GZ_surf = (G_sing(:,2) == max(G_sing(ind_GX,2)));
        G_sing_depth(i,:) = G_sing(ind_GX & ind_GZ_depth,:);
        G_sing_surf(i,:) = G_sing(ind_GX & ind_GZ_surf,:);
    end
    
    ind_nG = find(nG==1);
    for i = 1:length(ind_nG)
        Gx_tmp = G_sing_X(ind_nG(i));
        Gz_tmp = G_sing(G_sing(:,1) == Gx_tmp,2);
        diff_GZ_depth = abs(Gz_tmp - G_sing_depth(:,2));
        diff_GZ_surf = abs(Gz_tmp - G_sing_surf(:,2));
        if min(diff_GZ_depth) < min(diff_GZ_surf)
            G_sing_depth = [G_sing_depth; G_sing(G_sing(:,1) == Gx_tmp,:)];
        else
            G_sing_surf = [G_sing_surf; G_sing(G_sing(:,1) == Gx_tmp,:)];
        end
    end
    G_sing_depth = sortrows(G_sing_depth,1);
    G_sing_surf = sortrows(G_sing_surf,1);
else
    G_sing_depth = [];
    %     IS_depth = [];
    G_sing_surf = G_sing;
    %     IG_surf = ismember([Gx Gz],G_sing_surf,'rows');
end

if any(nS>1)
    ind_nS = find(nS>1);
    S_sing_depth = []; S_sing_surf = [];
    for i = 1:length(ind_nS)
        ind_SX = (S_sing(:,1) == S_sing_X(ind_nS(i)));
        ind_SZ_depth = (S_sing(:,2) == min(S_sing(ind_SX,2)));
        ind_SZ_surf = (S_sing(:,2) == max(S_sing(ind_SX,2)));
        S_sing_depth(i,:) = S_sing(ind_SX & ind_SZ_depth,:);
        S_sing_surf(i,:) = S_sing(ind_SX & ind_SZ_surf,:);
    end
    
    ind_nS = find(nS==1);
    for i = 1:length(ind_nS)
        Sx_tmp = S_sing_X(ind_nS(i));
        Sz_tmp = S_sing(S_sing(:,1) == Sx_tmp,2);
        diff_SZ_depth = abs(Sz_tmp - G_sing_depth(:,2));
        diff_SZ_surf = abs(Sz_tmp - G_sing_surf(:,2));
        if min(diff_SZ_depth) < min(diff_SZ_surf)
            S_sing_depth = [S_sing_depth; S_sing(S_sing(:,1) == Sx_tmp,:)];
        else
            S_sing_surf = [S_sing_surf; S_sing(S_sing(:,1) == Sx_tmp,:)];
        end
    end
    S_sing_depth = sortrows(S_sing_depth,1);
    S_sing_surf = sortrows(S_sing_surf,1);
else
    S_sing_depth = [];
    %     IS_depth = [];
    S_sing_surf = S_sing;
    %     IG_surf = ismember([Gx Gz],G_sing_surf,'rows');
end

for i = 1:4
    if i == 1
        if ~isempty(G_sing_surf) && ~isempty(S_sing_surf)
            array{i}.G_sing = G_sing_surf;
            array{i}.S_sing = S_sing_surf;
        else
            array{i} = [];
        end
    elseif i == 2
        if ~isempty(G_sing_depth) && ~isempty(S_sing_surf)
            array{i}.G_sing = G_sing_depth;
            array{i}.S_sing = S_sing_surf;
        else
            array{i} = [];
        end
    elseif i == 3
        if ~isempty(G_sing_depth) && ~isempty(S_sing_depth)
            array{i}.G_sing = G_sing_depth;
            array{i}.S_sing = S_sing_depth;
        else
            array{i} = [];
        end
        
    elseif i == 4
        if ~isempty(G_sing_surf) && ~isempty(S_sing_depth)
            array{i}.G_sing = G_sing_surf;
            array{i}.S_sing = S_sing_depth;
        else
            array{i} = [];
        end
    end
end