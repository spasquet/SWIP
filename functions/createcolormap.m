function [col,youp,ind]=createcolormap(data,map,colorlogscale,data_min,data_max)

%%% S. Pasquet - V17.02.17
% Create colormap
% [col,youp,ind]=createcolormap(data,map,colorlogscale)

if exist('map','var')==0 || isempty(map)==1
    map=colormap;
end
if exist('colorlogscale','var')==0 || isempty(colorlogscale)==1
    colorlogscale=0;
end
if exist('data_min','var')==0 || isempty(data_min)==1 || exist('data_max','var')==0 || isempty(data_max)==1
    diverg = 0;
else
    diverg = 1;
end
% Creation of colormaps
if isempty(data)==1
    col=[];
    youp=[];
    ind=[];
elseif length(data)==1
    col=map(1,:);
    ind=1;
else
    if diverg == 0
        data_min = min(data);
        data_max = max(data);
    end
    if colorlogscale==0
        youp=linspace(data_min,data_max,size(map,1));
    else
        youp=logspace(log10(data_min),log10(data_max),size(map,1));
    end
    youp(1)=data_min;
    youp(end)=data_max;
    col=ones(length(data),3);
    ind=ones(length(data),1);
    for e=1:size(map,1)-1;
        col(data>=youp(e) & data<youp(e+1),:)=...
            repmat(map(e,:),...
            size(col(data>=youp(e) & data<youp(e+1),:),1),1);
        ind(data>=youp(e) & data<youp(e+1))=e;
    end
    col(data==youp(e+1),:)=repmat(map(e+1,:),...
        size(col(data==youp(e+1),:),1),1);
    ind(data==youp(e+1))=e;
end