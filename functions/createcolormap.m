function [col,youp]=createcolormap(data,map,colorlogscale)
%%% S. Pasquet - V16.4.14
% Create colormap

if exist('map','var')==0 || isempty(map)==1
    map=colormap;
end
if exist('colorlogscale','var')==0 || isempty(colorlogscale)==1
    colorlogscale=0;
end
% Creation of colormaps
if isempty(data)==1
    col=[];
    youp=[];
elseif length(data)==1
    col=map(1,:);
else
    if colorlogscale==0
        youp=linspace(min(data),max(data),size(map,1));
    else
        youp=logspace(log10(min(data)),log10(max(data)),size(map,1));
    end
    youp(1)=min(data);
    youp(end)=max(data);
    col=ones(length(data),3);
    for e=1:size(map,1)-1;
        col(data>=youp(e) & data<youp(e+1),:)=...
            repmat(map(e,:),...
            size(col(data>=youp(e) & data<youp(e+1),:),1),1);
    end
    col(data==youp(e+1),:)=repmat(map(e+1,:),...
        size(col(data==youp(e+1),:),1),1);
end