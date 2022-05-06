function [VI,XI,ZI]=readtomo_new(Velfile,mask,X,Z,xsca,vaverage,nW,dx)

%%% S. Pasquet - V16.11.18
% Read tomography velocity file with x,z,v columns and resample it along
% input X and Z (with window averaging or not)
% [VelI,XI,ZI]=readtomo(Velfile,mask,X,Z,xsca,vaverage,nW,dx)

if exist('mask','var')==0 || isempty(mask)==1
    mask=0;
end
if exist('xsca','var')==0 || isempty(xsca)==1
    xsca=100;
end
if exist('vaverage','var')==0 || isempty(vaverage)==1
    vaverage=0; nW=0; dx=0;
end

% Read refraction velocity model
bb=load(Velfile);
bb(:,1)=round(bb(:,1)*xsca)/xsca;
bb(:,2)=round(bb(:,2)*xsca)/xsca;
if mask>=4 && size(bb,2)>=4
    bb(:,3)=bb(:,3).*bb(:,mask);
end
x = bb(:,1);
z = bb(:,2);
v = bb(:,3);

if exist('X','var')==0 || isempty(X)==1
    Xi=unique(bb(:,1));
    Zi=unique(bb(:,2));
    dxv=min(unique(diff(Xi)));
    dzv=min(unique(diff(Zi)));
    minx=min(Xi);
    minz=min(Zi);
    maxx=max(Xi);
    maxz=max(Zi);
    X=minx:dxv:maxx;
    Z=minz:dzv:maxz;
end
[XI,ZI] = meshgrid(X,Z);
F = scatteredInterpolant(x,z,v,'natural');
VI=F(XI,ZI);

end
