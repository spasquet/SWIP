function [VelI,XI,ZI]=readtomo(Velfile,mask,X,Z,xsca,vaverage,nW,dx)

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
if mask==1 && size(bb,2)>=4
    bb(:,3)=bb(:,3).*bb(:,4);
elseif mask==2 && size(bb,2)>=5
    bb(:,3)=bb(:,3).*bb(:,5);
elseif mask==2 && size(bb,2)<5 && size(bb,2)>3
    bb(:,3)=bb(:,3).*bb(:,4);
end
bb=sortrows(bb,1);

Xi=unique(bb(:,1));
Zi=unique(bb(:,2));
nxi=length(Xi);
% nzi=length(Zi);
dxv=min(unique(diff(Xi)));
dzv=min(unique(diff(Zi)));
dxz=round(mean([dxv,dzv])*0.5*xsca)/xsca;
% minx=min(Xi);
minz=min(Zi);
% maxx=max(Xi);
maxz=max(Zi);

% Xv=minx:dxv:maxx;
Xv=Xi';
Zv=minz:dxz:maxz;
% XV=repmat(Xv,length(Zv),1);
% ZV=repmat(Zv',1,length(Xv));
% nxv=length(Xv);
nzv=length(Zv);

Vel=zeros(nxi,nzv)*NaN;
for i=1:nxi
    [Ztmp,II]=sort(bb(bb(:,1)==Xi(i),2));
    Veltmp=bb(bb(:,1)==Xi(i),3);
    Veltmp=Veltmp(II);
    if length(find(Veltmp>0))>1
        Vel(i,:)=interp1qr(Ztmp(Veltmp>0),Veltmp(Veltmp>0),Zv');
    end
end
% Vel=zeros(nzv,nxv);
% for i=1:nzv
%     Veltmp=Veli(isnan(Veli(:,i))==0,i);
%     Xtmp=Xi(isnan(Veli(:,i))==0);
%     Vel(i,:)=interp1(Xtmp,Veltmp,Xv,'linear');
% end
% Vel=Vel';

% Xv=reshape(bb(:,1),nz,nx);
% Zv=reshape(bb(:,2),nz,nx);
% dxv=(max(bb(:,1))-min(bb(:,1)))/(nx-1);
% dzv=(max(bb(:,2))-min(bb(:,2)))/(nz-1);
% Vel=reshape(bb(:,3),nz,nx);

if exist('X','var')==1 && isempty(X)==0 && exist('Z','var')==1 && isempty(Z)==0
    [XI,ZI]=meshgrid(X,Z);
    if vaverage==1
        % Average velocity below extraction window
        [XI2,ZI2]=meshgrid(X,Zv);
        VelI2=zeros(size(Vel,2),length(X));
        for i=1:length(X)
            vtemp2=zeros(size(Vel,2),1);
            for j=1:size(Vel,2)
                vtemp=Vel(Xv>=X(i)-dx*(nW-1)*0.5 ...
                & Xv<=X(i)+dx*(nW-1)*0.5,j);
                vtemp=vtemp(~isnan(vtemp));
                vtemp2(j)=mean(vtemp);
            end
            VelI2(:,i)=vtemp2;
        end
        VelI=interp2(XI2,ZI2,VelI2,XI,ZI);
    else
        VelI=interp2(Xv',Zv',Vel',XI,ZI);
    end
else
    VelI=Vel';
    [XI,ZI]=meshgrid(Xi,Zv);
end

end
