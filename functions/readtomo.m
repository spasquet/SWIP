function [VelI,XI,ZI,VelI_std]=readtomo(Velfile,mask,X,Z,xsca,vaverage,nW,dx)

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

VelI_std = [];

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
nzi=length(Zi);

test_sizeX = length(unique(hist(bb(:,1),Xi)));
test_sizeZ = length(unique(hist(bb(:,2),Zi)));

if test_sizeX == 1 && (~exist('X','var') || isempty(X) || ~exist('Z','var') || isempty(Z))
    XI = reshape(bb(:,1),length(bb(:,1))/nxi,nxi);
    ZI = reshape(bb(:,2),length(bb(:,2))/nxi,nxi);
    VelI = reshape(bb(:,3),length(bb(:,3))/nxi,nxi);
else
    dxv=min(unique(diff(Xi)));
    dzv=min(unique(diff(Zi)));
    dxz=round(mean([dxv,dzv])*0.5*xsca)/xsca;
    minz=min(Zi);
    maxz=max(Zi);
    
    Xv=Xi';
    Zv=minz:dxz:maxz;
    nzv=length(Zv);
    
    Vel=zeros(nxi,nzv)*NaN;
    topo=zeros(nxi,1)*NaN; 
    
    for i=1:nxi
        [Ztmp,II]=sort(bb(bb(:,1)==Xi(i),2));
        Veltmp=bb(bb(:,1)==Xi(i),3);
        Veltmp=Veltmp(II);
        if length(find(Veltmp>0))>1
            Vel(i,:)=interp1qr(Ztmp(Veltmp>0),Veltmp(Veltmp>0),Zv');
        end
        ind_topo = find(~isnan(Vel(i,:)),1,'last');
        if ~isempty(ind_topo)
            topo(i) = Zv(ind_topo);
        end
    end
    DEPTH = bsxfun(@minus,topo,Zv);
    
end

if exist('X','var')==1 && isempty(X)==0 && exist('Z','var')==1 && isempty(Z)==0
    if min(size(X))==1
        [XI,ZI]=meshgrid(X,Z);
    else
        XI = X; ZI = Z;
    end
    if vaverage==1
        % Average velocity below extraction window
        [XI2,ZI2]=meshgrid(Xv,Zv);
        VelI2 = zeros(size(XI2))*NaN;
        VelI2_std = VelI2;
        
        for i = 1:size(Vel,1)
            vtemp2 = zeros(size(Vel,2),1)*NaN;
            vtemp2_std = vtemp2;
            if length(unique(nW)) == 2
                max_depth_win = max(nW).*10.^(1./sqrt(0.5*max(nW)));
                nW_tmp = interp1([0 max_depth_win],nW,DEPTH(i,:),'linear',max(nW));
                nW_tmp(find(nW_tmp==min(nW))+1:end) = NaN;
%                 nW_tmp = interp1([0 max(DEPTH(:))],nW,DEPTH(i,:));
                nW_tmp = flipud(nW_tmp);
            else
                nW_tmp = unique(nW)*ones(size(Vel,2));
            end
            for j = 1:size(Vel,2)
                if DEPTH(i,j) < 0
                    continue
                end
                vtemp = Vel(Xv>=Xv(i)-dx*(nW_tmp(j))*0.5 & Xv<=Xv(i)+dx*(nW_tmp(j))*0.5,DEPTH(i,:) == DEPTH(i,j));
                vtemp = vtemp(~isnan(vtemp));
                vtemp2(j) = mean(vtemp);
                vtemp2_std(j) = std(vtemp);
            end
            VelI2(:,i)= vtemp2;
            VelI2_std(:,i)= vtemp2_std;
        end
        VelI = interp2(XI2,ZI2,VelI2,XI,ZI);
        VelI_std = interp2(XI2,ZI2,VelI2_std,XI,ZI);
%         plot_img([],XI2,ZI2,Vel');
%         plot_img([],XI2,ZI2,VelI2);
%         plot_img([],XI,ZI,VelI);
    else
        VelI = interp2(Xv',Zv',Vel',XI,ZI);
        VelI_std = 0.1*VelI;
    end
else
    if test_sizeX ~=1
        VelI = Vel';
        VelI_std = 0.1*VelI;
        [XI,ZI] = meshgrid(Xi,Zv);
    end
end