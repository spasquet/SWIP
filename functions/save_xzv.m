function [XX,ZZ,VV,MASK1,MASK2]=save_xzv(filexzv,xi,zi,vi,mask1,mask2,convert)

%%% S. Pasquet - V16.12.02
% Save velocity model and mask in .xzv column ascii file 
% [XX,ZZ,VV,MASK1,MASK2]=save_xzv(filexzv,xi,zi,vi,mask1,mask2,convert)

if exist('mask1','var')==0 || isempty(mask1)==1
    mask1=[];
end
if exist('mask2','var')==0 || isempty(mask2)==1
    mask2=[];
end
if exist('convert','var')==0 || isempty(convert)==1
   convert=0;
end


% Save velocity model for Vp/Vs calculation
sizevi=size(vi);
if min(size(xi))>1
    XX=reshape(xi,size(xi,1)*size(xi,2),1);
else
    if size(xi,1)==1
        xi=xi';
    end
    sizeok=sizevi(sizevi~=length(xi));
    XX=sort(repmat(xi,sizeok,1));
end
if min(size(zi))>1
    ZZ=reshape(zi,size(zi,1)*size(zi,2),1);
else
    if size(zi,1)==1
        zi=zi';
    end
    sizeok=sizevi(sizevi~=length(zi));
    ZZ=repmat(zi,sizeok,1);
end
VV=reshape(vi,size(vi,1)*size(vi,2),1);
MASK1=reshape(mask1,size(mask1,1)*size(mask1,2),1);
MASK2=reshape(mask2,size(mask2,1)*size(mask2,2),1);

XX(isnan(VV)==1)=[];
ZZ(isnan(VV)==1)=[];
if isempty(mask1)==0
    MASK1(isnan(VV)==1)=[];
end
if isempty(mask2)==0
    MASK2(isnan(VV)==1)=[];
end
VV(isnan(VV)==1)=[];
if convert==1;
    VV=VV*1000;
end

dlmwrite(filexzv,[XX,ZZ,VV,MASK1,MASK2],'delimiter','\t','precision','%6.2f');
