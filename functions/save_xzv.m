function [XX,ZZ,VV,MASK]=save_xzv(filexzv,xi,zi,vi,clearnan,mask,convert,delimiter)

%%% S. Pasquet - V17.02.13
% Save velocity model and mask in .xzv column ascii file 
% [XX,ZZ,VV,MASK]=save_xzv(filexzv,xi,zi,vi,clearnan,mask,convert)

if exist('mask','var')==0 || isempty(mask)==1
    mask=[];
end
if exist('clearnan','var')==0 || isempty(clearnan)==1
   clearnan=1;
end
if exist('convert','var')==0 || isempty(convert)==1
   convert=0;
end
if exist('delimiter','var')==0 || isempty(delimiter)==1
   delimiter='\t';
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
MASK=reshape(mask,size(mask,1)*size(mask,2),1);

if clearnan==1
    XX(isnan(VV)==1)=[];
    ZZ(isnan(VV)==1)=[];
end
if isempty(mask)==0 && clearnan==1
    MASK(isnan(VV)==1)=[];
end
if clearnan==1
    VV(isnan(VV)==1)=[];
end
if convert==1
    VV=VV*1000;
end

dlmwrite(filexzv,[XX,ZZ,VV,MASK],'delimiter',delimiter,'precision','%6.4f');
