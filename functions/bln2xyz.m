function [ZI,XI,YI]=bln2xyz(datafile,blnfile,xsca)
%%% S. Pasquet - V16.3.14
% Blank X,Y,Z file
% jhjhjh
if nargin<3
    xsca=100;
end

[ZI,XI,YI]=readtomo(datafile,[],[],xsca); % Read data file
if isnumeric(blnfile)==0
    blndat=load(blnfile); % Read bln file
    if blndat(1,1)==length(blndat)-1
        blndat=blndat(2:end,:);
    end
    Xbln=blndat(:,1);
Ybln=blndat(:,2);
mask=inpolygon(XI,YI,Xbln,Ybln);
else
    mask=ZI<blnfile;
end

XI(mask==0)=NaN;
YI(mask==0)=NaN;
ZI(mask==0)=NaN;
