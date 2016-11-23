function [specmat,f,x]=matspecfx(seismofile,xsca,specfile,datsave,norm)

%%% S. Pasquet - V16.11.18
% SUSPECFX for matlab
% [specmat,f,x]=matspecfx(seismofile,xsca,specfile,datsave,norm)

if exist('norm','var')==0 || isempty(norm)==1
    norm=1;
end

if norm==1
    com1=sprintf('suspecfx < %s | suop op=norm > %s',...
        seismofile,specfile);
else
    com1=sprintf('suspecfx < %s > %s',...
        seismofile,specfile);
end
unix(com1);

% Save in ASCII .dat file
% Nb of frequency samples
[~,nf]=unix(['sugethw < ',specfile,' ns output=geom | uniq']);
nf=str2double(nf);
% Frequency sampling
[~,df]=unix(['sugethw < ',specfile,' d1 output=geom | uniq']);
df=str2double(df);
% Geophones positions
[~,Gx]=unix(['sugethw < ',specfile,' key=gx output=geom']);
Gx=str2num(Gx)/xsca;

if isunix==1
    com1=sprintf('sustrip < %s | b2a n1=%d > %s',specfile,nf,[specfile,'.dat']);
else
    com1=sprintf('sustrip < %s head=head outpar=outpar | b2a n1=%d outpar=outpar > %s',...
        specfile,nf,[specfile,'.dat']);
end
[~,~]=unix(com1);

% X table
x=unique(Gx);
% Time table
f=0:df:(nf-1)*df;
% Load image file
specmat=load([specfile,'.dat']);

if datsave==0
    delete([specfile,'.dat']);
end
end