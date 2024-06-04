function [specmat,f,x]=matspecfx(seismofile,xsca,specfile,datsave,norm)

%%% S. Pasquet - V22.05.04
% SUSPECFX for matlab
% [specmat,f,x]=matspecfx(seismofile,xsca,specfile,datsave,norm)

wsl = ispc_wsl;

seismofile_unix = unix_wsl_path(seismofile,wsl);
specfile_unix = unix_wsl_path(specfile,wsl);

if exist('norm','var')==0 || isempty(norm)==1
    norm=1;
end

if norm>0
    com1=sprintf('suspecfx < %s | suop op=norm > %s',...
        seismofile_unix,specfile_unix);
else
    com1=sprintf('suspecfx < %s > %s',...
        seismofile_unix,specfile_unix);
end
unix_cmd(com1,wsl);

% Save in ASCII .dat file
% Nb of frequency samples
[~,nf]=unix_cmd(['sugethw < ',specfile_unix,' ns output=geom | uniq'],wsl);
nf=str2double(nf);
% Frequency sampling
[~,df]=unix_cmd(['sugethw < ',specfile_unix,' d1 output=geom | uniq'],wsl);
df=str2double(df);
% Geophones positions
[~,Gx]=unix_cmd(['sugethw < ',specfile_unix,' key=gx output=geom'],wsl);
Gx=str2num(Gx)/xsca;

if isunix==1 || ispc_wsl==1
    com1=sprintf('sustrip < %s | b2a n1=%d > %s',specfile_unix,nf,[specfile_unix,'.dat']);
else
    com1=sprintf('sustrip < %s head=head outpar=outpar | b2a n1=%d outpar=outpar > %s',...
        specfile_unix,nf,[specfile_unix,'.dat']);
end
[~,~]=unix_cmd(com1,wsl);

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
