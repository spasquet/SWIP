function [dspmat,f,v]=matpomegal(sufile,coord,nray,fmin,fmax,vmin,vmax,flip,xsca,tsca,norm,dspfile,datsave)

%%% S. Pasquet - V22.05.04
% SUPOMEGAL for matlab
% [dspmat,f,v]=matpomegal(sufile,coord,nray,fmin,fmax,vmin,vmax,flip,xsca,tsca,norm,dspfile,datsave)

wsl = ispc_wsl;

sufile_unix = unix_wsl_path(sufile,wsl);
dspfile_unix = unix_wsl_path(dspfile,wsl);

com1=sprintf(['supomegal < %s coord=%d nray=%d fmin=%d fmax=%d',...
    ' vmin=%d vmax=%d flip=%d xsca=%d tsca=%d > %s'],...
    sufile_unix,coord,nray,fmin,fmax,vmin,vmax,flip,xsca,tsca,dspfile_unix);
unix_cmd(com1,wsl);

if norm>0
    matop(dspfile_unix,'norm',flip);
end

% Save in ASCII .dat file
if  flip==0
    % Nb of frequency samples
    [~,nf]=unix_cmd(['sugethw < ',dspfile_unix,' key=ns output=geom | uniq'],wsl);
    nf=str2double(nf);
    n1=nf;
else
    % Nb of frequency samples
    [~,nf]=unix_cmd(['sugethw < ',dspfile_unix,' key=ntr output=geom | uniq'],wsl);
    nf=str2double(nf);
    n1=nray;
end
if isunix==1 || ispc_wsl==1
    com1=sprintf('sustrip < %s | b2a n1=%d > %s',dspfile_unix,n1,[dspfile_unix,'.dat']);
else
    com1=sprintf('sustrip < %s head=head outpar=outpar | b2a n1=%d outpar=outpar > %s',...
        dspfile_unix,n1,[dspfile_unix,'.dat']); 
end
[~,~]=unix_cmd(com1,wsl);

% Velocity step
dv=(vmax-vmin)/(nray-1);
% Velocity table
v=vmin:dv:vmax;
% Frequency step
df=(fmax-fmin)/(nf-1);
% Frequency table
f=fmin:df:fmax;
% Load image file
dspmat=load([dspfile,'.dat']);

if datsave==0
    delete([dspfile,'.dat']);
end
end