function [dspmat,f,v]=matpomegal(sufile,coord,nray,fmin,fmax,vmin,vmax,...
    flip,xsca,tsca,norm,dspfile,datsave)

%%% S. Pasquet - V16.5.6
% SUPOMEGAL for matlab

com1=sprintf(['supomegal < %s coord=%d nray=%d fmin=%d fmax=%d',...
    ' vmin=%d vmax=%d flip=%d xsca=%d tsca=%d > %s'],...
    sufile,coord,nray,fmin,fmax,vmin,vmax,flip,xsca,tsca,dspfile);
unix(com1);

if norm==1
    matop(dspfile,'norm',flip);
end

% Save in ASCII .dat file
if  flip==0
    % Nb of frequency samples
    [~,nf]=unix(['sugethw < ',dspfile,' key=ns output=geom | uniq']);
    nf=str2double(nf);
    n1=nf;
else
    % Nb of frequency samples
    [~,nf]=unix(['sugethw < ',dspfile,' key=ntr output=geom | uniq']);
    nf=str2double(nf);
    n1=nray;
end
if isunix==1
    com1=sprintf('sustrip < %s | b2a n1=%d > %s',dspfile,n1,[dspfile,'.dat']);
else
    com1=sprintf('sustrip < %s head=head outpar=outpar | b2a n1=%d outpar=outpar > %s',...
        dspfile,n1,[dspfile,'.dat']); 
end
[~,~]=unix(com1);

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