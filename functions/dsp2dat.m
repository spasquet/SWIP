function [dspmat,f,v]=dsp2dat(dspfile,flip,datsave)

%%% S. Pasquet - V22.05.04
% Convert .dsp su file in ASCII file .dat for matlab
% [dspmat,f,v]=dsp2dat(dspfile,flip,datsave)

wsl = ispc_wsl;

dspfile_unix = unix_wsl_path(dspfile,wsl);

if  flip==0
    % Frequency minimum
    [~,fmin]=unix_cmd(['sugethw < ',dspfile_unix,' key=f1 output=geom | uniq'],wsl);
    fmin=str2double(fmin);
    % Velocity minimum
    [~,vmin]=unix_cmd(['sugethw < ',dspfile_unix,' key=f2 output=geom | uniq'],wsl);
    vmin=str2double(vmin);
    
    % Nb of frequency samples
    [~,nf]=unix_cmd(['sugethw < ',dspfile_unix,' key=ns output=geom | uniq'],wsl);
    nf=str2double(nf);
    % Nb of velocity samples
    [~,nv]=unix_cmd(['sugethw < ',dspfile_unix,' key=ntr output=geom | uniq'],wsl);
    nv=str2double(nv);
    
    % Frequency sampling
    [~,df]=unix_cmd(['sugethw < ',dspfile_unix,' key=d1 output=geom | uniq'],wsl);
    df=str2double(df);
    % Velocity sampling
    [~,dv]=unix_cmd(['sugethw < ',dspfile_unix,' key=d2 output=geom | uniq'],wsl);
    dv=str2double(dv);
    
    n1=nf;
else
    % Frequency minimum
    [~,fmin]=unix_cmd(['sugethw < ',dspfile_unix,' key=f2 output=geom | uniq'],wsl);
    fmin=str2double(fmin);
    % Velocity minimum
    [~,vmin]=unix_cmd(['sugethw < ',dspfile_unix,' key=delrt output=geom | uniq'],wsl);
    vmin=str2double(vmin);
    
    % Nb of frequency samples
    [~,nf]=unix_cmd(['sugethw < ',dspfile_unix,' key=ntr output=geom | uniq'],wsl);
    nf=str2double(nf);
    % Nb of velocity samples
    [~,nv]=unix_cmd(['sugethw < ',dspfile_unix,' key=ns output=geom | uniq'],wsl);
    nv=str2double(nv);
    
     % Frequency sampling
    [~,df]=unix_cmd(['sugethw < ',dspfile_unix,' key=d2 output=geom | uniq'],wsl);
    df=str2double(df);
    % Velocity sampling
    [~,dv]=unix_cmd(['sugethw < ',dspfile_unix,' key=dt output=geom | uniq'],wsl);
    dv=str2double(dv)/1000;
    
    n1=nv;
end
if isunix==1 || ispc_wsl==1
    com1=sprintf('sustrip < %s | b2a n1=%d > %s',dspfile_unix,n1,[dspfile_unix,'.dat']);
else
    com1=sprintf('sustrip < %s head=head outpar=outpar | b2a n1=%d outpar=outpar > %s',...
        dspfile_unix,n1,[dspfile_unix,'.dat']);
end
[~,~]=unix_cmd(com1,wsl);

% Velocity table
v=vmin:dv:vmin+dv*(nv-1);
% Frequency table
f=fmin:df:fmin+df*(nf-1);

dspmat=load([dspfile,'.dat']);
if datsave~=1
    delete([dspfile,'.dat']);
end
if ~isunix && ~ispc_wsl
    delete('head','outpar');
end
end