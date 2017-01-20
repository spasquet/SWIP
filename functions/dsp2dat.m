function [dspmat,f,v]=dsp2dat(dspfile,flip,datsave)

%%% S. Pasquet - V17.01.13
% Convert .dsp su file in ASCII file .dat for matlab
% [dspmat,f,v]=dsp2dat(dspfile,flip,datsave)

if  flip==0
    % Frequency minimum
    [~,fmin]=unix(['sugethw < ',dspfile,' key=f1 output=geom | uniq']);
    fmin=str2double(fmin);
    % Velocity minimum
    [~,vmin]=unix(['sugethw < ',dspfile,' key=f2 output=geom | uniq']);
    vmin=str2double(vmin);
    
    % Nb of frequency samples
    [~,nf]=unix(['sugethw < ',dspfile,' key=ns output=geom | uniq']);
    nf=str2double(nf);
    % Nb of velocity samples
    [~,nv]=unix(['sugethw < ',dspfile,' key=ntr output=geom | uniq']);
    nv=str2double(nv);
    
    % Frequency sampling
    [~,df]=unix(['sugethw < ',dspfile,' key=d1 output=geom | uniq']);
    df=str2double(df);
    % Velocity sampling
    [~,dv]=unix(['sugethw < ',dspfile,' key=d2 output=geom | uniq']);
    dv=str2double(dv);
    
    n1=nf;
else
    % Frequency minimum
    [~,fmin]=unix(['sugethw < ',dspfile,' key=f2 output=geom | uniq']);
    fmin=str2double(fmin);
    % Velocity minimum
    [~,vmin]=unix(['sugethw < ',dspfile,' key=delrt output=geom | uniq']);
    vmin=str2double(vmin);
    
    % Nb of frequency samples
    [~,nf]=unix(['sugethw < ',dspfile,' key=ntr output=geom | uniq']);
    nf=str2double(nf);
    % Nb of velocity samples
    [~,nv]=unix(['sugethw < ',dspfile,' key=ns output=geom | uniq']);
    nv=str2double(nv);
    
     % Frequency sampling
    [~,df]=unix(['sugethw < ',dspfile,' key=d2 output=geom | uniq']);
    df=str2double(df);
    % Velocity sampling
    [~,dv]=unix(['sugethw < ',dspfile,' key=dt output=geom | uniq']);
    dv=str2double(dv)/1000;
    
    n1=nv;
end
if isunix==1
    com1=sprintf('sustrip < %s | b2a n1=%d > %s',dspfile,n1,[dspfile,'.dat']);
else
    com1=sprintf('sustrip < %s head=head outpar=outpar | b2a n1=%d outpar=outpar > %s',...
        dspfile,n1,[dspfile,'.dat']);
end
[~,~]=unix(com1);

% Velocity table
v=vmin:dv:vmin+dv*(nv-1);
% Frequency table
f=fmin:df:fmin+df*(nf-1);

dspmat=load([dspfile,'.dat']);
if datsave~=1
    delete([dspfile,'.dat']);
end
if isunix==0
    delete('head','outpar');
end
end