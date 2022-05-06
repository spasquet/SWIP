function [seismomat,t,x]=seismo2dat(seismofile,datsave)

%%% S. Pasquet - V22.05.04
% Convert .su seismo file in ASCII file .dat for matlab
% [seismomat,t,x]=seismo2dat(seismofile,datsave)

wsl = ispc_wsl;

% Time minimum
[~,tmin]=unix_cmd(['sugethw < ',seismofile,' key=delrt output=geom | uniq'],wsl);
tmin=str2double(tmin)/1000;
% Nb of time samples
[~,nt]=unix_cmd(['sugethw < ',seismofile,' key=ns output=geom | uniq'],wsl);
nt=str2double(nt);
% Time sampling
[~,dt]=unix_cmd(['sugethw < ',seismofile,' key=dt output=geom | uniq'],wsl);
dt=str2double(dt)/1e6;
% Frequency table
t=tmin:dt:tmin+dt*(nt-1);

% % X scaling
% [~,xsca]=unix_cmd(['sugethw < ',seismofile,' key=scalco output=geom | uniq']);
% xsca=abs(str2double(xsca));

% Gx table
[~,x]=unix_cmd(['sugethw < ',seismofile,' key=gx output=geom | uniq'],wsl);
pause(1);
% x=str2num(x)/xsca;
x=str2num(x);
x=x';

n1=nt;

if isunix==1
    com1=sprintf('sustrip < %s | b2a n1=%d > %s',seismofile,n1,[seismofile,'.dat']);
else
    com1=sprintf('sustrip < %s head=head outpar=outpar | b2a n1=%d outpar=outpar > %s',...
        seismofile,n1,[seismofile,'.dat']);
end
[~,~]=unix_cmd(com1,wsl);

seismomat=load([seismofile,'.dat']);
if datsave~=1
    delete([seismofile,'.dat']);
end
if isunix==0
    unix_cmd('rm -rf head outpar');
end
end