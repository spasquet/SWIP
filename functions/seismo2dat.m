function [seismomat,t,x]=seismo2dat(seismofile,datsave)
%%% S. Pasquet - V16.6.14
% Convert .su seismo file in ASCII file .dat for matlab

% Time minimum
[~,tmin]=unix(['sugethw < ',seismofile,' key=delrt output=geom | uniq']);
tmin=str2double(tmin)/1000;
% Nb of time samples
[~,nt]=unix(['sugethw < ',seismofile,' key=ns output=geom | uniq']);
nt=str2double(nt);
% Time sampling
[~,dt]=unix(['sugethw < ',seismofile,' key=dt output=geom | uniq']);
dt=str2double(dt)/1e6;
% Frequency table
t=tmin:dt:tmin+dt*(nt-1);

% % X scaling
% [~,xsca]=unix(['sugethw < ',seismofile,' key=scalco output=geom | uniq']);
% xsca=abs(str2double(xsca));

% Gx table
[~,x]=unix(['sugethw < ',seismofile,' key=gx output=geom | uniq']);
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
[~,~]=unix(com1);

seismomat=load([seismofile,'.dat']);
if datsave~=1
    delete([seismofile,'.dat']);
end
if isunix==0
    unix('rm -rf head outpar');
end
end