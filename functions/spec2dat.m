function [specmat,f,x]=spec2dat(specfile,datsave)

%%% S. Pasquet - V22.05.04
% Convert .spec su file in ASCII file .dat for matlab
% [specmat,f,x]=spec2dat(specfile,datsave)

wsl = ispc_wsl;

% Frequency minimum
fmin=0;
% Nb of frequency samples
[~,nf]=unix_cmd(['sugethw < ',specfile,' key=ns output=geom | uniq'],wsl);
nf=str2double(nf);
% Frequency sampling
[~,df]=unix_cmd(['sugethw < ',specfile,' key=d1 output=geom | uniq'],wsl);
df=str2double(df);
% Frequency table
f=fmin:df:fmin+df*(nf-1);

% % X scaling
% [~,xsca]=unix_cmd(['sugethw < ',specfile,' key=scalco output=geom | uniq']);
% xsca=abs(str2double(xsca));

% Gx table
[~,x]=unix_cmd(['sugethw < ',specfile,' key=gx output=geom | uniq'],wsl);
% x=str2num(x)/xsca;
x=str2num(x);
x=x';

n1=nf;

if isunix==1
    com1=sprintf('sustrip < %s | b2a n1=%d > %s',specfile,n1,[specfile,'.dat']);
else
    com1=sprintf('sustrip < %s head=head outpar=outpar | b2a n1=%d outpar=outpar > %s',...
        specfile,n1,[specfile,'.dat']);
end
[~,~]=unix_cmd(com1,wsl);

specmat=load([specfile,'.dat']);
if datsave~=1
    delete([specfile,'.dat']);
end
if isunix==0
    unix_cmd('rm -rf head outpar');
end
end