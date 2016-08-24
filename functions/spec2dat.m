function [specmat,f,x]=spec2dat(specfile,datsave)
%%% S. Pasquet - V16.6.14
% Convert .spec su file in ASCII file .dat for matlab

% Frequency minimum
fmin=0;
% Nb of frequency samples
[~,nf]=unix(['sugethw < ',specfile,' key=ns output=geom | uniq']);
nf=str2double(nf);
% Frequency sampling
[~,df]=unix(['sugethw < ',specfile,' key=d1 output=geom | uniq']);
df=str2double(df);
% Frequency table
f=fmin:df:fmin+df*(nf-1);

% % X scaling
% [~,xsca]=unix(['sugethw < ',specfile,' key=scalco output=geom | uniq']);
% xsca=abs(str2double(xsca));

% Gx table
[~,x]=unix(['sugethw < ',specfile,' key=gx output=geom | uniq']);
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
[~,~]=unix(com1);

specmat=load([specfile,'.dat']);
if datsave~=1
    delete([specfile,'.dat']);
end
if isunix==0
    unix('rm -rf head outpar');
end
end