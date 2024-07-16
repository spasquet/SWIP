function [seismomat,x,z,t,Sx,Sz]=shotselec(sufile,fldr,xsca,loadmat,keepsu,normalize)
% S. Pasquet - V16.3.28
% Select shot with suwind in matlab

if exist('keepsu','var')==0 || isempty(keepsu)==1
   keepsu = 0; 
end

if exist('normalize','var')==0 || isempty(normalize)==1
   normalize = 0; 
end

seismofile=strcat(num2str(fldr),'.su');
% Select shot then window between gmin and gmax
if normalize == 0
    com1=sprintf('suwind < %s key=fldr min=%d max=%d | suop op=despike nw=3 > %s',...
        sufile,fldr,fldr,seismofile);
else
    com1=sprintf('suwind < %s key=fldr min=%d max=%d | suop op=despike nw=5 | suop op=norm > %s',...
        sufile,fldr,fldr,seismofile); %
end
unix(com1);

% Check number of trace
% Geophones positions
[~,Gx]=unix(['sugethw < ',seismofile,' key=gx output=geom']);
Gx=str2num(Gx)/xsca;
% Geophones elevation
[~,Gz]=unix(['sugethw < ',seismofile,' key=gelev output=geom']);
Gz=str2num(Gz)/xsca;
% Source position
[~,Sx]=unix(['sugethw < ',seismofile,' key=sx output=geom | uniq']);
Sx=str2double(Sx)/xsca;
% Source elevation
[~,Sz]=unix(['sugethw < ',seismofile,' key=selev output=geom | uniq']);
Sz=str2double(Sz)/xsca;

% Save in ASCII .dat file
% Time delay
[~,delay]=unix(['sugethw < ',seismofile,' key=delrt output=geom | uniq']);
delay=str2double(delay)/1000;
% Time sampling
[~,dt]=unix(['sugethw < ',seismofile,' key=dt output=geom | uniq']);
dt=str2double(dt)/1000000;
% Nb of time samples
[~,ns]=unix(['sugethw < ',seismofile,' key=ns output=geom | uniq']);
ns=str2double(ns);

% X table
x=Gx;
% Z table
z=Gz;
% Time table
t=delay:dt:delay+(ns-1)*dt;

if loadmat==1
    com1=sprintf('sustrip < %s | b2a n1=%d > %s',seismofile,ns,[seismofile,'.dat']);
    [~,~]=unix(com1);
    % Load image file
    seismomat=load([seismofile,'.dat']);
    delete([seismofile,'.dat']);
else
    seismomat=[];
end
if keepsu == 0
    delete(seismofile);
end

end