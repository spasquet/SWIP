function [seismomat,x,t,ntr]=matwind(sufile,sx,gmin,gmax,xsca,winsize,...
    seismofile,datsave,mute,tmin1,tmin2,tmax1,tmax2)

% S. Pasquet - V17.01.13
% SUWIND for matlab
% [seismomat,x,t,ntr]=matwind(sufile,sx,gmin,gmax,xsca,winsize,seismofile,datsave,...
% mute,tmin1,tmin2,tmax1,tmax2)

supath=fileparts(seismofile);
% Select shot then window between gmin and gmax
com1=sprintf(['suwind < %s key=sx min=%d max=%d |',...
    ' suwind key=gx min=%d max=%d > %s'],...
    sufile,round(sx*xsca),round(sx*xsca),...
    round(gmin*xsca),round(gmax*xsca),fullfile(supath,'tmp.su'));
unix(com1);

% Get fldr in case two shots at the same position
com1=sprintf('sugethw < %s key=fldr output=geom | uniq',fullfile(supath,'tmp.su'));
[~,fldr]=unix(com1);
fldr=str2num(fldr);

GxOK=[]; flag=0;
for i=1:length(fldr)
    [~,Gx]=unix(sprintf('suwind < %s key=fldr min=%d max=%d | sugethw key=gx output=geom',fullfile(supath,'tmp.su'),fldr(i),fldr(i)));
    Gx=str2num(Gx)/xsca;
    ntr=length(Gx);
    if ntr==winsize
        com1=sprintf('suwind < %s key=fldr min=%d max=%d > %s',fullfile(supath,'tmp.su'),fldr(i),fldr(i),fullfile(supath,'tmp2.su'));
        unix(com1);
        flag=1; break
    elseif ntr<winsize
        if exist(fullfile(supath,'tmp2.su'),'file')~=2
            com1=sprintf('suwind < %s key=fldr min=%d max=%d > %s',fullfile(supath,'tmp.su'),fldr(i),fldr(i),fullfile(supath,'tmp2.su'));
            unix(com1);
            GxOK=[GxOK;Gx];
        else
            GxKO=Gx(ismember(Gx,GxOK)==0);
            GxOK=[GxOK;Gx];
            for j=1:length(GxKO)
                com1=sprintf('suwind < %s key=fldr min=%d max=%d | suwind key=gx min=%d max=%d > %s',...
                    fullfile(supath,'tmp.su'),fldr(i),fldr(i),round(GxKO(j)*xsca),round(GxKO(j)*xsca),fullfile(supath,'tmp4.su'));
                unix(com1);
                movefile(fullfile(supath,'tmp2.su'),fullfile(supath,'tmp3.su'));
                unix(sprintf('cat %s %s > %s',fullfile(supath,'tmp3.su'),fullfile(supath,'tmp4.su'),fullfile(supath,'tmp2.su')));
            end
            [~,Gx]=unix(sprintf('sugethw < %s key=gx output=geom',fullfile(supath,'tmp2.su')));
            Gx=str2num(Gx)/xsca;
            ntr=length(Gx);
            if ntr==winsize
                unix(sprintf('susort < %s gx > %s',fullfile(supath,'tmp2.su'),fullfile(supath,'tmp3.su')));
                movefile(fullfile(supath,'tmp3.su'),fullfile(supath,'tmp2.su'));
                flag=1; break
            end
        end
    end
end

if isempty(fldr)==1
    seismomat=[]; x=[]; t=[]; ntr=0;
    return
end

if flag==0
    com1=sprintf('suwind < %s key=fldr min=%d max=%d > %s',fullfile(supath,'tmp.su'),fldr(1),fldr(1),fullfile(supath,'tmp2.su'));
    unix(com1);
end

% Mute if enabled
if mute==1
    if sx<=gmin
        com1=sprintf(['sumute < %s key=gx xmute=%d,%d tmute=%d,%d mode=0 ntaper=50 |',...
            ' sumute key=gx xmute=%d,%d tmute=%d,%d mode=1 ntaper=50 > %s'],...
            fullfile(supath,'tmp2.su'),round(gmin*xsca),round(gmax*xsca),tmin1,tmin2,...
            round(gmin*xsca),round(gmax*xsca),tmax1,tmax2,fullfile(supath,'tmp.su'));
    else
        com1=sprintf(['sumute < %s key=gx xmute=%d,%d tmute=%d,%d mode=0 ntaper=50 |',...
            ' sumute key=gx xmute=%d,%d tmute=%d,%d mode=1 ntaper=50 > %s'],...
            fullfile(supath,'tmp2.su'),round(gmin*xsca),round(gmax*xsca),tmin2,tmin1,...
            round(gmin*xsca),round(gmax*xsca),tmax2,tmax1,fullfile(supath,'tmp.su'));
    end
    unix(com1);
    delete(fullfile(supath,'tmp2.su'))
else
    movefile(fullfile(supath,'tmp2.su'),fullfile(supath,'tmp.su'));
end

% Offset
if sx<=gmin
    com1=sprintf('suchw < %s key1=offset key2=gx key3=sx a=0 b=1 c=-1 d=1 e=1 f=1 > %s',fullfile(supath,'tmp.su'),fullfile(supath,'tmp2.su'));
    unix(com1);
    movefile(fullfile(supath,'tmp2.su'),fullfile(supath,'tmp.su'));
else
    com1=sprintf('suchw < %s key1=offset key2=gx key3=sx a=0 b=-1 c=1 d=1 e=1 f=1 > %s',fullfile(supath,'tmp.su'),fullfile(supath,'tmp2.su'));
    unix(com1);
    movefile(fullfile(supath,'tmp2.su'),fullfile(supath,'tmp.su'));
end

% If first shot overlap first trace
if sx==gmin
    com1=sprintf('sushw < %s key=sx a=%d > %s',fullfile(supath,'tmp.su'),xsca*(gmin-median(diff(Gx))/2),fullfile(supath,'tmp2.su'));
    unix(com1);
    movefile(fullfile(supath,'tmp2.su'),fullfile(supath,'tmp.su'));
elseif sx==gmax
    com1=sprintf('sushw < %s key=sx a=%d > %s',fullfile(supath,'tmp.su'),xsca*(gmax+median(diff(Gx))/2),fullfile(supath,'tmp2.su'));
    unix(com1);
    movefile(fullfile(supath,'tmp2.su'),fullfile(supath,'tmp.su'));
end

% Remove continuous part and save in .SU file
com1=sprintf('suop < %s op=avg > %s',fullfile(supath,'tmp.su'),seismofile);
unix(com1);
delete(fullfile(supath,'tmp*.su'));

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

if isunix==1
    com1=sprintf('sustrip < %s | b2a n1=%d > %s',seismofile,ns,[seismofile,'.dat']);
else
     com1=sprintf('sustrip < %s head=head outpar=outpar | b2a n1=%d outpar=outpar > %s',...
        seismofile,ns,[seismofile,'.dat']);
end
[~,~]=unix(com1);

% X table
x=unique(Gx);
% Time table
t=delay:dt:delay+(ns-1)*dt;
% Load image file
seismomat=load([seismofile,'.dat']);

if datsave==0
    delete([seismofile,'.dat']);
end
end