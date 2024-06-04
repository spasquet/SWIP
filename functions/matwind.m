function [seismomat,x,t,ntr]=matwind(sufile,sx,gmin,gmax,xsca,winsize,...
    seismofile,datsave,mute,tmin1,tmin2,tmax1,tmax2)

% S. Pasquet - V22.05.04
% SUWIND for matlab
% [seismomat,x,t,ntr]=matwind(sufile,sx,gmin,gmax,xsca,winsize,seismofile,datsave,...
% mute,tmin1,tmin2,tmax1,tmax2)

wsl = ispc_wsl;

supath=fileparts(seismofile);
supath_unix = unix_wsl_path(supath,wsl);
seismofile_unix = unix_wsl_path(seismofile,wsl);

% Select shot then window between gmin and gmax
com1=sprintf(['suwind < %s key=sx min=%d max=%d |',...
    ' suwind key=gx min=%d max=%d > %s'],...
    sufile,round(sx*xsca),round(sx*xsca),...
    round(gmin*xsca),round(gmax*xsca),strcat(supath_unix,'/','tmp.su'));
unix_cmd(com1,wsl);

% Get fldr in case two shots at the same position
com1=sprintf('sugethw < %s key=fldr output=geom | uniq',strcat(supath_unix,'/','tmp.su'));
[~,fldr]=unix_cmd(com1,wsl);
fldr=str2num(fldr);

GxOK=[]; flag=0;
for i=1:length(fldr)
    [~,Gx]=unix_cmd(sprintf('suwind < %s key=fldr min=%d max=%d | sugethw key=gx output=geom',strcat(supath_unix,'/','tmp.su'),fldr(i),fldr(i)),wsl);
    Gx=str2num(Gx)/xsca;
    ntr=length(Gx);
    if ntr==winsize
        com1=sprintf('suwind < %s key=fldr min=%d max=%d > %s',strcat(supath_unix,'/','tmp.su'),fldr(i),fldr(i),strcat(supath_unix,'/','tmp2.su'));
        unix_cmd(com1,wsl);
        flag=1; break
    elseif ntr<winsize
        if exist(strcat(supath_unix,'/','tmp2.su'),'file')~=2
            com1=sprintf('suwind < %s key=fldr min=%d max=%d > %s',strcat(supath_unix,'/','tmp.su'),fldr(i),fldr(i),strcat(supath_unix,'/','tmp2.su'));
            unix_cmd(com1,wsl);
            GxOK=[GxOK;Gx];
        else
            GxKO=Gx(ismember(Gx,GxOK)==0);
            GxOK=[GxOK;Gx];
            for j=1:length(GxKO)
                com1=sprintf('suwind < %s key=fldr min=%d max=%d | suwind key=gx min=%d max=%d > %s',...
                    strcat(supath_unix,'/','tmp.su'),fldr(i),fldr(i),round(GxKO(j)*xsca),round(GxKO(j)*xsca),strcat(supath_unix,'/','tmp4.su'));
                unix_cmd(com1,wsl);
                movefile(fullfile(supath,'tmp2.su'),fullfile(supath,'tmp3.su'));
                unix_cmd(sprintf('cat %s %s > %s',strcat(supath_unix,'/','tmp3.su'),strcat(supath_unix,'/','tmp4.su'),strcat(supath_unix,'/','tmp2.su')),wsl);
            end
            [~,Gx]=unix_cmd(sprintf('sugethw < %s key=gx output=geom',strcat(supath_unix,'/','tmp2.su')),wsl);
            Gx=str2num(Gx)/xsca;
            ntr=length(Gx);
            if ntr==winsize
                unix_cmd(sprintf('susort < %s gx > %s',strcat(supath_unix,'/','tmp2.su'),strcat(supath_unix,'/','tmp3.su')),wsl);
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
    com1=sprintf('suwind < %s key=fldr min=%d max=%d > %s',strcat(supath_unix,'/','tmp.su'),fldr(1),fldr(1),strcat(supath_unix,'/','tmp2.su'));
    unix_cmd(com1,wsl);
end

% Mute if enabled
if mute==1
    if sx<=gmin
        com1=sprintf(['sumute < %s key=gx xmute=%d,%d tmute=%d,%d mode=0 ntaper=50 |',...
            ' sumute key=gx xmute=%d,%d tmute=%d,%d mode=1 ntaper=50 > %s'],...
            strcat(supath_unix,'/','tmp2.su'),round(gmin*xsca),round(gmax*xsca),tmin1,tmin2,...
            round(gmin*xsca),round(gmax*xsca),tmax1,tmax2,strcat(supath_unix,'/','tmp.su'));
    else
        com1=sprintf(['sumute < %s key=gx xmute=%d,%d tmute=%d,%d mode=0 ntaper=50 |',...
            ' sumute key=gx xmute=%d,%d tmute=%d,%d mode=1 ntaper=50 > %s'],...
            strcat(supath_unix,'/','tmp2.su'),round(gmin*xsca),round(gmax*xsca),tmin2,tmin1,...
            round(gmin*xsca),round(gmax*xsca),tmax2,tmax1,strcat(supath_unix,'/','tmp.su'));
    end
    unix_cmd(com1,wsl);
    delete(strcat(supath_unix,'/','tmp2.su'))
else
    movefile(fullfile(supath,'tmp2.su'),fullfile(supath,'tmp.su'));
end

% Offset
if sx<=gmin
    com1=sprintf('suchw < %s key1=offset key2=gx key3=sx a=0 b=1 c=-1 d=1 e=1 f=1 > %s',strcat(supath_unix,'/','tmp.su'),strcat(supath_unix,'/','tmp2.su'));
    unix_cmd(com1,wsl);
    movefile(fullfile(supath,'tmp2.su'),fullfile(supath,'tmp.su'));
else
    com1=sprintf('suchw < %s key1=offset key2=gx key3=sx a=0 b=-1 c=1 d=1 e=1 f=1 > %s',strcat(supath_unix,'/','tmp.su'),strcat(supath_unix,'/','tmp2.su'));
    unix_cmd(com1,wsl);
    movefile(fullfile(supath,'tmp2.su'),fullfile(supath,'tmp.su'));
end

% If first shot overlap first trace
if sx==gmin
    com1=sprintf('sushw < %s key=sx a=%d > %s',strcat(supath_unix,'/','tmp.su'),xsca*(gmin-median(diff(Gx))/2),strcat(supath_unix,'/','tmp2.su'));
    unix_cmd(com1,wsl);
    movefile(fullfile(supath,'tmp2.su'),fullfile(supath,'tmp.su'));
elseif sx==gmax
    com1=sprintf('sushw < %s key=sx a=%d > %s',strcat(supath_unix,'/','tmp.su'),xsca*(gmax+median(diff(Gx))/2),strcat(supath_unix,'/','tmp2.su'));
    unix_cmd(com1,wsl);
    movefile(fullfile(supath,'tmp2.su'),fullfile(supath,'tmp.su'));
end

% Remove continuous part and save in .SU file
com1=sprintf('suop < %s op=avg > %s',strcat(supath_unix,'/','tmp.su'),seismofile_unix);
unix_cmd(com1,wsl);
delete(fullfile(supath,'tmp*.su'));

% Save in ASCII .dat file
% Time delay
[~,delay]=unix_cmd(['sugethw < ',seismofile_unix,' key=delrt output=geom | uniq'],wsl);
delay=str2double(delay)/1000;
% Time sampling
[~,dt]=unix_cmd(['sugethw < ',seismofile_unix,' key=dt output=geom | uniq'],wsl);
dt=str2double(dt)/1000000;
% Nb of time samples
[~,ns]=unix_cmd(['sugethw < ',seismofile_unix,' key=ns output=geom | uniq'],wsl);
ns=str2double(ns);

if isunix==1 || ispc_wsl==1
    com1=sprintf('sustrip < %s | b2a n1=%d > %s',seismofile_unix,ns,[seismofile_unix,'.dat']);
else
     com1=sprintf('sustrip < %s head=head outpar=outpar | b2a n1=%d outpar=outpar > %s',...
        seismofile_unix,ns,[seismofile_unix,'.dat']);
end
[~,~]=unix_cmd(com1,wsl);

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