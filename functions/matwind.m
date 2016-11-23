function [seismomat,x,t,ntr]=matwind(sufile,sx,gmin,gmax,xsca,winsize,...
    seismofile,datsave,mute,tmin1,tmin2,tmax1,tmax2)

% S. Pasquet - V16.11.18
% SUWIND for matlab
% [seismomat,x,t,ntr]=matwind(sufile,sx,gmin,gmax,xsca,winsize,seismofile,datsave,...
% mute,tmin1,tmin2,tmax1,tmax2)

% Select shot then window between gmin and gmax
com1=sprintf(['suwind < %s key=sx min=%d max=%d |',...
    ' suwind key=gx min=%d max=%d > tmp.su'],...
    sufile,round(sx*xsca),round(sx*xsca),...
    round(gmin*xsca),round(gmax*xsca));
unix(com1);

% Get fldr in case two shots at the same position
[~,fldr]=unix('sugethw < tmp.su key=fldr output=geom | uniq');
fldr=str2num(fldr);

GxOK=[]; flag=0;
for i=1:length(fldr)
    [~,Gx]=unix(sprintf('suwind < tmp.su key=fldr min=%d max=%d | sugethw key=gx output=geom',fldr(i),fldr(i)));
    Gx=str2num(Gx)/xsca;
    ntr=length(Gx);
    if ntr==winsize
        com1=sprintf('suwind < tmp.su key=fldr min=%d max=%d > tmp2.su',fldr(i),fldr(i));
        unix(com1);
        flag=1; break
    elseif ntr<winsize
        if exist('tmp2.su','file')~=2
            com1=sprintf('suwind < tmp.su key=fldr min=%d max=%d > tmp2.su',fldr(i),fldr(i));
            unix(com1);
            GxOK=[GxOK;Gx];
        else
            GxKO=Gx(ismember(Gx,GxOK)==0);
            GxOK=[GxOK;Gx];
            for j=1:length(GxKO)
                com1=sprintf(['suwind < tmp.su key=fldr min=%d max=%d |',...
                    ' suwind key=gx min=%d max=%d > tmp4.su'],...
                    fldr(i),fldr(i),round(GxKO(j)*xsca),round(GxKO(j)*xsca));
                unix(com1);
                movefile('tmp2.su','tmp3.su');
                unix('cat tmp3.su tmp4.su > tmp2.su');
            end
            [~,Gx]=unix(sprintf('sugethw < tmp2.su key=gx output=geom'));
            Gx=str2num(Gx)/xsca;
            ntr=length(Gx);
            if ntr==winsize
                unix(sprintf('susort < tmp2.su gx > tmp3.su'));
                movefile('tmp3.su','tmp2.su');
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
    com1=sprintf('suwind < tmp.su key=fldr min=%d max=%d > tmp2.su',fldr(1),fldr(1));
    unix(com1);
end

% Mute if enabled
if mute==1
    if sx<gmin
        com1=sprintf(['sumute < tmp2.su key=gx xmute=%d,%d',...
            ' tmute=%d,%d mode=0 ntaper=50 |',...
            ' sumute key=gx xmute=%d,%d tmute=%d,%d mode=1 ntaper=50',...
            ' > tmp.su'],round(gmin*xsca),round(gmax*xsca),tmin1,tmin2,...
            round(gmin*xsca),round(gmax*xsca),tmax1,tmax2);
    else
        com1=sprintf(['sumute < tmp2.su key=gx xmute=%d,%d',...
            ' tmute=%d,%d mode=0 ntaper=50 |',...
            ' sumute key=gx xmute=%d,%d tmute=%d,%d mode=1 ntaper=50',...
            ' > tmp.su'],round(gmin*xsca),round(gmax*xsca),tmin2,tmin1,...
            round(gmin*xsca),round(gmax*xsca),tmax2,tmax1);
    end
    unix(com1);
    unix('rm tmp2.su');
else
    unix('mv tmp2.su tmp.su');
end

if sx==gmin
    com1=sprintf('sushw < tmp.su key=sx a=%d > tmp2.su',xsca*(gmin-median(diff(Gx))/2));
    unix(com1);
    unix('mv tmp2.su tmp.su');
elseif sx==gmax
    com1=sprintf('sushw < tmp.su key=sx a=%d > tmp2.su',xsca*(gmax+median(diff(Gx))/2));
    unix(com1);
    unix('mv tmp2.su tmp.su');
end

% Remove continuous part and save in .SU file
com1=sprintf('suop < tmp.su op=avg > %s',seismofile);
unix(com1);
unix('rm -f tmp*.su');

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