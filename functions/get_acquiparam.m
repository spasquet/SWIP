function acquiparam=get_acquiparam(sufile,xsca,istomo)

%%% S. Pasquet - V16.9.14
% get_acquiparam.m retrieves .SU file acquisition parameters
%
% acquiparam=get_acquiparam(sufile,xsca)
%
% sufile (required) => name of the .SU file
% xsca (optionnal) => scaling factor for spatial coordinates
%
% acquiparam.dt => sampling interval in sec
% acquiparam.dx => mean inter-geophone spacing in m
% acquiparam.Gx => geophones positions in m
% acquiparam.Sx => sources positions in m
% acquiparam.Gz => geophones elevation in m
% acquiparam.Sz => sources elevation in m
% acquiparam.NGx => Nb of geophones positions
% acquiparam.NSx => Nb of sources positions
% acquiparam.Gxsing => Single geophones positions
% acquiparam.Sxsing => Single sources positions
% acquiparam.Gzsing => Single geophones elevations
% acquiparam.Szsing => Single sources elevations
% acquiparam.fldr => Source shot number
% acquiparam.topo => Profile topography (X,Z)

if exist('xsca','var')==0 || isempty(xsca)==1
    [~,xsca]=unix(['sugethw < ',sufile,' scalco output=geom | uniq']);
    xsca=abs(str2double(xsca));
    if xsca==0
       answer1=inputdlg({'Scaling factor'},'',1);
    if isempty(answer1)==1
        xsca=0;
    else
        xsca=str2double(answer1(1));
    end
    end
end

if nargin<3
    istomo=0;
end

% Sampling interval (in second)
[~,dt]=unix(['sugethw < ',sufile,' dt output=geom | uniq']);
dt=str2double(dt)/1000000;
% Mean inter-geophone spacing
[~,dx]=unix(['sugethw < ',sufile,' gdel output=geom | uniq']);
dx=str2double(dx)/xsca;
if dx==0 && istomo==0
    answer1=inputdlg({'Mean geophone spacing'},'',1);
    if isempty(answer1)==1
        dx=0;
    else
        dx=str2double(answer1(1));
    end
end
% Geophones positions
[~,Gx]=unix(['sugethw < ',sufile,' key=gx output=geom']);
Gx=str2num(Gx)/xsca;
% Sources positions
[~,Sx]=unix(['sugethw < ',sufile,' key=sx output=geom']);
Sx=str2num(Sx)/xsca;
% Geophones elevations
[~,Gz]=unix(['sugethw < ',sufile,' key=gelev output=geom']);
Gz=str2num(Gz)/xsca;
% Sources elevations
[~,Sz]=unix(['sugethw < ',sufile,' key=selev output=geom']);
Sz=str2num(Sz)/xsca;
% Sources number
[~,fldr]=unix(['sugethw < ',sufile,' key=fldr output=geom | uniq']);
fldr=str2num(fldr);

[Gxsing,I]=unique(Gx,'first'); % Single geophone positions
Gzsing=Gz(I); % Single geophone altitude
[Sxsing,I]=unique(Sx,'first'); % Single source positions
Szsing=Sz(I); % Single source altitude
NGx=length(Gxsing); % Nb of geophones positions
NSx=length(Sxsing); % Nb of sources positions
[Xsing,I]=unique([Gxsing;Sxsing],'first'); % X single positions
Zsing=[Gzsing;Szsing]; % Z single positions
topo=sortrows([Xsing,Zsing(I)],1); % Topo

acquiparam=struct('dt',dt,'dx',dx,'Gx',Gx,'Gz',Gz,'Sx',Sx,'Sz',Sz,...
    'NGx',NGx,'NSx',NSx,'Gxsing',Gxsing,'Gzsing',Gzsing,...
    'Sxsing',Sxsing,'Szsing',Szsing,'fldr',fldr,'topo',topo,'xsca',xsca);
end