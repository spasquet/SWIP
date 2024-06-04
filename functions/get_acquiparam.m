function acquiparam=get_acquiparam(sufile,xsca,istomo)

%%% S. Pasquet - V22.05.04
% get_acquiparam.m retrieves .SU file acquisition parameters
%
% acquiparam=get_acquiparam(sufile,xsca,istomo)
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

wsl = ispc_wsl;

if exist('xsca','var')==0 || isempty(xsca)==1
    [~,xsca]=unix_cmd(['sugethw < ',sufile,' scalco output=geom | uniq'],wsl);
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
[~,dt]=unix_cmd(['sugethw < ',sufile,' dt output=geom | uniq'],wsl);
dt=str2double(dt)/1000000;
% Mean inter-geophone spacing
[~,dx]=unix_cmd(['sugethw < ',sufile,' gdel output=geom | uniq'],wsl);
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
[~,Gx]=unix_cmd(['sugethw < ',sufile,' key=gx output=geom'],wsl);
Gx=str2num(Gx)/xsca;
% Sources positions
[~,Sx]=unix_cmd(['sugethw < ',sufile,' key=sx output=geom'],wsl);
Sx=str2num(Sx)/xsca;
% Geophones elevations
[~,Gz]=unix_cmd(['sugethw < ',sufile,' key=gelev output=geom'],wsl);
Gz=str2num(Gz)/xsca;
% Sources elevations
[~,Sz]=unix_cmd(['sugethw < ',sufile,' key=selev output=geom'],wsl);
Sz=str2num(Sz)/xsca;
% Sources number
[~,fldr]=unix_cmd(['sugethw < ',sufile,' key=fldr output=geom'],wsl);
[fldr,I]=unique(str2num(fldr),'first');
% Sort source numbers with shot positions
[~,J]=sort(Sx(I));
fldr=fldr(J);

[Gxsing,I]=unique(Gx,'first'); % Single geophone positions
Gzsing=Gz(I); % Single geophone altitude
[Sxsing,I]=unique(Sx,'first'); % Single source positions
Szsing=Sz(I); % Single source altitude
NGx=length(Gxsing); % Nb of geophones positions
NSx=length(Sxsing); % Nb of sources positions
[Xsing,I]=unique([Gxsing;Sxsing],'first'); % X single positions
Zsing=[Gzsing;Szsing]; % Z single positions

array  = check_depth_array(Gx,Gz,Sx,Sz);
if ~isempty(array{1})
    topo = unique(sortrows([array{1}.G_sing; array{1}.S_sing],1),'rows');
else
    topo=sortrows([Xsing,Zsing(I)],1); % Topo
end

acquiparam=struct('dt',dt,'dx',dx,'Gx',Gx,'Gz',Gz,'Sx',Sx,'Sz',Sz,...
    'NGx',NGx,'NSx',NSx,'Gxsing',Gxsing,'Gzsing',Gzsing,...
    'Sxsing',Sxsing,'Szsing',Szsing,'fldr',fldr,'topo',topo,'xsca',xsca);
end