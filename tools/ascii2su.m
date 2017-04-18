clear all; clc; close all;

%%% ASCII2SU - Convert ASCII files to SU
%%% S. Pasquet - V17.01.16

%-----------------------%
% START OF INITIALIZATION

profilename='swip_profile'; % Name of the profile
% newrectime=1; % New recording length (add zeros at the end of each trace)

% END OF INITIALIZATION
%-----------------------%

%% FILE CONVERSION

dir_start=pwd;
dir_asc=uigetdir('./','Select folder containing ASCII files');
if dir_asc==0
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n   Please select a folder');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!\n\n');
    return
end

ascstruct=[dir(fullfile(dir_asc,'*txt'));dir(fullfile(dir_asc,'*dat'))];
nasc=length(ascstruct);

if nasc==0
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n   No ASCII files in the selected folder');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
    return
end

aa=load(fullfile(dir_asc,ascstruct(1).name));
ngeo=min(size(aa));
nsamp=max(size(aa));
% keyboard
if nsamp>65535 && nsamp==size(aa,1)
    aa=aa(1:65535,:);
    fprintf('\n   Too many samples - Reducing number of samples to 65535\n');
    movefile(fullfile(dir_asc,ascstruct(1).name),fullfile(dir_asc,[ascstruct(1).name,'.old']));
    dlmwrite(fullfile(dir_asc,ascstruct(1).name),aa,'delimiter',' ','precision','%10.5e');
end

linetest=size(aa,2)==max(size(aa));
if linetest==1
    n1=ngeo;
    n2=nsamp;
else
    n1=nsamp;
    n2=ngeo;
end

fprintf('\n   Acquisition settings\n');
fprintf(['\n   ',num2str(nasc),' shot(s)']);
fprintf(['\n   ',num2str(ngeo),' traces per shot']);
fprintf(['\n   ',num2str(nsamp),' samples per trace\n']);
answer1=inputdlg({'First trace pos. (m)','Trace spacing (m)','First shot pos. (m)','Shot spacing (m)',...
    'Starting time (s)','Sampling time (ms)','Scaling factor'},'Acquisition settings',1,{'0','1','0','10','0','0.125','100'});
GX0=str2num(answer1{1});
DGX=str2num(answer1{2});
SX0=str2num(answer1{3});
DSX=str2num(answer1{4});
time0=str2num(answer1{5});
Dtime=str2num(answer1{6});
xsca=str2num(answer1{7});

fprintf(['\n   dt = ',num2str(Dtime),' ms']);
fprintf(['\n   Record length = ',num2str(nsamp*Dtime/1000),' s\n']);

fprintf('\n   Creating SEGY files...\n');
for i=1:nasc
    com1=sprintf('a2b < %s > %s n1=%d',fullfile(dir_asc,ascstruct(i).name),fullfile(dir_asc,'tmp.bin'),n1);
    [~,~]=unix(com1);
    if linetest~=1
        com1=sprintf('suaddhead ns=%d < %s | suflip flip=0 | sushw key=dt,gx,sx,tracf,gdel,fldr,scalco,scalel a=%d,%d,%d,%d,%d,%d,%d,%d b=0,%d,0,%d,0,0,0,0 > %s',...
            n2,fullfile(dir_asc,'tmp.bin'),Dtime*1000,GX0*xsca,SX0*xsca+(i-1)*DSX*xsca,1,DGX*xsca,1000+i,xsca,xsca,DGX*xsca,1,fullfile(dir_asc,[num2str(1000+i),'.su']));
        [~,~]=unix(com1);
    else
        com1=sprintf('suaddhead ns=%d < %s | sushw key=dt,gx,sx,tracf,gdel,fldr,scalco,scalel a=%d,%d,%d,%d,%d,%d,%d,%d b=0,%d,0,%d,0,0,0,0 > %s',...
            n2,fullfile(dir_asc,'tmp.bin'),Dtime*1000,GX0*xsca,SX0*xsca+(i-1)*DSX*xsca,1,DGX*xsca,1000+i,xsca,xsca,DGX*xsca,1,fullfile(dir_asc,[num2str(1000+i),'.su']));
        unix(com1);
    end
    [~,~]=unix(sprintf('segyhdrs < %s',fullfile(dir_asc,[num2str(1000+i),'.su'])));
    [~,~]=unix(sprintf('segywrite < %s tape=%s endian=0',fullfile(dir_asc,[num2str(1000+i),'.su']),fullfile(dir_asc,[num2str(1000+i),'.sgy'])));
    delete(fullfile(dir_asc,[num2str(1000+i),'.su']),fullfile(dir_asc,'tmp.bin'));
end

segystruct=[dir(fullfile(dir_asc,'*segy'));dir(fullfile(dir_asc,'*sgy'))];
nsegy=length(segystruct);

xshift=0; zshift=0;
rawsufile=fullfile(dir_start,[profilename,'.raw.su']);
sufile=fullfile(dir_start,[profilename,'.su']);

cd(dir_asc);
nfile=length(segystruct);
for i=1:nfile
    shotfile=segystruct(i).name;
    [test,fldrcell{i},extension]=fileparts(segystruct(i).name);
    if strcmp(extension,'.segy')==1
        unix(['mv ',shotfile,' ',shotfile(1:end-5),'.sgy']);
    end
    if nfile>1
        finddot=strfind(fldrcell{i},'.');
        if isempty(finddot)==1
            fldr(i)=str2num(fldrcell{i});
        else
            fldr(i)=str2num(fldrcell{i}(1:finddot(1)));
        end
    end
    com1=['segyread tape=',fldrcell{i},'.sgy endian=0 | segyclean >',...
        fldrcell{i},'.su'];
    [~,~]=unix(com1);
end
unix(['cat *.su > ',rawsufile]);
delete('*.su','*.sgy',fullfile(dir_asc,'binary'),fullfile(dir_asc,'header'),fullfile(dir_start,'binary'),fullfile(dir_start,'header'));

cd(dir_start);

fprintf('\n   Checking acquisition settings in SEGY file(s)\n');
% Get acquisition settings
[~,hdrs]=unix(['sugethw < ',rawsufile,' key=fldr,delrt,dt,ns,sx,gx,tracf output=geom']);
hdrs=str2num(hdrs);
fldrall=hdrs(:,1);
fldr_raw=unique(fldrall);
if length(fldr_raw)==nfile
    fldr=fldr_raw;
end
NSX=length(fldr);
sxall=hdrs(:,5)/xsca;
gxall=hdrs(:,6)/xsca;
tracfall=hdrs(:,7);
gzall=gxall*0; szall=sxall*0;

% Time delay
delay=unique(hdrs(:,2))/1000;
if length(delay)>1
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n   Different delays found in headers');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
    delay=max(delay);
    unix(['sushw < ',rawsufile,' key=delrt a=',num2str(delay*1000),'>tmp.su']);
    movefile('tmp.su',rawsufile);
end

% Time sampling
dt=unique(hdrs(:,3))/1000000;
if length(dt)>1
    fprintf('\n\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n   Different dt found in headers');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
    dt=max(dt);
    unix(['sushw < ',rawsufile,' key=dt a=',num2str(dt*1000000),'>tmp.su']);
    movefile('tmp.su',rawsufile);
end

% Nb of time samples
ns=unique(hdrs(:,4));
if length(ns)>1
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n   Different ns found in headers');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
    ns=max(ns);
    unix(['sushw < ',rawsufile,' key=ns a=',num2str(ns),'>tmp.su']);
    movefile('tmp.su',rawsufile);
end

fprintf(['\n   ',num2str(NSX),' shot(s)']);
fprintf(['\n   delay = ',num2str(delay),' s']);
fprintf(['\n   dt = ',num2str(dt*1000),' ms']);
fprintf(['\n   ns = ',num2str(ns)]);
fprintf(['\n   Record length = ',num2str(ns*dt),' s\n']);

if length(fldr_raw)==nfile
    % Deployments settings
    fprintf('\n   Original deployment settings\n');
    nroll=0;
    for i=1:NSX
        sx(i)=unique(sxall(fldrall==fldr(i)));
        gx{i}=gxall(fldrall==fldr(i));
        
        if (i>1 && isequal(gx{i},gx{i-1})==0) || i==NSX
            nroll=nroll+1;
            fprintf(['\n   Deployment ',num2str(nroll),' settings']);
            if nroll==1
                if i==NSX
                    nsx(nroll)=i;
                else
                    nsx(nroll)=i-1;
                end
                sx1(nroll)=sx(1);
                dsx(nroll)=median(diff(sx(1:nsx(nroll))));
            elseif i==NSX
                nsx(nroll)=i-sum(nsx(1:nroll-1));
                sx1(nroll)=sx(sum(nsx(1:nroll-1))+1);
                dsx(nroll)=median(diff(sx(1+nsx(nroll-1):nsx(nroll-1)+nsx(nroll))));
            else
                nsx(nroll)=i-1-sum(nsx(1:nroll-1));
                sx1(nroll)=sx(sum(nsx(1:nroll-1))+1);
                dsx(nroll)=median(diff(sx(1+nsx(nroll-1):nsx(nroll-1)+nsx(nroll))));
            end
            fprintf(['\n   ',num2str(nsx(nroll)),' shot(s)']);
            fprintf(['\n   First shot at ',num2str(sx1(nroll)),' m']);
            fprintf(['\n   dsx = ',num2str(dsx(nroll)),' m']);
            if i==NSX
                ngx(nroll)=length(gx{i});
                gx1(nroll)=gx{i}(1);
                dgx(nroll)=median(diff(gx{i}));
            else
                ngx(nroll)=length(gx{i-1});
                gx1(nroll)=gx{i-1}(1);
                dgx(nroll)=median(diff(gx{i-1}));
            end
            fprintf(['\n   ',num2str(ngx(nroll)),' trace(s)']);
            fprintf(['\n   First trace at ',num2str(gx1(nroll)),' m']);
            fprintf(['\n   dgx = ',num2str(dgx(nroll)),' m\n']);
        end
    end
    
    cmap=hsv(nroll); RGB=[];
    for i=1:nroll
        RGB=[RGB;repmat(cmap(i,:),ngx(i)*nsx(i),1)];
    end
    
    figure(1);
    scatter(sxall,gxall,50,RGB,'.');
    axis equal; grid on; box on;
    title('Original acquisition settings');
    xlabel('Source position (m)');
    ylabel('Geophone position (m)');
    h=get(gca,'title'); set(h,'FontSize',14);
    h=get(gca,'xlabel'); set(h,'FontSize',14);
    h=get(gca,'ylabel'); set(h,'FontSize',14);
    set(gca,'TickDir','out','linewidth',1,'XMinorTick','on','YMinorTick','on');
    h=findall(gcf,'Type','Axes'); set(h,'FontSize',14);
    set(gcf,'Units','centimeters');
    set(gcf,'Position',[2 2 24 18]);
else
    fprintf('\n   Unable to read original deployment settings\n');
    nroll=1;
    gx1=0; dgx=1; ngx=length(hdrs)/nfile;
    sx1=0; dsx=1; nsx=nfile;
end

answer1=inputdlg({'Scaling factor'},...
    '',1,{'100'});
if isempty(answer1)==1
    return
end
xsca=str2double(answer1(1)); % Scaling factor

% Accept current configuration
flagOK=1; flaghdrs=0;
while flagOK~=0
    if length(fldr_raw)<nfile && flaghdrs==0
        flagOK=2;
    else
        choice = questdlg('Do you want to keep the current settings?', ...
            '', 'Yes','No, reset all headers','Yes');
        % Handle response
        switch choice
            case 'Yes'
                flagOK = 0;
            case 'No, reset all headers'
                flagOK = 2;
        end
    end
    
    if flagOK==1
        fprintf('\n   This option is not yet available\n');
        flaghdrs=1;
        
    elseif flagOK==2
        fprintf('\n   Updated deployment settings\n');
        flaghdrs=1;
        answer1=inputdlg({'Number of deployments'},...
            '',1,{num2str(nroll)});
        if isempty(answer1)==1
            return
        end
        nroll=str2double(answer1(1)); % Number of deployments
        
        gxall=[]; sxall=[]; RGB=[]; cmap=hsv(nroll);
        fldrall=[]; tracfall=[];
        for i=1:nroll
            fprintf(['\n   Deployment ',num2str(nroll),' settings']);
            answer2=inputdlg({['Vector of shot positions (in m) for deployment ',num2str(i)],...
                ['Vector of traces positions (in m) for deployment ',num2str(i)]},...
                ['Deployment ',num2str(i),' settings'],2,...
                {[num2str(sx1(i)),':',num2str(dsx(i)),':',num2str(sx1(i)+(nsx(i)-1)*dsx(i))],...
                [num2str(gx1(i)),':',num2str(dgx(i)),':',num2str(gx1(i)+(ngx(i)-1)*dgx(i))]});
            sxroll{i}=str2num(answer2{1});
            gxroll{i}=str2num(answer2{2});
            nsx(i)=length(sxroll{i});
            ngx(i)=length(gxroll{i});
            dsx(i)=mean(unique(diff(sxroll{i})));
            dgx(i)=mean(unique(diff(gxroll{i})));
            for j=1:nsx(i)
                sxall=[sxall;repmat(sxroll{i}(j),ngx(i),1)];
                gxall=[gxall;gxroll{i}'];
                RGB=[RGB;repmat(cmap(i,:),ngx(i),1)];
                fldrall=[fldrall;repmat(fldr(j),ngx(i),1)];
                tracfall=[tracfall;(1:ngx(i))'];
            end
            fprintf(['\n   ',num2str(nsx(nroll)),' shot(s)']);
            fprintf(['\n   First shot at ',num2str(sx1(nroll)),' m']);
            fprintf(['\n   dsx = ',num2str(dsx(nroll)),' m']);
            fprintf(['\n   ',num2str(ngx(nroll)),' trace(s)']);
            fprintf(['\n   First trace at ',num2str(gx1(nroll)),' m']);
            fprintf(['\n   dgx = ',num2str(dgx(nroll)),' m\n']);
        end
        
        figure(3);
        scatter(sxall,gxall,50,RGB,'.');
        axis equal; grid on; box on;
        title('Updated acquisition settings');
        xlabel('Source position (m)');
        ylabel('Geophone position (m)');
        h=get(gca,'title'); set(h,'FontSize',14);
        h=get(gca,'xlabel'); set(h,'FontSize',14);
        h=get(gca,'ylabel'); set(h,'FontSize',14);
        set(gca,'TickDir','out','linewidth',1,'XMinorTick','on','YMinorTick','on');
        h=findall(gcf,'Type','Axes'); set(h,'FontSize',14);
        set(gcf,'Units','centimeters');
        set(gcf,'Position',[28 2 24 18]);
    end
end

% Import topography
choice = questdlg('Import topography?', ...
    '', 'No','Yes','No');
% Handle response
switch choice
    case 'No'
        flagtopo = 0;
%         if flaghdrs==0
%             while fix(xsca*dgx/2)~=xsca*dgx/2
%                 xsca=xsca*10;
%             end
%             fprintf(['\n   X scaling factor (xsca) set to ',num2str(xsca),' to enable Xmid position with decimals\n']);
%         end
    case 'Yes'
        flagtopo = 1;
end

if flagtopo==1
    [filetopo,pathtopo]=uigetfile('*','Select topo ASCII file with 2 columns (X,Z)');
    if filetopo==0
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        fprintf('\n   No topo file selected - altitude set as 0');
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
        X=unique([gxall;sxall]);
        Z=zeros(size(X));
        flagtopo=0;
    else
        try
            topo=load(fullfile(pathtopo,filetopo));
            X=unique([gxall;sxall]); Z=interp1(topo(:,1),topo(:,2),X,'linear','extrap');
%             answer1=inputdlg({'Z shift (m)'},'',1,{'0'});
%             zshift=str2double(answer1(1)); % Z shift
%             xsca=str2double(answer1(2)); % Scaling factor
%             Z=Z+zshift;
        catch
            fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
            fprintf('\n   Invalid format - altitude set as 0');
            fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
            X=unique([gxall;sxall]);
            Z=zeros(size(X));
            flagtopo=0;
        end
    end
else
    X=unique([gxall;sxall]);
    Z=zeros(size(X));
end
figure(2);
plot(X,Z,'k-','linewidth',1);
grid on; box on;
title('Imported topography profile');
xlabel('X (m)');
ylabel('Altitude (m)');
h=get(gca,'title'); set(h,'FontSize',14);
h=get(gca,'xlabel'); set(h,'FontSize',14);
h=get(gca,'ylabel'); set(h,'FontSize',14);
set(gca,'TickDir','out','linewidth',1,'XMinorTick','on','YMinorTick','on');
h=findall(gcf,'Type','Axes'); set(h,'FontSize',14);
set(gcf,'Units','centimeters');
set(gcf,'Position',[2 24 24 6]);

for i=1:length(Z)
    szall(sxall==X(i))=Z(i);
    gzall(gxall==X(i))=Z(i);
end

dx=mean(diff(unique(gxall)))*xsca;

dlmwrite('infile.txt',[fldrall,tracfall,(gxall+xshift)*xsca,(gzall+zshift)*xsca,...
    (sxall+xshift)*xsca,(szall+zshift)*xsca],'delimiter',' ','precision',6);
unix('a2b < infile.txt > infile.bin');
unix(['sushw < ',rawsufile,' key=scalel,scalco,gdel a=',...
    num2str(-xsca),',',num2str(-xsca),',',num2str(dx),...
    ' | sushw key=fldr,tracf,gx,gelev,sx,selev infile=infile.bin > ',sufile]);
if exist('newrectime','var')==1 && isempty(newrectime)==0
    nsnew=newrectime/dt;
    unix(['suvlength < ',sufile,' ns=',num2str(nsnew),' > tmp.su']);
    unix(['mv tmp.su ',sufile]);
end
unix(['rm -rf infile.txt infile.bin ',rawsufile]);

sufile=strrep(sufile,'\','\\');
fprintf(['\n   Shots saved in ',sufile,'\n']);
