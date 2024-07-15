clear all; clc; close all;

%%% SEG2SU - Convert SEG2 or SEGY files to SU
%%% S. Pasquet - V19.03.11
%%% Clean version with no stacking option

%-----------------------%
% START OF INITIALIZATION

profilename = 'swip_profile'; % Name of the profile
newrectime  = 1;           % New recording length (add zeros at the end of each trace)
colx        = 1;              % Column for X position (in topofile) -- Leave empty or comment for default
colz        = 4;              % Column for Z position (in topofile) -- Leave empty or comment for default
% trace_selec = [1 48];        % Traces to select in file Leave empty or comment for default
% trace_shift = 48;
% xshift      = -192;          % Shift along horizontal axis
% zshift      = 0;             % Shift along vertical axis
% delay = -0.05;

% END OF INITIALIZATION
%-----------------------%

%% FILE CONVERSION

if (~exist('colx','var') || isempty(colx)) || (~exist('colz','var') || isempty(colz))
    colx = 1;
    colz = 2;
end

if ~exist('stack','var') || isempty(stack)
    stack = 0;
end

if ~exist('xshift','var') || isempty(xshift)
    xshift = 0;
end

if ~exist('zshift','var') || isempty(zshift)
    zshift = 0;
end

dir_start=pwd;
dir_seg=uigetdir('./','Select folder containing SEG2 or SEGY files');
if dir_seg==0
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n   Please select a folder');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!\n\n');
    return
end

seg2struct=[dir(fullfile(dir_seg,'*sg2'));dir(fullfile(dir_seg,'*seg2'));dir(fullfile(dir_seg,'*dat'))];
nseg2=length(seg2struct);
segystruct=[dir(fullfile(dir_seg,'*segy'));dir(fullfile(dir_seg,'*sgy'))];
nsegy=length(segystruct);

if nseg2==0 && nsegy==0
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n   No SEG2 or SEGY files in the selected folder');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
    return
elseif nseg2>0 && nsegy>0
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n   Both SEG2 and SEGY files in the selected folder');
    fprintf('\n   Keep only one format in the folder and re-run script');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
    return
end

xsca=1;
rawsufile=fullfile(dir_start,[profilename,'.raw.su']);
sufile=fullfile(dir_start,[profilename,'.su']);

if exist(sufile,'file')==2
    answer1=inputdlg({'New file name'},'File already exists',1,{profilename});
    profilename=answer1{1};
    rawsufile=fullfile(dir_start,[profilename,'.raw.su']);
    sufile=fullfile(dir_start,[profilename,'.su']);
end

cd(dir_seg);
if nseg2>0
    nfile=length(seg2struct);
    for i=1:nfile
        [~,fldrcell{i},extension]=fileparts(seg2struct(i).name);
        fldr(i)=str2num(fldrcell{i});
        if fldr(i) < 2000
%             movefile(seg2struct(i).name,sprintf('%i.dat',1000+nfile-i+1));
            movefile(seg2struct(i).name,sprintf('%i.dat',2000+i));
            [~,fldrcell{i},extension]=fileparts(seg2struct(i).name);
            fldr(i)=str2num(fldrcell{i});
        end
    end
    firstshot=min(fldr);
    lastshot=max(fldr);
    fprintf('\n   Reading SEG2 files...\n');
    com1=['seg2segy ',num2str(firstshot),extension,' ',num2str(lastshot-firstshot+1)];
    [~,~]=unix(com1);
        
    movefile([num2str(firstshot),'.sgy'],[profilename,'.sgy']);
    com1=['segyread tape=',profilename,'.sgy endian=0 | segyclean > ',rawsufile];
    [~,~]=unix(com1);
    delete([profilename,'.sgy'],'binary','header');
else
    nfile=length(segystruct);
    for i=1:nfile
        shotfile=segystruct(i).name;
        [~,fldrcell{i},extension]=fileparts(segystruct(i).name);
        if strcmp(extension,'.segy')==1
            unix(['mv ',shotfile,' ',shotfile(1:end-5),'.sgy']);
        end
        if nfile>1
            findstring = isstrprop(fldrcell{i},'digit');
            fldr(i)=str2num(fldrcell{i}(findstring));
            if fldr(i) < 3000
                movefile(segystruct(i).name,sprintf('%i.sgy',3000+i));
                [~,fldrcell{i},extension]=fileparts(segystruct(i).name);
                fldr(i)=str2num(fldrcell{i});
            end
        end
        com1=['segyread tape=',fldrcell{i},'.sgy endian=0 | segyclean > tmp.su'];
        [~,~]=unix(com1);
        [~,~]=unix(sprintf('sushw < %s key=fldr a=%i > %s','tmp.su', fldr(i), [fldrcell{i},'.su']));
        delete('tmp.su');
    end
    unix(['cat *.su > ',rawsufile]);
    delete('*.su');
end
cd(dir_start);

unix(sprintf('susort < %s +fldr +tracf > tmp.su',rawsufile));
movefile('tmp.su',rawsufile);

if exist('trace_shift','var')==1 && ~isempty(trace_shift)
    unix(['suchw < ',rawsufile,' key1=tracf key2=tracf a=',num2str(-trace_shift-1),' d=-1 > tmp.su']);
    movefile('tmp.su',rawsufile);
else
    trace_shift = 0;
end

if exist('trace_selec','var')==1 && ~isempty(trace_selec)
    unix(['suwind < ',rawsufile,' key=tracf min=',num2str(trace_selec(1)),' max=',num2str(trace_selec(2)),' > tmp.su']);
    movefile('tmp.su',rawsufile);
end

%%
fprintf('\n   Extract original acquisition settings from SEG2(Y) file(s)\n');
% Get acquisition settings
[~,hdrs] = unix(['sugethw < ',rawsufile,' key=fldr,delrt,dt,ns,sx,gx,tracf,scalco output=geom']);
hdrs = str2num(hdrs);
fldrall = hdrs(:,1);
if ~any(fldrall)
    [~,~,fldrall] = unique(hdrs(:,5));
end
fldr_raw = unique(fldrall);
if length(fldr_raw) == nfile
    fldr = fldr_raw;
else
    fldr_pb = fldrall(find(ismember(fldrall,fldr),1,'last'));
    fprintf(['\n   Acquisition parameters changed at shot #' num2str(fldr_pb)]);
    fprintf('\n   Run seg2su separately for shots recorded with identical parameters\n');
    fprintf('\n   Run sumerge to merge SU files\n');
    return
end
NSX = length(fldr);
xsca = hdrs(:,8);
xsca = abs(xsca).^-sign(xsca);
sxall = (hdrs(:,5)./xsca) - xshift;
gxall = (hdrs(:,6)./xsca) - xshift;
tracfall = hdrs(:,7);

% Time delay
if ~exist('delay','var') || isempty(delay)
    delay=unique(hdrs(:,2))/1000;
    if length(delay)>1
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        fprintf('\n   Different delays found in headers');
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
        delay=max(delay);
        unix(['sushw < ',rawsufile,' key=delrt a=',num2str(delay*1000),'>tmp.su']);
        unix(['mv -f tmp.su ',rawsufile]);
    end
end

% Time sampling
dt=unique(hdrs(:,3))/1000000;
if length(dt)>1
    fprintf('\n\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n   Different dt found in headers');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
    dt=max(dt);
    unix(['sushw < ',rawsufile,' key=dt a=',num2str(dt*1000000),'>tmp.su']);
    unix(['mv -f tmp.su ',rawsufile]);
end

% Nb of time samples
ns=unique(hdrs(:,4));
if length(ns)>1
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n   Different ns found in headers');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
    ns=max(ns);
    unix(['sushw < ',rawsufile,' key=ns a=',num2str(ns),'>tmp.su']);
    unix(['mv -f tmp.su ',rawsufile]);
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
        
        if any(abs(gx{i})>1e6)
            gx{i} = gx{i}*0;
        end
        
%         if i>1 && length(gx{i-1}) - length(gx{i}) == 1 %% Fix if one trace has been deactivated
%             fprintf('\n   One trace missing from original acquistion setup - add empty trace\n');
%             new_tracf = find(~ismember(gx{i-1},gx{i}));
%             new_gx = gx{i-1}(new_tracf);
%             
%             unix(sprintf('suwind < %s count=1 > tmp.su',rawsufile));
%             unix(sprintf('suop2 tmp.su tmp.su op=diff > tmp2.su'));
%             unix(sprintf('sushw < tmp2.su key=tracf,sx,gx a=%i,%i,%i > tmp.su',new_tracf,new_gx,sx(i)));
%             unix(sprintf('cat %s tmp.su > tmp2.su ',rawsufile));
%             unix(sprintf('susort < tmp2.su +fldr +tracf > %s',rawsufile));
%             
%             gx{i} = gx{i-1};
%             gxall = [gxall; new_gx];
%             sxall = [sxall; sx(i)];
%             fldrall = [fldrall; fldr(i)];
%             tracfall = [tracfall; new_tracf];
%             all = sortrows([tracfall fldrall gxall sxall],[2 1]);
%         end
        
        if (i>1 && isequal(gx{i},gx{i-1})==0) || i==NSX
            nroll=nroll+1;
            fprintf(['\n   Deployment ',num2str(nroll),' settings']);
            if nroll==1
                if i==NSX
                    nsx(nroll)=i;
                    ngx(nroll)=length(gx{i});
                    gx1(nroll)=gx{i}(1);
                    dgx(nroll)=nanmedian(diff(gx{i}));
                else
                    nsx(nroll)=i-1;
                    ngx(nroll)=length(gx{i-1});
                    gx1(nroll)=gx{i-1}(1);
                    dgx(nroll)=nanmedian(diff(gx{i-1}));
                end
                sx1(nroll)=sx(1);
                dsx(nroll)=nanmedian(diff(sx(1:nsx(nroll))));
            else
                if i==NSX
                    nsx(nroll)=i-sum(nsx(1:nroll-1));
                    ngx(nroll)=length(gx{i});
                    gx1(nroll)=gx{i}(1);
                    dgx(nroll)=nanmedian(diff(gx{i}));
                else
                    nsx(nroll)=i-1-sum(nsx(1:nroll-1));
                    ngx(nroll)=length(gx{i-1});
                    gx1(nroll)=gx{i-1}(1);
                    dgx(nroll)=nanmedian(diff(gx{i-1}));
                end
                sx1(nroll)=sx(sum(nsx(1:nroll-1))+1);
                dsx(nroll)=nanmedian(diff(sx(1+nsx(nroll-1):nsx(nroll-1)+nsx(nroll))));
            end
            fprintf(['\n   ',num2str(nsx(nroll)),' shot(s)']);
            all_shots = sx(1+sum(nsx(1:nroll-1)):nsx(nroll)+sum(nsx(1:nroll-1)));
            vec_shots = 1:10:length(all_shots);
            for iii = 1:length(vec_shots)
                if iii == 1
                    ind_shot = min([length(all_shots) 10]);
                    fprintf(['\n   Shots at: ',num2str(all_shots(1:ind_shot))]);
                elseif iii == length(vec_shots)
                    fprintf(['\n   ',num2str(all_shots(vec_shots(iii):end))]);
                else
                    fprintf(['\n   ',num2str(all_shots(vec_shots(iii):vec_shots(iii)+9))]);
                end
            end
            fprintf(['\n   dsx = ',num2str(dsx(nroll)),' m']);
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

gzall = gxall*0; szall=sxall*0;

answer1=inputdlg({'Scaling factor'},'',1,{'100'});
if isempty(answer1)==1
    xsca=100; % Scaling factor
else
    xsca=str2double(answer1(1)); % Scaling factor
end

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
        answer1=inputdlg({'Number of deployments'},'',1,{num2str(nroll)});
        if isempty(answer1)==1
            return
        end
        nroll=str2double(answer1(1)); % Number of deployments
        
        gxall=[]; sxall=[]; RGB=[]; cmap=hsv(nroll);
        fldrall=[]; tracfall=[]; k=0;
        if length(dsx)<nroll
            sx1 = sx1(1)*ones(1,nroll);
            dsx = dsx(1)*ones(1,nroll);
            nsx = (NSX/nroll)*ones(1,nroll);
        end
        if length(dgx)<nroll
            gx1 = gx1(1)*ones(1,nroll);
            dgx = dgx(1)*ones(1,nroll);
            ngx = ngx*ones(1,nroll);
        end
        for i=1:nroll
            fprintf(['\n   Deployment ',num2str(i),' settings']);
            if dsx(i)~=0 && isnan(dsx(i))==0
                def_shots_pos = [num2str(sx1(i)),':',num2str(dsx(i)),':',num2str(sx1(i)+(nsx(i)-1)*dsx(i))];
            else
                def_shots_pos = num2str(sx1(i));
            end
            if dgx(i)~=0 && isnan(dgx(i))==0
                def_traces_pos = [num2str(gx1(i)),':',num2str(dgx(i)),':',num2str(gx1(i)+(ngx(i)-1)*dgx(i))];
            else
                def_traces_pos = num2str(gx1(i));
            end
            
            answer2=inputdlg({['Vector of shot positions (in m) for deployment ',num2str(i)],...
                ['Vector of traces positions (in m) for deployment ',num2str(i)]},...
                ['Deployment ',num2str(i),' settings'],2,{def_shots_pos,def_traces_pos});
            
            sxroll{i}=str2num(answer2{1});
            gxroll{i}=str2num(answer2{2});
            sx1(i)=sxroll{i}(1);
            gx1(i)=gxroll{i}(1);
            nsx(i)=length(sxroll{i});
            ngx(i)=length(gxroll{i});
            dsx(i)=mean(unique(diff(sxroll{i})));
            dgx(i)=mean(unique(diff(gxroll{i})));
            for j=1:nsx(i)
                k=k+1;
                sxall=[sxall;repmat(sxroll{i}(j),ngx(i),1)];
                gxall=[gxall;gxroll{i}'];
                RGB=[RGB;repmat(cmap(i,:),ngx(i),1)];
                fldrall=[fldrall;repmat(fldr(k),ngx(i),1)];
%                 tracfall=[tracfall;-((1:ngx(i))'-trace_shift-1)+trace_shift];
                tracfall=[tracfall;(1:ngx(i))'+trace_shift];
            end
            fprintf(['\n   ',num2str(nsx(i)),' shot(s)']);
            all_shots = sxroll{i};
            vec_shots = 1:10:length(all_shots);
            for iii = 1:length(vec_shots)
                if iii == 1
                    ind_shot = min([length(all_shots) 10]);
                    fprintf(['\n   Shots at: ',num2str(all_shots(1:ind_shot))]);
                elseif iii == length(vec_shots)
                    fprintf(['\n   ',num2str(all_shots(vec_shots(iii):end))]);
                else
                    fprintf(['\n   ',num2str(all_shots(vec_shots(iii):vec_shots(iii)+9))]);
                end
            end
            fprintf(['\n   dsx = ',num2str(dsx(i)),' m']);
            fprintf(['\n   ',num2str(ngx(i)),' trace(s)']);
            fprintf(['\n   First trace at ',num2str(gx1(i)),' m']);
            fprintf(['\n   dgx = ',num2str(dgx(i)),' m\n']);
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
    '','Yes (linear)','Yes (curvilinear)','No','Yes (linear)');
% Handle response
switch choice
    case 'No'
        flagtopo = 0;
    case 'Yes (linear)'
        flagtopo = 1;
    case 'Yes (curvilinear)'
        flagtopo = 2;
end
%%
if flagtopo==1 || flagtopo==2
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
            [~,II] = unique(topo(:,1));
            topo = topo(II,:);
            
            if flagtopo==1
                X=unique([gxall;sxall]);
                Z=interp1(topo(~isnan(topo(:,colz)),colx),topo(~isnan(topo(:,colz)),colz),X,'linear','extrap');
            elseif flagtopo==2
                X_raw = unique([gxall;sxall]);
                X_lin=topo(:,colx);
                X_lin=[min(X_raw):min(dgx)/2:ceil(max(topo(:,colx)))];
                Z_lin=interp1(topo(:,colx),topo(:,colz),X_lin,'linear','extrap');
                [dist_curv,X_curv]=arclength(X_lin,Z_lin);
                X_curv=[0;cumsum(X_curv)];
                X = round(xsca*interp1(X_curv,X_lin,X_raw,'linear','extrap'))/xsca;
                Z = round(xsca*interp1(X_lin,Z_lin,X,'linear','extrap'))/xsca;
                [dist_curv2,X_curv2]=arclength(X,Z);
                X_curv2=[0;cumsum(X_curv2)];
                fprintf(['\n   Curvilinear distance = ',num2str(round(dist_curv2*100)/100),' m\n']);
                
                for i=1:length(X)
                    sxall(sxall==X_raw(i))=X(i);
                    gxall(gxall==X_raw(i))=X(i);
                end
            end
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
if flagtopo>0
    plot_curv(2,X,Z,[],'.-','k',1,1,0,0,14,'X (m)','Altitude (m)',[],[min(X) max(X)],[min(Z)-0.01*min(Z) max(Z)+0.01*max(Z)],...
        [],[],[],[],[],[],[2 20 24 7],[],1);
    hold on;
    % plot(X(72),Z(72),'ro');
end
%%
for i=1:length(Z)
    szall(sxall==X(i))=Z(i);
    gzall(gxall==X(i))=Z(i);
end

dx=mean(diff(unique(gxall)))*xsca;

dlmwrite('infile.txt',[fldrall,tracfall,gxall*xsca,(gzall+zshift)*xsca,...
    sxall*xsca,(szall+zshift)*xsca],'delimiter',' ','precision',6);
unix('a2b < infile.txt > infile.bin');
unix(['sushw < ',rawsufile,' key=scalel,scalco,gdel a=',...
    num2str(-xsca),',',num2str(-xsca),',',num2str(dx),...
    ' | sushw key=fldr,tracf,gx,gelev,sx,selev infile=infile.bin > ',sufile]);

if exist('newrectime','var')==1 && isempty(newrectime)==0
    nsnew=newrectime/dt;
    unix(['suvlength < ',sufile,' ns=',num2str(nsnew),' > tmp.su']);
    unix(['mv tmp.su ',sufile]);
end

delete('infile.txt','infile.bin',rawsufile);

if exist('delay','var')==1 && ~isempty(delay)
    [~,~]=unix(sprintf('sushw < %s key=delrt a=%1.3f > %s',sufile,delay*1000,'tmp.su'));
    unix(['mv tmp.su ',sufile]);
end

sufile=strrep(sufile,'\','\\');
fprintf(['\n   Shots saved in ',sufile,'\n']);
