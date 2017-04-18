%%% SURFACE-WAVE dispersion INVERSION & PROFILING (SWIP)
%%% MODULE A : SWIPdisp.m
%%% S. Pasquet - V17.04.17
%%% SWIPdisp.m performs windowing and stacking of surface-wave dispersion
%%% It allows to pick dispersion curves and save figures of dispersion, 
%%% spectrograms and shot gathers

%%% This module calls the following Seismic Unix native functions:
%%% suwind - sugethw - sushw - suchw - sufilter - suop2 - suop - susort -
%%% sustack - suspecfx - sustrip - sumute - suflip
%%% It also calls the following third party SU function:
%%% supomegal
%%% Some SU functions are called through the following MATLAB functions:
%%% get_acquiparam - matpomegal - matwind - matspecfx - matop
%%% dsp2dat - spec2dat - seismo2dat
%%% It also calls the Linux code pdfcrop if correctly installed

%%%-------------------------%%%
%%% START OF INITIALIZATION %%%

run('SWIP_defaultsettings') % Read default settings

% Check main option
if (exist('calc','var')==1 && isempty(calc)==1) || exist('calc','var')==0
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n   Please choose a calculation option');
    fprintf('\n   calc = 0 to select an existing SWIP folder');
    fprintf('\n   calc = 1 to extract dispersion images from SU file');
    fprintf('\n   calc = 2 to import dispersion curves');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
    return
end

% Initialize folders
if calc==1
    if (exist('nWmin','var')==1 && isempty(nWmin)==0) && (exist('nWmax','var')==1 && isempty(nWmax)==0)
        nWvec=nWmin:2:nWmax;
    end
    nWmin=min(nWvec);
    nWmax=max(nWvec);
    dir_all=dir_create(1,nWmin,nWmax,dW,dSmin,dSmax,side); % Create subproject folder
    sizeax=[]; f=[]; v=[]; tseis=[]; fspec=[];
    if plotspec==1
        plotdisp=1; % To avoid sizeax problems
    end
elseif calc==0
    dir_all=dir_create(0); % Select existing subproject folder
    if dir_all.dir_main==0
        return
    end
elseif calc==2
    % Get imported dispersion curves extraction settings to create subproject folder
    answer1=inputdlg({'Window size (nb of traces)','Inter-geophone spacing (m)',...
        'Window lateral shift (nb of traces)','Spatial scaling xsca'},...
        'Acquisition Settings',1,{'31','1','1','100'});
    if isempty(answer1)==1 || isnan(str2double(answer1(1)))==1 || ...
            isnan(str2double(answer1(1)))==1 || isnan(str2double(answer1(3)))==1
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        fprintf('\n   Please provide all requested parameters');
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
        return
    end
    dx=str2double(answer1(2));
    acquiparam.dx=dx;
    nWmin=str2double(answer1(1));
    nWmax=nWmin;
    dW=str2double(answer1(3));
    xsca=str2double(answer1(4));
    dSmin=1; dSmax=1; side='imported';
    dir_all=dir_create(1,nWmin,nWmax,dW,dSmin,dSmax,side);
    sizeax=[]; f=fmin:1:fmax; v=0:10:vmax; fspec=[]; tseis=[];
    sufile='imported';
else
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n   Please select a valid option for calc');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
    return
end
% Read subfolders in subproject folder
dir_start=dir_all.dir_start;
dir_dat=dir_all.dir_dat;
dir_pick=dir_all.dir_pick;
dir_targ=dir_all.dir_targ;
dir_img=dir_all.dir_img;
dir_img_xmid=dir_all.dir_img_xmid;

% Create them if not existing
dir_img_disp=fullfile(dir_img_xmid,'disp');
if exist(dir_img_disp,'dir')~=7 && plotdisp==1 && strcmp(side,'imported')==0
    mkdir(dir_img_disp);
end
dir_img_pick=fullfile(dir_img_xmid,'disp_pick');
if exist(dir_img_pick,'dir')~=7 && plotpckdisp==1
    mkdir(dir_img_pick);
end
dir_img_spec=fullfile(dir_img_xmid,'spectro');
if exist(dir_img_spec,'dir')~=7 && plotspec==1 && strcmp(side,'imported')==0
    mkdir(dir_img_spec);
end
dir_img_seismo=fullfile(dir_img_xmid,'seismo');
if exist(dir_img_seismo,'dir')~=7 && plotseismo==1 && strcmp(side,'imported')==0
    mkdir(dir_img_seismo);
end
dir_img_single=fullfile(dir_img_xmid,'prestack');
if exist(dir_img_single,'dir')~=7 && plotsingle==1 && strcmp(side,'imported')==0
    mkdir(dir_img_single);
end
dir_img_stkdisp=fullfile(dir_img_xmid,'synstack');
if exist(dir_img_stkdisp,'dir')~=7 && plotstkdisp==1 && strcmp(side,'imported')==0
    mkdir(dir_img_stkdisp);
end

% Check if targets already exist when dispersion curves plots are active
if (plotpckdisp==1 || plot1dobs==1 || plot2dobs==1) && target==0
    targstruct=dir(fullfile(dir_targ,'*.target'));
    if isempty(targstruct)==0
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        fprintf('\n       No target files in file.targ');
        fprintf('\n   No dispersion curves will be plotted');
        fprintf('\n       Re-run script with target = 1');
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
    end
end

% Check acquisition parameters
if calc==1
    sustruct=dir(fullfile(dir_start,'*.su')); % Get all .SU file in main folder
    if length(sustruct)>1
        fprintf('\n  Multiple SU files in working folder - Select correct SU file\n');
        sufile=uigetfile('*.su','Select SU file to use');
        if sufile==0
            fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!');
            fprintf('\n   Please provide SU file');
            fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!\n\n');
            return
        end
    elseif isempty(sustruct)==1
        fprintf('\n  !!!!!!!!!!!!!!!!!!');
        fprintf('\n   No SU file found');
        fprintf('\n  !!!!!!!!!!!!!!!!!!\n\n');
        return
    else
        sufile=sustruct.name;
    end
    matfile=fullfile(dir_dat,[sufile,'.param.mat']); % .mat file containing subproject parameters
    if exist(matfile,'file')==2
        load(matfile); % Read existing parameters from .mat file
        dir_all=dir_create(1,nWmin,nWmax,dW,dSmin,dSmax,side);
    end
    acquiparam=get_acquiparam(sufile,xsca); % Get acquisition parameters from .SU file
    xsca=acquiparam.xsca; % Scaling factor
    dx=acquiparam.dx; % Mean inter-geophone spacing (m)
    if dx==0
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        fprintf('\n   Please provide inter-geophone spacing');
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
        return
    end
    topo=acquiparam.topo; % Get topography from SU file (X,Z in m)
elseif calc==0
    % Read stack and p-w parameters from .mat file
    matstruct=dir(fullfile(dir_dat,'*.param.mat'));
    matfile=fullfile(dir_dat,matstruct.name); % .mat file containing subproject parameters
    try
        load(matfile); % Read existing parameters from .mat file
    catch
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        fprintf('\n   Missing .mat file in file.dat folder');
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
        return
    end
    % Read parameters
    dx=acquiparam.dx; % Mean inter-geophone spacing (m)
    nWmin=stackdisp.nWmin; % Min. window size (no. of traces)
    nWmax=stackdisp.nWmax; % Max. window size (no. of traces)
    dW=stackdisp.dW; % Shift between two successive windows (no. of traces)
    dSmin=stackdisp.dSmin; % Min. offset between window and shot (no. of traces)
    dSmax=stackdisp.dSmax; % Max. offset between window and shot (no. of traces)
    side=stackdisp.side; % Sources side (B=both, L=left, R=right)
    nray=pomega.nray; % Number of velocity samples for p-omega transform
    fmin=pomega.fmin; % Min. frequency for p-omega transform (Hz)
    fmax=pomega.fmax; % Max. frequency for p-omega transform (Hz)
    vmin=pomega.vmin; % Min. velocity for p-omega transform (m/s)
    vmax=pomega.vmax; % Max. velocity for p-omega transform (m/s)
    xsca=pomega.xsca; % Scaling factor
    flip=pomega.flip; % Flip .SU file
    sizeax=plotopt.sizeax; % Previous axis size (to keep same X axis dimension between plots)
    flim=xmidparam.flim; % Minimum frequency defined from amplitude threshold on spectrogram
    if strcmp(side,'imported')==1
        freqlim=0;
    end
    if freqlim==0
        flim=flim*NaN;
    end
    if plotspec==1 && isempty(sizeax)==1
        plotdisp=1; % To avoid sizeax problems
    end
elseif calc==2
    matfile=fullfile(dir_dat,'imported.param.mat'); % .mat file containing subproject parameters
    if exist(matfile,'file')==2
        load(matfile); % Read existing parameters from .mat file
    end
    acquiparam.dx=dx;
    plotdisp=0; plotspec=0; plotseismo=0;
    plotsingle=0; plotstkdisp=0;
    freqlim=0; pick=0; plotflim=0;
    acquiparam.dt=[]; acquiparam.Gx=[]; acquiparam.Gz=[];
    acquiparam.Sx=[]; acquiparam.Sz=[]; acquiparam.NGx=[];
    acquiparam.NSx=[]; acquiparam.Gxsing=[]; acquiparam.Gzsing=[];
    acquiparam.Sxsing=[]; acquiparam.Szsing=[];
end

% XmidT position initialization
xmidformat=['%12.',num2str(log(xsca)/log(10)),'f']; % Precision
if calc~=2 % From SU file or existing subproject
    dt=acquiparam.dt; % Sampling interval (s)
    Gxsing=acquiparam.Gxsing; % Single geophones positions
    Sxsing=acquiparam.Sxsing; % Single sources positions
    if calc==1 % From SU file
        xmin=min(Gxsing); % Get starting X coordinate (m)
        xmax=max(Gxsing); % Get ending X coordinate (m)
        winsize=nWvec; % Window sizes
        maxwinsize=(winsize-1)*dx; % Max. window size
        nwin=length(winsize); % Number of windows
        if mod(nWmin,2)==1 % Uneven number of traces
            XmidT=Gxsing(1+(nWmin-1)/2:dW:end-(nWmin-1)/2);
        else % Even number of traces
            XmidT=mean([Gxsing((nWmin)/2:dW:end-(nWmin)/2),...
                Gxsing(1+(nWmin)/2:dW:1+end-(nWmin)/2)],2);
        end
        XmidT=round(XmidT'*xsca)/xsca; % Xmid position vector
        Xlength=length(XmidT); % Number of Xmids
        if exist(matfile,'file')==2
            nshot=xmidparam.nshot; % Number of shots per Xmid
            Gmin=xmidparam.Gmin; Gmax=xmidparam.Gmax; % First and last geophone for each Xmid
            lmaxpick=targopt.lmaxpick; % Max pick wavelength
        else
            nshot=zeros(Xlength,nwin);
            Gmin=NaN*nshot; Gmax=Gmin;
            lmaxpick=zeros(Xlength,1)*NaN;
        end
    else % From existing subproject
        winsize=xmidparam.winsize; maxwinsize=xmidparam.maxwinsize;
        nwin=xmidparam.nwin; nshot=xmidparam.nshot;
        XmidT=xmidparam.XmidT; Xlength=xmidparam.Xlength;
        Gmin=xmidparam.Gmin; Gmax=xmidparam.Gmax;
        flim=xmidparam.flim; zround=xmidparam.zround;
        lmaxpick=targopt.lmaxpick;
    end
    if isempty(dt)==1
        plotdisp=0; plotspec=0; plotseismo=0;
        plotsingle=0; plotstkdisp=0;
        freqlim=0; pick=0; plotflim=0;
    end
else % From imported dispersion curves
    dt=[];
    % Get folder containing dispersion curves
    dir_pick_ext=uigetdir('./','Select folder containing dispersion curves');
    if dir_pick_ext==0
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        fprintf('\n   Please select a folder containing dispersion curves');
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
        return
    end
    dispstruct=dir(fullfile(dir_pick_ext));
    xmidlocal=ones(length(dispstruct)-2,1);
    dspfile_sum='fake';
    if length(dispstruct)==2
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        fprintf('\n   Empty folder - Please select a folder containing dispersion curves');
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
        return
    end
    ii=0;
    % Read dispersion curves
    for ip=1:length(dispstruct)
        if strcmp(dispstruct(ip).name,'.')==1 || strcmp(dispstruct(ip).name,'..')==1
            continue
        else
            ii=ii+1;
        end
        dispfile=dispstruct(ip).name;
        [~,~,extension]=fileparts(dispfile);
        if strcmp(extension,'.pvc')==1 % Read dispersion curves in SWIP .pvc format
            try
                xmidlocal(ii)=str2double(dispfile(1:strfind(dispfile,'.M')-1));
                mi=str2double(dispfile(strfind(dispfile,'.M')+2:strfind(dispfile,'.pvc')-1));
                load(fullfile(dir_pick_ext,dispfile));
                if strcmp([dir_pick_ext,'/'],fullfile(dir_start,dir_pick))==0
                    newfilename=[num2str(xmidlocal(ii),xmidformat),'.M',num2str(mi),'.pvc'];
                    copyfile(fullfile(dir_pick_ext,dispfile),fullfile(dir_pick,newfilename));
                end
            catch
                fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
                fprintf('\n   Wrong format for dispersion curve - Go to next file');
                fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
                continue
            end
        else % Read dispersion curves in any 3-(or 2-)column ASCII file
            try
                load(fullfile(dir_pick_ext,dispfile));
                answer1=inputdlg({'Xmid position (m)','Mode number (fundamental=0)'},...
                    [dispfile,' settings'],1);
                if isempty(answer1)==1 || isnan(str2double(answer1(1)))==1 || ...
                        isnan(str2double(answer1(1)))==1
                    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
                    fprintf('\n   Please provide Xmid position');
                    fprintf('\n         and mode number');
                    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
                    return
                end
                xmidlocal(ii)=str2double(answer1(1));
                mi=str2double(answer1(2));
                newfilename=[num2str(xmidlocal(ii),xmidformat),'.M',num2str(mi),'.pvc'];
                copyfile(fullfile(dir_pick_ext,dispfile),fullfile(dir_pick,newfilename));
            catch
                fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
                fprintf('\n   Wrong format for dispersion curve - Go to next file');
                fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
                continue
            end
        end
    end
    XmidT_import=unique(xmidlocal(isnan(xmidlocal)==0)); % Xmid position vector
    Xlength_import=length(unique(xmidlocal(isnan(xmidlocal)==0))); % Number of Xmids
    if Xlength_import==0
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        fprintf('\n   No valid dispersion curves in the folder');
        fprintf('\n         Check file format and retry');
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
        return
    end
    XmidT=[];
    for ii=2:Xlength_import
        XmidT=[XmidT,XmidT_import(ii-1):dW*dx:XmidT_import(ii)-0.5*dW*dx];
        if ii==Xlength_import
            XmidT=[XmidT,XmidT_import(end)];
        end
    end
    XmidT=unique(XmidT);
    Xlength=length(XmidT);
    xmin=min(XmidT)-dx; xmax=max(XmidT)+dx;
    nshot=ones(Xlength,1);
    Gmin=nshot; Gmax=nshot; flim=nshot-1;
    nshot(ismember(XmidT,XmidT_import)==0)=-1;
    winsize=nWmin; maxwinsize=(winsize-1)*dx;
    nwin=length(winsize);
    if exist(matfile,'file')==2
        lmaxpick=targopt.lmaxpick;
    else
        lmaxpick=zeros(Xlength,1)*NaN;
    end
end
% Select Xmids
if exist('Xmidselec','var')~=1 || isempty(Xmidselec)==1
    Xmidselec=1:Xlength;
end
if max(Xmidselec)>Xlength
    Xmidselec=Xmidselec(Xmidselec<=Xlength);
end

% Bandpass filter with SU
if filt==1 && calc==1
    sufilefilt=[sufile,'.filt'];
    com1=sprintf('sufilter < %s f=%d,%d,%d,%d amps=0,1,1,0 > %s',...
        sufile,fcutlow,fcutlow+taper,fcuthigh-taper,fcuthigh,sufilefilt);
    unix(com1);
    sufileOK=sufilefilt;
else
    if calc~=2
        sufileOK=sufile;
    end
end

% Topography and cut-off frequency initialization
if calc~=0
    if calc==1
        if exist(matfile,'file')==2
            flim=xmidparam.flim;
        else
            flim=zeros(Xlength,1);
        end
    elseif calc==2
        [topofile,topopath]=uigetfile('*','Select topo file (X,Z) or cancel for flat topo');
        if sum(topofile)~=0
            topo=load(fullfile(topopath,topofile));
        else
            topo(:,1)=XmidT;
            topo(:,2)=zeros(Xlength,1);
        end
        acquiparam.topo=topo;
    end
    z=interp1(topo(:,1),topo(:,2),XmidT,'linear','extrap'); % Interpolate topo for all Xmid
    zround=fix(z*100)/100; % Keep only 2 decimals for topo
    xmidparam.zround=zround;
end

% Colormap and colorbar initialization for curve plot
if plot1dobs==1
    cticks=XmidT(1:ceil(Xlength/5):end); % Colorbar ticks
    ccurve=flipud(autumn(Xlength)); % Colormap for dispersion curve position
end

% Initialization of the maximum number of modes
if isempty(maxmodeinv)==1
    maxmodeinv=10;
end

% Initialization of phase velocity pseudo-section
if plot2dobs==1
    vph2dobs=cell(maxmodeinv+1,1);
    for ip=1:maxmodeinv+1
        vph2dobs{ip}=zeros(length(resampvec),Xlength)*NaN;
    end
end

% Save settings in .mat file
if calc~=0
    % Store parameters in structure
    stackdisp=struct('nWmin',nWmin,'nWmax',nWmax,'dW',dW,...
        'dSmin',dSmin,'dSmax',dSmax,'side',side,'xmidformat',xmidformat);
    pomega=struct('nray',nray,'fmin',fmin,'fmax',fmax,...
        'vmin',vmin,'vmax',vmax,'xsca',xsca,'tsca',tsca,'flip',flip);
    xmidparam=struct('winsize',winsize,'maxwinsize',maxwinsize,...
        'nwin',nwin,'XmidT',XmidT,'Xlength',Xlength,'nshot',nshot,...
        'Gmin',Gmin,'Gmax',Gmax,'flim',flim,'zround',zround);
    targopt=struct('sampling',sampling,'resampvec',resampvec,'wave',wave,'lmaxpick',lmaxpick);
    filtmute=struct('filt',filt,'fcutlow',fcutlow,'fcuthigh',fcuthigh,...
        'taper',taper,'mute',mute,'tmin1',tmin1,'tmin2',tmin2,...
        'tmax1',tmax1,'tmax2',tmax2);
    % Save param in .mat file
    if exist(matfile,'file')==2
        save(matfile,'-append','acquiparam','stackdisp','pomega','sufile','xmidparam','targopt','filtmute');
    else
        save(matfile,'acquiparam','stackdisp','pomega','sufile','xmidparam','targopt','filtmute');
    end
end

% Save target settings
if target==1
    targopt.sampling=sampling;
    targopt.resampvec=resampvec;
    targopt.wave=wave;
    % Save param in .mat file
    save(matfile,'-append','targopt');
end

% Colorbar position
if cb_disp==1 && cbpos==2
    cb_disp=2;
end

fprintf('\n  **********************************************************');
fprintf('\n  **********************************************************\n');

%%% END OF INITIALIZATION %%%
%%%-----------------------%%%

%% CALCULATIONS FOR ALL XMIDS

%%%%%% Loop over all Xmids %%%%%%

Tstart=tic; % Start clock
i=0; % Xmid number flag
while i<length(Xmidselec)
    i=i+1; ix=Xmidselec(i);
    
    %% %% %%
    if sum(nshot(ix,:))>=0
        fprintf(['\n  Xmid',num2str(ix),' = ',num2str(XmidT(ix),xmidformat),' m \n']);
    end
    
    if calc~=2
        
        %%%%% Folder and file name initialization %%%%%
        
        % Create folder to store intermediate files for each Xmid
        dir_dat_xmid=fullfile(dir_dat,['Xmid_',num2str(XmidT(ix),xmidformat)]);
        if exist(dir_dat_xmid,'dir')~=7 && (calc==1 || plotsingle==1 || plotstkdisp==1)
            mkdir(dir_dat_xmid);
        end
        % Stacked dispersion file name (delete if exists and calc=1)
        dspfile_sum=fullfile(dir_dat,[num2str(XmidT(ix),xmidformat),'.sum.dsp']);
        if exist(dspfile_sum,'file')==2 && (calc==1 || plotstkdisp==1)
            delete(dspfile_sum);
        end
        % Stacked spectrogram file name
        specfile_sum=fullfile(dir_dat,[num2str(XmidT(ix),xmidformat),'.sum.spec']);
        % Example of shot gather file name
        seisfile_sum=fullfile(dir_dat,[num2str(XmidT(ix),xmidformat),'.sum.su']);
        % String containing all spectrogram file names
        specfile_all=[];
        
        if stack==2  % Alternative method for extracting dispersion (experimental)
            % Stacked shot gather file name (delete if exists and calc=1)
            seisfile_sum_new=fullfile(dir_dat,[num2str(XmidT(ix),xmidformat),'.sum_new.su']);
            dspfile_sum_new=fullfile(dir_dat,[num2str(XmidT(ix),xmidformat),'.sum_new.dsp']);
            specfile_sum_new=fullfile(dir_dat,[num2str(XmidT(ix),xmidformat),'.sum_new.spec']);
            if exist(seisfile_sum_new,'file')==2 && (calc==1 || plotstkdisp==1)
                delete(seisfile_sum_new);
            end
            seismofile_left=[];
            seismofile_right=[];
        end
        
        %%%%% Loop over window sizes %%%%%
        
        if calc==1 || plotsingle==1 || plotstkdisp==1
            j=0; % Stack flag
            for jw=1:nwin
                % Retrieve first and last geophone position for the current window
                if mod(winsize(jw),2)==1 % Non-even number of traces
                    Gleft=Gxsing(find(Gxsing<XmidT(ix),(winsize(jw)-1)/2,'last'));
                    Gright=Gxsing(find(Gxsing>XmidT(ix),(winsize(jw)-1)/2));
                    ntr=length(Gleft)+length(Gright)+1;
                    if ntr~=winsize(jw) && (calc==1 || plotsingle==1 || plotstkdisp==1) % Check number of extracted traces
                        fprintf(['\n  Not enough traces with nW = ',num2str(winsize(jw)),...
                            ' - Go to next shot or Xmid\n']);
                        continue
                    end
                    Gmin(ix,jw)=min(Gleft);
                    Gmax(ix,jw)=max(Gright);
                else % Even number of traces
                    Gleft=Gxsing(find(Gxsing<XmidT(ix),(winsize(jw))/2,'last'));
                    Gright=Gxsing(find(Gxsing>XmidT(ix),(winsize(jw))/2));
                    ntr=length(Gleft)+length(Gright);
                    if ntr~=winsize(jw) && (calc==1 || plotsingle==1 || plotstkdisp==1) % Check number of extracted traces
                        fprintf(['\n  Not enough traces with nW = ',num2str(winsize(jw)),...
                            ' - Go to next shot or Xmid\n']);
                        continue
                    end
                    Gmin(ix,jw)=min(Gleft);
                    Gmax(ix,jw)=max(Gright);
                end
                % Retrieve min and max sources position on both sides of the window
                Smin=XmidT(ix)-(maxwinsize(jw)/2)-(dx*dSmax+dx);
                Smax=XmidT(ix)+(maxwinsize(jw)/2)+(dx*dSmax+dx);
                Smed1=XmidT(ix)-(maxwinsize(jw)/2)-(dx*dSmin);
                Smed2=XmidT(ix)+(maxwinsize(jw)/2)+(dx*dSmin);
                % Select existing sources according to the specified side
                if strcmp(side,'L')==1
                    Sselec=Sxsing(Sxsing>Smin & Sxsing<=Smed1);
                elseif strcmp(side,'R')==1
                    Sselec=Sxsing(Sxsing<Smax & Sxsing>=Smed2);
                else
                    Sselec=Sxsing((Sxsing>Smin & Sxsing<=Smed1) | ...
                        (Sxsing<Smax & Sxsing>=Smed2));
                end
                leftOK=0;
                % Get nb of selected shot for the current window
                nshot(ix,jw)=length(Sselec);
                if mod(winsize(jw),2)~=mod(winsize(1),2) % Check number of extracted traces
                    fprintf(['\n  Window sizes should all be even or uneven - Discard nW = ',num2str(winsize(jw)),'\n']);
                    nshot(ix,jw)=0;
                    continue
                end
                if (calc==1 || plotsingle==1 || plotstkdisp==1)
                    fprintf(['\n  ',num2str(nshot(ix,jw)),...
                        ' shot(s) with nW = ',num2str(winsize(jw)),'\n']);
                end
                
                %%%% Loop over all selected shots %%%%
                
                for ks=1:nshot(ix,jw)
                    % Windowing, muting with SU and saving shot gather in .su file
                    seismofile=fullfile(dir_dat_xmid,[num2str(XmidT(ix),xmidformat),'.',...
                        num2str(winsize(jw)),'.',num2str(Sselec(ks)),'.su']);                    
                    [seismomat,xseis,tseis,ntr]=matwind(sufileOK,Sselec(ks),Gmin(ix,jw),Gmax(ix,jw),xsca,...
                        winsize(jw),seismofile,0,mute,tmin1,tmin2,tmax1,tmax2);
                    
                    if stack==2 % Alt. method (xp)
                        if ((strcmp(side,'L')==1 || strcmp(side,'B')==1) && XmidT(ix)>Sselec(ks))
                            seismofile_left=[seismofile_left,' ',seismofile];
                        elseif ((strcmp(side,'R')==1 || strcmp(side,'B')==1) && XmidT(ix)<Sselec(ks))
                            seismofile_right=[seismofile_right,' ',seismofile];
                        end
                    end
                    
                    if ntr~=winsize(jw) % Check number of extracted traces
                        fprintf('\n  Not enough traces - Go to next shot\n');
                        nshot(ix,jw)=nshot(ix,jw)-1;
                        continue
                    end
                    j=j+1; % Stack flag
                    
                    % P-Omega transform on shot gather and saving in .dsp file
                    dspfile=fullfile(dir_dat_xmid,[num2str(XmidT(ix),xmidformat),'.',...
                        num2str(winsize(jw)),'.',num2str(Sselec(ks)),'.dsp']);
                    [dspmat,f,v]=matpomegal(seismofile,1,nray,fmin,fmax,vmin,vmax,...
                        flip,xsca,tsca,1,dspfile,0);
                    % Spectrogram calculation on shot gather and saving in .spec file
                    specfile=fullfile(dir_dat_xmid,[num2str(XmidT(ix),xmidformat),'.',...
                        num2str(winsize(jw)),'.',num2str(Sselec(ks)),'.spec']);
                    [specmat,fspec,xspec]=matspecfx(seismofile,xsca,specfile,0);
                    specfile_all=[specfile_all,' ',specfile]; % Concat. file names
                    % Get example seismogram (shortest offset available, left side in priority)
                    if strcmp(side,'L')==1
                        seismofileOK = seismofile;
                    elseif strcmp(side,'R')==1 && j==1
                        seismofileOK = seismofile;
                    elseif strcmp(side,'B')==1
                        if XmidT(ix)>Sselec(ks)
                            seismofileOK = seismofile;
                            leftOK=1;
                        elseif XmidT(ix)<Sselec(ks) && j==1 && leftOK==0
                            seismofileOK = seismofile;
                        end
                    end
                                       
                    %%% Plot and save single images for each selected shot %%%
                    
                    if plotsingle==1 && exist(dir_dat_xmid,'dir')==7
                        % Initialize and create folder to store these images
                        dir_img_xmid_single=fullfile(dir_img_single,['Xmid_',num2str(XmidT(ix),xmidformat)]);
                        if exist(dir_img_xmid_single,'dir')~=7
                            mkdir(dir_img_xmid_single);
                        end
                        fprintf(['\n  Plot and save single images for shot ' num2str(j) '\n']);
                        if plotflim==1
                            % Search for cut-off frequency in single spectrogram
                            specstruct=dir(fullfile(dir_dat_xmid,[num2str(XmidT(ix),xmidformat),'.',...
                                num2str(winsize(jw)),'.',num2str(Sselec(ks)),'.spec']));
                            flimsing=fmin_search(specstruct,dir_dat_xmid,dt,specampmin,fminpick);
                        else
                            flimsing=[];
                        end
                        % Plot single dispersion image
                        if plotdisp==1
                            if Dlogscale==0
                                fig1=plot_img(showplot,f,v,dspmat',flipud(map0),axetop,axerev,cb_disp,fs,...
                                    freqtitle_long,'Phase velocity (m/s)',...
                                    'Norm. ampli.',[fMIN fMAX],[VphMIN VphMAX],...
                                    [],fticks,Vphticks,[],[],flimsing,[],[0 0 24 18],[]);
                            else
                                dspmatinv=1./(1-dspmat);
                                dspmatinv(isinf(dspmatinv))=max(max(dspmatinv(isinf(dspmatinv)==0)));
                                fig1=plot_img_log(showplot,f,v,dspmatinv',flipud(map0),axetop,axerev,cb_disp,fs,...
                                    freqtitle_long,'Phase velocity (m/s)','1/(1-Norm. ampli.)',...
                                    [fMIN fMAX],[VphMIN VphMAX],[1 length(map0)],fticks,Vphticks,...
                                    [],[],flimsing,[],[0 0 24 18],[]);
                            end
                            hold on
                            if Flogscale==1
                                set(gca,'xscale','log');
                            end
                            sizeax=get(findobj(fig1,'Type','Axes'),'Position');
                            if cb_disp==1
                                sizeax=sizeax{2};
                            end
                            file1=fullfile(dir_img_xmid_single,[num2str(XmidT(ix),xmidformat),...
                                '.',num2str(winsize(jw)),'.',num2str(Sselec(ks)),...
                                '.disp.',imgform]);
                            save_fig(fig1,file1,imgform,imgres,1);
                            close(fig1)
                        end
                        % Plot single spectrogram image
                        if plotspec==1
                            fig2=plot_img(showplot,fspec,xspec,specmat,flipud(map0),axetop,axerev,cb_disp,fs,...
                                freqtitle_long,'Gx (m)','Norm. ampli.',[fMIN fMAX],[min(xspec) max(xspec)],...
                                [],[],[],[],[],[],[],[0 0 24 18],[],[],0);
                            hold on
                            if Flogscale==1
                                set(gca,'xscale','log');
                            end
                            set(findobj(fig2,'Type','Axes'),'ActivePositionProperty','Position');
                            if cb_disp==1
                                axeok=findobj(fig2,'Type','Axes');
                                set(axeok(2),'position',[sizeax(1),sizeax(2),sizeax(3),sizeax(4)/3]);
                            else
                                set(findobj(fig2,'Type','Axes'),'position',...
                                    [sizeax(1),sizeax(2),sizeax(3),sizeax(4)/3]);
                            end
                            if isempty(flimsing)==0
                                yL=get(gca,'YLim');
                                han3=dashline([flimsing flimsing],yL,3,3,3,3,'color',[1 0 0],'linewidth',5);
                            end
                            file2=fullfile(dir_img_xmid_single,[num2str(XmidT(ix),xmidformat),...
                                '.',num2str(winsize(jw)),'.',num2str(Sselec(ks)),...
                                '.spec.',imgform]);
                            save_fig(fig2,file2,imgform,imgres,1);
                            close(fig2)
                        end
                        % Plot single shot gather image
                        if plotseismo==1
                            fig3=plot_wiggle(showplot,-seismomat',xseis,tseis*1000,...
                                1,1,99,fs,'Gx (m)','Time (ms)',[],[tMIN tMAX],[],tticks,[0 0 18 24],[]);
                            file3=fullfile(dir_img_xmid_single,[num2str(XmidT(ix),xmidformat),...
                                '.',num2str(winsize(jw)),'.',num2str(Sselec(ks)),...
                                '.seismo.',imgform]);
                            save_fig(fig3,file3,imgform,imgres,1);
                            close(fig3)
                        end
                    end
                    
                    %%% Dispersion stacking calculation %%%
                    
                    % Create zeros .dsp file for first iteration stacking
                    if exist(dspfile_sum,'file')~=2 && (calc==1 || plotstkdisp==1)
                        com1=sprintf('suop2 %s %s op=diff > %s',dspfile,dspfile,dspfile_sum);
                        unix(com1);
                    end
                    
                    if calc==1 || plotstkdisp==1
                        % Stack current dispersion image with previous stack file
                        dspfile_sum_inter=[dspfile_sum,'.new'];
                        com1=sprintf('suop2 %s %s op=sum > %s',dspfile_sum,dspfile,dspfile_sum_inter);
                        unix(com1);
                        movefile(dspfile_sum_inter,dspfile_sum)
                        
                        % Plot and save intermediate stack
                        if  plotstkdisp==1
                            dspfile_sum_stack=[dspfile_sum,'.stack'];
                            com1=sprintf('suop < %s op=norm > %s',dspfile_sum,dspfile_sum_stack);
                            unix(com1);
                            [dspmat,f,v]=dsp2dat(dspfile_sum_stack,flip,0);
                            delete(dspfile_sum_stack);
                            dir_img_xmid_stack=fullfile(dir_img_stkdisp,['Xmid_',num2str(XmidT(ix),xmidformat)]);
                            if exist(dir_img_xmid_stack,'dir')~=7
                                mkdir(dir_img_xmid_stack);
                            end
                            if plotflim==1
                                % Search for cut-off frequency all stacked spectrograms
                                specstruct=dir(fullfile(dir_dat_xmid,...
                                    [num2str(XmidT(ix),xmidformat),'.*.spec']));
                                flimsing=fmin_search(specstruct,dir_dat_xmid,dt,specampmin,fminpick);
                            else
                                flimsing=[];
                            end
                            fprintf(['\n  Plot and save intermediate stack ',num2str(j),'\n']);
                            if Dlogscale==0
                                fig1=plot_img(showplot,f,v,dspmat',flipud(map0),axetop,axerev,cb_disp,fs,...
                                    freqtitle_long,'Phase velocity (m/s)',...
                                    'Norm. ampli.',[fMIN fMAX],[VphMIN VphMAX],...
                                    [],fticks,Vphticks,[],[],flimsing,[],[0 0 24 18],[]);
                            else
                                dspmatinv=1./(1-dspmat);
                                dspmatinv(isinf(dspmatinv))=max(max(dspmatinv(isinf(dspmatinv)==0)));
                                fig1=plot_img_log(showplot,f,v,dspmatinv',flipud(map0),axetop,axerev,cb_disp,fs,...
                                    freqtitle_long,'Phase velocity (m/s)',...
                                    '1/(1-Norm. ampli.)',[fMIN fMAX],[VphMIN VphMAX],...
                                    [1 length(map0)],fticks,Vphticks,[],[],flimsing,[],[0 0 24 18],[]);
                            end
                            hold on
                            if Flogscale==1
                                set(gca,'xscale','log');
                            end
                            sizeax=get(findobj(fig1,'Type','Axes'),'Position');
                            if cb_disp==1
                                sizeax=sizeax{2};
                            end
                            file1=fullfile(dir_img_xmid_stack,[num2str(XmidT(ix),xmidformat),...
                                '.stack',num2str(j),'.disp.',imgform]);
                            save_fig(fig1,file1,imgform,imgres,1);
                            close(fig1)
                        end
                    end
                end
            end
            
            if sum(nshot(ix,:))>0 % Check if there is enough shots to extract dispersion for current Xmid
                com1=sprintf('cat %s > %s',specfile_all,specfile_sum);
                unix(com1); % Concatenate all spectrograms into composite spectrogram
                com1=sprintf('susort < %s +gx | sustack key=gx > %s',specfile_sum,fullfile(dir_dat_xmid,'sort.spec'));
                unix(com1); % Sort composite spectrogram by offset and stack common offset traces
                movefile(fullfile(dir_dat_xmid,'sort.spec'),specfile_sum); % Rename file
                
                if stack==2 % Alt. method (xp)
                    % Alternative method consists in concatenating all traces from all selected shots, 
                    % then rearrange them by offset, to finally transform the wavefield of the
                    % composite shot gather
                    if isempty(seismofile_left)==0
                        com1=sprintf('cat %s > %s',seismofile_left,fullfile(dir_dat_xmid,'cat_left.su'));
                        unix(com1); % Concatenate left shots
                    end
                    if isempty(seismofile_right)==0
                        com1=sprintf('cat %s > %s',seismofile_right,fullfile(dir_dat_xmid,'cat_right.su'));
                        unix(com1); % Concatenate right shots
                        com1=sprintf('suop < %s op=neg > %s',fullfile(dir_dat_xmid,'cat_right.su'),fullfile(dir_dat_xmid,'cat_right_neg.su'));
                        unix(com1); % Reverse polarity
                    end
                    if isempty(seismofile_left)==0 && isempty(seismofile_right)==0
                        com1=sprintf('cat %s %s > %s',fullfile(dir_dat_xmid,'cat_left.su'),fullfile(dir_dat_xmid,'cat_right_neg.su'),fullfile(dir_dat_xmid,'cat.su'));
                        unix(com1); % Concatenate all shots
                        delete(fullfile(dir_dat_xmid,'cat_left.su'),fullfile(dir_dat_xmid,'cat_right.su'),fullfile(dir_dat_xmid,'cat_right_neg.su'));
                    elseif isempty(seismofile_left)==0 && isempty(seismofile_right)==1
                        movefile(fullfile(dir_dat_xmid,'cat_left.su'),fullfile(dir_dat_xmid,'cat.su'));
                    elseif isempty(seismofile_left)==1 && isempty(seismofile_right)==0
                        movefile(fullfile(dir_dat_xmid,'cat_right_neg.su'),fullfile(dir_dat_xmid,'cat.su'));
                        delete(fullfile(dir_dat_xmid,'cat_right.su'));
                    end
                    com1=sprintf('susort < %s +offset | sustack key=offset | sushw key=sx,fldr,gelev,selev a=-1,1,0,0 | suchw key1=gx key2=offset > %s',...
                        fullfile(dir_dat_xmid,'cat.su'),seisfile_sum_new);
                    unix(com1); % Sort and stack common traces + add headers
                    delete(fullfile(dir_dat_xmid,'cat.su'));
                    
                    % P-Omega transform on composite shot gather and saving in .dsp file
                    [dspmat_new,f_new,v_new]=matpomegal(seisfile_sum_new,1,nray,fmin,fmax,vmin,vmax,...
                        flip,xsca,tsca,1,dspfile_sum_new,0);
                    % Spectrogram calculation on composite shot gather and saving in .spec file
                    [specmat_new,fspec_new,xspec_new]=matspecfx(seisfile_sum_new,xsca,specfile_sum_new,0);
                end
            end
        end
    end
    
    %% %% %%
    
    %%%%% End of main loops %%%%%
    
    if (sum(nshot(ix,:))>0 && exist(dspfile_sum,'file')==2) || isempty(dt)==1
        % Normalize final dispersion image
        if calc==1 || plotstkdisp==1
            matop(dspfile_sum,'norm',flip);
            if stack==2
                matop(dspfile_sum_new,'norm',flip);
            end
            copyfile(seismofileOK,seisfile_sum); % Example of shot gather
        end
        
        % Search for cut-off frequency on composite spectrogram
        if (calc==1 || target==1 || pick==1 || plotdisp==1 || plotpckdisp==1 || plotspec==1) && freqlim==1
            if stack==2 % Alt. method (xp)
                specstruct=dir(specfile_sum_new);
            else
                specstruct=dir(specfile_sum);
            end
            [flim(ix),~]=fmin_search(specstruct,dir_dat,dt,specampmin,fminpick);
        end
        
        if plotflim==1
            flimsing=flim(ix);
        else
            flimsing=[];
        end
        
        if stack==2 % Alt. method (xp)
            [~,offsets]=unix(sprintf('sugethw < %s key=offset output=geom',seisfile_sum_new));
            offsets=str2num(offsets)/xsca;
            noffsets=length(offsets);
            nWmin=noffsets; nWmax=noffsets;
        end
        
        %% %% %%
        
        %%%%% Plot and save final images %%%%%
        
        % Read dispersion .dsp file
        if pick>0 || plotdisp==1 || plotpckdisp==1
            if exist(dspfile_sum,'file')==2
                [dspmat,f,v]=dsp2dat(dspfile_sum,flip,0);
                if pick==1 || pick==2
                    % Downsample dispersion image to dvmin m/s in velocity to speed up display when picking
                    % (also decrease min. resolution in velocity of dispersion curve)
                    dspmat2=dspmat; v2=v;
                    while min(diff(v2))<dvmin
                        v2=v2(1:2:end);
                        dspmat2=dspmat2(:,1:2:end);
                    end
                end
            end
            if stack==2 % Alt. method (xp)
                if exist(dspfile_sum_new,'file')==2
                    [dspmat_new,f_new,v_new]=dsp2dat(dspfile_sum_new,flip,0);
                    if pick==1 || pick==2
                        % Downsample dispersion image to dvmin m/s in velocity to speed up display when picking
                        % (also decrease min. resolution in velocity of dispersion curve)
                        dspmat2_new=dspmat_new; v2_new=v_new;
                        while min(diff(v2_new))<dvmin
                            v2_new=v2_new(1:2:end);
                            dspmat2_new=dspmat2_new(:,1:2:end);
                        end
                    end
                end
            end  
        end
        
        % Plot and save stacked dispersion image
        if  plotdisp==1
            fprintf('\n  Plot and save stacked dispersion image\n');
            if Dlogscale==0
                fig1=plot_img(showplot,f,v,dspmat',flipud(map0),axetop,axerev,cb_disp,fs,...
                    freqtitle_long,'Phase velocity (m/s)',...
                    'Norm. ampli.',[fMIN fMAX],[VphMIN VphMAX],...
                    [],fticks,Vphticks,[],[],flimsing,[],[0 0 24 18],[]);
            else
                dspmatinv=1./(1-dspmat);
                dspmatinv(isinf(dspmatinv))=max(max(dspmatinv(isinf(dspmatinv)==0)));
                fig1=plot_img_log(showplot,f,v,dspmatinv',flipud(map0),axetop,axerev,cb_disp,fs,...
                    freqtitle_long,'Phase velocity (m/s)',...
                    '1/(1-Norm. ampli.)',[fMIN fMAX],[VphMIN VphMAX],...
                    [1 length(map0)],fticks,Vphticks,[],[],flimsing,[],[0 0 24 18],[]);
            end
            hold on
            if Flogscale==1
                set(gca,'xscale','log');
            end
            sizeax=get(findobj(fig1,'Type','Axes'),'Position');
            if cb_disp==1
                sizeax=sizeax{2};
            end
            file1=fullfile(dir_img_disp,[num2str(XmidT(ix),xmidformat),'.disp.',imgform]);
            save_fig(fig1,file1,imgform,imgres,1);
            close(fig1)
            
            if stack==2 % Alt. method (xp)
                if Dlogscale==0
                    fig1=plot_img(showplot,f_new,v_new,dspmat_new',flipud(map0),axetop,axerev,cb_disp,fs,...
                        freqtitle_long,'Phase velocity (m/s)',...
                        'Norm. ampli.',[fMIN fMAX],[VphMIN VphMAX],...
                        [],fticks,Vphticks,[],[],flimsing,[],[0 0 24 18],[]);
                else
                    dspmatinv_new=1./(1-dspmat_new);
                    dspmatinv_new(isinf(dspmatinv_new))=max(max(dspmatinv_new(isinf(dspmatinv_new)==0)));
                    fig1=plot_img_log(showplot,f_new,v_new,dspmatinv_new',flipud(map0),axetop,axerev,cb_disp,fs,...
                        freqtitle_long,'Phase velocity (m/s)',...
                        '1/(1-Norm. ampli.)',[fMIN fMAX],[VphMIN VphMAX],...
                        [1 length(map0)],fticks,Vphticks,[],[],flimsing,[],[0 0 24 18],[]);
                end
                hold on
                if Flogscale==1
                    set(gca,'xscale','log');
                end
                sizeax=get(findobj(fig1,'Type','Axes'),'Position');
                if cb_disp==1
                    sizeax=sizeax{2};
                end
                file1=fullfile(dir_img_disp,[num2str(XmidT(ix),xmidformat),'.disp_new.',imgform]);
                save_fig(fig1,file1,imgform,imgres,1);
                close(fig1)
            end
            
        end
        
        % Plot and save spectrogram image
        if plotspec==1
            fprintf('\n  Plot and save final spectrogram image\n');
            [specmat,fspec,xspec]=spec2dat(specfile_sum,0);
            xspec=xspec/xsca;
            fig2=plot_img(showplot,fspec,xspec,specmat,flipud(map0),axetop,axerev,cb_disp,fs,...
                freqtitle_long,'Gx (m)','Norm. ampli.',...
                [fMIN fMAX],[min(xspec) max(xspec)],[],[],[],...
                [],[],[],[],[0 0 24 18],[],[],0);
            hold on
            if Flogscale==1
                set(gca,'xscale','log');
            end
            set(findobj(fig2,'Type','Axes'),'ActivePositionProperty','Position');
            if cb_disp==1
                axeok=findobj(fig2,'Type','Axes');
                set(axeok(2),'position',[sizeax(1),sizeax(2),sizeax(3),sizeax(4)/3]);
            else
                set(findobj(fig2,'Type','Axes'),'position',...
                    [sizeax(1),sizeax(2),sizeax(3),sizeax(4)/3]);
            end
            if isempty(flimsing)==0
                yL=get(gca,'YLim');
                han3=dashline([flimsing flimsing],yL,3,3,3,3,'color',[1 0 0],'linewidth',5);
            end
            file2=fullfile(dir_img_spec,[num2str(XmidT(ix),xmidformat),'.spec.',imgform]);
            save_fig(fig2,file2,imgform,imgres,1);
            close(fig2)
            
            if stack==2 % Alt. method (xp)
                [specmat_new,fspec_new,xspec_new]=spec2dat(specfile_sum_new,0);
                xspec_new=xspec_new/xsca;
                fig2=plot_img(showplot,fspec_new,xspec_new,specmat_new,flipud(map0),axetop,axerev,cb_disp,fs,...
                    freqtitle_long,'Offset (m)','Norm. ampli.',[fMIN fMAX],[min(xspec_new) max(xspec_new)],[],[],[],...
                    [],[],[],[],[0 0 24 18],[],[],0);
                hold on
                if Flogscale==1
                    set(gca,'xscale','log');
                end
                set(findobj(fig2,'Type','Axes'),'ActivePositionProperty','Position');
                if cb_disp==1
                    axeok=findobj(fig2,'Type','Axes');
                    set(axeok(2),'position',[sizeax(1),sizeax(2),sizeax(3),sizeax(4)/3]);
                else
                    set(findobj(fig2,'Type','Axes'),'position',...
                        [sizeax(1),sizeax(2),sizeax(3),sizeax(4)/3]);
                end
                if isempty(flimsing)==0
                    yL=get(gca,'YLim');
                    han3=dashline([flimsing flimsing],yL,3,3,3,3,'color',[1 0 0],'linewidth',5);
                end
                file2=fullfile(dir_img_spec,[num2str(XmidT(ix),xmidformat),'.spec_new.',imgform]);
                save_fig(fig2,file2,imgform,imgres,1);
                close(fig2)
            end
            
        end
        
        % Plot and save shot gather image
        if plotseismo==1
            fprintf('\n  Plot and save final shot gather image\n');
            [seismomat,tseis,xseis]=seismo2dat(seisfile_sum,0);
            xseis=xseis/xsca;
            fig3=plot_wiggle(showplot,-seismomat',xseis,tseis*1000,1,1,99,...
                fs,'Gx (m)','Time (ms)',[],[tMIN tMAX],[],tticks,[0 0 18 24],[]);
            file3=fullfile(dir_img_seismo,[num2str(XmidT(ix),xmidformat),'.seismo.',imgform]);
            save_fig(fig3,file3,imgform,imgres,1);
            close(fig3)
            
            if stack==2 % Alt. method (xp)
                [seismomat_new,tseis_new,xseis_new]=seismo2dat(seisfile_sum_new,0);
                xseis_new=xseis_new/xsca;
                fig3=plot_wiggle(showplot,-seismomat_new',xseis_new,tseis_new*1000,1,1,99,...
                    fs,'Offset (m)','Time (ms)',[],[tMIN tMAX],[],tticks,[0 0 18 24],[]);
                file3=fullfile(dir_img_seismo,[num2str(XmidT(ix),xmidformat),'.seismo_new.',imgform]);
                save_fig(fig3,file3,imgform,imgres,1);
                close(fig3)
            end
        end
        
        %% %% %%
        
        %%%%% Pick dispersion curves %%%%%
        
        if pick==1 % Manual picking
            if stack==2 % Alt. method (xp)
               dspmat2=dspmat2_new;
               v2=v2_new; f=f_new;
            end
            % Plot previous Xmid to help picking
            if exist('dspmatprev','var')==1
                % Plot previous dispersion image if existing
                if xmidprev==0
                    flimplot=flim(Xmidselec(i-1));
                    figname=['Xmid = ',num2str(XmidT(Xmidselec(i-1))),' m'];
                else
                    flimplot=flim(Xmidselec(i+1));
                    figname=['Xmid = ',num2str(XmidT(Xmidselec(i+1))),' m'];
                end
                if mappicklog==0
                    fig2=plot_img(2,f,v2,-dspmatprev',mappick,axetop,axerev,0,12,...
                        freqtitle_long,'Phase velocity (m/s)',...
                        'Norm. ampli.',[fMIN fMAX],[VphMIN VphMAX],...
                        [],fticks,Vphticks,[],[],flimplot,[],[],[],[],0);
                else
                    dspmatinvprev=1./(1-dspmatprev);
                    dspmatinvprev(isinf(dspmatinvprev))=max(max(dspmatinvprev(isinf(dspmatinvprev)==0)));
                    fig2=plot_img_log(2,f,v2,-dspmatinvprev',flipud(mappick),axetop,axerev,0,12,...
                        freqtitle_long,'Phase velocity (m/s)',...
                        '1/(1-Norm. ampli.)',[fMIN fMAX],[VphMIN VphMAX],...
                        [1 length(mappick)],fticks,Vphticks,[],[],flimplot,[],[],[],[],0);
                end
                hold on
                if plotlamlim==1 && sampling==1
                    dashline(f,f*max(resampvec),3,3,3,3,'color',[0 0 0],'linewidth',5);
                end
                % Window positions during picking
                set(fig2,'name',figname,'numbertitle','off');
                % set(gcf,'units','normalized','outerposition',[0 0 1 1]); % Full screen left
                % set(gcf,'units','normalized','outerposition',[-0.02 0.5 0.5 0.5]); % 1/4 screen left bottom
                set(gcf,'units','normalized','outerposition',[0.65 0.5 0.35 0.5]); % 1/4 screen right top
                if Flogscale==1
                    set(gca,'xscale','log');
                end
                % Plot previous dispersion curves if existing
                if xmidprev==0
                    pvcstruct=dir(fullfile(dir_pick,...
                        [num2str(XmidT(Xmidselec(i-1)),xmidformat),'.*.pvc']));
                else
                    pvcstruct=dir(fullfile(dir_pick,...
                        [num2str(XmidT(Xmidselec(i+1)),xmidformat),'.*.pvc']));
                end
                for ip=1:length(pvcstruct)
                    pvcfile=pvcstruct(ip).name;
                    Vprev=load(fullfile(dir_pick,pvcfile));
                    hold on;
                    plot(Vprev(:,1),Vprev(:,2),'c.');
                end
            end
            modenext=modeinit;
            while isempty(modenext)==0
                % Plot current dispersion image
                filepick=fullfile(dir_pick,[num2str(XmidT(ix),xmidformat),...
                    '.M',num2str(modenext),'.pvc']);
                if mappicklog==0
                    [fig1,h1,~,h0]=plot_img(1,f,v2,-dspmat2',mappick,axetop,axerev,0,12,...
                        freqtitle_long,'Phase velocity (m/s)',...
                        'Norm. ampli.',[fMIN fMAX],[VphMIN VphMAX],...
                        [],fticks,Vphticks,[],[],flim(ix),[],[],[],[],0);
                else
                    dspmatinv2=1./(1-dspmat2);
                    dspmatinv2(isinf(dspmatinv2))=max(max(dspmatinv2(isinf(dspmatinv2)==0)));
                    [fig1,h1,~,h0]=plot_img_log(1,f,v2,-dspmatinv2',flipud(mappick),axetop,axerev,0,12,...
                        freqtitle_long,'Phase velocity (m/s)',...
                        '1/(1-Norm. ampli.)',[fMIN fMAX],[VphMIN VphMAX],...
                        [1 length(mappick)],fticks,Vphticks,[],[],flim(ix),[],[],[],[],0);
                end
                hold on
                if plotlamlim==1 && sampling==1
                    dashline(f,f*max(resampvec),3,3,3,3,'color',[0 0 0],'linewidth',5);
                end
                % Window positions during picking
                set(fig1,'name',['Xmid = ',num2str(XmidT(ix)),' m'],'numbertitle','off');
                % set(gcf,'units','normalized','outerposition',[0.3 0 1 1]); % Full screen right
                % set(gcf,'units','normalized','outerposition',[-0.02 -0.02 0.5 0.5]); % % 1/4 screen left top
                set(gcf,'units','normalized','outerposition',[-0.02 -0.02 0.6 1]); % % 1/2 screen left
                if Flogscale==1
                    set(gca,'xscale','log');
                end
                pvcstruct=dir(fullfile(dir_pick,[num2str(XmidT(ix),xmidformat),'.*.pvc']));
                % Plot current dispersion curves if existing
                for ip=1:length(pvcstruct)
                    pvcfile=pvcstruct(ip).name;
                    m=str2double(pvcfile(strfind(pvcfile,'.M')+2:strfind(pvcfile,'.pvc')-1)); % Mode number
                    if m~=modenext
                        Vprev=load(fullfile(dir_pick,pvcfile));
                        hold on;
                        plot(Vprev(:,1),Vprev(:,2),'m.');
                    end
                end
                % Pick dispersion curves
                [fi,Vi,deltac,modenext,closefig,xmidprev]=matpickamp(dspmat2,f,v2,filepick,pickstyle,...
                    modenext,err,smoothpick,mean([nWmin,nWmax]),dx,nWfac,maxerrrat,minerrvel,sigma);
                if closefig==0
                    close(fig1);
                end
                if xmidprev==1 && i==1
                    fprintf('\n  No previous Xmid - Stay on first Xmid\n');
                end
                if isempty(Vi)==0 && length(Vi)>1 && sum(isnan(Vi))~=length(Vi)
                    dlmwrite(filepick,[fi(isnan(Vi)==0);Vi(isnan(Vi)==0);...
                        deltac(isnan(Vi)==0)]','delimiter','\t','precision','%.6f');
                    apvcfile=[filepick(1:end-4),'.apvc'];
                    if exist(apvcfile,'file')==2
                        delete(apvcfile);
                    end
                elseif exist(filepick,'file')==2 && (isempty(Vi)==1 || (length(Vi)>1 && sum(isnan(Vi))~=length(Vi)))
                    delete(filepick);
                end
            end
            if exist('dspmatprev','var')==1
                close(fig2);
            end
            dspmatprev=dspmat2; % Store current dispersion matrix to show along next one
            
        elseif pick==2 % Automatic picking (requires at least one manually picked file - experimental)
            % Get previously picked dispersion curves
            filepick=fullfile(dir_pick,[num2str(XmidT(ix),xmidformat),...
                '.M',num2str(modeinit),'.pvc']);
            pvcstructauto=dir(fullfile(dir_pick,['*.M',num2str(modeinit),'.pvc']));
            if length(pvcstructauto)<1
                fprintf('\n  Autopick requires at least one manually picked file\n');
                continue
            end
            % Read existing dispersion curves
            if i>1 && exist('pvcfileauto','var')==1
                filepickprev=fullfile(dir_pick,[num2str(XmidT(Xmidselec(i-1)),...
                    xmidformat),'.M',num2str(modeinit),'.pvc']);
                if exist(filepick,'file')==2
                    pvcfileauto=filepick;
                elseif exist(filepickprev,'file')==2
                    pvcfileauto=filepickprev;
                end
            else
                if exist(filepick,'file')==0
                    [~,idx]=sort([pvcstructauto.datenum]);
                    pvcfileauto=fullfile(dir_pick,pvcstructauto(idx(end)).name);
                else
                    pvcfileauto=filepick;
                end
            end
            Vprevauto=load(pvcfileauto);
            fpick=Vprevauto(:,1);
            vpick=Vprevauto(:,2);
            wl=Vprevauto(:,3);
            % Perform autopick, median filter and moving average
            [vpickauto,fpickauto]=findpeak(dspmat2,f,v2,fpick,vpick,1.5*wl);
            vpickauto=median_filt(vpickauto,5,1,length(vpickauto));
            vpickauto=mov_aver(vpickauto',3,1,length(vpickauto));
            deltacauto=lorentzerr(vpickauto',vpickauto'./fpickauto,mean([nWmin,nWmax]),dx,...
                nWfac,maxerrrat,minerrvel);
            dlmwrite(filepick,[fpickauto;vpickauto';deltacauto']','delimiter','\t','precision','%.6f');
            fprintf(['\n  Automatic pick for mode ',num2str(modeinit),'\n']);
        end
        
        %% %% %%
        
        %%%%% Dinver target creation from .pvc dispersion curves %%%%%
        
        % Check existence of all dispersion curves for this Xmid
        pvcstruct=dir(fullfile(dir_pick,[num2str(XmidT(ix),xmidformat),'.*.*pvc']));
        npvc=length(pvcstruct);
        
        % Convert .pvc files to dinver .target
        if target==1
            nametarg=fullfile(dir_targ,[num2str(XmidT(ix),xmidformat),'.target']);
            if npvc>0 && maxmodeinv>=0
                for ip=1:npvc % Read all picked modes
                    pvcfile=pvcstruct(ip).name;
                    [~,pvcname,extension]=fileparts(pvcfile);
                    if strcmp(extension,'.apvc')==1 % Check if previously unused
                        movefile(fullfile(dir_pick,pvcfile),fullfile(dir_pick,[pvcname,'.pvc']));
                        pvcfile=[pvcname,'.pvc'];
                    end
                    m=str2double(pvcfile(strfind(pvcfile,'.M')+2:strfind(pvcfile,'.pvc')-1)); % Mode number
                    if m>maxmodeinv
                        % Rename if mode is unused (cf maxmodeinv parameter)
                        apvcfile=[pvcfile(1:end-4),'.apvc'];
                        movefile(fullfile(dir_pick,pvcfile),fullfile(dir_pick,apvcfile))
                        continue
                    end
                    Vprev=load(fullfile(dir_pick,pvcfile));
                    % Uncertainty range calculation
                    if err==1
                        Vprev(:,3)=lorentzerr(Vprev(:,2)',Vprev(:,2)'./Vprev(:,1)',...
                            mean([nWmin,nWmax]),dx,nWfac,maxerrrat,minerrvel);
                    elseif err==2
                        Vprev(:,3)=Vprev(:,2)*0.01*sigma;
                    else
                        Vprev(:,3)=0;
                    end
                    % Write new file with uncertainty
                    dlmwrite(fullfile(dir_pick,pvcfile),Vprev,'delimiter','\t','precision','%.6f');
                end
                % Get all .pvc files
                pvcstruct=dir(fullfile(dir_pick,[num2str(XmidT(ix),xmidformat),'.*.pvc']));
                if isempty(pvcstruct)==0
                    % Convert .pvc in .target
                    if freqlim==1
                        lmpick=pvc2targ(pvcstruct,dir_pick,nametarg,wave,...
                            sampling,resampvec,flim(ix),maxerrrat);
                    else
                        lmpick=pvc2targ(pvcstruct,dir_pick,nametarg,wave,...
                            sampling,resampvec,fminpick,maxerrrat);
                    end
                    lmaxpick(ix)=max(lmpick);
                else
                    lmaxpick(ix)=NaN;
                    if exist(nametarg,'file')==2
                        delete(nametarg); % Delete old target if no .pvc exist anymore for this Xmid
                    end
                end
            else
                if exist(nametarg,'file')==2
                    delete(nametarg);
                end
                lmaxpick(ix)=NaN;
            end
        else
            lmaxpick(ix)=NaN;
        end
        
        %% %% %%
        
        %%%%% Read all dispersion curves %%%%%
        
        nametarg=fullfile(dir_targ,[num2str(XmidT(ix),xmidformat),'.target']);
        if exist(nametarg,'file')==2
            % Read target file to get picked dispersion curves
            [freqresamp,vresamp,deltaresamp,modes]=targ2pvc(nametarg);
            npvc=length(modes);
            for ip=1:npvc
                % Resample in lambda or frequency
                if length(freqresamp{modes(ip)+1})>1
                    [freqresamp{modes(ip)+1},vresamp{modes(ip)+1},deltaresamp{modes(ip)+1}]=...
                        resampvel(freqresamp{modes(ip)+1},vresamp{modes(ip)+1},...
                        deltaresamp{modes(ip)+1},resampvec,sampling,1);
                end
            end
        else
            npvc=0;
        end
        if sum(nshot(ix,:))>=0 && exist(nametarg,'file')~=2 && (plotpckdisp==1 || plot1dobs==1 || plot2dobs==1)
            fprintf('\n  No target file for this Xmid\n');
        end
        
        %% %% %%
        
        %%%%% Plot all dispersion curves together %%%%%
        
        if (plot1dobs==1 || plot2dobs==1) && npvc>0
            if exist('f','var')~=1
                if isempty(fMIN)==0
                    f=fMIN:1:fMAX;
                else
                    f=0:1:1000;
                end
            end
            for ip=1:npvc % Loop on number of modes
                % Plot all dispersion curves in 1D
                if plot1dobs==1
                    fprintf(['\n  Plot dispersion curve for mode ',num2str(modes(ip)),'\n']);
                    
                    % All modes in one figure
                    if length(Xmidselec)==1 % Only one Xmid
                        if eb==1
                            plot_curv(4,freqresamp{modes(ip)+1},vresamp{modes(ip)+1},...
                                deltaresamp{modes(ip)+1},'.-',[1 0 0],[],axetop,axerev,...
                                0,fs,freqtitle_long,'Phase velocity (m/s)',...
                                'X (m)',[fMIN fMAX],[VphMIN VphMAX],[],fticks,Vphticks,[],...
                                [],[],[1 1 24 18],[],[],Flogscale);
                        else
                            plot_curv(4,freqresamp{modes(ip)+1},vresamp{modes(ip)+1},[],'.-',[1 0 0],[],axetop,axerev,...
                                0,fs,freqtitle_long,'Phase velocity (m/s)',...
                                'X (m)',[fMIN fMAX],[VphMIN VphMAX],[],fticks,Vphticks,[],...
                                [],[],[1 1 24 18],[],[],Flogscale);
                        end
                        hold on
                        if plotlamlim==1 && sampling==1
                            dashline(f,f*max(resampvec),3,3,3,3,'color','k','linewidth',5);
                        end
                    else % Multiple Xmids
                        if ishandle(4)==0 % First Xmid
                            if eb==1
                                plot_curv(4,freqresamp{modes(ip)+1},vresamp{modes(ip)+1},...
                                    deltaresamp{modes(ip)+1},'.-',ccurve(ix,:),[],axetop,axerev,...
                                    cbpos,fs,freqtitle_long,'Phase velocity (m/s)',...
                                    'X (m)',[fMIN fMAX],[VphMIN VphMAX],[min(XmidT) max(XmidT)],...
                                    fticks,Vphticks,[],[],[],[1 1 24 18],[],[],Flogscale);
                                colormap(ccurve);
                            else
                                plot_curv(4,freqresamp{modes(ip)+1},vresamp{modes(ip)+1},[],'.-',ccurve(ix,:),[],axetop,axerev,...
                                    cbpos,fs,freqtitle_long,'Phase velocity (m/s)',...
                                    'X (m)',[fMIN fMAX],[VphMIN VphMAX],[min(XmidT) max(XmidT)],...
                                    fticks,Vphticks,[],[],[],[1 1 24 18],[],[],Flogscale);
                                colormap(ccurve);
                            end
                            hold on
                            if plotlamlim==1 && sampling==1
                                dashline(f,f*max(resampvec),3,3,3,3,'color','k','linewidth',5);
                            end
                        else % Other Xmids
                            figure(4);
                            plot(freqresamp{modes(ip)+1},vresamp{modes(ip)+1},'.-','Color',ccurve(ix,:),...
                                'linewidth',2,'markersize',10);
                            if eb==1
                                if str2double(matrelease(1:4))>2014
                                    han=terrorbar(freqresamp{modes(ip)+1},vresamp{modes(ip)+1},deltaresamp{modes(ip)+1},1,'units');
                                    set(han,'LineWidth',1.5,'Color',ccurve(ix,:))
                                else
                                    han=errorbar(freqresamp{modes(ip)+1},vresamp{modes(ip)+1},deltaresamp{modes(ip)+1},...
                                        '.-','Color',ccurve(ix,:),'linewidth',2,'markersize',10);
                                    xlimits=xlim;
                                    tick_length=diff(xlimits)/100;
                                    errorbar_tick(han,tick_length,'units');
                                end
                            end
                        end
                    end
                    
                    % Single modes on separate figures
                    if length(Xmidselec)==1 % Only one Xmid
                        if eb==1
                            plot_curv(modes(ip)+5,freqresamp{modes(ip)+1},vresamp{modes(ip)+1},...
                                deltaresamp{modes(ip)+1},'.-',[1 0 0],[],axetop,axerev,...
                                0,fs,freqtitle_long,'Phase velocity (m/s)',...
                                'X (m)',[fMIN fMAX],[VphMIN VphMAX],[],fticks,Vphticks,[],...
                                [],[],[30 1 24 18],[],[],Flogscale);
                        else
                            plot_curv(modes(ip)+5,freqresamp{modes(ip)+1},vresamp{modes(ip)+1},[],'.-',[1 0 0],[],axetop,axerev,...
                                0,fs,freqtitle_long,'Phase velocity (m/s)',...
                                'X (m)',[fMIN fMAX],[VphMIN VphMAX],[],fticks,Vphticks,[],...
                                [],[],[30 1 24 18],[],[],Flogscale);
                        end
                        hold on
                        if plotlamlim==1 && sampling==1
                            dashline(f,f*max(resampvec),3,3,3,3,'color','k','linewidth',5);
                        end
                    else % Multiple Xmids
                        if ishandle(modes(ip)+5)==0 % First Xmid
                            if eb==1
                                plot_curv(modes(ip)+5,freqresamp{modes(ip)+1},vresamp{modes(ip)+1},...
                                    deltaresamp{modes(ip)+1},'.-',ccurve(ix,:),[],axetop,axerev,...
                                    cbpos,fs,freqtitle_long,'Phase velocity (m/s)',...
                                    'X (m)',[fMIN fMAX],[VphMIN VphMAX],[min(XmidT) max(XmidT)],...
                                    fticks,Vphticks,[],[],[],[30 1 24 18],[],[],Flogscale);
                                colormap(ccurve);
                            else
                                plot_curv(modes(ip)+5,freqresamp{modes(ip)+1},vresamp{modes(ip)+1},[],'.-',ccurve(ix,:),[],axetop,axerev,...
                                    cbpos,fs,freqtitle_long,'Phase velocity (m/s)',...
                                    'X (m)',[fMIN fMAX],[VphMIN VphMAX],[min(XmidT) max(XmidT)],...
                                    fticks,Vphticks,[],[],[],[30 1 24 18],[],[],Flogscale);
                                colormap(ccurve);
                            end
                            hold on
                            if plotlamlim==1 && sampling==1
                                dashline(f,f*max(resampvec),3,3,3,3,'color','k','linewidth',5);
                            end
                        else % Other Xmids
                            figure(modes(ip)+5);
                            plot(freqresamp{modes(ip)+1},vresamp{modes(ip)+1},'.-','Color',ccurve(ix,:),...
                                'linewidth',2,'markersize',10);
                            if eb==1
                                if str2double(matrelease(1:4))>2014
                                    han=terrorbar(freqresamp{modes(ip)+1},vresamp{modes(ip)+1},deltaresamp{modes(ip)+1},1,'units');
                                    set(han,'LineWidth',1.5,'Color',ccurve(ix,:))
                                else
                                    han=errorbar(freqresamp{modes(ip)+1},vresamp{modes(ip)+1},deltaresamp{modes(ip)+1},...
                                        '.-','Color',ccurve(ix,:),'linewidth',2,'markersize',10);
                                    xlimits=xlim;
                                    tick_length=diff(xlimits)/100;
                                    errorbar_tick(han,tick_length,'units');
                                end
                            end
                        end
                    end
                    drawnow;
                end
                % Save all dispersion curves in matrix
                if plot2dobs==1
                    vph2dobs{modes(ip)+1}(:,ix)=vresamp{modes(ip)+1}';
                end
            end
            
            % Save 1D image with dispersion curves at the end of main loop
            if plot1dobs==1 && ishandle(4)==1 && (ix==Xmidselec(find(Xmidselec(sum(nshot(Xmidselec,:)>0)),1,'last')) || (pick==1 && xmidprev==-1))
                fprintf('\n  Save picked dispersion curves\n');
                file1=fullfile(dir_img,['Dispcurve.allmode.',imgform]);
                save_fig(4,file1,imgform,imgres,1);
                close(4)
                figHandles = findall(0,'Type','figure');
                if str2double(matrelease(1:4))>2014
                    figHandles=get(figHandles,'Number');
                end
                for ifig=figHandles'
                    if str2double(matrelease(1:4))>2014 && length(figHandles)>1
                        ifigok=ifig{1};
                    else
                        ifigok=ifig;
                    end
                    file1=fullfile(dir_img,['Dispcurve.M',num2str(ifigok-5),'.',imgform]);
                    save_fig(ifigok,file1,imgform,imgres,1);
                    close(ifigok)
                end
            end
        end
        
        %% %% %%
        
        %%%%% Plot and save picked dispersion image %%%%%
        
        if  plotpckdisp==1 && npvc>0
            fprintf('\n  Plot and save picked dispersion images\n');
            if exist(dspfile_sum,'file')==2
                if Dlogscale==0
                    fig1=plot_img(showplot,f,v,dspmat',flipud(map0),axetop,axerev,cb_disp,fs,...
                        freqtitle_long,'Phase velocity (m/s)',...
                        'Norm. ampli.',[fMIN fMAX],[VphMIN VphMAX],...
                        [],fticks,Vphticks,[],[],flimsing,[],[0 0 24 18],[]);
                else
                    dspmatinv=1./(1-dspmat);
                    dspmatinv(isinf(dspmatinv))=max(max(dspmatinv(isinf(dspmatinv)==0)));
                    fig1=plot_img_log(showplot,f,v,dspmatinv',flipud(map0),axetop,axerev,cb_disp,fs,...
                        freqtitle_long,'Phase velocity (m/s)',...
                        '1/(1-Norm. ampli.)',[fMIN fMAX],[VphMIN VphMAX],...
                        [1 length(map0)],fticks,Vphticks,[],[],flimsing,[],[0 0 24 18],[]);
                end
            else
                fig1=plot_curv(showplot,NaN,NaN,[],'.',[0 0 0],[],axetop,axerev,...
                    0,fs,freqtitle_long,'Phase velocity (m/s)',[],...
                    [fMIN fMAX],[VphMIN VphMAX],[],fticks,Vphticks,[],...
                    [],[],[0 0 24 18],[]);
            end
            hold on
            hh=dashline(f,f,3,3,3,3,'color',[0 0 1],'linewidth',5);
            set(hh,'visible','off');
            if plotlamlim==1 && sampling==1
                dashline(f,f*max(resampvec),3,3,3,3,'color',[0 0 1],'linewidth',5);
            end
            if Flogscale==1
                set(gca,'xscale','log');
            end
            for ip=1:npvc
                hold on
                if mod(modes(ip),2)==0
                    col=pickcol1;
                else
                    col=pickcol2;
                end
                han=plot(freqresamp{modes(ip)+1},vresamp{modes(ip)+1},'.','Color',col,...
                    'linewidth',2,'markersize',9);
                if eb==1
                    if str2double(matrelease(1:4))>2014
                        han=terrorbar(freqresamp{modes(ip)+1},vresamp{modes(ip)+1},deltaresamp{modes(ip)+1},1,'units');
                        set(han,'LineWidth',1.5,'Color',col)
                    else
                        han=errorbar(freqresamp{modes(ip)+1},vresamp{modes(ip)+1},deltaresamp{modes(ip)+1},...
                            '.','Color',col,'linewidth',2,'markersize',9);
                        xlimits=xlim;
                        tick_length=diff(xlimits)/100;
                        errorbar_tick(han,tick_length,'units');
                    end
                end
            end
            hold off
            file1=fullfile(dir_img_pick,[num2str(XmidT(ix),xmidformat),'.disp.pick.',imgform]);
            save_fig(fig1,file1,imgform,imgres,1);
            close(fig1)
            
            if stack==2 % Alt. method (xp)
                if exist(dspfile_sum_new,'file')==2
                   if Dlogscale==0
                       fig1=plot_img(showplot,f_new,v_new,dspmat_new',flipud(map0),axetop,axerev,cb_disp,fs,...
                           freqtitle_long,'Phase velocity (m/s)',...
                           'Norm. ampli.',[fMIN fMAX],[VphMIN VphMAX],...
                           [],fticks,Vphticks,[],[],flimsing,[],[0 0 24 18],[]);
                   else
                       dspmatinv=1./(1-dspmat_new);
                       dspmatinv(isinf(dspmatinv))=max(max(dspmatinv(isinf(dspmatinv)==0)));
                       fig1=plot_img_log(showplot,f_new,v_new,dspmatinv',flipud(map0),axetop,axerev,cb_disp,fs,...
                           freqtitle_long,'Phase velocity (m/s)',...
                           '1/(1-Norm. ampli.)',[fMIN fMAX],[VphMIN VphMAX],...
                           [1 length(map0)],fticks,Vphticks,[],[],flimsing,[],[0 0 24 18],[]);
                   end
               else
                   fig1=plot_curv(showplot,NaN,NaN,[],'.',[0 0 0],[],axetop,axerev,...
                       0,fs,freqtitle_long,'Phase velocity (m/s)',[],...
                       [fMIN fMAX],[VphMIN VphMAX],[],fticks,Vphticks,[],...
                       [],[],[0 0 24 18],[]);
               end
               hold on
               hh=dashline(f_new,f_new,3,3,3,3,'color',[0 0 1],'linewidth',5);
               set(hh,'visible','off');
               if plotlamlim==1 && sampling==1
                   dashline(f_new,f_new*max(resampvec),3,3,3,3,'color',[0 0 1],'linewidth',5);
               end
               if Flogscale==1
                   set(gca,'xscale','log');
               end
               for ip=1:npvc
                   hold on
                   if mod(modes(ip),2)==0
                       col=pickcol1;
                   else
                       col=pickcol2;
                   end
                   han=plot(freqresamp{modes(ip)+1},vresamp{modes(ip)+1},'.','Color',col,...
                       'linewidth',2,'markersize',9);
                   if eb==1
                       if str2double(matrelease(1:4))>2014
                           han=terrorbar(freqresamp{modes(ip)+1},vresamp{modes(ip)+1},deltaresamp{modes(ip)+1},1,'units');
                           set(han,'LineWidth',1.5,'Color',col)
                       else
                           han=errorbar(freqresamp{modes(ip)+1},vresamp{modes(ip)+1},deltaresamp{modes(ip)+1},...
                               '.','Color',col,'linewidth',2,'markersize',9);
                           xlimits=xlim;
                           tick_length=diff(xlimits)/100;
                           errorbar_tick(han,tick_length,'units');
                       end
                   end
               end
               hold off
               file1=fullfile(dir_img_pick,[num2str(XmidT(ix),xmidformat),'.disp_new.pick.',imgform]);
               save_fig(fig1,file1,imgform,imgres,1);
               close(fig1)
            end
        end
    else
        fprintf('\n  No dispersion data for this Xmid\n');
    end
    
    if nshot(ix)>=0
        fprintf('\n  **********************************************************');
        fprintf('\n  **********************************************************\n');
    end
    
    %% %% %%
    
    %%%%% Save settings and remove temp files %%%%%
    
    if calc~=2
        if clearmem==1 && exist(dir_dat_xmid,'dir')==7
            rmdir(dir_dat_xmid,'s');
        end
    end
    if calc~=0
        xmidparam.nshot(ix,:)=nshot(ix,:);
        xmidparam.Gmin(ix,:)=Gmin(ix,:);
        xmidparam.Gmax(ix,:)=Gmax(ix,:);
        save(matfile,'-append','xmidparam');
    end
    if freqlim==1
        xmidparam.flim(ix,:)=flim(ix,:);
        save(matfile,'-append','xmidparam');
    end
    if exist('plotopt','var')==0
        plotopt=struct('f',f,'v',v,'fspec',fspec,'tseis',tseis,...
            'sizeax',sizeax,'xmin',xmin,'xmax',xmax);
        save(matfile,'-append','plotopt');
    end
    if exist('plotopt','var')==1 && isempty(plotopt.sizeax)==1 && isempty(sizeax)==0
        plotopt.sizeax=sizeax;
        save(matfile,'-append','plotopt');
    end
    if exist('plotopt','var')==1 && calc==1
        if exist('f','var')==1 && isempty(f)==0
            plotopt.f=f; plotopt.v=v;
        end
        if exist('fspec','var')==1 && isempty(fspec)==0
            plotopt.fspec=fspec;
        end
        if  exist('tseis','var')==1 && isempty(tseis)==0
            plotopt.tseis=tseis;
        end
        save(matfile,'-append','plotopt');
    end
    if target==1
        targopt.lmaxpick(ix)=lmaxpick(ix);
        % Save param in .mat file
        save(matfile,'-append','targopt');
    end
    
    % Initialize next Xmid position
    if pick==1 && exist('xmidprev','var')==1 && xmidprev==1 && i~=1
        ix=ix-2;
        i=i-2;
    elseif pick==1 && exist('xmidprev','var')==1 && xmidprev==1 && i==1
        ix=ix-1;
        i=i-1;
        clear('dspmatprev');
    elseif pick==1 && exist('xmidprev','var')==1 && xmidprev==-1
        break
    end
end

%% %% %%

%%%%%% Plot and save picked dispersion 2D pseudo-section %%%%%%

if plot2dobs==1 && Xlength>1
    flagprint=0;
    for ip=1:maxmodeinv+1
        if sum(sum(isnan(vph2dobs{ip})))==numel(vph2dobs{ip})
            continue
        end
        if flagprint==0
            fprintf('\n  Saving observed phase velocity sections\n');
            flagprint=1;
        end
        fprintf(['\n      Mode ',num2str(ip-1),'\n']);
        if sampling==0
            f1=plot_img(showplot,XmidT,resampvec,vph2dobs{ip},map1,0,0,cbpos,fs,'X (m)',...
                freqtitle_short,'Vph_{obs} (m/s)',[xMIN xMAX],[0 max(resampvec)],...
                [vphMIN vphMAX],xticks,fticks,vphticks,[],[],vphISO,[0 0 24 12],[],1,0);
        else
            f1=plot_img(showplot,XmidT,resampvec,vph2dobs{ip},map1,1,1,cbpos,fs,'X (m)',...
                lamtitle,'Vph_{obs} (m/s)',[xMIN xMAX],[lamMIN lamMAX],...
                [vphMIN vphMAX],xticks,lticks,vphticks,[],[],vphISO,[0 0 24 12],[],1,0);
        end
        szfig1=get(f1,'Position');
        file1=fullfile(dir_img,['Vphobs.M',num2str(ip-1),'.',imgform]);
        save_fig(f1,file1,imgform,imgres,1);
        if showplot==0
            close(f1);
        else
            showplot=showplot+1;
        end
    end
    if flagprint==0
        fprintf('\n  No dispersion to plot - Re-run with target=1\n');
    end
end

%% %% %%

%%%%%% Remove temp files and check time %%%%%%

% Remove filtered SU file
if filt==1 && calc==1
    delete(sufilefilt);
end

% Check elapsed time
Tend=toc(Tstart);
[time_string]=secs2hms(Tend);
fprintf(['\n  Total elapsed time: ',time_string,'\n']);