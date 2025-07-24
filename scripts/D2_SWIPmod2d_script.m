%%% SURFACE-WAVE dispersion INVERSION & PROFILING (SWIP)
%%% MODULE D2 : SWIPmod2d.m
%%% S. Pasquet - V22.05.04
%%% SWIPmod2d.m plots observed, calculated and residual pseudo-sections
%%% It also plots Vp, Vs, Vp/Vs, Poisson's ratio and auxiliary data 2D sections

%%% This module calls the following Geopsy native function:
%%% gpdc
%%% Geopsy function are called through the following MATLAB function:
%%% matgpdc
%%% The following Linux codes are also called if correctly installed:
%%% ImageMagick (convert, montage) - pdfjam - pdfcrop

%%%-------------------------%%%
%%% START OF INITIALIZATION %%%

dpMIN=0; % Min depth (m)

run('SWIP_defaultsettings') % Read default settings

% Check user input
if input_vel==0 || input_vel>2
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n      Select at least one model');
    fprintf('\n   Set either "input_vel" to 1 or 2');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
    return
end

if plot2dcal==0 && plot2dmod==0 && plothisto==0 && savexzv==0
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n   Select at least one plot/save option');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
    return
end

% Initialization (same for D1 and D2)
if input_vel==1
    % Get SWIP subproject and inversion folders
    [dir_all,dir_inv_img]=dir_create(2);
    if dir_all.dir_main==0
        return
    end
    if dir_inv_img.dir_rep_inv==0
        return
    end
    % Read previous inversion settings
    dir_rep_inv=dir_inv_img.dir_rep_inv;
    dir_img_inv=dir_inv_img.dir_img_inv;
    dir_xzv_inv=dir_inv_img.dir_xzv_inv; % Only for D2
    matstruct=dir(fullfile(dir_rep_inv,'*.invparam.mat'));
    matfileinv=fullfile(dir_rep_inv,matstruct.name);
    try
        load(matfileinv);
    catch
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        fprintf('\n   Missing .mat file in inversion folder');
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
        return
    end
    if isfield(inv_set,'maxmodeinv')==1
        maxmodeinv=inv_set.maxmodeinv;
    else
        nmodeinv=inv_set.nmodeinv;
        maxmodeinv=nmodeinv-1;
    end
    paramtype=inv_set.paramtype;
    nmod = inv_set.nmod;
else
    % Get SWIP suproject folder
    dir_all=dir_create(0);
    if dir_all.dir_main==0
        return
    end
    maxmodeinv=[];
end
dir_dat=dir_all.dir_dat;
dir_img=dir_all.dir_img;
dir_pick=dir_all.dir_pick;
dir_param=dir_all.dir_param;
dir_targ=dir_all.dir_targ;
dir_inv=dir_all.dir_inv;
dir_xzv=dir_all.dir_xzv;
matstruct=dir(fullfile(dir_dat,'*.param.mat'));
matfile=fullfile(dir_dat,matstruct.name);
try
    load(matfile);
catch
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n   Missing .mat file in file.dat folder');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
    return
end
dx=acquiparam.dx;
topo=acquiparam.topo;
nWmin=stackdisp.nWmin;
nWmax=stackdisp.nWmax;
xsca=pomega.xsca;
XmidT=xmidparam.XmidT; % Get Xmids
Xlength=length(XmidT); % Number of Xmids
xmidformat=stackdisp.xmidformat;

% Select Xmids
if exist('Xmidselec','var')~=1 || isempty(Xmidselec)==1
    Xmidselec=1:Xlength;
end
if max(Xmidselec)>Xlength
    Xmidselec=Xmidselec(Xmidselec<=Xlength);
end
nshot=xmidparam.nshot;
fmin=pomega.fmin;
fmax=pomega.fmax;
f=plotopt.f;
nf=length(f);
if exist('targopt_inv','var')==0 || input_vel == 2
    wave=targopt.wave(1);
    resampvec=targopt.resampvec;
    sampling=targopt.sampling;
    lmaxpick=targopt.lmaxpick;
else
    wave=targopt_inv.wave(1);
    resampvec=targopt_inv.resampvec;
    sampling=targopt_inv.sampling;
    lmaxpick=targopt_inv.lmaxpick;
end

% Initialize depth vector and topography
zround=xmidparam.zround; % Get topography
if isempty(dpMAX)==1 % Get maximum depth from parameterization if not setup in launcher
    if input_vel==1
        fprintf('\n  Looking for maximum DOI...\n');
        dpMAX=zeros(size(Xmidselec));
        for ix=Xmidselec
            dir_rep_ind=fullfile(dir_rep_inv,[num2str(XmidT(ix),xmidformat),'_reports']);
            if exist(dir_rep_ind,'dir')==7
                paramstruct=dir(fullfile(dir_rep_ind,'*.param'));
                for ip=1:length(paramstruct)
                    paramfile=fullfile(dir_rep_ind,paramstruct(ip).name);
                    [~,dpMAX(ix)]=param2mod(paramfile);
                end
            end
        end
        dpMAX=max(dpMAX);
        if dpMAX == 0
            if paramtype>0
                paramstruct=dir(fullfile(dir_targ,['*.type',num2str(paramtype),'.param']));
                dpMAX=zeros(size(paramstruct));
                for ip=1:length(paramstruct)
                    paramfile=fullfile(dir_targ,paramstruct(ip).name);
                    [~,dpMAX(ip)]=param2mod(paramfile);
                end
                dpMAX=max(dpMAX);
            else
                if isempty(strfind(dir_rep_inv,dir_inv))==1
                    dir_rep_inv=fullfile(dir_all.dir_start,dir_rep_inv);
                end
                paramname=fullfile(dir_param,[dir_rep_inv(length(dir_inv)+1:end),'.param']);
                if exist(paramname,'file')==2
                    [~,dpMAX]=param2mod(paramname);
                else
                    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!');
                    fprintf('\n   Missing .param file');
                    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!\n\n');
                    return
                end
            end
        end
    else
        dpMAX=max(resampvec);
    end
    dpMIN=0;
end
% Create depth vector
maxdepth=ceil(dpMAX/dz)*dz;
depth = bsxfun(@plus,zround,(0:-dz:-maxdepth)');
ZZ=0:dz:maxdepth;
nZ=length(ZZ);

% File and folder names initialization
if nbest==0
    extens=['.bweb',num2str(outpoints)]; % Best within error bars
else
    extens=['.best',num2str(nbest)]; % Arbitrary nb
end
if input_vel==1
    dir_img_inv_mod=fullfile(dir_img_inv,['models',extens]);
    dir_img_inv_2d=fullfile(dir_img_inv_mod,'2dmodels');
    dir_xzv_inv_mod=fullfile(dir_xzv_inv,['models',extens]);
else
    dir_img_ind=fullfile(dir_img,'Usermodels');
    dir_img_inv_mod=dir_img_ind;
    dir_img_inv_2d=fullfile(dir_img_ind,'2dmodels');
    dir_xzv_inv_mod=fullfile(dir_xzv,'Usermodels');
end
if exist(dir_img_inv_2d,'dir')~=7
    mkdir(dir_img_inv_2d);
end
if exist(dir_xzv_inv_mod,'dir')~=7
    mkdir(dir_xzv_inv_mod);
end
% Final average model types
if modeltype==1
    modeltype='best'; avertype='Vms';
elseif modeltype==2
    modeltype='layered'; avertype='Vms';
elseif modeltype==3
    modeltype='smooth'; avertype='Vms';
elseif modeltype==4
    modeltype='layered'; avertype='Vws';
elseif modeltype==5
    modeltype='smooth'; avertype='Vws';
elseif modeltype==6
    modeltype='ridge'; avertype='Vms';
else
    modeltype='smooth'; avertype='Vws';
    fprintf('\n  Weighted smooth model selected by default\n');
end
XmidT_vp=XmidT; depth_vp=depth;

fprintf('\n  ---------------------\n');

% Select refraction velocity models from file
if ((input_vel==1 && usevptomo==1) || input_vel==2) && (plot2dcal==1 || plothisto==1 || plot2dmod==1 || savexzv>0)
    
    % Select VP
    modstruct = [dir('*.model');dir('*.dat');dir('*.xzv');dir('*.txt')];
    if length(modstruct) ~= 1
        fprintf('\n  Select Vp model file');
        [filevel,pathvel] = uigetfile({'*.model;*.dat;*.xzv;*.txt'},'Select Vp model');
    else
        fprintf('\n  Vp model file');
        filevel = modstruct.name;
        pathvel = pwd;
    end
    
    [~,selected_folder,ext] = fileparts(filevel);
    if length(modstruct) ~= 1
        fprintf('\n  Manual selection: %s\n',strcat(selected_folder,ext));
    else
        fprintf('\n  Auto selection: %s\n',strcat(selected_folder,ext));
    end
    
    Vpfile=fullfile(pathvel,filevel); % File with velocity (3 columns X,Z,Vp)
    if pathvel==0
        VpI=[]; VpItomo=[]; VsItomo=[]; usevptomo=0; zinc=0:dz:maxdepth;
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        fprintf('\n   No Vp model file selected - Use Vp from SWIP');
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
        if input_vel==2
            return
        end
    else
        try
            X_plot = repmat(XmidT,size(depth,1),1);
            X_plot_vp = repmat(XmidT_vp,size(depth_vp,1),1);
            
            if input_vel==1 && usevptomo==1
                VpI = readtomo(Vpfile,0,X_plot,depth,xsca,vpaver,[nWmin,nWmax],dx); % Read Vp tomo file
                max_vp_depth = zeros(size(XmidT));
                if vpmask==1
                    X_plot_vp=X_plot; depth_vp=depth;
                else
                    [VpItomo,X_plot_vp,depth_vp]=readtomo(Vpfile,0,[],[],xsca); % Read Vp tomo file
                end
            elseif input_vel==2
                [VpItomo,Xi,Zi]=readtomo(Vpfile,0,X_plot,depth,xsca,vpaver,[nWmin,nWmax],dx); % Read Vp tomo file
                zround=interp1(XmidT,zround,unique(Xi),'linear','extrap');
                XmidT=unique(Xi)';
                Xlength=length(XmidT);
                Xmidselec=1:Xlength;
                X_plot_vp=X_plot; depth_vp=depth;
            end
        catch
            fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
            fprintf('\n   Invalid Vp model file - Use Vp from SWIP');
            fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
            VpI=[]; VpItomo=[]; VsItomo=[]; usevptomo=0; zinc=0:dz:maxdepth;
            if input_vel==2
                return
            end
        end
    end
    
    if input_vel==2 && isempty(VpItomo)==0
        % Select VS
        fprintf('\n  Select Vs model file (cancel to skip Vs)\n');
        [filevel,pathvel]=uigetfile({'*.model;*.dat;*.xzv;*.txt'},'Select Vs model (cancel if no Vs model available)');
        if pathvel==0
            VsItomo=sqrt(VpItomo.^2/((1/(1-2*pois_test))+1));
            fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
            fprintf('\n   No Vs model file selected - Use Poisson''s ratio of %1.2f',pois_test);
            fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
            
        else
            Vsfile=fullfile(pathvel,filevel); % File with velocity (3 columns X,Z,Vs)
            try
                VsItomo=readtomo(Vsfile,0,X_plot,depth,xsca,vpaver,[nWmin,nWmax],dx); % Read Vp tomo file
                ind_nan = find(isnan(VsItomo) | isnan(VpItomo));
                VpItomo(ind_nan)=NaN; %VpItomo = reshape(VpItomo,length(depth),length(XmidT));
                VsItomo(ind_nan)=NaN; %VsItomo = reshape(VsItomo,length(depth),length(XmidT));
            catch
                fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
                fprintf('\n   Invalid Vs model file - Ignore Vstomo');
                fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
                VsItomo=[]; VpItomo=[];
            end
        end
    end
else
    zinc=0:dz:maxdepth;
    X_plot_vp = repmat(XmidT_vp,size(depth_vp,1),1);
end

% Specific settings for D2
indf=ones(Xlength,1);
DOI=indf*NaN;

% Select auxiliary data models from file
if input_aux==1
    fprintf('\n  Select auxiliary data file\n');
    [fileaux,pathaux]=uigetfile({'*.model;*.dat;*.xzv;*.txt'},'Select auxiliary data file');
    Auxfile=fullfile(pathaux,fileaux); % File with auxiliary data (3 columns X,Z,Aux)
    if pathaux==0
        AuxI=[];
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        fprintf('\n   No auxiliary data file selected - Ignore auxiliary data');
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
        input_aux=0;
    else
        if input_vel==1 && auxmask==1
            try
                AuxI=readtomo(Auxfile,0,X_plot,depth,xsca); % Read auxiliary file
                XmidT_aux=X_plot; depth_aux=depth;
                auxmat=AuxI;
            catch
                fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
                fprintf('\n   Invalid auxiliary data file - Ignore auxiliary data');
                fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
                AuxI=[];
            end
        else
            try
                [AuxI,XmidT_aux,depth_aux]=readtomo(Auxfile,0,[],[],xsca); % Read auxiliary file
                auxmat=AuxI;
            catch
                fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
                fprintf('\n   Invalid auxiliary datafile - Ignore auxiliary data');
                fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
                return
            end
        end
    end
end

% Initialization of the maximum number of modes
if isempty(maxmodeinv)==1
    pvcstruct=dir(fullfile(dir_pick,'*.pvc'));
    npvc=length(pvcstruct);
    M=[];
    for ip=1:npvc
        pvcfile=pvcstruct(ip).name;
        m=str2double(pvcfile(end-4));
        if ismember(m,M)==0
            M=[M,m];
        end
    end
    if isempty(M)==0
        maxmodeinv=ones(1,Xlength)*max(M);
    else
        maxmodeinv=zeros(1,Xlength);
    end
    npvc=length(M);
    modes=M;
else
    npvc=nmodemax;
    modes=(0:npvc-1);
end

% Initialization of phase velocity matrices
vph2dobs=cell(max(maxmodeinv)+1,1);
for ip=1:max(maxmodeinv)+1
    vph2dobs{ip}=zeros(length(resampvec),Xlength)*NaN;
end
delta2dobs=vph2dobs;
vph2dobsALL=[]; delta2dobsALL=[];
if plot2dcal==1 || plothisto==1
    vph2dcal=vph2dobs;
    vph2dres=vph2dobs;
    misfitall=zeros(1,Xlength)*NaN;
    vph2dcalALL=[];
end

% Initialization of velocity matrices
vpmat=zeros(size(depth)).*NaN;
vpmatok=vpmat; vsmat=vpmat;
vpvsmat=vpmat; vptvsmat=vpmat;
poismat=vpmat; rhomat=vpmat;
vsstdmat=vpmat; maskmat=vpmat;
maskmatvp=vpmat;

% Check if image concatenation functions are installed
[testimgmgck,~]=unix_cmd('which montage');
[testpdfjam,~]=unix_cmd('which pdfjam');
testplot=((testpdfjam==0 && strcmp(imgform,'pdf')==1) || (testimgmgck==0 && strcmp(imgform,'pdf')==0 && strcmp(imgform,'fig')==0));
if concat == 0
    testplot = 0;
end

% Check length of selected Xmid
if Xlength==1
    fprintf('\n  Only one Xmid - cannot plot 2D sections\n');
    return
end

fprintf('\n  **********************************************************');
fprintf('\n  **********************************************************\n');

%%% END OF INITIALIZATION %%%
%%%-----------------------%%%

%% CALCULATIONS FOR ALL XMIDS

%%%%%% Loop over all Xmids %%%%%%

Tstart=tic; % Start clock
modexist=zeros(1,Xlength);
fprintf('\n  Reading models...\n');
for ix=Xmidselec
    
    % Initialization
    vph2dobsall=[]; delta2dobsall=[]; vph2dcalall=[];
    if input_vel==1
        dir_rep_ind=[dir_rep_inv,'/',num2str(XmidT(ix),xmidformat),'_reports'];
        if exist(dir_rep_ind,'dir')~=7 && sum(nshot(ix,:))>=0
            fprintf(['\n  No data/inversion for Xmid',num2str(ix),' = ',...
                num2str(XmidT(ix),xmidformat),' m\n']);
            continue
        end
    end
    
    %% %% %%
    
    %%%%% Read picked dispersion curves %%%%%
    
    if plot2dcal==1 || plothisto==1
        if input_vel==1
            nametarg=fullfile(dir_rep_ind,[num2str(XmidT(ix),xmidformat),'.target']);
            if exist(nametarg,'file')==2
                % Read target file to get picked dispersion curves
                [freqresamp,vresamp,deltaresamp,modes]=targ2pvc(nametarg);
                npvc=length(modes);
                lmaxpicktmp=zeros(npvc,1);
                for ip=1:npvc
                    % Resample in lambda or frequency
                    [freqresamp{modes(ip)+1},vresamp{modes(ip)+1},deltaresamp{modes(ip)+1}]=...
                        resampvel(freqresamp{modes(ip)+1},vresamp{modes(ip)+1},...
                        deltaresamp{modes(ip)+1},resampvec,sampling,1);
                    lmaxpicktmp(ip)=max(vresamp{modes(ip)+1}./freqresamp{modes(ip)+1});
                end
                lmaxpick(ix)=max(lmaxpicktmp);
            else
                if sum(nshot(ix,:))>=0
                    fprintf(['\n  No dispersion picked for Xmid',num2str(ix),' = ',...
                        num2str(XmidT(ix),xmidformat),' m\n']);
                end
            end
        end
    end
    
    %% %% %%
    
    %%%%% Read velocity models %%%%%
    
    %%%% From SWIP inversion results %%%%
    
    if input_vel==1
        % Velocity file name
        filevel=fullfile(dir_rep_ind,[num2str(XmidT(ix),xmidformat),extens,'.',...
            avertype,'.',modeltype]);
        % Standard deviation file name
        filestd=fullfile(dir_rep_ind,[num2str(XmidT(ix),xmidformat),extens,'.',...
            'VmsStd.',modeltype]);
        % Get min and max std from ridge search
        filemin=fullfile(dir_rep_ind,[num2str(XmidT(ix),xmidformat),extens,'.',...
            'VmsMin.ridge']);
        filemax=fullfile(dir_rep_ind,[num2str(XmidT(ix),xmidformat),extens,'.',...
            'VmsMax.ridge']);
        if strcmp(modeltype,'best')==1
            filestd=fullfile(dir_rep_ind,[num2str(XmidT(ix),xmidformat),extens,'.',...
                'VmsStd.layered']);
            filemin=fullfile(dir_rep_ind,[num2str(XmidT(ix),xmidformat),extens,'.',...
                'VmsMin.layered']);
            filemax=fullfile(dir_rep_ind,[num2str(XmidT(ix),xmidformat),extens,'.',...
                'VmsMax.layered']);
        end
        if strcmp(modeltype,'ridge')==1
            filestd=fullfile(dir_rep_ind,[num2str(XmidT(ix),xmidformat),extens,'.',...
                'VmsStd.smooth']);
        end
        
        %%%% Create velocity file in gpdc format %%%%
        
        % Read velocity file
        if exist(filevel,'file')==0
            if exist(dir_rep_ind,'dir')==7 && sum(nshot(ix,:))>=0
                fprintf(['\n  No SWIP model for Xmid',num2str(ix),' = ',...
                    num2str(XmidT(ix),xmidformat),' m\n']);
            end
            modvel=[];
        else
            modvel=dlmread(filevel,'',1,0);
            moddepth=[0;cumsum(modvel(:,1))];
            if exist(filemin,'file')==2 && exist(filemax,'file')==2
                modmin=dlmread(filemin,'',1,0);
                modmax=dlmread(filemax,'',1,0);
                modstd=modmax;
                %%% A checker !
%                 modstd(:,3) = max([abs(modmax(:,3)-modvel(:,3)),abs(modmin(:,3)-modvel(:,3))],[],2);
                modstd(:,3) = (modmax(:,3)-modmin(:,3))/4; % Check division by two for one standard deviation
                
                % test
%                 modstd=dlmread(filestd,'',1,0);

            elseif exist(filestd,'file')==2
                modstd=dlmread(filestd,'',1,0);
            else
                modstd = modvel;
                modstd(:,2:4) = modstd(:,2:4).*0.2;
            end
            depthstd=[0;modstd(:,1)];
            
            if maxdepth>moddepth(end)
                moddepth(end)=maxdepth;
                zinc=0:dz:maxdepth;
            else
                modvel=modvel(moddepth<maxdepth,:);
                modstd=modstd(moddepth<maxdepth,:);
                moddepth=[0;cumsum(modvel(:,1))];
                depthstd=depthstd(1:length(moddepth));
            end
            thick=modvel(:,1);
            vpsw=modvel(:,2);
            vssw=modvel(:,3);
            rhosw=modvel(:,4);
            
            if strcmp(modeltype,'ridge') == 1
                vssw = median_filt(vssw',9,1,length(vssw));
                vssw = mov_aver(vssw',5,1,length(vssw));
            end
            
            depthstd(end)=0;
            vpstd=modstd(:,2);
            vsstd=modstd(:,3);
            rhostd=modstd(:,4);
            if ~any(vsstd)
                vsstd = vssw*0.15;
            end
            
            vsstd_perc = 100*(vsstd./vssw);
            vsstd_perc_log_up = abs(100*((log10(vssw+vsstd) - log10(vssw))./log10(vssw)));
            vsstd_perc_log_low = abs(100*((log10(vssw-vsstd) - log10(vssw))./log10(vssw)));
            vsstd_perc_log = (vsstd_perc_log_up + vsstd_perc_log_low)/2;

            vsstd_test = 2;
            if vsstd_test == 1
                vsstd_ok = vsstd;
            elseif vsstd_test == 2
                vsstd_ok = vsstd_perc;
            end
            
            %%% Replace VP from SWIP with VP from tomo file if required %%%
            
            if usevptomo==1
                filevel=[filevel,'_vptomo'];
                filedisp=[filevel,'.disp'];
                vptomo=VpI(VpI(:,ix)>0,ix);
                if isempty(vptomo)~=1
                    vptomo=[vptomo(1);vptomo;vptomo(end)];
                    ztmp=dz.*ones(1,size(VpI(VpI(:,ix)>0,ix),1));
                    zinc=[0,cumsum(ztmp)];
                    if max(zinc)>maxdepth && abs(max(zinc)-maxdepth)>1e-10
                        max_vp_depth(ix) = maxdepth;
                        zinc=0:dz:maxdepth;
                        vptomo=vptomo(1:length(zinc)+1);
                    elseif max(zinc)<maxdepth && abs(max(zinc)-maxdepth)>1e-10
                        max_vp_depth(ix) = max(zinc);
                        zinc=0:dz:maxdepth;
                        vptomo2=vptomo(end)*ones(length(zinc)+1,1);
                        vptomo2(1:length(vptomo))=vptomo;
                        vptomo=vptomo2;
                    else
                        zinc=0:dz:maxdepth;
                        max_vp_depth(ix) = max(zinc);
                    end
                    [vpsw,~,~,vssw]=velresamp(zinc,vptomo,moddepth,vssw,0.1,0,0);
                else
                    modvel=[];
                end
            else
                filedisp=[filevel,'.disp'];
            end
        end
        
        %%%% From refraction tomography models (Vp and Vs) %%%%
        
    elseif input_vel==2
        filevel=fullfile(dir_dat,[num2str(XmidT(ix),xmidformat),'.tomo']);
        filedisp=fullfile(dir_dat,[num2str(XmidT(ix),xmidformat),'.tomo.disp']);
        nlay=size(VpItomo(VpItomo(:,ix)>0,ix),1);
        ztomo=dz.*ones(nlay,1);
        zinc=[0;cumsum(ztomo)];
        vptomo=VpItomo(VpItomo(:,ix)>0,ix);
        if isempty(VsItomo)==0
            vstomo=VsItomo(VsItomo(:,ix)>0,ix);
        else
            vstomo=[];
        end
        if max(zinc)<maxdepth
            zinc=0:dz:maxdepth;
            vptomo2=NaN*ones(length(zinc),1);
            vptomo2(1:length(vptomo))=vptomo;
            vptomo=vptomo2;
            if isempty(VsItomo)==0
                if isempty(vstomo)==0
                    vstomo2=NaN*ones(length(zinc),1);
                    vstomo2(1:length(vstomo))=vstomo;
                    vstomo=vstomo2;
                else
                    vstomo=vptomo*NaN;
                end
            end
        end
        if isempty(vptomo)~=1
            if isempty(vstomo)~=1
                for ll=1:nlay
                    while poisson(vptomo(ll),vstomo(ll))<=0.1
                        vptomo(ll)=vptomo(ll)+1;
                        vstomo(ll)=vstomo(ll)-1;
                    end
                end
            end
        else
            if sum(nshot(ix,:))>=0
                fprintf(['\n  No Vp or Vs from tomo for Xmid',num2str(ix),' = ',...
                    num2str(XmidT(ix),xmidformat),' m\n']);
            end
            vptomo=[]; vstomo=[];
        end
    else
        filevel=[];
    end
    
    %% %% %%
    
    %%%%% Compute theoretical dispersion %%%%%
    
    D=[];
    if (input_vel==1 && isempty(modvel)==0) || (input_vel==2 && isempty(vptomo)~=1)
        modexist(ix)=1;
        
        if plot2dcal==1 || plothisto==1
            % Save in gpdc format to run forward model
            if input_vel==1 && usevptomo==1 && isempty(vptomo)~=1
                dinsave(filevel,thick,vpsw,vssw,rhosw);
            elseif input_vel==2 && isempty(vptomo)~=1 && isempty(vstomo)~=1
                if isempty(rhoMIN)==1
                    flagrho=1;
                    rhoMIN=Rhomin; rhoMAX=Rhomax;
                else
                    flagrho=0;
                end
                dinsave(filevel,ztomo,vptomo,vstomo,mean([rhoMIN,rhoMAX]));
                if flagrho==1
                    rhoMIN=[]; rhoMAX=[];
                end
            end
            
            % Run forward modeling
            nftest=nf;
            while nftest>10
                matgpdc(filevel,maxmodeinv(ix)+1,wave,nftest-1,fmin+(fmax-fmin)/nftest,fmax,sampling,filedisp);
                D=readdisp(filedisp,maxmodeinv(ix)+1);
                if isempty(D)==1 % Check if it worked, otherwise try with less frequency samples
                    if nftest==nf
                        if input_vel==1
                            fprintf(['\n  Re-run forward calculation with less frequency samples for Xmid',...
                                num2str(ix),' = ',num2str(XmidT(ix),xmidformat),' m\n']);
                        else
                            fprintf(['\n  Re-run forward calculation with less frequency samples for X = ',...
                                num2str(XmidT(ix),xmidformat),' m\n']);
                        end
                    end
                    nftest=nftest-10;
                else
                    break
                end
            end
            
            % Delete temp files
            if usevptomo==1 && input_vel==1
                delete(filevel);
            elseif input_vel==2
                delete(filevel);
            end
            delete(filedisp);
        end
        
        %%%%% Get Depth Of Investigation (DOI) %%%%%
        
        if plot2dmod==1 || savexzv>0
            if (plotDOI==0 || maskDOI==0) && input_vel==1
                indf(ix)=round(dpMAX/dz);
            end
            % Empirical DOI (Lmax*doifact)
            if (plotDOI==1 || maskDOI==1) && input_vel==1
                DOI(ix)=zround(ix)-lmaxpick(ix)*doifact;
                indf(ix)=round((lmaxpick(ix)*doifact)/dz);
            end
            
            % DOI from VS standard deviation threshold
            if (plotDOI==2 || maskDOI==2) && input_vel==1
                if exist('vsstd_ok','var')==1 && isempty(vsstd_ok)==0
                    flipmoddepth = flipud(moddepth);
                    flipvsstd = flipud([vsstd_ok;vsstd_ok(end)]);
                    indhsd = find(flipvsstd<std_mask,1,'first');
%                     keyboard
                    if indhsd <= 2
                            indhsd = find(flipvsstd<median(flipvsstd)-0.5*std(flipvsstd),1,'first');
                            
%                             indhsd = find(flipvsstd<flipvsstd(1),1,'first');
%                             indmax = find(flipvsstd==max(flipvsstd),1);
%                             indhsd_all = find(flipvsstd<std_mask);
%                             indhsd = find(indhsd_all>indmax,1,'first');

                    elseif isempty(indhsd)
                        if ix>1 && ~isnan(DOI(ix-1))
                            indhsd = find(abs(flipmoddepth-(zround(ix)-DOI(ix-1))) == min(abs(flipmoddepth-(zround(ix)-DOI(ix-1)))))+2;
                        else
                            indhsd = find(flipvsstd<median(flipvsstd)-0.75*std(flipvsstd),1,'first');
                        end
                    end
                    if isempty(indhsd) || indhsd <= 2
                        indhsd = find(flipvsstd~=flipvsstd(1),1,'first');
                    end
                    
                    hsdtmp=flipmoddepth(indhsd-2);
                    if plotDOI==2 || maskDOI==2
                        DOI(ix)=zround(ix)-hsdtmp;
                    end
                    if maskDOI==2
                        indf(ix)=round((hsdtmp)/dz);
                    end
                else
                    if plotDOI==2 || maskDOI==2
                        DOI(ix)=zround(ix)-maxdepth;
                    end
                    if maskDOI==2
                        indf(ix)=round((maxdepth)/dz);
                    end
                end
            end
            
            % DOI from VS standard deviation threshold (experimental)
            if (plotDOI==3 || maskDOI==3) && input_vel==1
                                
                lam_rat = 0.1;
                depth_rat = 0.5;
                first_vs_rat = 0.95;
                
                flipmoddepth = flipud(moddepth);
                flipvssw = flipud([vssw;vssw(end)]);
                flipvsstd = (flipud([vsstd_ok;vsstd_ok(end)]));
                ind_hsdtmp = find(flipvssw/flipvssw(1)>first_vs_rat,1,'last');
                hsdtmp = flipmoddepth(ind_hsdtmp);
                ind_std = find(flipvsstd < std_mask);
                
                if ~isempty(ind_std > ind_hsdtmp) && any(flipvsstd(flipmoddepth> lmaxpick(ix)*lam_rat) > std_mask)
                    ind_max = find(flipvsstd == max(flipvsstd(flipmoddepth > lmaxpick(ix)*lam_rat)),1,'last');
                    ind_hsdtmp = ind_std(find(ind_std > ind_max,1,'first'))-1;
                    if flipmoddepth(ind_hsdtmp) > lmaxpick(ix)*lam_rat || ind_hsdtmp < depth_rat*length(flipmoddepth)
                        hsdtmp = flipmoddepth(ind_hsdtmp);
                    end
                else
                    ind_hsdtmp = find(flipvsstd == max(flipvsstd(flipmoddepth > lmaxpick(ix)*lam_rat)),1,'last');
                    if (flipmoddepth(ind_hsdtmp) > lmaxpick(ix)*lam_rat || ind_hsdtmp <= depth_rat*length(flipmoddepth)) && flipmoddepth(ind_hsdtmp) < 2*lam_rat*lmaxpick(ix)
                        if flipvssw(ind_hsdtmp) < flipvssw(find(flipvssw/flipvssw(1)>first_vs_rat,1,'last'))
                            hsdtmp = flipmoddepth(ind_hsdtmp);
                        end
                    else
                        hsdtmp = 2*lam_rat*lmaxpick(ix);
                    end
                end
                
                if isempty(hsdtmp)==1
                    hsdtmp=0;
                end
                if plotDOI==3 || maskDOI==3
                    DOI(ix)=zround(ix)-hsdtmp;
                end
                if maskDOI==3
                    indf(ix)=round((hsdtmp)/dz);
                end
                flipvsstd_filt = flipvsstd;
            end
            
            % DOI from VS standard deviation threshold
            if (plotDOI==4 || maskDOI==4) && input_vel==1
                if exist('vsstd_ok','var')==1 && isempty(vsstd_ok)==0
                    flipmoddepth = flipud(moddepth);
                    flipvsstd = flipud([vsstd_ok;vsstd_ok(end)]);
                    indhsd = find(flipvsstd<std_mask,1,'first');
                    if isempty(indhsd) || indhsd <= 2
                        if ix>1 && ~isnan(DOI(ix-1))
                            indhsd = find(abs(flipmoddepth-(zround(ix)-DOI(ix-1))) == min(abs(flipmoddepth-(zround(ix)-DOI(ix-1)))))+2;
                        else
                            indhsd = find(flipvsstd<flipvsstd(1),1,'first');
%                             indmax = find(flipvsstd==max(flipvsstd),1);
%                             indhsd_all = find(flipvsstd<std_mask);
%                             indhsd = find(indhsd_all>indmax,1,'first');
                            if isempty(indhsd) || indhsd <= 2
                                indhsd = find(flipvsstd~=flipvsstd(1),1,'first');
                            end
                        end
                    end
                    hsdtmp=flipmoddepth(indhsd-2);
                    if plotDOI==2 || maskDOI==2
                        DOI(ix)=zround(ix)-hsdtmp;
                    end
                    if maskDOI==2
                        indf(ix)=round((hsdtmp)/dz);
                    end
                else
                    if plotDOI==2 || maskDOI==2
                        DOI(ix)=zround(ix)-maxdepth;
                    end
                    if maskDOI==2
                        indf(ix)=round((maxdepth)/dz);
                    end
                end
            end
            
            % DOI from VS standard deviation threshold (xp)
            if (plotDOI==5 || maskDOI==5) && input_vel==1
                
                flipmoddepth = flipud(moddepth);
                flipvssw = flipud([vssw;vssw(end)]);
                flipvsstd = flipud([vsstd_ok;vsstd_ok(end)]);
                flipvsstd_filt = mov_aver(flipvsstd',5,1,length(flipvsstd));
                gradvsstd = mov_aver(gradient(flipvsstd_filt)',5,1,length(gradient(flipvsstd_filt)));
                ind_nochange = find(flipvssw/flipvssw(1)>0.95,1,'last');
                
                if abs(mean(gradvsstd)) < 0.01
                    gradvsstd = gradvsstd - mean(gradvsstd);
                    gradvsstd(abs(gradvsstd)<0.025) = 0;
                end
                
                ind_maxgrad = find(abs(gradvsstd)>0.01,1,'first');
                ind_first_sign = find(sign(gradvsstd)~=0,1,'first');
                first_sign = sign(gradvsstd(ind_first_sign));
                ind_sign = find(sign(gradvsstd) ~= first_sign);
                if isempty(ind_sign)
                    ind_sign = 1;
                end
                
                ind_first_max = ind_sign(find(ind_sign > ind_maxgrad,1,'first'));
                if isempty(ind_first_max)
                    ind_first_max = 1;
                end
                
                ind_depth_min = length(flipmoddepth)-round(0.2*length(flipmoddepth));
                ind_depth_min = max([ind_depth_min find(flipmoddepth > min(lmaxpick/10),1,'last')]);
                ind_max = find(flipvsstd_filt == max(flipvsstd_filt(ind_nochange:end)));
                
                if isempty(find(flipvsstd_filt(ind_nochange:ind_depth_min) >= std_mask, 1)) % no std higher than std_mask
                    
                    if max(flipvsstd_filt(ind_nochange:ind_depth_min)) >= 0.5*std_mask && flipvsstd_filt(1) < max(flipvsstd_filt(ind_nochange:ind_depth_min))
                        ind_hsdtmp_max = find(flipvsstd_filt <= 0.5*std_mask);
                        ind_hsdtmp_max = ind_hsdtmp_max(find(ind_hsdtmp_max > ind_first_max & ind_hsdtmp_max < ind_depth_min,1,'first'));
                        hsdtmp = flipmoddepth(ind_hsdtmp_max);
                    else
                        ind_hsdtmp_max = round(0.5*(ind_max + ind_nochange));
                        hsdtmp = flipmoddepth(ind_hsdtmp_max);
                    end
                    
                    if isempty(hsdtmp)
                        hsdtmp = flipmoddepth(ind_first_max);
                    end
                    
                else % some std higher than std_mask
                    
                    if flipvsstd_filt(1) < std_mask && first_sign==-1
                        ind_sign2 = find(sign(gradvsstd) == first_sign);
                        ind_second_max = ind_sign2(find(ind_sign2 > ind_first_max,1,'first'));
                        hsdtmp = flipmoddepth(ind_second_max);
                    else
                        ind_hsdtmp_max = find(flipvsstd_filt <= std_mask);
                        ind_hsdtmp_max = ind_hsdtmp_max(find(ind_hsdtmp_max > ind_first_max,1,'first'));
                        hsdtmp = flipmoddepth(ind_hsdtmp_max);
                    end
                end
                                
                if isempty(hsdtmp)==1
                    hsdtmp = flipmoddepth(ind_nochange);
                end
                if plotDOI==4 || maskDOI==4
                    DOI(ix)=zround(ix)-hsdtmp;
                end
                if maskDOI==4
                    indf(ix)=round((hsdtmp)/dz);
                end  
            end
            
            % DOI from max VS standard deviation
            if (plotDOI==6 || maskDOI==6) && input_vel==1
                if exist('vsstd_ok','var')==1 && isempty(vsstd_ok)==0
                    flipmoddepth = flipud(moddepth);
                    flipvsstd = flipud([vsstd_ok;vsstd_ok(end)]);
                    indhsd = find(flipvsstd==max(flipvsstd),1,'first');
                    if indhsd <= 2
                        indhsd = find(flipvsstd<median(flipvsstd)-0.5*std(flipvsstd),1,'first');
                    elseif isempty(indhsd)
                        if ix>1 && ~isnan(DOI(ix-1))
                            indhsd = find(abs(flipmoddepth-(zround(ix)-DOI(ix-1))) == min(abs(flipmoddepth-(zround(ix)-DOI(ix-1)))))+2;
                        else
                            indhsd = find(flipvsstd<median(flipvsstd)-0.75*std(flipvsstd),1,'first');
                        end
                    end
                    if isempty(indhsd) || indhsd <= 2
                        indhsd = find(flipvsstd~=flipvsstd(1),1,'first');
                    end
                    
                    hsdtmp=flipmoddepth(indhsd-2);
                    if plotDOI==6 || maskDOI==6
                        DOI(ix)=zround(ix)-hsdtmp;
                    end
                    if maskDOI==6
                        indf(ix)=round((hsdtmp)/dz);
                    end
                else
                    if plotDOI==6 || maskDOI==6
                        DOI(ix)=zround(ix)-maxdepth;
                    end
                    if maskDOI==6
                        indf(ix)=round((maxdepth)/dz);
                    end
                end
            end
            
            % DOI from max VS standard deviation
            if (plotDOI==7 || maskDOI==7) && input_vel==1
                if exist('vsstd_ok','var')==1 && isempty(vsstd_ok)==0
                    flipmoddepth = flipud(moddepth);
                    flipvsstd = flipud([vsstd_ok;vsstd_ok(end)]);
                    indhsd = find(flipvsstd<=0.9*median(flipvsstd),1,'first');
                    
                    
%                     if indhsd <= 2
%                         indhsd = find(flipvsstd<median(flipvsstd)-0.5*std(flipvsstd),1,'first');
%                     elseif isempty(indhsd)
%                         if ix>1 && ~isnan(DOI(ix-1))
%                             indhsd = find(abs(flipmoddepth-(zround(ix)-DOI(ix-1))) == min(abs(flipmoddepth-(zround(ix)-DOI(ix-1)))))+2;
%                         else
%                             indhsd = find(flipvsstd<median(flipvsstd)-0.75*std(flipvsstd),1,'first');
%                         end
%                     end
                    if isempty(indhsd) || indhsd <= 2
                        indhsd = find(flipvsstd~=flipvsstd(1),1,'first');
                    end
                    
                    hsdtmp=flipmoddepth(indhsd-2);
                    if plotDOI==7 || maskDOI==7
                        DOI(ix)=zround(ix)-hsdtmp;
                    end
                    if maskDOI==7
                        indf(ix)=round((hsdtmp)/dz);
                    end
                else
                    if plotDOI==7 || maskDOI==7
                        DOI(ix)=zround(ix)-maxdepth;
                    end
                    if maskDOI==7
                        indf(ix)=round((maxdepth)/dz);
                    end
                end
            end
            
            if DOI(ix)<zround(ix)-maxdepth
                DOI(ix)=zround(ix)-maxdepth;
            end
            
            % Look for topo index
%             crit=abs(zround(ix)-depth);
%             indi(ix)=find(crit==min(crit),1);
            
            % Fill matrices with velocities
            if input_vel==1
                vp0=velresamp(moddepth,[vpsw;vpsw(end)],zinc);
                vs0=velresamp(moddepth,[vssw;vssw(end)],zinc);
                rho0=velresamp(moddepth,[rhosw;rhosw(end)],zinc);
                vsstd0=velresamp(moddepth,[vsstd;vsstd(end)],zinc);
            elseif input_vel==2
                vp0=vptomo;
                vs0=vstomo;
            end
            if indf(ix)>length(depth) || isnan(indf(ix))==1
                indf(ix)=length(depth);
            end
            if indf(ix)>length(vp0)
                indf(ix)=length(vp0);
            end
            vpmat(1:length(vp0),ix)=vp0;
            vpmatok(1:length(vp0),ix)=vp0;
            if input_vel==1 || (input_vel==2 && isempty(VsItomo)==0)
                vsmat(1:length(vs0),ix)=vs0;
                vpvsmat(1:length(vs0),ix)=vp0./vs0;
                vptvsmat(1:length(vs0),ix)=vp0.*vs0;
                poismat(1:length(vs0),ix)=poisson(vp0,vs0);
            end
            if input_vel==1
                rhomat(1:length(rho0),ix)=rho0;
                vsstdmat(1:length(vsstd0),ix)=vsstd0;
            end
            
%             if hsdtmp>25
%                 figure(1);plot(flipmoddepth,flipvsstd_filt,'r-x');hold on; drawnow;
%                 line([hsdtmp hsdtmp],[min(flipvsstd_filt) max(flipvsstd_filt)]);
%                 line([min(flipmoddepth) max(flipmoddepth)],[std_mask std_mask]);
%                 set(gcf,'position',[53   412   720   550]);
%                 hold off;
%                 
%                 if (plotDOI==4 || maskDOI==4) && input_vel==1
%                     figure(2);plot(flipmoddepth,gradvsstd,'r-x'); hold on; drawnow;
%                     line([hsdtmp hsdtmp],[min(gradvsstd) max(gradvsstd)]);
%                     line([min(flipmoddepth) max(flipmoddepth)],[0 0]);
%                     set(gcf,'position',[773   412   720   550]);
%                     hold off;
%                 end
%                 
%                 figure(3); close(3);
%                 plot_img(3,XmidT,depth,vsstdmat,map7,axetop,0,cbpos,fs,'Distance (m)',...
%                     'Elevation (m)','Vs STD (%)',[xMIN xMAX],[zMIN zMAX],...
%                     [stdMIN stdMAX],xticks,zticks,stdticks,[],[],stdISO,[25 0 24 10],[],vertex,0);
%                 hold on; plot(XmidT,DOI,'-w'); hold off; drawnow;
%                 
%                 keyboard
%                 close all;
%             end
        end
    end
    
    %% %% %%
    
    %%%%% Compute dispersion residuals %%%%%
    if plot2dcal==1 || plothisto==1
        for ip=1:npvc
            if input_vel==1
                if exist(nametarg,'file')==2 && exist('vresamp','var')==1 && isempty(vresamp{modes(ip)+1})==0
                    vph2dobs{modes(ip)+1}(:,ix)=vresamp{modes(ip)+1}';
                    delta2dobs{modes(ip)+1}(:,ix)=deltaresamp{modes(ip)+1}';
                else
                    vph2dobs{modes(ip)+1}(:,ix)=NaN;
                    delta2dobs{modes(ip)+1}(:,ix)=NaN;
                end
            end
            if isempty(D)==0
                freqcal=D{modes(ip)+1,1}; % Frequency
                vcal=1./D{modes(ip)+1,2}; % Velocity
                if isnan(freqcal)==0
                    % Resample in lambda or frequency
                    [~,vresamp{modes(ip)+1}]=resampvel(freqcal,vcal,vcal,resampvec,sampling,1);
                    vph2dcal{modes(ip)+1}(:,ix)=vresamp{modes(ip)+1};
                else
                    vph2dcal{modes(ip)+1}(:,ix)=NaN;
                end
                if input_vel==1
                    vph2dcal{modes(ip)+1}(isnan(vph2dobs{modes(ip)+1}(:,ix)),ix)=NaN;
                    vph2dobs{modes(ip)+1}(isnan(vph2dcal{modes(ip)+1}(:,ix)),ix)=NaN;
                end
                vph2dobsall=[vph2dobsall;vph2dobs{modes(ip)+1}(:,ix)];
                delta2dobsall=[delta2dobsall;delta2dobs{modes(ip)+1}(:,ix)];
                vph2dcalall=[vph2dcalall;vph2dcal{modes(ip)+1}(:,ix)];
            else
                vph2dobsall=[vph2dobsall;vph2dobs{modes(ip)+1}(:,ix)*NaN];
                delta2dobsall=[delta2dobsall;vph2dobs{modes(ip)+1}(:,ix)*NaN];
                vph2dcalall=[vph2dcalall;vph2dobs{modes(ip)+1}(:,ix)*NaN];
            end
        end
        if sum(isnan(vph2dcalall))~=length(vph2dcalall)
            % Single Xmid misfit
            misfitall(ix)=sqrt(sum(((vph2dcalall(vph2dobsall>0)-vph2dobsall(vph2dobsall>0)).^2./...
                (length(vph2dobsall)*(delta2dobsall(vph2dobsall>0).^2)))));
            vph2dcalALL=[vph2dcalALL;vph2dcalall];
            vph2dobsALL=[vph2dobsALL;vph2dobsall];
            delta2dobsALL=[delta2dobsALL;delta2dobsall];
        else
            misfitall(ix)=NaN;
        end
    end
end

if input_vel==1
    %     keyboard
    %     DOI_old = DOI; indf_old=indf;
    %     DOI = DOI_old; ind = indf_old;
    
    if any(~isnan(DOI))
        npoints_med = 9;
        npoints_aver = 3;
        
        %         DOI(~isnan(DOI)) = mov_aver(DOI(~isnan(DOI))',npoints-2,1,length(DOI(~isnan(DOI))));
        DOI(~isnan(DOI)) = median_filt(DOI(~isnan(DOI))',npoints_med,1,length(DOI(~isnan(DOI))),1);
        DOI(~isnan(DOI)) = mov_aver(DOI(~isnan(DOI))',npoints_aver,1,length(DOI(~isnan(DOI))));
        
        %         indf(~isnan(DOI)) = round(mov_aver(indf(~isnan(DOI))',npoints-2,1,length(indf(~isnan(DOI)))));
        indf(~isnan(DOI)) = median_filt(indf(~isnan(DOI))',npoints_med,1,length(indf(~isnan(DOI))),1);
        indf(~isnan(DOI)) = ceil(mov_aver(indf(~isnan(DOI))',npoints_aver,1,length(indf(~isnan(DOI)))));
    end
        
    for ix = Xmidselec
        maskmat(1:indf(ix),ix) = ones(size(vs0(1:indf(ix))));
        maskmatvp(1:indf(ix),ix) = ones(size(vs0(1:indf(ix))));
    end
    
%     f1=plot_img([],XmidT,depth,vsstdmat.*maskmat,map7,axetop,0,cbpos,fs,'Distance (m)',...
%         'Elevation (m)','Vs STD (%)',[],[],[],[],[],[],[],[],[],[25 0 30 15],[],vertex,1);
%     if plottopo==1
%         hold on
%         plot(topo(:,1),topo(:,2),'k-','linewidth',2);
%     end
%     if input_vel==1 && length(find(isnan(DOI)==0))>1
%         hold on
%         %             dashline(XmidT,DOI_old,2,2,2,2,'color','k','linewidth',1.5);
%         dashline(XmidT,DOI,2,2,2,2,'color','r','linewidth',1.5);
%     end
end

if sum(modexist)==0 && (plot2dcal==1 || plot2dmod==1 || plothisto==1 || savexzv>0)
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n                No model existing with these settings');
    fprintf('\n         Check settings => "modeltype", "nbest" and "outpoints"');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
end

fprintf('\n  **********************************************************');
fprintf('\n  **********************************************************\n');

%% %% %%

%%%%% Plot 2D models %%%%%

% Initialization
if usevptomo==1 && input_vel==1
    avertype=[avertype,'_Vptomo'];
    for ix = 1:length(XmidT)
        poismat(depth(:,ix)<(zround(ix)-max_vp_depth(ix)),ix) = NaN;
        vpvsmat(depth(:,ix)<(zround(ix)-max_vp_depth(ix)),ix) = NaN;
        vptvsmat(depth(:,ix)<(zround(ix)-max_vp_depth(ix)),ix) = NaN;
        vpmat(depth(:,ix)<(zround(ix)-max_vp_depth(ix)),ix) = NaN;
    end
    if vpmask==0
        vpmatok=VpItomo;
        maskmatvp=ones(size(vpmatok));
    else
        vpmatok=vpmat;
    end
elseif input_vel==2
    vpmatok=VpItomo;
    maskmat=ones(size(maskmat));
    maskmatvp=ones(size(vpmatok));
    if isempty(VsItomo)==0
        vsmat=VsItomo;
        vpvsmat=vpmatok./vsmat;
        vptvsmat=vpmatok.*vsmat;
        poismat=poisson(vpmatok,vsmat);
        %     else
        %         vsmat=vpmatok*NaN;
    end
    avertype='user';
    modeltype='tomo';
end

if sum(modexist)>0
    
    %% %% %%
    
    %%%%% Plot phase velocity pseudo-sections and residuals %%%%%
    
    if plot2dcal==1 || plothisto==1
        fprintf('\n  Saving observed, calculated and residual phase velocity sections\n');
        nmodeinv=0;
        res_tot = []; vph_obs_tot = []; vph_cal_tot = [];
        
        %%%% Loop over all existing modes %%%%
        for ip=1:max(maxmodeinv)+1
            fprintf(['\n      Mode ',num2str(ip-1),'\n']);
            if sum(sum(isnan(vph2dcal{ip})))==numel(vph2dcal{ip})
                continue
            end
            nmodeinv=nmodeinv+1;
            if plot2dcal==1
                if input_vel==1
                    
                    %%% Plot observed phase velocity pseudo-section %%%
                    
                    if sampling==0
                        f1=plot_img(showplot,XmidT,resampvec,vph2dobs{ip},map1,axetop,0,cbpos,fs,'Distance (m)',...
                            freqtitle_short,'V_\phi obs. (m/s)',[xMIN xMAX],[0 max(resampvec)],...
                            [Vphmin Vphmax],xticks,lticks,vphticks,[],[],vphISO,[12 0 30 15],[],vertex,3);
                    else
                        f1=plot_img(showplot,XmidT,resampvec,vph2dobs{ip},map1,axetop,1,cbpos,fs,'Distance (m)',...
                            lamtitle,'V_\phi obs. (m/s)',[xMIN xMAX],[lamMIN lamMAX],...
                            [vphMIN vphMAX],xticks,lticks,vphticks,[],[],vphISO,[12 0 30 15],[],vertex,3);
                    end
                    sizeax=get(gca,'position');
                    fileobs=fullfile(dir_img_inv_2d,['Vphobs','.M',num2str(ip-1),'.',avertype,...
                        '.',modeltype,'.',imgform]);
%                     save_fig(f1,fileobs,imgform,imgres,1,1-testplot);
                    export_fig(fileobs,strcat('-r',num2str(imgres)));
                    if showplot==0
                        close(f1);
                    else
                        showplot=showplot+1;
                    end
                end
                
                %%% Plot calculated phase velocity pseudo-section %%%
                
                if sampling==0
                    f1=plot_img(showplot,XmidT,resampvec,vph2dcal{ip},map1,axetop/home/pasquet/Documents/Sites_OZCAR/Guadeloupe/SWIPprocess/Line5,0,cbpos,fs,'Distance (m)',...
                        freqtitle_short,'V_\phi calc. (m/s)',[xMIN xMAX],[0 max(resampvec)],...
                        [Vphmin Vphmax],xticks,lticks,vphticks,[],[],vphISO,[12 0 30 15],[],vertex,3);
                else
                    f1=plot_img(showplot,XmidT,resampvec,vph2dcal{ip},map1,axetop,1,cbpos,fs,'Distance (m)',...
                        lamtitle,'V_\phi calc. (m/s)',[xMIN xMAX],[lamMIN lamMAX],...
                        [vphMIN vphMAX],xticks,lticks,vphticks,[],[],vphISO,[12 0 30 15],[],vertex,3);
                end
                filecal=fullfile(dir_img_inv_2d,['Vphcalc','.M',num2str(ip-1),'.',avertype,...
                    '.',modeltype,'.',imgform]);
%                 save_fig(f1,filecal,imgform,imgres,1,1-testplot);
                export_fig(filecal,strcat('-r',num2str(imgres)));
                if showplot==0
                    close(f1);
                else
                    showplot=showplot+1;
                end
            end
            
            if input_vel==1
                
                %%% Plot phase velocity residuals pseudo-section %%%
                
                vph2dres{ip}=100*(vph2dobs{ip}-vph2dcal{ip})./vph2dobs{ip}; % Residuals in percent
                if plot2dcal==1
                    if sampling==0
                        f1=plot_img(showplot,XmidT,resampvec,vph2dres{ip},map4,axetop,0,cbpos,fs,'Distance (m)',...
                            freqtitle_short,'Residuals (%)',[xMIN xMAX],[0 max(resampvec)],...
                            [residMIN residMAX],xticks,lticks,residticks,[],[],[],[12 8.5 30 15],[],vertex,3);
                    else
                        f1=plot_img(showplot,XmidT,resampvec,vph2dres{ip},map4,axetop,1,cbpos,fs,'Distance (m)',...
                            lamtitle,'Residuals (%)',[xMIN xMAX],[lamMIN lamMAX],...
                            [residMIN residMAX],xticks,lticks,residticks,[],[],[],[12 8.5 30 15],[],vertex,3);
                    end
                    fileaux=fullfile(dir_img_inv_2d,['Vphres','.M',num2str(ip-1),'.',avertype,...
                        '.',modeltype,'.',imgform]);
%                     save_fig(f1,fileaux,imgform,imgres,1,1-testplot);
                    export_fig(fileaux,strcat('-r',num2str(imgres)));
                    if showplot==0
                        close(f1);
                    else
                        showplot=showplot+1;
                    end
                end
                
                % Compute RMS
                res=reshape(vph2dres{ip},size(vph2dres{ip},1)*size(vph2dres{ip},2),1);
                res=res(~isnan(res));
                res_tot=[res_tot;res];
                
                stdRMS=std(res);
                meanRMS=mean(res);
                goodRMS=res(res<meanRMS+2*stdRMS(end) & res>meanRMS-2*stdRMS(end));
                percent=fix(1000*length(goodRMS)/length(res))/10;
                
                vph_obs = reshape(vph2dobs{ip},size(vph2dobs{ip},1)*size(vph2dobs{ip},2),1);
                vph_cal = reshape(vph2dcal{ip},size(vph2dcal{ip},1)*size(vph2dcal{ip},2),1);
                vph_obs = vph_obs(~isnan(vph_obs));
                vph_cal = vph_cal(~isnan(vph_cal));
                vph_obs_tot=[vph_obs_tot;vph_obs];
                vph_cal_tot=[vph_cal_tot;vph_cal];
                
                % Gaussian curve
                edges=2*max(abs(min(goodRMS)),max(goodRMS));
                xgaus=linspace(-edges,edges,50); % Plotting range
                ygaus=0.5*length(res)*exp(-0.5*((xgaus-mean(res))/std(res)).^ 2)/(std(res)*sqrt(2*pi));
                
                %%% Plot residual histograms %%%
                
                if plothisto==1
                    if showplot==0
                        f2=figure(1);
                        set(f2,'visible','off')
                    else
                        f2=figure(showplot);
                    end
                    N=hist(res,xgaus);
                    hist(res(abs(res)<edges),xgaus,'linewidth',3);
                    hi=findobj(gca,'Type','patch');
                    set(hi,'FaceColor','r','EdgeColor','k');
                    xlabel('Residuals (%)'); ylabel('Nb of samples');
                    h=get(gca,'xlabel'); set(h,'FontSize',fs*1.75);
                    h=get(gca,'ylabel'); set(h,'FontSize',fs*1.75);
                    sizetick=get(gca,'ticklength');
                    set(gca,'TickDir','out','linewidth',3,'XMinorTick','on','YMinorTick','on',...
                        'ticklength',[sizetick(1)*3 sizetick(2)]);
                    h=findall(gcf,'Type','Axes'); set(h,'FontSize',fs*1.75);
                    xlim([-edges edges]);
                    set(f2,'Units','centimeters');
                    set(gcf,'Position',[25 0 18 18]);
                    hold on;
                    axis square
                    titr=sprintf(['mu = ',num2str(round(meanRMS*100)/100),' %%  |  sigma = ',num2str(round(stdRMS*10)/10),...
                        ' %%\n',num2str(length(goodRMS)),' / ',num2str(length(res)),' samples (',num2str(percent),' %%) < 2*sigma']);
                    title(titr,'fontsize',0.75*fs);
                    file1=fullfile(dir_img_inv_2d,['HistRes','.M',num2str(ip-1),'.',avertype,...
                        '.',modeltype,'.',imgform]);
%                     save_fig(f2,file1,imgform,imgres,1,1-testplot);
                    export_fig(file1,strcat('-r',num2str(imgres)));
                    if showplot==0
                        close(f2);
                    else
                        showplot=showplot+1;
                    end
                end
            end
            
            %%% Concatenate figures %%%
            
            if testplot==1 && plot2dcal==1 && input_vel~=2
                fileobs_unix = unix_wsl_path(fileobs,wsl);
                filecal_unix = unix_wsl_path(filecal,wsl);
                fileaux_unix = unix_wsl_path(fileaux,wsl);
                panel1=fullfile(dir_img_inv_2d,['Vph_Obs_Cal_Res','.M',num2str(ip-1),'.',avertype,...
                    '.',modeltype,'.',imgform]);
                panel1_unix = unix_wsl_path(panel1,wsl);
                cat_img([fileobs_unix,' ',filecal_unix,' ',fileaux_unix],imgform,1,[],panel1_unix,1);
                delete(fileobs,filecal,fileaux);
            end
            
            %%% Save phase velocities in .xzv files %%%
            
            if savexzv>0
                if input_vel==1
                    save_xzv(fullfile(dir_xzv_inv_mod,['Vphobs','.M',num2str(ip-1),'.',avertype,'.',modeltype,'.xzv']),...
                        XmidT,resampvec,vph2dobs{ip});
                end
                save_xzv(fullfile(dir_xzv_inv_mod,['Vphcalc','.M',num2str(ip-1),'.',avertype,'.',modeltype,'.xzv']),...
                    XmidT,resampvec,vph2dcal{ip},0);
            end
        end
        
        %%%% Save residual histograms figure %%%%
        
        if testplot==1 && plothisto==1 && input_vel~=2
            fprintf('\n  Saving residual histograms\n');
            if nmodeinv<=3
                columns=nmodeinv;
            else
                columns=ceil(nmodeinv/2);
            end
            panel2=fullfile(dir_img_inv_2d,['HistRes.',avertype,'.',modeltype,'.',imgform]);
            panel2_unix = unix_wsl_path(panel2,wsl);
            histres_file = fullfile(dir_img_inv_2d,['HistRes','.M*.',avertype,'.',modeltype,'.',imgform]);
            histres_file_unix = unix_wsl_path(histres_file,wsl);
            cat_img(histres_file_unix,imgform,columns,[],panel2_unix,1);
            delete(histres_file);
        end
        
        % Residuals phase velocities
        f1 = plot_curv(showplot,[vphMIN vphMAX],[vphMIN vphMAX],[],'-'); hold on;
        plot_scat(f1,vph_obs_tot,vph_cal_tot,res_tot,'o',[],1,map4,0,0,1,1.3*fs,'V_\phi obs. (m/s)','V_\phi calc. (m/s)','V_\phi residuals (%)',...
            [vphMIN vphMAX],[vphMIN vphMAX],[residMIN residMAX],vphticks,vphticks,residticks,[],[],[0 30 26 20],[],1);
        hold on;
        
        file1 = fullfile(dir_img_inv_2d,['PlotRes.all.',avertype,'.',modeltype,'.',imgform]);
        save_fig(f1,file1,imgform,imgres,1);
        if showplot==0
            close(f1);
        else
            showplot=showplot+1;
        end
        
        %%%% Save and display inversion QC %%%%
        
        if input_vel==1 && plot2dcal==1
            stdRMS=std(res_tot);
            meanRMS=mean(res_tot);
            goodRMS=res_tot(res_tot<meanRMS+2*stdRMS(end) & res_tot>meanRMS-2*stdRMS(end));
            percent=fix(1000*length(goodRMS)/length(res_tot))/10;
            fprintf('\n  Residual statistics (all modes included)\n');
            fprintf(['\n      mu = ',num2str(round(meanRMS*100)/100),' %%  |  sigma = ',num2str(round(stdRMS*10)/10),' %%']);
            fprintf(['\n      ',num2str(length(goodRMS)),' / ',num2str(length(res_tot)),' samples (',num2str(percent),' %%) < 2*sigma\n']);
            
            fprintf('\n  Saving inversion misfit graph\n');
            % Plot QC figure
            vph2dobsALL=vph2dobsALL(isnan(vph2dobsALL)==0);
            vph2dcalALL=vph2dcalALL(isnan(vph2dcalALL)==0);
            RMSfinal = sqrt(sum((vph2dcalALL-vph2dobsALL).^2)/length(vph2dcalALL));
            %             RMSfinal = meanabs(vph2dcalALL-vph2dobsALL);
            
            % Misfit for each Xmid
            f4=plot_curv(showplot,XmidT,misfitall,[],'.-',[0 0 0],[],axetop,0,0,fs,'Distance (m)',...
                'Misfit',[],[xMIN xMAX],[0 max([misfitall+0.05 1])],[],xticks,[],[],[],[],...
                [26 0 30 15],[],[]);
            sizeax2=get(gca,'position');
            set(findobj(f4,'Type','Axes'),'ActivePositionProperty','Position');
            set(findobj(f4,'Type','Axes'),'position',...
                [sizeax2(1),sizeax2(2),sizeax(3),sizeax2(4)/3]);
            xlabh=get(get(gca,'XLabel'),'extent');
            text(sizeax(1)+5,xlabh(2)/1.5,['Final RMS = ',...
                num2str(round(RMSfinal*10)/10),' m/s'],'FontSize',fs);
            
            % Save figure
            file1=fullfile(dir_img_inv_2d,['Misfit.',avertype,'.',modeltype,'.',imgform]);
            %             save_fig(f4,file1,imgform,imgres,1);
            export_fig(file1,strcat('-r',num2str(imgres)));
            if showplot==0
                close(f4);
            else
                showplot=showplot+1;
            end
        end
    end
    
    %% %% %%
    
    %%%%% Plot 2D models (Vs, Vp, VsSTD, Vp/Vs, Poisson, auxiliary) %%%%%
    
    if plot2dmod==1
        fprintf('\n  Saving 2D sections\n');
        
        % Initialization
        if plotiso==1
            specmat=vpmatok;
        elseif plotiso==2
            specmat=vsmat;
        elseif plotiso==3
            specmat=vsstdmat;
        elseif plotiso==4
            specmat=vpvsmat;
        elseif plotiso==5
            specmat=poismat;
        elseif plotiso==6
            specmat=auxmat;
        else
            specmat=[];
        end
        %         if input_vel==2 && plotiso>0
        %             specmat=flipud(specmat);
        %         end
        
        X_plot = repmat(XmidT,size(depth,1),1);
        if usevptomo == 0
            X_plot_vp = repmat(XmidT,size(depth,1),1);
        end
        
        %%
        %%%% Saving Vs section %%%%
        
        if isempty(zMIN) || isempty(zMAX)
            zMIN=floor(min(min(depth))/10)*10;
            zMAX=ceil(max(max(depth))/10)*10;
        end
        if length(find(isnan(DOI)==0))==1
            blocky=0;
        end
        if transpa == 0
            [f1,han1,~,~,c]=plot_img(showplot,X_plot,depth,vsmat.*maskmat,map5,axetop,0,cbpos,fs,'Distance (m)',...
                'Elevation (m)','Vs (m/s)',[xMIN xMAX],[zMIN zMAX],...
                [vsMIN vsMAX],xticks,zticks,vsticks,[],[],vsISO,[12 0 30 15],[],vertex,blocky);
        else
            [f1,han1,~,~,c]=plot_img(showplot,X_plot,depth,vsmat,map5,axetop,0,cbpos,fs,'Distance (m)',...
                'Elevation (m)','Vs (m/s)',[xMIN xMAX],[zMIN zMAX],...
                [vsMIN vsMAX],xticks,zticks,vsticks,[],[],vsISO,[12 0 30 15],[],vertex,blocky);
            
            limvs = get(c,'Ylim'); tickV = get(c,'Ytick');
            limx = get(gca,'Xlim'); tickX = get(gca,'Xtick');
            limz = get(gca,'Ylim'); tickZ = get(gca,'Ytick');
        end
        
        if transpa == 1
            mask_transp2 = maskmat;
            mask_transp2(isnan(mask_transp2)) = 10;
            mask_transp2(mask_transp2 == 1) = NaN;
            [f3,han3,~,~,c] = plot_img(0,X_plot,depth,mask_transp2,map5*0+1,axetop,0,cbpos,fs,'Distance (m)',...
                'Elevation (m)','Vs (m/s)',limx,limz,limvs,tickX,tickZ,tickV,[],[],[],...
                [12 0 30 15],[],vertex,3);
            set(gca, 'color', 'none');
            set(c,'visible','off');
%             set(cbhandle,'color','none');
%             set(get(cbhandle,'Children'),'visible','off');
%             set(get(cbhandle,'ylabel'),'visible','on')
%             mask_transp2 = maskmat;
%             mask_transp2(isnan(mask_transp2)) = 0.75;
%             set(han3,'AlphaDataMapping','none','AlphaData',mask_transp2,'facealpha','flat','edgealpha','flat');
            if plottopo==1
                hold on
                plot(topo(:,1),topo(:,2),'k-','linewidth',2);
            end
            if plotDOI>0 && input_vel==1 && length(find(isnan(DOI)==0))>1
                hold on
                plot(XmidT,DOI,'w-','linewidth',1.5);
            end
            if str2double(matrelease(1:4))>2014
                update_minortick(4);
            end
            
            filevsmask=fullfile(dir_img_inv_2d,['VSmask.',avertype,'.',modeltype,'.',imgform]);
            export_fig(f3,filevsmask,'-transparent',strcat('-r',num2str(imgres)));
            close(f3);
        end
        
        sizeax=get(gca,'Position');
        if plotiso>0 && isempty(specISO)==0
            hold on;
            if length(specISO)==1
                isoline=[specISO specISO];
            else
                isoline=specISO;
            end
            [cs,hc]=contour(X_plot,depth,specmat,isoline,'color',[0 0 0],'linewidth',1);
            clabel(cs, hc,'Color', 'k', 'Rotation', 0,'fontsize',12,'labelspacing', 500);
            hold off;
        end
        if plottopo==1
            hold on
            plot(topo(:,1),topo(:,2),'k-','linewidth',2);
        end
        if plotDOI>0 && input_vel==1 && length(find(isnan(DOI)==0))>1
            hold on
            plot(XmidT,DOI,'w-','linewidth',1.5);
        end
        if str2double(matrelease(1:4))>2014
            update_minortick(4);
        end
        
        filevs=fullfile(dir_img_inv_2d,['VS.',avertype,...
            '.',modeltype,'.',imgform]);
        filevs_unix = unix_wsl_path(filevs,wsl);

        if input_vel==1 || (input_vel==2 && isempty(VsItomo)==0)
            %             save_fig(f1,filevs,imgform,imgres,1,1-testplot);
            export_fig(filevs,strcat('-r',num2str(imgres)));
        end
        if showplot==0
            close(f1);
        else
            showplot=showplot+1;
        end
        if testimgmgck==0 && transpa == 1
            filevsmask_unix = unix_wsl_path(filevsmask,wsl);
            if input_vel ~= 2
                unix_cmd(sprintf('convert %s -alpha set -background none -channel A -evaluate multiply 0.75 +channel %s',filevsmask_unix,filevsmask_unix));
                if cbpos == 1
                    unix_cmd(sprintf('composite %s %s -gravity West %s',filevsmask_unix,filevs_unix,filevs_unix));
                else
                    unix_cmd(sprintf('composite %s %s -gravity NorthWest %s',filevsmask_unix,filevs_unix,filevs_unix));
                end
            delete(filevsmask);
            end
        end
        fprintf(['\n      Saved as ',filevs_unix,'\n']);
        
        %%
        %%%% Saving Vp section %%%%
        
        f1=plot_img(showplot,X_plot_vp,depth_vp,vpmatok.*maskmatvp,map5,axetop,0,cbpos,fs,'Distance (m)',...
            'Elevation (m)','Vp (m/s)',[xMIN xMAX],[zMIN zMAX],...
            [vpMIN vpMAX],xticks,zticks,vpticks,[],[],vpISO,[12 8.5 30 15],sizeax,vertex,blocky);
        if plotiso>0 && isempty(specISO)==0
            hold on;
            if length(specISO)==1
                isoline=[specISO specISO];
            else
                isoline=specISO;
            end
            [cs,hc]=contour(X_plot,depth,specmat,isoline,'color',[0 0 0],'linewidth',1);
            clabel(cs, hc,'Color', 'k', 'Rotation', 0,'fontsize',12,'labelspacing', 500);
            hold off;
        end
        if plottopo==1
            hold on
            plot(topo(:,1),topo(:,2),'k-','linewidth',2);
        end
        if plotDOI>0 && input_vel==1 && length(find(isnan(DOI)==0))>1
            hold on
%             dashline(XmidT,DOI,2,2,2,2,'color','w','linewidth',1.5);
        end
        if str2double(matrelease(1:4))>2014
            update_minortick(4);
        end
        
        filevp=fullfile(dir_img_inv_2d,['VP.',avertype,...
            '.',modeltype,'.',imgform]);
        filevp_unix = unix_wsl_path(filevp,wsl);
        %         save_fig(f1,filevp,imgform,imgres,1,1-testplot);
        export_fig(filevp,strcat('-r',num2str(imgres)));
        if showplot==0
            close(f1);
        else
            showplot=showplot+1;
        end
        fprintf(['\n      Saved as ',filevp_unix,'\n']);
        
        %%
        %%%% Saving VsStd section %%%%
        
        if input_vel==1
            if vsstd_test == 1
                vsstd_plot = vsstdmat;
                vsstd_title = 'Vs standard deviation (m/s)';
            elseif vsstd_test == 2
                vsstd_plot = 100*(vsstdmat./vsmat);
                vsstd_title = 'Vs standard deviation (%)';
            end
            f1=plot_img(showplot,X_plot,depth,vsstd_plot,map7,axetop,0,cbpos,fs,'Distance (m)',...
                'Elevation (m)',vsstd_title,[xMIN xMAX],[zMIN zMAX],...
                [stdMIN stdMAX],xticks,zticks,stdticks,[],[],stdISO,[25 0 30 15],sizeax,vertex,blocky);
            if plottopo==1
                hold on
                plot(topo(:,1),topo(:,2),'k-','linewidth',2);
            end
            if input_vel==1 && length(find(isnan(DOI)==0))>1
                hold on
                dashline(XmidT,DOI,2,2,2,2,'color','w','linewidth',1.5);
            end
            filestd=fullfile(dir_img_inv_2d,['VSstd.',avertype,...
                '.',modeltype,'.',imgform]);
            filestd_unix = unix_wsl_path(filestd,wsl);
            %             save_fig(f1,filestd,imgform,imgres,1,1-testplot);
            export_fig(f1,filestd,strcat('-r',num2str(imgres)));
            if showplot==0
                close(f1);
            else
                showplot=showplot+1;
            end
            fprintf(['\n      Saved as ',filestd_unix,'\n']);
            
            
            f1=plot_img(showplot,X_plot,depth,vsmat,map5,axetop,0,cbpos,fs,'Distance (m)',...
                'Elevation (m)','Vs (m/s)',[xMIN xMAX],[zMIN zMAX],...
                [vsMIN vsMAX],xticks,zticks,vsticks,[],[],vsISO,[12 0 30 15],[],vertex,blocky);
            if plottopo==1
                hold on
                plot(topo(:,1),topo(:,2),'k-','linewidth',2);
            end
            if input_vel==1 && length(find(isnan(DOI)==0))>1
                hold on
                dashline(XmidT,DOI,2,2,2,2,'color','w','linewidth',1.5);
            end
            filevs2=fullfile(dir_img_inv_2d,['VS2.',avertype,...
                '.',modeltype,'.',imgform]);
            filevs2_unix = unix_wsl_path(filevs2,wsl);
            if input_vel==1 || (input_vel==2 && isempty(VsItomo)==0)
                %                 save_fig(f1,filevs2,imgform,imgres,1,1-testplot);
                export_fig(filevs2,strcat('-r',num2str(imgres)));
                
            end
            if showplot==0
                close(f1);
            else
                showplot=showplot+1;
            end
            fprintf(['\n      Saved as ',filevs2_unix,'\n']);
        end
        
        %%
        %%%% Saving Vp/Vs section %%%%
        
        if input_vel==1 || (input_vel==2 && isempty(VsItomo)==0)
            if transpa == 0
                [f1,han1,~,~,c]=plot_img(showplot,X_plot,depth,vpvsmat.*maskmat,map6,axetop,0,cbpos,fs,'Distance (m)',...
                    'Elevation (m)','Vp/Vs',[xMIN xMAX],[zMIN zMAX],...
                    [vpvsMIN vpvsMAX],xticks,zticks,vpvsticks,[],[],vpvsISO,[25 8.5 30 15],sizeax,vertex,blocky);
            else
                [f1,han1,~,~,c]=plot_img(showplot,X_plot,depth,vpvsmat,map6,axetop,0,cbpos,fs,'Distance (m)',...
                    'Elevation (m)','Vp/Vs',[xMIN xMAX],[zMIN zMAX],...
                    [vpvsMIN vpvsMAX],xticks,zticks,vpvsticks,[],[],vpvsISO,[25 8.5 30 15],sizeax,vertex,blocky);
                limvpvs = get(c,'Ylim'); tickVPVS = get(c,'Ytick');
                limx = get(gca,'Xlim'); tickX = get(gca,'Xtick');
                limz = get(gca,'Ylim'); tickZ = get(gca,'Ytick');
            end
            
            if transpa == 1
                mask_transp2 = maskmat;
                mask_transp2(isnan(mask_transp2)) = 10;
                mask_transp2(mask_transp2 == 1) = NaN;
                [f3,han3,~,~,c] = plot_img(0,X_plot,depth,mask_transp2,map6*0+1,axetop,0,cbpos,fs,'Distance (m)',...
                    'Elevation (m)','Vp/Vs',limx,limz,limvpvs,tickX,tickZ,tickVPVS,[],[],[],...
                    [25 8.5 30 15],sizeax,vertex,3);
                set(gca, 'color', 'none');
                set(c,'visible','off');
%                 set(cbhandle,'color','none');
%                 set(get(cbhandle,'Children'),'visible','off');
%                 set(get(cbhandle,'ylabel'),'visible','on')
%                 mask_transp2 = maskmat;
%                 mask_transp2(isnan(mask_transp2)) = 0.75;
%                 set(han3,'AlphaDataMapping','none','AlphaData',mask_transp2,'facealpha','flat','edgealpha','flat');
                if plottopo==1
                    hold on
                    plot(topo(:,1),topo(:,2),'k-','linewidth',2);
                end
                if plotDOI>0 && input_vel==1 && length(find(isnan(DOI)==0))>1
                    hold on
                    plot(XmidT,DOI,'w-','linewidth',1.5);
                end
                
                filevsmask=fullfile(dir_img_inv_2d,['VSmask.',avertype,'.',modeltype,'.',imgform]);
                export_fig(f3,filevsmask,'-transparent',strcat('-r',num2str(imgres)));
                close(f3);
            end
            
            if plotiso>0 && isempty(specISO)==0
                hold on;
                if length(specISO)==1
                    isoline=[specISO specISO];
                else
                    isoline=specISO;
                end
                [cs,hc]=contour(X_plot,depth,specmat,isoline,'color',[0 0 0],'linewidth',1);
                clabel(cs, hc,'Color', 'k', 'Rotation', 0,'fontsize',12,'labelspacing', 500);
                hold off;
            end
            if plottopo==1
                hold on
                plot(topo(:,1),topo(:,2),'k-','linewidth',2);
            end
            if plotDOI>0 && input_vel==1 && length(find(isnan(DOI)==0))>1
                hold on
                plot(XmidT,DOI,'w-','linewidth',1.5);
            end
            filevpvs=fullfile(dir_img_inv_2d,['VPVS.',avertype,...
                '.',modeltype,'.',imgform]);
            filevpvs_unix = unix_wsl_path(filevpvs,wsl);
            %             save_fig(f1,filevpvs,imgform,imgres,1,1-testplot);
            export_fig(filevpvs,strcat('-r',num2str(imgres)));
            if showplot==0
                close(f1);
            else
                showplot=showplot+1;
            end
            if testimgmgck==0 && transpa == 1
                filevsmask_unix = unix_wsl_path(filevsmask,wsl);
                if input_vel ~= 2
                    unix_cmd(sprintf('convert %s -alpha set -background none -channel A -evaluate multiply 0.75 +channel %s',filevsmask_unix,filevsmask_unix));
                    if cbpos == 1
                        unix_cmd(sprintf('composite %s %s -gravity West %s',filevsmask_unix,filevpvs_unix,filevpvs_unix));
                    else
                        unix_cmd(sprintf('composite %s %s -gravity NorthWest %s',filevsmask_unix,filevpvs_unix,filevpvs_unix));
                    end
                end
                delete(filevsmask);
            end
            fprintf(['\n      Saved as ',filevpvs_unix,'\n']);
        end
        
        %%
        %%%% Saving Vp * Vs section %%%%
        
        if input_vel==1 || (input_vel==2 && isempty(VsItomo)==0)
            if transpa == 0
                [f1,han1,~,~,c]=plot_img(showplot,X_plot,depth,log10(vptvsmat).*maskmat,map6,axetop,0,cbpos,fs,'Distance (m)',...
                    'Elevation (m)','log10(Vp*Vs)',[xMIN xMAX],[zMIN zMAX],...
                    [vptvsMIN vptvsMAX],xticks,zticks,vptvsticks,[],[],vptvsISO,[25 8.5 30 15],sizeax,vertex,blocky);
            else
                [f1,han1,~,~,c]=plot_img(showplot,X_plot,depth,log10(vptvsmat),map6,axetop,0,cbpos,fs,'Distance (m)',...
                    'Elevation (m)','log10(Vp*Vs)',[xMIN xMAX],[zMIN zMAX],...
                    [vptvsMIN vptvsMAX],xticks,zticks,vptvsticks,[],[],vptvsISO,[25 8.5 30 15],sizeax,vertex,blocky);
                pause(0.1);
                limx = get(gca,'Xlim'); tickX = get(gca,'Xtick');
                limz = get(gca,'Ylim'); tickZ = get(gca,'Ytick');
                limvptvs = get(c,'Ylim'); tickVPtVS = get(c,'Ytick');
            end
            
            if transpa == 1
                mask_transp2 = maskmat;
                mask_transp2(isnan(mask_transp2)) = 10;
                mask_transp2(mask_transp2 == 1) = NaN;
                [f3,han3,~,~,c] = plot_img(0,X_plot,depth,mask_transp2,map6*0+1,axetop,0,cbpos,fs,'Distance (m)',...
                    'Elevation (m)','log10(Vp*Vs)',limx,limz,limvptvs,tickX,tickZ,tickVPtVS,[],[],[],...
                    [25 8.5 30 15],sizeax,vertex,3);
                set(gca, 'color', 'none');
                set(c,'visible','off');
%                 set(cbhandle,'color','none');
%                 set(get(cbhandle,'Children'),'visible','off');
%                 set(get(cbhandle,'ylabel'),'visible','on')
%                 mask_transp2 = maskmat;
%                 mask_transp2(isnan(mask_transp2)) = 0.75;
%                 set(han3,'AlphaDataMapping','none','AlphaData',mask_transp2,'facealpha','flat','edgealpha','flat');
                if plottopo==1
                    hold on
                    plot(topo(:,1),topo(:,2),'k-','linewidth',2);
                end
                if plotDOI>0 && input_vel==1 && length(find(isnan(DOI)==0))>1
                    hold on
                    plot(XmidT,DOI,'w-','linewidth',1.5);
                end
                
                filevsmask=fullfile(dir_img_inv_2d,['VSmask.',avertype,'.',modeltype,'.',imgform]);
                export_fig(f3,filevsmask,'-transparent',strcat('-r',num2str(imgres)));
                close(f3);
            end
            
            if plotiso>0 && isempty(specISO)==0
                hold on;
                if length(specISO)==1
                    isoline=[specISO specISO];
                else
                    isoline=specISO;
                end
                [cs,hc]=contour(X_plot,depth,specmat,isoline,'color',[0 0 0],'linewidth',1);
                clabel(cs, hc,'Color', 'k', 'Rotation', 0,'fontsize',12,'labelspacing', 500);
                hold off;
            end
            if plottopo==1
                hold on
                plot(topo(:,1),topo(:,2),'k-','linewidth',2);
            end
            if plotDOI>0 && input_vel==1 && length(find(isnan(DOI)==0))>1
                hold on
                plot(XmidT,DOI,'w-','linewidth',1.5);
            end
            filevptvs=fullfile(dir_img_inv_2d,['VPtVS.',avertype,...
                '.',modeltype,'.',imgform]);
            filevptvs_unix = unix_wsl_path(filevptvs,wsl);
            %             save_fig(f1,filevpvs,imgform,imgres,1,1-testplot);
            export_fig(filevptvs,strcat('-r',num2str(imgres)));
            if showplot==0
                close(f1);
            else
                showplot=showplot+1;
            end
            if testimgmgck==0 && transpa == 1
                if input_vel ~= 2
                    filevsmask_unix = unix_wsl_path(filevsmask,wsl);
                    unix_cmd(sprintf('convert %s -alpha set -background none -channel A -evaluate multiply 0.75 +channel %s',filevsmask,filevsmask));
                    if cbpos == 1
                        unix_cmd(sprintf('composite %s %s -gravity West %s',filevsmask_unix,filevptvs_unix,filevptvs_unix));
                    else
                        unix_cmd(sprintf('composite %s %s -gravity NorthWest %s',filevsmask_unix,filevptvs_unix,filevptvs_unix));
                    end
                end
                delete(filevsmask);
            end
            fprintf(['\n      Saved as ',filevptvs_unix,'\n']);
        end
        
        %%
        %%%% Saving Poisson's ratio section %%%%
        
        if input_vel==1 || (input_vel==2 && isempty(VsItomo)==0)
            if transpa == 0
                [f1,han1,~,~,c]=plot_img(showplot,X_plot,depth,poismat.*maskmat,map6,axetop,0,cbpos,fs,'Distance (m)',...
                    'Elevation (m)','Poisson''s ratio',[xMIN xMAX],[zMIN zMAX],...
                    [poisMIN poisMAX],xticks,zticks,poisticks,[],[],poisISO,[25 16 30 15],sizeax,vertex,blocky);
            else
                [f1,han1,~,~,c]=plot_img(showplot,X_plot,depth,poismat,map6,axetop,0,cbpos,fs,'Distance (m)',...
                    'Elevation (m)','Poisson''s ratio',[xMIN xMAX],[zMIN zMAX],...
                    [poisMIN poisMAX],xticks,zticks,poisticks,[],[],poisISO,[25 16 30 15],sizeax,vertex,blocky);
                limpois = get(c,'Ylim'); tickPOIS = get(c,'Ytick');
                limx = get(gca,'Xlim'); tickX = get(gca,'Xtick');
                limz = get(gca,'Ylim'); tickZ = get(gca,'Ytick');
            end
            
            if transpa == 1
                mask_transp2 = maskmat;
                mask_transp2(isnan(mask_transp2)) = 10;
                mask_transp2(mask_transp2 == 1) = NaN;
                [f3,han3,~,~,c] = plot_img(0,X_plot,depth,mask_transp2,map6*0+1,axetop,0,cbpos,fs,'Distance (m)',...
                    'Elevation (m)','Poisson''s ratio',limx,limz,limpois,tickX,tickZ,tickPOIS,[],[],[],...
                    [25 16 30 15],sizeax,vertex,3);
                set(gca, 'color', 'none');
                set(c,'visible','off');
%                 set(cbhandle,'color','none');
%                 set(get(cbhandle,'Children'),'visible','off');
%                 set(get(cbhandle,'ylabel'),'visible','on')
%                 mask_transp2 = maskmat;
%                 mask_transp2(isnan(mask_transp2)) = 0.75;
%                 set(han3,'AlphaDataMapping','none','AlphaData',mask_transp2,'facealpha','flat','edgealpha','flat');
                if plottopo==1
                    hold on
                    plot(topo(:,1),topo(:,2),'k-','linewidth',2);
                end
                if plotDOI>0 && input_vel==1 && length(find(isnan(DOI)==0))>1
                    hold on
                    plot(XmidT,DOI,'w-','linewidth',1.5);
                end
                
                filevsmask=fullfile(dir_img_inv_2d,['VSmask.',avertype,'.',modeltype,'.',imgform]);
                export_fig(f3,filevsmask,'-transparent',strcat('-r',num2str(imgres)));
                close(f3);
            end
            
            if plotiso>0 && isempty(specISO)==0
                hold on;
                if length(specISO)==1
                    isoline=[specISO specISO];
                else
                    isoline=specISO;
                end
                [cs,hc]=contour(X_plot,depth,specmat,isoline,'color',[0 0 0],'linewidth',1);
                clabel(cs, hc,'Color', 'k', 'Rotation', 0,'fontsize',12,'labelspacing', 500);
                hold off;
            end
            if plottopo==1
                hold on
                plot(topo(:,1),topo(:,2),'k-','linewidth',2);
            end
            if plotDOI>0 && input_vel==1 && length(find(isnan(DOI)==0))>1
                hold on
                plot(XmidT,DOI,'w-','linewidth',1.5);
            end
            filepois=fullfile(dir_img_inv_2d,['Poisson.',avertype,...
                '.',modeltype,'.',imgform]);
            filepois_unix = unix_wsl_path(filepois,wsl);
            %             save_fig(f1,filepois,imgform,imgres,1,1-testplot);
            export_fig(filepois,strcat('-r',num2str(imgres)));
            if showplot==0
                close(f1);
            else
                showplot=showplot+1;
            end
            if testimgmgck==0 && transpa == 1
                filevsmask_unix = unix_wsl_path(filevsmask,wsl);
                if input_vel ~= 2
                    unix_cmd(sprintf('convert %s -alpha set -background none -channel A -evaluate multiply 0.75 +channel %s',filevsmask_unix,filevsmask_unix));
                    if cbpos == 1
                        unix_cmd(sprintf('composite %s %s -gravity West %s',filevsmask_unix,filepois_unix,filepois_unix));
                    else
                        unix_cmd(sprintf('composite %s %s -gravity NorthWest %s',filevsmask_unix,filepois_unix,filepois_unix));
                    end
                end
                delete(filevsmask);
            end
            fprintf(['\n      Saved as ',filepois_unix,'\n']);
        end
        %%
        
        %%%% Concatenate figures %%%%
        
        if testplot==1
            if input_vel==1 || (input_vel==2 && isempty(VsItomo)==0)
                if input_aux==1 && exist('AuxI','var')==1
                    dispmsg=0;
                else
                    dispmsg=1;
                end
                panel3=fullfile(dir_img_inv_2d,['VP_VS_Pois.',avertype,'.',modeltype,'.',imgform]);
                panel3_unix = unix_wsl_path(panel3,wsl);
                cat_img([filevp_unix,' ',filevs_unix,' ',filepois_unix],imgform,1,[],panel3_unix,dispmsg);
                fprintf(['\n      Saved as ',panel3_unix,'\n']);
                
                panel4=fullfile(dir_img_inv_2d,['VP_VS_VPVS.',avertype,'.',modeltype,'.',imgform]);
                panel4_unix = unix_wsl_path(panel4,wsl);
                cat_img([filevp_unix,' ',filevs_unix,' ',filevpvs_unix],imgform,1,[],panel4_unix,dispmsg);
                fprintf(['\n      Saved as ',panel4_unix,'\n']);
            end
            if input_vel==1
                panel5=fullfile(dir_img_inv_2d,['VSstd_VS.',avertype,'.',modeltype,'.',imgform]);
                panel5_unix = unix_wsl_path(panel5,wsl);
                cat_img([filestd_unix,' ',filevs2_unix],imgform,1,[],panel5_unix,1);
                fprintf(['\n      Saved as ',panel5_unix,'\n']);
            end
        end
    end
end

%% %% %%

%%%%% Saving auxiliary data section %%%%%

if plot2dmod==1 && exist('auxmat','var')==1 && input_aux==1 && isempty(auxmat)==0
    if input_vel==2
        fileAUX=fullfile(dir_img_inv_2d,['AUXI.user.tomo.',imgform]);
    else
        fileAUX=fullfile(dir_img_inv_2d,['AUXI.',avertype,'.',modeltype,'.',imgform]);
        if maskDOI>0 && auxmask==1
            auxmat(isnan(maskmat)==1)=NaN;
        end
    end
    fileAUX_unix = unix_wsl_path(fileAUX,wsl);
    if input_vel==2 || (input_vel==1 && sum(modexist)>0)
        % Plot auxiliary data section
        if exist('sizeax','var')~=1
            sizeax=[];
        end
        if auxlogscal==1
            f1=plot_img_log(showplot,XmidT_aux,depth_aux,auxmat,map8,axetop,0,cbpos,fs,'Distance (m)',...
                'Elevation (m)',auxtitle,[xMIN xMAX],[zMIN zMAX],...
                [auxMIN auxMAX],xticks,zticks,auxticks,[],[],auxISO,[25 16 30 15],sizeax,vertex,blocky);
        else
            f1=plot_img(showplot,XmidT_aux,depth_aux,auxmat,map8,axetop,0,cbpos,fs,'Distance (m)',...
                'Elevation (m)',auxtitle,[xMIN xMAX],[zMIN zMAX],...
                [auxMIN auxMAX],xticks,zticks,auxticks,[],[],auxISO,[25 16 30 15],sizeax,vertex,blocky);
        end
        if plotiso>0 && isempty(specISO)==0
            hold on;
            if length(specISO)==1
                isoline=[specISO specISO];
            else
                isoline=specISO;
            end
            [cs,hc]=contour(X_plot,depth,specmat,isoline,'color',[0 0 0],'linewidth',1);
            clabel(cs, hc,'Color', 'k', 'Rotation', 0,'fontsize',12,'labelspacing', 500);
            hold off;
        end
        if plottopo==1
            hold on
            plot(topo(:,1),topo(:,2),'k-','linewidth',2);
        end
%         if plotDOI>0 && input_vel==1 && length(find(isnan(DOI)==0))>1
%             hold on
%             plot(XmidT,DOI,'w-','linewidth',1.5);
%         end
        %         save_fig(f1,fileAUX,imgform,imgres,1,1-testplot);
        export_fig(fileAUX,strcat('-r',num2str(imgres)));
        if showplot==0
            close(f1);
        else
            showplot=showplot+1;
        end
    end
    
    %%%% Concatenate figures %%%%
    
    if testplot==1
        if (input_vel==1 && sum(modexist)>0) || (input_vel==2 && isempty(VsItomo)==0)
            panel6=fullfile(dir_img_inv_2d,['VP_VS_Pois_Aux.',avertype,'.',modeltype,'.',imgform]);
            panel6_unix = unix_wsl_path(panel6,wsl);
            cat_img([filevp_unix,' ',filevs_unix,' ',filepois_unix,' ',fileAUX_unix],imgform,1,[],panel6_unix,1);
            delete(panel3);
            panel7=fullfile(dir_img_inv_2d,['VP_VS_VPVS_Aux.',avertype,'.',modeltype,'.',imgform]);
            panel7_unix = unix_wsl_path(panel7,wsl);
            cat_img([filevp_unix,' ',filevs_unix,' ',filevpvs_unix,' ',fileAUX_unix],imgform,1,[],panel7_unix,1);
            delete(panel4);
        elseif input_vel==2 && isempty(VpItomo)==0 && isempty(VsItomo)==1
            panel8=fullfile(dir_img_inv_2d,['Aux_VP.',avertype,'.',modeltype,'.',imgform]);
            panel8_unix = unix_wsl_path(panel8,wsl);
            cat_img([fileAUX_unix,' ',filevp_unix],imgform,1,[],panel8_unix,1);
        elseif input_vel==2 && isempty(VpItomo)==0 && isempty(VsItomo)==0
            panel9=fullfile(dir_img_inv_2d,['Aux_VP_VS.',avertype,'.',modeltype,'.',imgform]);
            panel9_unix = unix_wsl_path(panel9,wsl);
            cat_img([fileAUX_unix,' ',filevp_unix,' ',filevs_unix],imgform,1,[],panel9_unix,1);
        end
    end
end

%%%% Delete temp files %%%%

if testplot==1 && plot2dmod==1 && concat==1
    if exist('filevp','var')==1 && exist(filevp,'file')==2
        if input_vel==1 || (input_vel==2 && isempty(VsItomo)==0)
            delete(filevp);
        end
    end
    if exist('filevs','var')==1 && exist(filevs,'file')==2
        delete(filevs);
    end
    if exist('filevs2','var')==1 && exist(filevs2,'file')==2
        delete(filevs2);
    end
    if exist('filepois','var')==1 && exist(filepois,'file')==2
        delete(filepois);
    end
    if exist('filevpvs','var')==1 && exist(filevpvs,'file')==2
        delete(filevpvs);
    end
    if exist('filestd','var')==1 && exist(filestd,'file')==2
        delete(filestd);
    end
    if exist('fileAUX','var')==1 && exist(fileAUX,'file')==2
        delete(fileAUX);
    end
end

%% %% %%

%%%% Save models in .xzv files %%%%

if savexzv>0 && input_vel==1
    if savexzv==1
        save_xzv(fullfile(dir_xzv_inv_mod,['VS.',avertype,'.',modeltype,'.xzv']),XmidT,depth,vsmat.*maskmat);
        save_xzv(fullfile(dir_xzv_inv_mod,['VP.',avertype,'.',modeltype,'.xzv']),XmidT,depth,vpmat.*maskmat);
        save_xzv(fullfile(dir_xzv_inv_mod,['VPVS.',avertype,'.',modeltype,'.xzv']),XmidT,depth,vpvsmat.*maskmat);
        save_xzv(fullfile(dir_xzv_inv_mod,['VPtVS.',avertype,'.',modeltype,'.xzv']),XmidT,depth,vptvsmat.*maskmat);
        save_xzv(fullfile(dir_xzv_inv_mod,['Poisson.',avertype,'.',modeltype,'.xzv']),XmidT,depth,poismat.*maskmat);
        save_xzv(fullfile(dir_xzv_inv_mod,['VSstd.',avertype,'.',modeltype,'.xzv']),XmidT,depth,vsstdmat.*maskmat);
        if input_aux==1 && auxmask==1
            auxmat(isnan(auxmat)==1 & isnan(vpmat)==0)=min(auxmat(isnan(auxmat)==0));
            save_xzv(fullfile(dir_xzv_inv_mod,['Aux.',avertype,'.',modeltype,'.xzv']),XmidT,depth,auxmat.*maskmat);
        end
    elseif savexzv==2
        save_xzv(fullfile(dir_xzv_inv_mod,['VS.',avertype,'.',modeltype,'.xzv']),XmidT,depth,vsmat.*maskmat,0);
        save_xzv(fullfile(dir_xzv_inv_mod,['VP.',avertype,'.',modeltype,'.xzv']),XmidT,depth,vpmat.*maskmat,0);
        save_xzv(fullfile(dir_xzv_inv_mod,['VPVS.',avertype,'.',modeltype,'.xzv']),XmidT,depth,vpvsmat.*maskmat,0);
        save_xzv(fullfile(dir_xzv_inv_mod,['VPtVS.',avertype,'.',modeltype,'.xzv']),XmidT,depth,vptvsmat.*maskmat,0);
        save_xzv(fullfile(dir_xzv_inv_mod,['Poisson.',avertype,'.',modeltype,'.xzv']),XmidT,depth,poismat.*maskmat,0);
        save_xzv(fullfile(dir_xzv_inv_mod,['VSstd.',avertype,'.',modeltype,'.xzv']),XmidT,depth,vsstdmat.*maskmat,0);
        if input_aux==1 && auxmask==1
            auxmat(isnan(auxmat)==1 & isnan(vpmat)==0)=min(auxmat(isnan(auxmat)==0));
            save_xzv(fullfile(dir_xzv_inv_mod,['Aux.',avertype,'.',modeltype,'.xzv']),XmidT,depth,auxmat.*maskmat,0);
        end
    elseif savexzv==3
        save_xzv(fullfile(dir_xzv_inv_mod,['VS.',avertype,'.',modeltype,'.xzv']),XmidT,depth,vsmat,0,maskmat);
        save_xzv(fullfile(dir_xzv_inv_mod,['VP.',avertype,'.',modeltype,'.xzv']),XmidT,depth,vpmat,0,maskmat);
        save_xzv(fullfile(dir_xzv_inv_mod,['VPVS.',avertype,'.',modeltype,'.xzv']),XmidT,depth,vpvsmat,0,maskmat);
        save_xzv(fullfile(dir_xzv_inv_mod,['VPtVS.',avertype,'.',modeltype,'.xzv']),XmidT,depth,vptvsmat,0,maskmat);
        save_xzv(fullfile(dir_xzv_inv_mod,['Poisson.',avertype,'.',modeltype,'.xzv']),XmidT,depth,poismat,0,maskmat);
        save_xzv(fullfile(dir_xzv_inv_mod,['VSstd.',avertype,'.',modeltype,'.xzv']),XmidT,depth,vsstdmat,0,maskmat);
        if input_aux==1 && auxmask==1
            auxmat(isnan(auxmat)==1 & isnan(vpmat)==0)=min(auxmat(isnan(auxmat)==0));
            save_xzv(fullfile(dir_xzv_inv_mod,['Aux.',avertype,'.',modeltype,'.xzv']),XmidT,depth,auxmat,0,maskmat);
        end
    end
end

%% %% %%

% Check elapsed time
Tend=toc(Tstart);
[time_string]=secs2hms(Tend);
fprintf(['\n  Total elapsed time: ',time_string,'\n']);
