%%% SURFACE-WAVE dispersion INVERSION & PROFILING (SWIP)
%%% MODULE D1 : SWIPmod1d.m
%%% S. Pasquet - V20.04.03
%%% SWIPmod1d.m plots observed and calculated dispersion for each Xmid
%%% It also plots 1D Vp, Vs, Vp/Vs and Poisson's ratio models

%%% This module calls the following Geopsy native function:
%%% gpdc
%%% Geopsy function are called through the following MATLAB function:
%%% matgpdc
%%% The following Linux codes are also called if correctly installed:
%%% ImageMagick (convert, montage) - pdfjam - pdfcrop

%%%-------------------------%%%
%%% START OF INITIALIZATION %%%

run('SWIP_defaultsettings') % Read default settings

% Check user input
if swip==0 && tomo==0 && user==0
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n   Select at least one model for forward calculation');
    fprintf('\n       Set either "swip", "tomo" or "user" to 1');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
    return
end
if plot1dcal==0 && plot1dmod==0
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n       Select at least one plot option');
    fprintf('\n   Set either "plot1dcal" or "plot1dmod" to 1');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
    return
end

% Initialization (same for D1 and D2)
if swip==1
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
dx = acquiparam.dx;
topo = acquiparam.topo;
nWmin = stackdisp.nWmin;
nWmax = stackdisp.nWmax;
dW = stackdisp.dW;
if isfield(stackdisp,'nWvec')
    nWvec = stackdisp.nWvec;
else
    nWvec = nWmin:dW:nWmax;
end
xsca=pomega.xsca;
XmidT=xmidparam.XmidT; % Get Xmids
Xlength=length(XmidT); % Number of Xmids
xmidformat=stackdisp.xmidformat;
winsize=nWvec; % Window sizes (no. of geophones)
maxwinsize=(winsize-1)*dx;

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
if exist('targopt_inv','var')==0
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

% Default maximum resampling wavelength (power law depending on window size)
if length(maxwinsize)>1
    max_resamp_win = maxwinsize.*10.^(1./sqrt(0.5*maxwinsize));
else
    max_resamp_win = max(resampvec);
end

if ~exist('auto_resamp','var') || isempty(auto_resamp)
    auto_resamp = 1;
end

% Initialize depth vector and topography
zround=xmidparam.zround; % Get topography
if isempty(dpMAX)==1 % Get maximum depth from parameterization if not setup in launcher
    if swip==1
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

% Specific settings for D1
DOI=NaN*zeros(Xlength,1);
flip=pomega.flip;
flim=xmidparam.flim;

% File and folder names initialization
if nbest==0
    extens=['.bweb',num2str(outpoints)]; % Best within error bars
else
    extens=['.best',num2str(nbest)]; % Arbitrary nb
end
if swip==1
    dir_img_inv_mod = fullfile(dir_img_inv,['models',extens]);
    dir_img_inv_1d = fullfile(dir_img_inv_mod,'1dmodels');
    dir_img_inv_disp = fullfile(dir_img_inv_mod,'1ddisp');
else
    dir_img_inv_mod = fullfile(dir_img,'Usermodels');
    dir_img_inv_1d = fullfile(dir_img_inv_mod,'1dmodels');
    dir_img_inv_disp = fullfile(dir_img_inv_mod,'1ddisp');
end
if exist(dir_img_inv_1d,'dir') ~= 7 && plot1dmod == 1
    mkdir(dir_img_inv_1d);
end
if exist(dir_img_inv_disp,'dir') ~= 7 && plot1dcal == 1
    mkdir(dir_img_inv_disp);
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

fprintf('\n  ---------------------\n');

% Select refraction velocity models from file
if usevptomo==1 || tomo==1
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
        VpI=[]; VpItomo=[]; VsItomo=[];
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        fprintf('\n   No Vp model file selected - Ignore Vptomo');
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
    else
        try
            X_plot = repmat(XmidT,size(depth,1),1);
            VpI = readtomo(Vpfile,0,X_plot,depth,xsca,vpaver,[nWmin,nWmax],dx); % Read Vp tomo file
            VpItomo = VpI; % Read Vp tomo file
        catch
            fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
            fprintf('\n   Invalid Vp model file - Ignore Vptomo');
            fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
            VpI=[]; VpItomo=[]; VsItomo=[]; usevptomo=0;
        end
    end
    
    if tomo==1 && isempty(VpItomo)==0
        % Select VS
        fprintf('\n  Select Vs model file (cancel to skip Vs)\n');
        [filevel,pathvel]=uigetfile({'*.model;*.dat;*.xzv;*.txt'},'Select Vs model (cancel if no Vs model available)');
        if pathvel==0
            pois_test = 0.4;
            VsItomo=sqrt(VpItomo.^2/((1/(1-2*pois_test))+1));
            fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
            fprintf('\n   No Vs model file selected - Use Poisson''s ratio of %1.2f',pois_test);
            fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
            
        else
            Vsfile=fullfile(pathvel,filevel); % File with velocity (3 columns X,Z,Vs)
            try
                VsItomo=readtomo(Vsfile,0,X_plot,depth,xsca,vpaver,[nWmin,nWmax],dx); % Read Vp tomo file
                ind_nan = find(isnan(VsItomo) | isnan(VpItomo));
                VpItomo(ind_nan)=NaN;
                VsItomo(ind_nan)=NaN;
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
end

% Read user input gpdc format model
if user==2
    [fileman,pathman]=uigetfile('*','Velocity model with gpdc format');
    if fileman==0
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        fprintf('\n   No velocity model file selected - Ignore user model');
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
    end
    fileuser=fullfile(pathman,fileman);
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
end

% Check if image concatenation functions are installed
[testimgmgck,~]=unix('which montage');
[testpdfjam,~]=unix('which pdfjam');
testplot=((testpdfjam==0 && strcmp(imgform,'pdf')==1) || (testimgmgck==0 && strcmp(imgform,'pdf')==0 && strcmp(imgform,'fig')==0));
if concat == 0
    testplot = 0;
end
showplotinit=showplot;

fprintf('\n  **********************************************************');
fprintf('\n  **********************************************************\n');

%%% END OF INITIALIZATION %%%
%%%-----------------------%%%

%% CALCULATIONS FOR ALL XMIDS

%%%%%% Loop over all Xmids %%%%%%

Tstart=tic; % Start clock
for ix=Xmidselec
    
    % Initialization
    if Xlength>1
        close all
        showplot=showplotinit;
    end
    if sum(nshot(ix,:))>=0
        fprintf(['\n  Xmid',num2str(ix),' = ',num2str(XmidT(ix),xmidformat),' m\n']);
    end
    if swip==1
        dir_rep_ind = [dir_rep_inv,'/',num2str(XmidT(ix),xmidformat),'_reports'];
    end
    
    %% %% %%
    
    %%%%% Plot dispersion images and curves %%%%%
    
    if plot1dcal==1
        if plotflim==1
            flimsing=flim(ix);
        else
            flimsing=[];
        end
        
        %%%% Plot dispersion images %%%%
        
        % Check existence of dispersion image file
        dspfile_sum=fullfile(dir_dat,[num2str(XmidT(ix),xmidformat),'.sum.dsp']);
        dspfile_sum_new=fullfile(dir_dat,[num2str(XmidT(ix),xmidformat),'.weight.dsp']);
        if exist(dspfile_sum_new,'file')==2
            dspfile_sum = dspfile_sum_new;
        end
        
        if exist(dspfile_sum,'file')==2
            % Read dispersion image file
            [dspmat,f,v]=dsp2dat(dspfile_sum,flip,0);
            % Plot dispersion image
            if Dlogscale==0
                f1=plot_img(showplot,f,v,dspmat',flipud(map0),axetop,axerev,cb_disp,fs,...
                    freqtitle_long,'Phase velocity (m/s)',...
                    'Norm. ampli.',[fMIN fMAX],[VphMIN VphMAX],...
                    [],fticks,Vphticks,[],[],flimsing,[],[0 0 24 18],[]);
            else
                dspmatinv=1./(1-dspmat);
                dspmatinv(isinf(dspmatinv))=max(max(dspmatinv(isinf(dspmatinv)==0)));
                f1=plot_img_log(showplot,f,v,dspmatinv',flipud(map0),axetop,axerev,cb_disp,fs,...
                    freqtitle_long,'Phase velocity (m/s)',...
                    '1/(1-Norm. ampli.)',[fMIN fMAX],[VphMIN VphMAX],...
                    [1 length(map0)],fticks,Vphticks,[],[],flimsing,[],[0 0 24 18],[]);
            end
        else
            % Plot empty figure with correct axis
            f1=plot_curv(showplot,NaN,NaN,[],'.',[0 0 0],[],axetop,axerev,...
                0,fs,freqtitle_long,'Phase velocity (m/s)',[],...
                [fMIN fMAX],[VphMIN VphMAX],[],fticks,Vphticks,[],...
                [],[],[0 0 24 18],[]);
        end
        hold on; cm_saturation(map0sat);
        hh=dashline(f,f,2,2,2,2,'color',[0 0 1],'linewidth',5);
        set(hh,'visible','off');
        %         dashline(v./dx,1./v,2,2,2,2,'color',[0 0 1],'linewidth',5);
        if plotlamlim==1 && sampling==1
            if auto_resamp == 1
                max_resamp_xmid = ceil(max(max_resamp_win(nshot(ix,:)>0)));
            else
                max_resamp_xmid = max(resampvec);
            end
            dashline(f,f*max_resamp_xmid,2,2,2,2,'color',[0 0 1],'linewidth',3);
        end
        if Flogscale==1
            set(gca,'xscale','log');
        end
        
        %%%% Plot picked dispersion curves %%%%
        
        if plot1dobs==1
            % Check existence of picked dispersion curves
            if swip~=1
                nametarg=fullfile(dir_targ,[num2str(XmidT(ix),xmidformat),'.target']);
            else
                nametarg=fullfile(dir_rep_ind,[num2str(XmidT(ix),xmidformat),'.target']);
            end
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
                if swip==1 && sum(nshot(ix,:))>=0
                    fprintf(['\n  No dispersion picked for Xmid',num2str(ix),' = ',...
                        num2str(XmidT(ix),xmidformat),' m\n']);
                end
                npvc=0;
            end
            
            % Plot dispersion curves for each mode
            for ip=1:npvc
                if mod(modes(ip),2)==0
                    col=pickcol1;
                else
                    col=pickcol2;
                end
                
                plot(freqresamp{modes(ip)+1},vresamp{modes(ip)+1},'.','linewidth',1.5,'markersize',7,'color',col);
                if eb==1
                    if str2double(matrelease(1:4))>2014 && Flogscale~=1
                        han=terrorbar(freqresamp{modes(ip)+1},vresamp{modes(ip)+1},deltaresamp{modes(ip)+1},1,'units');
                        set(han,'LineWidth',1.5,'Color',col)
                    else
                        han=errorbar(freqresamp{modes(ip)+1},vresamp{modes(ip)+1},deltaresamp{modes(ip)+1},...
                            '.','Color',col,'linewidth',2.5,'markersize',7);
                        xlimits=xlim;
                        tick_length=diff(xlimits)/100;
                        errorbar_tick(han,tick_length,'units');
                    end
                end
            end
        end
    end
    
    %% %% %%
    
    %%%%% Get velocity models and calculate theoretical dispersion %%%%%
    
    %%%% From SWIP inversion results %%%%
    
    if swip==1
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
        
        D=[];
        if exist(filevel,'file')==0
            if exist(dir_rep_ind,'dir')==7 && sum(nshot(ix,:))>=0
                fprintf(['\n  No SWIP model for Xmid',num2str(ix),' = ',...
                    num2str(XmidT(ix),xmidformat),' m\n']);
            end
            vpsw=[]; vssw=[]; rhosw=[];
            vpstd=[]; vsstd_ok=[]; rhostd=[];
        else
            
            %%% Create velocity file in gpdc format %%%
            
            % Read velocity file
            modvel=dlmread(filevel,'',1,0);
            moddepth=[0;cumsum(modvel(:,1))];
            if exist(filemin,'file')==2 && exist(filemax,'file')==2
                modmin = dlmread(filemin,'',1,0);
                modmax = dlmread(filemax,'',1,0);
                modstd = modmax;
                %%% A checker !
                %                 modstd(:,3) = max([abs(modmax(:,3)-modvel(:,3)),abs(modmin(:,3)-modvel(:,3))],[],2);
                modstd(:,2) = (modmax(:,2)-modmin(:,2));  
                modstd(:,3) = (modmax(:,3)-modmin(:,3));
            elseif exist(filestd,'file')==2
                modstd=dlmread(filestd,'',1,0);
            else
                modstd = modvel;
                modstd(:,2:4) = modstd(:,2:4).*0.2;
            end
            depthstd=[0;cumsum(modstd(:,1))];
            
            if maxdepth>moddepth(end)
                moddepth(end)=maxdepth;
                depthstd(end)=maxdepth;
            else
                modvel=modvel(moddepth<maxdepth,:);
                modstd=modstd(depthstd<maxdepth,:);
                modmin=modmin(depthstd<maxdepth,:);
                modmax=modmax(depthstd<maxdepth,:);
                moddepth=[0;cumsum(modvel(:,1))];
                depthstd=[0;cumsum(modstd(:,1))];
            end
            thick=modvel(:,1);
            vpsw=modvel(:,2);
            vssw=modvel(:,3);
            rhosw=modvel(:,4);
            
            vpstd=modstd(:,2);
            vsstd=modstd(:,3);
            rhostd=modstd(:,4);
            vsstd_perc=100*(vsstd./vssw);
            
            if exist(filemin,'file')==2 && exist(filemax,'file')==2
                vsstd_min = modmin(:,3);
                vsstd_max = modmax(:,3);
            end
            vsstd_ok = vsstd;
            
            %%% Replace VP from SWIP with VP from tomo file if required %%%
            
            if usevptomo==0
                filedisp=[filevel,'.disp'];
            else
                filevel=[filevel,'_vptomo'];
                filedisp=[filevel,'.disp'];
                vptomo=VpI(VpI(:,ix)>0,ix);
                if isempty(vptomo)~=1
                    vptomo=[vptomo(1);vptomo;vptomo(end)];
                    ztmp=dz.*ones(1,size(VpI(VpI(:,ix)>0,ix),1));
                    zinc=[0,cumsum(ztmp)];
                    if max(zinc)>maxdepth && abs(max(zinc)-maxdepth)>1e-10
                        zinc=0:dz:maxdepth;
                        vptomo=vptomo(1:length(zinc)+1);
                    elseif max(zinc)<maxdepth && abs(max(zinc)-maxdepth)>1e-10
                        zinc=0:dz:maxdepth;
                        vptomo2=vptomo(end)*ones(length(zinc)+1,1);
                        vptomo2(1:length(vptomo))=vptomo;
                        vptomo=vptomo2;
                    else
                        zinc=0:dz:maxdepth;
                    end
                    nlay=size(moddepth,1)-1;
                    [vpsw,~,~,vssw]=velresamp(zinc,vptomo,moddepth,vssw,0.1,0,0); % Resample velocity model
                    dinsave(filevel,thick,vpsw,vssw,rhosw); % Save in gpdc format to run forward model
                else
                    if sum(nshot(ix,:))>=0
                        fprintf(['\n  No Vp from tomo for Xmid',num2str(ix),' = ',...
                            num2str(XmidT(ix),xmidformat),' m\n']);
                    end
                    vpsw=[]; vssw=[]; rhosw=[];
                    vpstd=[]; vsstd_ok=[]; rhostd=[];
                end
            end
            
            %%% Compute and plot theoretical dispersion %%%
            
            if plot1dcal==1 && isempty(vssw)==0
                nftest=nf;
                while nftest>10
                    % Run forward modeling
                    matgpdc(filevel,nmodemax,wave,nftest-1,fmin+(fmax-fmin)/nftest,fmax,sampling,filedisp);
                    D=readdisp(filedisp,nmodemax); % Get calculated dispersion curve
                    if isempty(D)==1 % Check if it worked, otherwise try with less frequency samples
                        if nftest==nf
                            fprintf('\n  Re-run forward calculation with less frequency samples\n');
                        end
                        nftest=nftest-10;
                    else
                        break
                    end
                end
                if usevptomo==1
                    delete(filevel); % Delete temp file
                end
                delete(filedisp);
                if isempty(D)==0
                    % Plot theoretical dispersion on dispersion image
                    for m=1:nmodemax
                        freqcal=D{m,1}; % Frequency
                        vcal=1./D{m,2}; % Frequency
                        if nmod(ix)<=25
                            col = [0 1 0];
                        else
                            col = [0 1 0];
                        end
                        plot(freqcal,vcal,'-','Color',col,'linewidth',2,'markersize',10);
                    end
                else
                    fprintf('\n  Failed to compute theoretical dispersion from SWIP model\n');
                end
            end
        end
    end
    
    %% %% %%
    
    %%%% From refraction tomography models (Vp and Vs) %%%%
    
    if tomo==1 && isempty(VpItomo)==0 && isempty(VsItomo)==0
        
        %%% Create velocity file in gpdc format %%%
        
        D=[];
        filevel=fullfile(dir_dat,[num2str(XmidT(ix),xmidformat),'.tomo']);
        filedisp=fullfile(dir_dat,[num2str(XmidT(ix),xmidformat),'.tomo.disp']);
        nlay=size(VpItomo(VpItomo(:,ix)>0,ix),1);
        ztomo=dz.*ones(nlay,1);
        
        vstomo=VsItomo(VsItomo(:,ix)>0,ix);
        vptomo=VpItomo(VpItomo(:,ix)>0,ix);
        
        if isempty(vptomo)~=1 && isempty(vstomo)~=1
            for ll=1:nlay
                while poisson(vptomo(ll),vstomo(ll))<=0.1
                    vptomo(ll)=vptomo(ll)+1;
                    vstomo(ll)=vstomo(ll)-1;
                end
            end
            if isempty(rhoMIN)==1
                flagrho=1;
                rhoMIN=1800; rhoMAX=1800;
            end
            dinsave(filevel,ztomo,vptomo,vstomo,mean([rhoMIN,rhoMAX]));
            if flagrho==1
                rhoMIN=[]; rhoMAX=[];
            end
        else
            if sum(nshot(ix,:))>=0
                fprintf('\n  No Vp or Vs from tomo for this Xmid\n');
            end
            vptomo=[]; vstomo=[];
            vpstd=[]; vsstd_ok=[];
        end
        
        %%% Compute and plot theoretical dispersion %%%
        
        if plot1dcal==1 && exist(filevel,'file')==2
            nftest=nf;
            while nftest>10
                matgpdc(filevel,nmodemax,wave,nftest-1,fmin+(fmax-fmin)/nftest,fmax,sampling,filedisp);
                D=readdisp(filedisp,nmodemax);
                if isempty(D)==1
                    if nftest==nf
                        fprintf('\n  Re-run forward calculation with less frequency samples\n');
                    end
                    nftest=nftest-10;
                else
                    break
                end
            end
            delete(filevel,filedisp);
            if isempty(D)==0
                for m=1:nmodemax
                    freqcal=D{m,1}; % Frequency
                    vcal=1./D{m,2}; % Frequency
                    plot(freqcal,vcal,'-','Color',[1 0 0],...
                        'linewidth',3,'markersize',10);
                end
            else
                fprintf('\n  Failed to compute theoretical dispersion from tomo model\n');
            end
        end
    end
    
    %% %% %%
    
    %%%% From user defined models %%%%
    
    if user>0
        
        %%% Create velocity file in gpdc format %%%
        
        D=[];
        if user==1 % Defined in launcher
            filevel=fullfile(dir_dat,[num2str(XmidT(ix),xmidformat),'.user']);
            filedisp=fullfile(dir_dat,[num2str(XmidT(ix),xmidformat),'.user.disp']);
            dinsave(filevel,[thkuser,thkuser(end)],vpuser,vsuser,rhouser);
        else % Defined in file
            filevel=fileuser;
            filedisp=[fileuser,'.disp'];
            try
                modveluser=dlmread(filevel,'',1,0); % Read model specified manually
            catch
                fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
                fprintf('\n   Invalid velocity model file - Ignore user model');
                fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
                user=0; filevel=[];
            end
        end
        
        %%% Compute and plot theoretical dispersion %%%
        
        if plot1dcal==1 && exist(filevel,'file')==2
            nftest=nf;
            while nftest>10
                matgpdc(filevel,nmodemax,wave,nftest-1,fmin+(fmax-fmin)/nftest,fmax,sampling,filedisp);
                D=readdisp(filedisp,nmodemax);
                if isempty(D)==1
                    if nftest==nf
                        fprintf('\n  Re-run forward calculation with less frequency samples\n');
                    end
                    nftest=nftest-10;
                else
                    break
                end
            end
            delete(filedisp);
            if user==1
                delete(filevel);
            end
            if isempty(D)==0
                for m=1:nmodemax
                    freqcal=D{m,1}; % Frequency
                    vcal=1./D{m,2}; % Frequency
                    plot(freqcal,vcal,'-','Color',[0 1 0],...
                        'linewidth',2,'markersize',10);
                end
            else
                fprintf('\n  Failed to compute theoretical dispersion from user model\n');
            end
        end
    end
    
    %% %% %%
    
    %%%%% Get Depth Of Investigation (DOI) %%%%%    
    
    if plotDOI==1 % Empirical DOI (Lmax*doifact)
        DOI(ix)=lmaxpick(ix)*doifact;
        
    elseif plotDOI==2 % DOI from VS standard deviation threshold
        if exist('vsstd_ok','var')==1 && isempty(vsstd_ok)==0
            flipmoddepth = flipud(moddepth);
            flipvsstd = flipud([vsstd_ok;vsstd_ok(end)]);
            indhsd = find(flipvsstd<std_mask,1,'first');
            if isempty(indhsd) || indhsd <= 2
                indhsd = find(flipvsstd<flipvsstd(1),1,'first');
%                 indmax = find(flipvsstd==max(flipvsstd),1);
%                 indhsd_all = find(flipvsstd<std_mask);
%                 indhsd = find(indhsd_all>indmax,1,'first');
                if isempty(indhsd) || indhsd <= 2
                    indhsd = find(flipvsstd~=flipvsstd(1),1,'first');
                end
            end
            hsdtmp=flipmoddepth(indhsd-2);
            DOI(ix)=hsdtmp;
        else
            DOI(ix)=maxdepth;
        end
        
    elseif plotDOI==3 % DOI from VS standard deviation median (experimental)
        if exist('vsstd_ok','var')==1 && isempty(vsstd_ok)==0
            flipmoddepth=flipud(moddepth);
            flipvssw=flipud([vssw;vssw(end)]);
            flipvsstd = flipud([vsstd_ok;vsstd_ok(end)]);
            ind_hsdtmp = find(flipvssw/flipvssw(1)>0.95,1,'last');
            hsdtmp = flipmoddepth(ind_hsdtmp);
            ind_std = find(flipvsstd < std_mask);
            
            if ~isempty(ind_std > ind_hsdtmp) && any(flipvsstd(flipmoddepth> lmaxpick(ix)/10) > std_mask)
                ind_max = find(flipvsstd == max(flipvsstd(flipmoddepth> lmaxpick(ix)/10)),1,'last');
                ind_hsdtmp = ind_std(find(ind_std > ind_max,1,'first'))-1;
                if flipmoddepth(ind_hsdtmp) > lmaxpick(ix)/10 || ind_hsdtmp < 0.66*length(flipmoddepth)
                    hsdtmp = flipmoddepth(ind_hsdtmp);
                end
            else
                ind_hsdtmp = find(flipvsstd == max(flipvsstd(flipmoddepth > lmaxpick(ix)/10)),1,'last');
                if (flipmoddepth(ind_hsdtmp) > lmaxpick(ix)/10 || ind_hsdtmp <= 0.66*length(flipmoddepth)) && flipmoddepth(ind_hsdtmp) < 0.2*lmaxpick(ix)
                    if flipvssw(ind_hsdtmp) < flipvssw(find(flipvssw/flipvssw(1)>0.95,1,'last'))
                        hsdtmp = flipmoddepth(ind_hsdtmp);
                    end
                else
                    hsdtmp = 0.2*lmaxpick(ix);
                end
            end
            
            if isempty(hsdtmp)==1
                hsdtmp=0;
            end
            DOI(ix)=hsdtmp;
        else
            DOI(ix)=maxdepth;
        end
        
    elseif plotDOI == 4 % DOI from VS standard deviation threshold (xp)
        if exist('vsstd_ok','var') == 1 && isempty(vsstd_ok) == 0
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
                hsdtmp=0;
            end
            DOI(ix)=hsdtmp;
        else
            DOI(ix)=maxdepth;
        end
    end
    
    %% %% %%
    
    %%%%% Save dispersion images with dispersion curves %%%%%
    
    if plot1dcal==1
        if (swip==1 && isempty(vssw)==0) || user>0 || (tomo==1 && isempty(VpItomo)==0 && isempty(VsItomo)==0)
            fprintf('\n  Saving calculated dispersion\n');
            file1=fullfile(dir_img_inv_disp,[num2str(XmidT(ix),xmidformat),'.xdisp.fwd.',imgform]);
            save_fig(f1,file1,imgform,imgres,1);
            if showplot==0
                close(f1);
            else
                showplot=showplot+1;
            end
        end
    end
    
    %% %% %%
    
    %%%%% Plot 1D models %%%%%
    
    if plot1dmod==1
        
        %% %% %%
        
        %%%% Plot 1D VS models %%%%
        
        if (swip==1 && isempty(vssw)==0) || user>0 || (tomo==1 && isempty(VpItomo)==0 && isempty(VsItomo)==0)
            fprintf('\n  Saving 1D models\n');
        end
        
        %%%% From SWIP inversion results %%%%
        
        if swip==1 && isempty(vssw)==0
            %%
            % Plot Vs
            [VSplot,Zplot]=stair2plot(vssw,moddepth);
            f1=plot_curv(showplot,VSplot,Zplot,[],'-',[1 0 0],[],1,1,...
                0,fs,'Vs (m/s)',depthtitle,[],[vsMIN vsMAX],[dpMIN dpMAX],[],...
                vsticks,dticks,[],[],[],[0 0 24 18],[],0);
            hold on
            
            % Plot VsSTD envelope on same plot
            if plot1dstd==1
                depthstdcalc=[0;diff(vssw)<0;0];
                depthstdcalc(depthstdcalc==0)=-1;
                if errstd>0
                    vsstd_plot=vssw*errstd/100;
                    depthstdup=sort(moddepth+depthstdcalc.*moddepth*errstd/100);
                    depthstddown=sort(moddepth-depthstdcalc.*moddepth*errstd/100);
                    [VSstdplotup,Zplotup]=stair2plot(vssw+vsstd_plot,depthstdup);
                    [VSstdplotdown,Zplotdown]=stair2plot(vssw-vsstd_plot,depthstddown);
                else
                    if exist(filemin,'file')==2 && exist(filemax,'file')==2
                        [VSstdplotup,Zplotup]=stair2plot(vsstd_max,depthstd);
                        [VSstdplotdown,Zplotdown]=stair2plot(vsstd_min,depthstd);
                    else
                        [VSstdplotup,Zplotup]=stair2plot(vsstd_ok,depthstd);
                        [VSstdplotdown,Zplotdown]=stair2plot(vsstd_ok,depthstd);
                        VSstdplotup = VSplot + VSstdplotup;
                        VSstdplotdown = VSplot - VSstdplotdown;
                    end
                end
                dashline(VSstdplotup,Zplotup,2,2,2,2,'color',[1,0,0],'linewidth',2);
                dashline(VSstdplotdown,Zplotdown,2,2,2,2,'color',[1,0,0],'linewidth',2);
                
                if isempty(vsMAX)
                    xlim([0 max(VSstdplotup)]);
                end
            end
            %%
            % Plot Vp on same plot
            if plot1dvp==1
                [VPplot,Zplot]=stair2plot(vpsw,moddepth);
                plot(VPplot,Zplot,'m-','linewidth',2);
                xlabel('V (m/s)');
                
                % Plot VpSTD envelope on same plot
                if plot1dstd==1
                    depthstdcalc=[0;diff(vpsw)<0;0];
                    depthstdcalc(depthstdcalc==0)=-1;
                    if errstd>0
                        vpstd=vpsw*errstd/100;
                        depthstdup=sort(moddepth+depthstdcalc.*moddepth*errstd/100);
                        depthstddown=sort(moddepth-depthstdcalc.*moddepth*errstd/100);
                        [VPstdplotup,Zplotup]=stair2plot(vpsw+vpstd,depthstdup);
                        [VPstdplotdown,Zplotdown]=stair2plot(vpsw-vpstd,depthstddown);
                    else
                        [VPstdplotup,Zplotup]=stair2plot(vpstd,depthstd);
                        [VPstdplotdown,Zplotdown]=stair2plot(vpstd,depthstd);
                        VPstdplotup = VPplot + VPstdplotup;
                        VPstdplotdown = VPplot - VPstdplotdown;
                    end
                    
                    dashline(VPstdplotup,Zplotup,2,2,2,2,'color','m','linewidth',2);
                    dashline(VPstdplotdown,Zplotdown,2,2,2,2,'color','m','linewidth',2);
                    
                    if isempty(vsMAX)
                        xlim([0 max(VPstdplotup)]);
                    end
                end
            end
        end
        
        %%%% From user defined models %%%%
        
        if user>0
            % Plot Vs
            if user==2
                thkuserok=modveluser(1:end-1,1);
                vsuserok=modveluser(:,3);
                vpuserok=modveluser(:,2);
            else
                thkuserok=thkuser';
                vsuserok=vsuser';
                vpuserok=vpuser';
            end
            if maxdepth<sum(thkuserok)
                moddepthuser=[0;cumsum(thkuserok);sum(thkuserok)];
            else
                moddepthuser=[0;cumsum(thkuserok);maxdepth];
            end
            [VSplotuser,Zplotuser]=stair2plot(vsuserok,moddepthuser);
            if swip==1 && isempty(vssw)==0
                plot(VSplotuser,Zplotuser,'b-','linewidth',2);
            else
                f1=plot_curv(showplot,VSplotuser,Zplotuser,[],'-',[0 0 1],[],1,1,...
                    0,fs,'Vs (m/s)',depthtitle,[],[vsMIN vsMAX],[dpMIN dpMAX],[],...
                    vsticks,dticks,[],[],[],[0 0 24 18],[],0);
                hold on
            end
            
            % Plot VsSTD envelope on same plot
            if plot1dstd==1 && errstd>0
                depthstdcalc=[0;diff(vsuserok)<0;0];
                depthstdcalc(depthstdcalc==0)=-1;
                
                vsstd_plot=vsuserok*errstd/100;
                depthstdup=sort(moddepthuser+depthstdcalc.*moddepthuser*errstd/100);
                depthstddown=sort(moddepthuser-depthstdcalc.*moddepthuser*errstd/100);
                [VSstdplotup,Zplotup]=stair2plot(vsuserok+vsstd_plot,depthstdup);
                [VSstdplotdown,Zplotdown]=stair2plot(vsuserok-vsstd_plot,depthstddown);
                
                dashline(VSstdplotup,Zplotup,2,2,2,2,'color',[0,0,1],'linewidth',2);
                dashline(VSstdplotdown,Zplotdown,2,2,2,2,'color',[0,0,1],'linewidth',2);  
            end
            
            % Plot Vp on same plot
            if plot1dvp==1
                [VPplotuser,Zplotuser]=stair2plot(vpuserok,moddepthuser);
                plot(VPplotuser,Zplotuser,'c-','linewidth',2);
                xlabel('V (m/s)');
                
                % Plot VpSTD envelope on same plot
                if plot1dstd==1 && errstd>0
                    depthstdcalc=[0;diff(vpuserok)<0;0];
                    depthstdcalc(depthstdcalc==0)=-1;
                    
                    vpstd=vpuserok*errstd/100;
                    depthstdup=sort(moddepthuser+depthstdcalc.*moddepthuser*errstd/100);
                    depthstddown=sort(moddepthuser-depthstdcalc.*moddepthuser*errstd/100);
                    [VPstdplotup,Zplotup]=stair2plot(vpuserok+vpstd,depthstdup);
                    [VPstdplotdown,Zplotdown]=stair2plot(vpuserok-vpstd,depthstddown);
                    
                    dashline(VPstdplotup,Zplotup,2,2,2,2,'color','c','linewidth',2);
                    dashline(VPstdplotdown,Zplotdown,2,2,2,2,'color','c','linewidth',2);
                    
                    if isempty(vsMAX)
                        xlim([0 max(VPstdplotup)]);
                    end
                end
            end
        end
        
        %%%% From refraction tomography models (Vp and Vs) %%%%
        
        if tomo==1 && isempty(VpItomo)==0 && isempty(VsItomo)==0
            % Plot Vs
            if (swip==1 && isempty(vssw)==0) || user==1
                plot(vstomo,cumsum(ztomo),'-','color',[0 0.75 0],'linewidth',2);
            else
                f1=plot_curv(showplot,vstomo,cumsum(ztomo),[],'-',[0 0.75 0],[],1,1,...
                    0,fs,'Vs (m/s)',depthtitle,[],[vsMIN vsMAX],[dpMIN dpMAX],[],...
                    vsticks,dticks,[],[],[],[0 0 24 18],[],0);
                hold on
            end
            
            % Plot VsSTD envelope on same plot
            if plot1dstd==1
                if errstd>0
                    dashline(vstomo+vstomo*errstd/100,cumsum(ztomo),2,2,2,2,'color',[0 0.75 0],'linewidth',2);
                    dashline(vstomo-vstomo*errstd/100,cumsum(ztomo),2,2,2,2,'color',[0 0.75 0],'linewidth',2);
                else
                    
                end
            end
            
            % Plot Vp on same plot
            if plot1dvp==1
                plot(vptomo,cumsum(ztomo),'g-','linewidth',2);
                xlabel('V (m/s)');
                
                % Plot VpSTD envelope on same plot
                if plot1dstd==1
                    if errstd>0
                        dashline(vptomo+vptomo*errstd/100,cumsum(ztomo),2,2,2,2,'color','g','linewidth',2);
                        dashline(vptomo-vptomo*errstd/100,cumsum(ztomo),2,2,2,2,'color','g','linewidth',2);
                    else
                        
                    end
                end
            end
        end
        
        %%%% Save figure %%%%
        
        if (swip==1 && isempty(vssw)==0) || user>0 || (tomo==1 && isempty(VpItomo)==0 && isempty(VsItomo)==0)
            % Plot DOI
            if plotDOI>0
                xL=get(gca,'XLim');
                han3=dashline(xL,[DOI(ix) DOI(ix)],2,2,2,2,'color','k','linewidth',2);
            end
            sizeax=get(gca,'Position');
            file_vs=fullfile(dir_img_inv_1d,[num2str(XmidT(ix),xmidformat),'.xVs1D.',imgform]);
            save_fig(f1,file_vs,imgform,imgres,1,0);
            if showplot==0
                close(f1);
            else
                showplot=showplot+1;
            end
        end
        
        %% %% %%
        
        %%%% Plot 1D VP models %%%%
        
        %%%% From SWIP inversion results %%%%
        
        if swip==1 && isempty(vssw)==0
            % Plot Vp
            [VPplot,Zplot]=stair2plot(vpsw,moddepth);
            f1=plot_curv(showplot,VPplot,Zplot,[],'-',[1 0 0],[],1,1,...
                0,fs,'Vp (m/s)',depthtitle,[],[vpMIN vpMAX],[dpMIN dpMAX],[],...
                vpticks,dticks,[],[],[],[0 0 24 18],sizeax,0);
            hold on

            % Plot VpSTD envelope on same plot
            if plot1dstd==1
                depthstdcalc=[0;diff(vpsw)<0;0];
                depthstdcalc(depthstdcalc==0)=-1;
                if errstd>0
                    vpstd=vpsw*errstd/100;
                    depthstdup=sort(moddepth+depthstdcalc.*moddepth*errstd/100);
                    depthstddown=sort(moddepth-depthstdcalc.*moddepth*errstd/100);
                    [VPstdplotup,Zplotup]=stair2plot(vpsw+vpstd,depthstdup);
                    [VPstdplotdown,Zplotdown]=stair2plot(vpsw-vpstd,depthstddown);
                else
                    [VPstdplotup,Zplotup]=stair2plot(vpstd,depthstd);
                    [VPstdplotdown,Zplotdown]=stair2plot(vpstd,depthstd);
                    VPstdplotup = VPplot + VPstdplotup;
                    VPstdplotdown = VPplot - VPstdplotdown;
                end
                
                dashline(VPstdplotup,Zplotup,2,2,2,2,'color','r','linewidth',2);
                dashline(VPstdplotdown,Zplotdown,2,2,2,2,'color','r','linewidth',2);
                
                if isempty(vpMAX)
                    xlim([0 max(VPstdplotup)]);
                end
            end  
        end
        
        %%%% From user defined models %%%%
        
        if user>0
            % Plot Vp
            [VPplotuser,Zplotuser]=stair2plot(vpuserok,moddepthuser);
            if swip==1 && isempty(vssw)==0
                plot(VPplotuser,Zplotuser,'b-','linewidth',2);
            else
                f1=plot_curv(showplot,VPplotuser,Zplotuser,[],'-',[0 0 1],[],1,1,...
                    0,fs,'Vp (m/s)',depthtitle,[],[vpMIN vpMAX],[dpMIN dpMAX],[],...
                    vpticks,dticks,[],[],[],[0 0 24 18],sizeax,0);
                hold on
            end
            
            % Plot VpSTD envelope on same plot
            if plot1dstd==1 && errstd>0
                depthstdcalc=[0;diff(vpuserok)<0;0];
                depthstdcalc(depthstdcalc==0)=-1;
                
                vpstd=vpuserok*errstd/100;
                depthstdup=sort(moddepthuser+depthstdcalc.*moddepthuser*errstd/100);
                depthstddown=sort(moddepthuser-depthstdcalc.*moddepthuser*errstd/100);
                [VPstdplotup,Zplotup]=stair2plot(vpuserok+vpstd,depthstdup);
                [VPstdplotdown,Zplotdown]=stair2plot(vpuserok-vpstd,depthstddown);
                
                dashline(VPstdplotup,Zplotup,2,2,2,2,'color','b','linewidth',2);
                dashline(VPstdplotdown,Zplotdown,2,2,2,2,'color','b','linewidth',2);
                
                if isempty(vpMAX)
                    xlim([0 max(VPstdplotup)]);
                end
            end
        end
        
        %%%% From refraction tomography models (Vp and Vs) %%%%
        
        if tomo==1 && isempty(VpItomo)==0 && isempty(VsItomo)==0
            % Plot Vp
            if (swip==1 && isempty(vssw)==0) || user>0
                plot(vptomo,cumsum(ztomo),'-','color',[0 0.75 0],'linewidth',2);
            else
                f1=plot_curv(showplot,vptomo,cumsum(ztomo),[],'-',[0 0.75 0],[],1,1,...
                    0,fs,'Vp (m/s)',depthtitle,[],[vpMIN vpMAX],[dpMIN dpMAX],[],...
                    vpticks,dticks,[],[],[],[0 0 24 18],sizeax,0);
                hold on
            end
            
            % Plot VpSTD envelope on same plot
            if plot1dstd==1
                if errstd>0
                    plot(vptomo+vptomo*errstd/100,cumsum(ztomo),'--','color',...
                        [0 0.75 0],'linewidth',2);
                    plot(vptomo-vptomo*errstd/100,cumsum(ztomo),'--','color',...
                        [0 0.75 0],'linewidth',2);
                else
                    
                end
            end
        end
        
        %%%% Save figure %%%%
        
        if (swip==1 && isempty(vssw)==0) || user>0 || (tomo==1 && isempty(VpItomo)==0 && isempty(VsItomo)==0)
            % Plot DOI
            if plotDOI>0
                xL=get(gca,'XLim');
                han3=dashline(xL,[DOI(ix) DOI(ix)],2,2,2,2,'color','k','linewidth',2);
            end
            file_vp=fullfile(dir_img_inv_1d,[num2str(XmidT(ix),xmidformat),'.xVp1D.',imgform]);
            save_fig(f1,file_vp,imgform,imgres,1,0);
            if showplot==0
                close(f1);
            else
                showplot=showplot+1;
            end
        end
        
        %% %% %%
        
        %%%% Plot 1D VP/VS models %%%%
        
        %%%% From SWIP inversion results %%%%
        
        if swip==1 && isempty(vssw)==0
            % Plot Vp/Vs
            f1=plot_curv(showplot,VPplot./VSplot,Zplot,[],'-',[1 0 0],[],1,1,...
                0,fs,'Vp/Vs',depthtitle,[],[vpvsMIN vpvsMAX],[dpMIN dpMAX],[],...
                vpvsticks,dticks,[],[],[],[0 0 24 18],sizeax,0);
            hold on
        end
        
        %%%% From user defined models %%%%
        
        if user>0
            % Plot Vp/Vs
            if swip==1 && isempty(vssw)==0
                plot(VPplotuser./VSplotuser,Zplotuser,'b-','linewidth',2);
            else
                f1=plot_curv(showplot,VPplotuser./VSplotuser,Zplotuser,[],'-',[0 0 1],[],1,1,...
                    0,fs,'Vp/Vs',depthtitle,[],[vpvsMIN vpvsMAX],[dpMIN dpMAX],[],...
                    vpvsticks,dticks,[],[],[],[0 0 24 18],sizeax,0);
                hold on
            end
        end
        
        %%%% From refraction tomography models (Vp and Vs) %%%%
        
        if tomo==1 && isempty(VpItomo)==0 && isempty(VsItomo)==0
            % Plot Vp/Vs
            if (swip==1 && isempty(vssw)==0) || user>0
                plot(vptomo./vstomo,cumsum(ztomo),'-','color',[0 0.75 0],'linewidth',2);
            else
                f1=plot_curv(showplot,vptomo./vstomo,cumsum(ztomo),[],'-',[0 0.75 0],[],1,1,...
                    0,fs,'Vp/Vs',depthtitle,[],[vpvsMIN vpvsMAX],[dpMIN dpMAX],[],...
                    vpvsticks,dticks,[],[],[],[0 0 24 18],sizeax,0);
                hold on
            end
        end
        
        %%%% Save figure %%%%
        
        if (swip==1 && isempty(vssw)==0) || user>0 || (tomo==1 && isempty(VpItomo)==0 && isempty(VsItomo)==0)
            % Plot DOI
            if plotDOI>0
                xL=get(gca,'XLim');
                han3=dashline(xL,[DOI(ix) DOI(ix)],2,2,2,2,'color','k','linewidth',2);
            end
            file_vpvs=fullfile(dir_img_inv_1d,[num2str(XmidT(ix),xmidformat),'.xVpVs1D.',imgform]);
            save_fig(f1,file_vpvs,imgform,imgres,1,0);
            if showplot==0
                close(f1);
            else
                showplot=showplot+1;
            end
        end
        
        %% %% %%
        
        %%%% Plot 1D Poisson's ratio models %%%%
        
        %%%% From SWIP inversion results %%%%
        
        if swip==1 && isempty(vssw)==0
            % Plot Poisson's ratio
            f1=plot_curv(showplot,poisson(VPplot,VSplot),Zplot,[],'-',[1 0 0],[],1,1,...
                0,fs,'Poisson',depthtitle,[],[poisMIN poisMAX],[dpMIN dpMAX],[],...
                poisticks,dticks,[],[],[],[0 0 24 18],sizeax,0);
            hold on
        end
        
        %%%% From user defined models %%%%
        
        if user>0
            % Plot Poisson's ratio
            if swip==1 && isempty(vssw)==0
                plot(poisson(VPplotuser,VSplotuser),Zplotuser,'b-','linewidth',2);
            else
                f1=plot_curv(showplot,poisson(VPplotuser,VSplotuser),Zplotuser,[],'-',[0 0 1],[],1,1,...
                    0,fs,'Poisson',depthtitle,[],[poisMIN poisMAX],[dpMIN dpMAX],[],...
                    poisticks,dticks,[],[],[],[0 0 24 18],sizeax,0);
                hold on
            end
        end
        
        %%%% From refraction tomography models (Vp and Vs) %%%%
        
        if tomo==1 && isempty(VpItomo)==0 && isempty(VsItomo)==0
            % Plot Poisson's ratio
            if (swip==1 && isempty(vssw)==0) || user>0
                plot(poisson(vptomo,vstomo),cumsum(ztomo),'-','color',[0 0.75 0],'linewidth',2);
            else
                f1=plot_curv(showplot,poisson(vptomo,vstomo),cumsum(ztomo),[],'-',[0 0.75 0],[],1,1,...
                    0,fs,'Poisson',depthtitle,[],[poisMIN poisMAX],[dpMIN dpMAX],[],...
                    poisticks,dticks,[],[],[],[0 0 24 18],sizeax,0);
                hold on
            end
        end
        
        %%%% Save figure %%%%
        
        if (swip==1 && isempty(vssw)==0) || user>0 || (tomo==1 && isempty(VpItomo)==0 && isempty(VsItomo)==0)
            % Plot DOI
            if plotDOI>0
                xL=get(gca,'XLim');
                han3=dashline(xL,[DOI(ix) DOI(ix)],2,2,2,2,'color','k','linewidth',2);
            end
            file_pois=fullfile(dir_img_inv_1d,[num2str(XmidT(ix),xmidformat),'.xPoisson1D.',imgform]);
            save_fig(f1,file_pois,imgform,imgres,1,0);
            if showplot==0
                close(f1);
            else
                showplot=showplot+1;
            end
            
            %%% Concatenate figures %%%
            
            if testplot==1
                filename_imgtmp=fullfile(dir_img_inv_1d,[num2str(XmidT(ix),xmidformat),...
                    '.mod1d.tmp.',imgform]);
                filename_imgtmp2=fullfile(dir_img_inv_1d,[num2str(XmidT(ix),xmidformat),...
                    '.mod1d.tmp2.',imgform]);
                filename_panel=fullfile(dir_img_inv_1d,[num2str(XmidT(ix),xmidformat),...
                    '.mod1d.final.',imgform]);
                cat_img([file_vs,' ',file_vpvs],imgform,1,'west',filename_imgtmp,0);
                cat_img([file_vp,' ',file_pois],imgform,1,'west',filename_imgtmp2,0);
                cat_img([filename_imgtmp,' ',filename_imgtmp2],imgform,2,'west',filename_panel,1);
                if concat==1
                    delete(file_vs,file_vp,file_vpvs,file_pois,filename_imgtmp,filename_imgtmp2);
                end
            end
        end
        
    end
    if sum(nshot(ix,:))>=0
        fprintf('\n  **********************************************************');
        fprintf('\n  **********************************************************\n');
    end
end

% Check elapsed time
Tend=toc(Tstart);
[time_string]=secs2hms(Tend);
fprintf(['\n  Total elapsed time: ',time_string,'\n']);