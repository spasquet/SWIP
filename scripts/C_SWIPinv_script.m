%%% SURFACE-WAVE dispersion INVERSION & PROFILING (SWIP)
%%% MODULE C : SWIPinv.m
%%% S. Pasquet - V20.04.03
%%% SWIPinv.m performs inversion of dispersion curves picked in module A
%%% and select best models for each Xmid to build a pseudo-2D Vs section
%%% It allows to plot all generated models and inversion parameters for 
%%% each Xmid with misfit based colorbars for selected and rejected models

%%% This module calls the following Geopsy native functions:
%%% dinver - gpdcreport - gpcurve
%%% Geopsy functions are called through the following MATLAB functions:
%%% matdinver - matgpdcreport
%%% The following Linux codes are also called if correctly installed:
%%% ImageMagick (convert, montage) - pdfjam - pdfcrop

%%%-------------------------%%%
%%% START OF INITIALIZATION %%%

run('SWIP_defaultsettings') % Read default settings

% Folders settings
if inversion==1
    dir_all=dir_create(0);
    if dir_all.dir_main==0
        return
    end
else
    [dir_all,dir_inv_img]=dir_create(2);
    if dir_all.dir_main==0
        return
    end
    if dir_inv_img.dir_rep_inv==0
        return
    end
end
dir_dat=dir_all.dir_dat;
dir_img=dir_all.dir_img;
dir_pick=dir_all.dir_pick;
dir_param=dir_all.dir_param;
dir_targ=dir_all.dir_targ;
dir_inv=dir_all.dir_inv;

% Read stack parameters from .mat file
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
XmidT=xmidparam.XmidT; % Get Xmids
Xlength=xmidparam.Xlength; % Get Number of Xmids
xmidformat=stackdisp.xmidformat;
% Select Xmids
if exist('Xmidselec','var')~=1 || isempty(Xmidselec)==1
    Xmidselec=1:Xlength;
end
if max(Xmidselec)>Xlength
    Xmidselec=Xmidselec(Xmidselec<=Xlength);
end
nshot=xmidparam.nshot;

if inversion==1 % Get parameterization file and create folders to store inversion results
    if paramtype==0
        % Select user-defined parameterization file
        fprintf('\n  Select parameterization file\n');
        [paramfile,parampath]=uigetfile(fullfile(dir_param,'*.param'),...
            'Select parameterization file');
        if paramfile==0
            fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
            fprintf('\n   Please select a parameterization file');
            fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
            return
        end
        paramname=fullfile(parampath,paramfile);
        % Inversion folder named after parameterization file
        dir_rep_inv=fullfile(dir_inv,paramfile(1:end-6));
        dir_img_inv=fullfile(dir_img,paramfile(1:end-6));
    else
        % Inversion folder named after parameterization type if paramtype>0
        dir_rep_inv=fullfile(dir_inv,['Type',num2str(paramtype),'_param']);
        dir_img_inv=fullfile(dir_img,['Type',num2str(paramtype),'_param']);
    end
    % Create folders if not existing
    if exist(dir_rep_inv,'dir')~=7
        mkdir(dir_rep_inv);
    end
    if exist(dir_img_inv,'dir')~=7
        mkdir(dir_img_inv);
    end
    % .mat file to store inversion settings
    matfileinv=fullfile(dir_rep_inv,[sufile,'.invparam.mat']);
    if exist(matfileinv,'file')==2 % Read .mat file if existing
        load(matfileinv);
        nmodeinv=inv_set.nmodeinv;
        maxmodeinv=inv_set.maxmodeinv;
        nC=inv_set.nC;
        nmod=inv_set.nmod;
        minmis=inv_set.minmis;
    else % Store inversion settings to be saved in .mat file later
        dir_inv_img.dir_rep_inv=dir_rep_inv;
        dir_inv_img.dir_img_inv=dir_img_inv;
        inv_set.paramtype=paramtype;
        inv_set.nrun=ones(Xlength,1)*nrun;
        inv_set.itmax=ones(Xlength,1)*itmax;
        inv_set.ns0=ones(Xlength,1)*ns0;
        inv_set.ns=ones(Xlength,1)*ns;
        inv_set.nr=ones(Xlength,1)*nr;
        inv_set.nmodeinv=zeros(Xlength,1)*NaN;
        inv_set.nC=inv_set.nmodeinv;
        inv_set.nmod=inv_set.nmodeinv;
        inv_set.minmis=inv_set.nmodeinv;
        inv_set.maxmodeinv=inv_set.nmodeinv;
        nmodeinv=inv_set.nmodeinv;
        nC=inv_set.nC;
        nmod=inv_set.nmod;
        minmis=inv_set.minmis;
    end
    targopt_inv.wave=targopt.wave(1);
    targopt_inv.resampvec=targopt.resampvec;
    targopt_inv.sampling=targopt.sampling;
    targopt_inv.lmaxpick=targopt.lmaxpick;
    
else % Read previous inversion settings
    dir_rep_inv=dir_inv_img.dir_rep_inv;
    dir_img_inv=dir_inv_img.dir_img_inv;
    % Read .mat file
    matfileinv=fullfile(dir_rep_inv,[sufile,'.invparam.mat']);
    try
        load(matfileinv);
    catch
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        fprintf('\n   Missing .mat file in inversion folder');
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
        return
    end
    % Read previous inversion settings
    paramtype=inv_set.paramtype;
    nrun=inv_set.nrun;
    itmax=inv_set.itmax;
    ns0=inv_set.ns0;
    ns=inv_set.ns;
    nr=inv_set.nr;
    nmodeinv=inv_set.nmodeinv;
    maxmodeinv=inv_set.maxmodeinv;
    nC=inv_set.nC;
    nmod=inv_set.nmod;
    minmis=inv_set.minmis;
    if exist('targopt_inv','var')==0
        targopt_inv.wave=targopt.wave(1);
        targopt_inv.resampvec=targopt.resampvec;
        targopt_inv.sampling=targopt.sampling;
        targopt_inv.lmaxpick=targopt.lmaxpick;
    end
end
nmaxmod=(inv_set.itmax.*inv_set.ns)+inv_set.ns0;

% Get previous settings
xmin=plotopt.xmin;
xmax=plotopt.xmax;
wave=targopt_inv.wave(1);
resampvec=targopt_inv.resampvec;
sampling=targopt_inv.sampling;
lmaxpick=targopt_inv.lmaxpick;

% Save inversion settings in .mat file
if inversion==1
    if exist(matfileinv,'file')==2
        save(matfileinv,'-append','targopt_inv');
    else
        save(matfileinv,'dir_inv_img','inv_set','targopt_inv');
    end
end

% Initialize depth vector and topography
zround=xmidparam.zround; % Get topography
if isempty(dpMAX)==1 % Get maximum depth from parameterization if not setup in launcher
    if inversion == 0
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
            if paramtype==0
                if isempty(strfind(dir_rep_inv,dir_inv))==1
                    dir_rep_inv=fullfile(dir_all.dir_start,dir_rep_inv);
                end
                paramname=fullfile(dir_param,[dir_rep_inv(length(dir_inv)+1:end),'.param']);
            end
        end
    end
    if inversion == 1 || (inversion == 0 && dpMAX == 0)
        if paramtype>0
            paramstruct=dir(fullfile(dir_targ,['*.type',num2str(paramtype),'.param']));
            dpMAX=zeros(size(paramstruct));
            for ip=1:length(paramstruct)
                paramfile=fullfile(dir_targ,paramstruct(ip).name);
                [~,dpMAX(ip)]=param2mod(paramfile);
            end
            dpMAX=max(dpMAX);
        else
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
    dpMIN=0;
end
% Create depth vector
maxdepth=ceil(dpMAX/dz)*dz;
depth=max(zround):-dz:min(zround)-maxdepth; % Depth vector with topo
ZZ=0:dz:maxdepth;
nZ=length(ZZ);

% Memory allocation
indf=nC;
indi=nC;
layered=zeros(length(depth),Xlength);

% File name extensions
if nbest==0
    extens=['.bweb',num2str(outpoints)]; % Best within error bars
else
    extens=['.best',num2str(nbest)]; % Arbitrary nb
end
dir_img_inv_mod = fullfile(dir_img_inv,['models',extens]);
dir_img_inv_1d = fullfile(dir_img_inv_mod,'1dinvres');
dir_img_inv_mean = fullfile(dir_img_inv_mod,'1dmean');
dir_img_inv_param = fullfile(dir_img_inv_mod,'1dparam');

if exist(dir_img_inv_1d,'dir') ~= 7
    mkdir(dir_img_inv_1d);
end
if exist(dir_img_inv_mean,'dir') ~= 7 && plot1dcal == 1
    mkdir(dir_img_inv_mean);
end
if exist(dir_img_inv_param,'dir') ~= 7 && plot1dcal == 1
    mkdir(dir_img_inv_param);
end

% Final average model type names
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

% Check if image concatenation functions are installed
[testimgmgck,~]=unix('which montage');
[testpdfjam,~]=unix('which pdfjam');
testplot=((testpdfjam==0 && strcmp(imgform,'pdf')==1) || (testimgmgck==0 && strcmp(imgform,'pdf')==0 && strcmp(imgform,'fig')==0));
if concat == 0
    testplot = 0;
end

% Check OS
if isunix==1
    zipmethod=5;
else
    zipmethod=1;
end

if calc == 0
    calcmod = 0;
end

fprintf('\n  **********************************************************');
fprintf('\n  **********************************************************\n');

%%% END OF INITIALIZATION %%%
%%%-----------------------%%%

%% CALCULATIONS FOR ALL XMIDS

%%%%%% Loop over all Xmids %%%%%%

Tstart=tic; % Start clock
for ix=Xmidselec
    tstart=tic; % Start single Xmid clock
    
    % Initialization
    if sum(nshot(ix,:))>=0 && (inversion==1 || calcmod==1 || plotinvres==1 || plotparam==1)
        fprintf(['\n  Xmid',num2str(ix),' = ',num2str(XmidT(ix),xmidformat),' m\n']);
    end
    % Create folder to store inversion results for each Xmid
    dir_rep_ind=fullfile(dir_rep_inv,[num2str(XmidT(ix),xmidformat),'_reports']);
    
    %% %% %%
    
    %%%%% Inversion %%%%%
    
    if inversion==1 % Run new inversion
        if sum(nshot(ix,:))>=0
            fprintf(['\n  Inversion with "type ',num2str(paramtype),'" parameterization\n']);
        end
        
        % Get target file
        nametarg=fullfile(dir_targ,[num2str(XmidT(ix),xmidformat),'.target']);
        if exist(nametarg,'file')~=2
            if sum(nshot(ix,:))>=0
                fprintf('\n  No target - Go to next Xmid\n\n');
                % Save settings in .mat file
                inv_set.nC(ix)=NaN;
                inv_set.nmod(ix)=NaN;
                inv_set.minmis(ix)=NaN;
                save(matfileinv,'-append','inv_set');
            end
            continue
        end
        
        % Get parameter file if autoparam
        if paramtype~=0
            paramname=fullfile(dir_targ,[num2str(XmidT(ix),xmidformat),'.type',...
                num2str(paramtype),'.param']);
            if exist(paramname,'file')~=2
                if sum(nshot(ix,:))>=0
                    fprintf('\n  No param - Go to next Xmid\n\n');
                end
                continue
            end
        end
        
        % Create report folder to store inversion results
        if exist(dir_rep_ind,'dir')~=7
            mkdir(dir_rep_ind);
        end
        
        % Read target file to get picked dispersion curves and test target
        % compatibility
        [~,~,~,modes,tst]=targ2pvc(nametarg);        
        if tst > 0
            fprintf('\n  Target file not compatible with dinver version!!\n\n');
            return
        end
        nmodeinv(ix)=length(modes);
        maxmodeinv(ix)=max(modes);
        
        % Run NA inversion with dinver
        status=matdinver(nametarg,paramname,nrun,itmax,ns0,ns,nr,dir_rep_ind,verbose);
        if status~=0 % Check if inversion ran correctly
            continue
        end
        % Copy target and param file in report folder
        copyfile(nametarg,dir_rep_ind);
        copyfile(paramname,dir_rep_ind);
        % Compress inversion results to save disk space
        matzip(1,fullfile(dir_rep_ind,'*.report'),zipmethod,1);
        
        % Save inversion settings        
        nmaxmod(ix)=(itmax*ns)+ns0;
        inv_set.nrun(ix)=nrun;
        inv_set.itmax(ix)=itmax;
        inv_set.ns0(ix)=ns0;
        inv_set.ns(ix)=ns;
        inv_set.nr(ix)=nr;
        inv_set.nmodeinv(ix)=nmodeinv(ix);
        inv_set.maxmodeinv(ix)=maxmodeinv(ix);
        if exist(matfileinv,'file')==2
            save(matfileinv,'-append','inv_set');
        end
    end
    
    % Get maximum wavelength from target file
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
    end
    
    %% %% %%
    
    %%%%% Extraction of inversion results / Building of average models %%%%%
    
    if calcmod==1 || plotinvres==1 || plotparam==1
        % Check if inversion exists for current Xmid
        if exist(dir_rep_ind,'dir')~=7
            if sum(nshot(ix,:))>=0
                fprintf('\n  No inversion - Go to next Xmid\n\n');
            end
            continue
        end
        if inversion==1 && verbose==1
            fprintf(['\n  Xmid',num2str(ix),' = ',num2str(XmidT(ix),xmidformat),' m\n']);
        end
        % Memory allocation
        OK=[]; KO=[]; misround=[]; ok=[];
        
        %% %% %%
        
        %%%% Looking for the best models %%%%
        
        % Get target file
        nametarg=fullfile(dir_rep_ind,[num2str(XmidT(ix),xmidformat),'.target']);
        if exist(nametarg,'file')~=2
            if sum(nshot(ix,:))>=0
                fprintf('\n  No target in inversion folder - Go to next Xmid\n\n');
            end
            continue
        end
        % Read target file to get picked dispersion curves
        [freqresamp,vresamp,deltaresamp,modes]=targ2pvc(nametarg);
        
        fprintf('\n  Looking for best models\n');
        
        % Export calculated dipsersion curves from binary report files
        if exist(fullfile(dir_rep_ind,'run_01.report'),'file')==2 % Compress if not previously zipped
            matzip(1,fullfile(dir_rep_ind,'*.report'),zipmethod,1);
        end
        matzip(0,fullfile(dir_rep_ind,'*.report.gz'),zipmethod,0); % Uncompress report
        
        export_gpdc=0;
        while export_gpdc<5
            try % Try reading reports (random crash probably due to memory lag)
                matgpdcreport(dir_rep_ind,0,inv_set.nrun(ix),maxmodeinv(ix)+1,nmaxmod(ix),wave);
                if nbest>nmaxmod(ix)
                    nbest=nmaxmod(ix);
                end
                D=cell(nmodeinv(ix),1);
                ok=cell(nmodeinv(ix),1);
                bestmis=zeros(inv_set.nrun(ix),nmodeinv(ix));
                misround=zeros(inv_set.nrun(ix)*nmaxmod(ix),nmodeinv(ix));
                
                % Loop over propagation modes
                for ip=1:nmodeinv(ix)
                    m=modes(ip);
                    fprintf(['\n      Mode ',num2str(m),'\n']);
                    % Get min and max velocity envelope from uncertainty range
                    minvelOK=vresamp{m+1}(isnan(vresamp{m+1})==0)-deltaresamp{m+1}(isnan(vresamp{m+1})==0);
                    maxvelOK=vresamp{m+1}(isnan(vresamp{m+1})==0)+deltaresamp{m+1}(isnan(vresamp{m+1})==0);
                    freqOK=freqresamp{m+1}(isnan(vresamp{m+1})==0);
                    D{ip}=cell(nmaxmod(ix)*inv_set.nrun(ix),2); % Frequency vs slowness cell array
                    dispselec=repmat(struct('modnum',[],'modok',[],...
                        'misround',[],'nfreqsample',[]),1,inv_set.nrun(ix)); % Models structure
                    % Loop over all runs
                    for n=1:inv_set.nrun(ix)
                        calcdisp=fullfile(dir_rep_ind,['best',num2str(n),'.M',num2str(m),'.txt']);
                        % Read dispersion curves and get the ones fitting the errorbar if nbest=0
                        [D{ip}(1+(n-1)*nmaxmod(ix):(n-1)*nmaxmod(ix)+nmaxmod(ix),:),dispselec(n)]=...
                            readdisp(calcdisp,nmaxmod(ix),n,nbest,outpoints,freqOK,minvelOK,maxvelOK);
                        bestmis(n,ip)=dispselec(n).misround(1); % Get best model misfit
                    end
                    names=fieldnames(dispselec); % Get field names
                    cellData=cellfun(@(f){vertcat(dispselec.(f))},names); % Collect field data into a cell array
                    dispselec=cell2struct(cellData,names); % Convert the cell array into a structure
                    ok{ip}=dispselec.modok; % Get selected models
                    misround(:,ip)=dispselec.misround; % Get model misfits
                end
                export_gpdc=10; % Success !
            catch
                export_gpdc=export_gpdc+1; % Try again
            end
            delete(fullfile(dir_rep_ind,'best*.txt')); % Delete temp file
        end
        if export_gpdc==5 % Only 5 tries, after that go to next Xmid
            fprintf('\n  Unable to read reports after 5 attempts - Go to next Xmid\n\n');
            delete(fullfile(dir_rep_ind,'*.report')); % Delete temp file
            continue
        end
        
        % Find best models indexes sorted by decreasing misfit
        bestrun=find(bestmis(:,1)==min(bestmis(:,1)),1,'first'); % Get best run nb to get best model
        mis=max(misround,[],2);
        [missort,I]=sortrows(mis,-1);
        if nbest==0 % Best models within errobars
            test=sum(cell2mat(ok'),2);
            OKsort=test==nmodeinv(ix);
            KOsort=test<nmodeinv(ix);
            OKsort=OKsort(I);
            KOsort=KOsort(I);
            if sum(OKsort) == 0
                OKsort=false(size(missort));
                OKsort(end-1:end)=true;
                KOsort=false(size(missort));
                KOsort(1:end-1)=true;
            end
        elseif nbest<0 % Test             
%             thres_misfit = 10^(log10(min(missort))+(nbest/100)*log10(min(missort)));
            thres_misfit = abs(nbest);
            OKsort=false(size(missort));
            OKsort(missort<=thres_misfit) = true;
            KOsort=false(size(missort));
            KOsort(missort>thres_misfit)=true;            
        else % Arbitrary number of models
            OKsort=false(size(missort));
            OKsort(end-nbest+1:end)=true;
            KOsort=false(size(missort));
            KOsort(1:end-nbest)=true;
        end
        % Sort by misfit for selected and rejected models
        misin=missort(OKsort);
        misout=missort(KOsort);
        % Get nb of selected models and misfit for the current Xmid
        nmod(ix)=size(misin,1);
        minmis(ix)=min(mis);
        
        % Display info
        fprintf(['\n  Minimum misfit = ',num2str(min(missort)),'\n']);
        if nbest==0
            fprintf(['\n  ',num2str(nmod(ix)),' models fit within the error bars\n']);
        elseif nbest<0
            fprintf(['\n  ',num2str(nmod(ix)),' models have misfit lower than %1.2f\n'],-nbest);
        else
            fprintf(['\n  The best ',num2str(nmod(ix)),' models have been selected\n']);
        end
        
        if nmod(ix)==0 % Go to next Xmid if no models selected
            export_gpdc=0;
            while export_gpdc<5
                try % Try reading reports (random crash probably due to memory lag)
                    fprintf('\n  Export best model only\n');
                    % Export calculated models from binary report files
                    matgpdcreport(dir_rep_ind,1,inv_set.nrun(ix),nmodeinv(ix),nmaxmod(ix),wave);
                    VSall=cell(nmaxmod(ix)*inv_set.nrun(ix),2); % Depth vs Vs cell array
                    VPall=cell(nmaxmod(ix)*inv_set.nrun(ix),2); % Depth vs Vp cell array
                    RHOall=cell(nmaxmod(ix)*inv_set.nrun(ix),2); % Depth vs Density cell array
                    
                    % Loop over all runs
                    for n=1:inv_set.nrun(ix)
                        fprintf(['\n      Run ',num2str(n),'\n']);
                        vsfile=fullfile(dir_rep_ind,['vs',num2str(n),'.txt']);
                        vpfile=fullfile(dir_rep_ind,['vp',num2str(n),'.txt']);
                        rhofile=fullfile(dir_rep_ind,['rho',num2str(n),'.txt']);
                        [VSall(1+(n-1)*nmaxmod(ix):(n-1)*nmaxmod(ix)+nmaxmod(ix),:)]=...
                            readmodel(vsfile,nmaxmod(ix),n); % Read file
                        [VPall(1+(n-1)*nmaxmod(ix):(n-1)*nmaxmod(ix)+nmaxmod(ix),:)]=...
                            readmodel(vpfile,nmaxmod(ix),n); % Read file
                        [RHOall(1+(n-1)*nmaxmod(ix):(n-1)*nmaxmod(ix)+nmaxmod(ix),:)]=...
                            readmodel(rhofile,nmaxmod(ix),n); % Read file
                    end
                    
                    % Read models
                    dim=ndims(VSall{1});
                    VS=cat(dim,VSall{:,1}); % All Vs
                    VSbest=VS(:,(nmaxmod(ix)*(bestrun-1))+1); % Best Vs
                    VS=VS(:,I);
                    VP=cat(dim,VPall{:,1}); % All Vp
                    VPbest=VP(:,(nmaxmod(ix)*(bestrun-1))+1); % Best Vp
                    VP=VP(:,I);
                    RHO=cat(dim,RHOall{:,1}); % All Rho
                    RHObest=RHO(:,(nmaxmod(ix)*(bestrun-1))+1); % Best Rho
                    RHO=RHO(:,I);
                    ZALL=cat(dim,VSall{:,2}); % All depth
                    ZALL(end,:)=maxdepth;
                    Zbest=ZALL(:,(nmaxmod(ix)*(bestrun-1))+1); % Best depth
                    THbest=diff(Zbest);
                    THbest=round(THbest(1:2:end)*1e6)/1e6; % Best thickness
                    ZALL=ZALL(:,I);
                    THALL=diff(ZALL);
                    THALL=round(THALL(1:2:end,:)*1e6)/1e6; % All thickness
                    
%                     % Read models
%                     dim=ndims(VSall{1});
%                     VS=cat(dim,VSall{:,1}); % All Vs
%                     VSbest=VS(:,(nmaxmod(ix)*(bestrun-1))+1); % Best Vs
%                     VP=cat(dim,VPall{:,1}); % All Vp
%                     VPbest=VP(:,(nmaxmod(ix)*(bestrun-1))+1); % Best Vp
%                     RHO=cat(dim,RHOall{:,1}); % All Rho
%                     RHObest=RHO(:,(nmaxmod(ix)*(bestrun-1))+1); % Best Rho
%                     ZALL=cat(dim,VSall{:,2}); % All depth
%                     ZALL(end,:)=maxdepth;
%                     Zbest=ZALL(:,(nmaxmod(ix)*(bestrun-1))+1); % Best depth
%                     THbest=diff(Zbest);
%                     THbest=round(THbest(1:2:end)*1e6)/1e6; % Best thickness
                    
                    % Best model
                    nameDINbest=fullfile(dir_rep_ind,[num2str(XmidT(ix),xmidformat),...
                        extens,'.Vms.best']); % Best model
                    % Save best model
                    dinsave(nameDINbest,THbest,VPbest(1:2:end),VSbest(1:2:end),RHObest(1:2:end));
                    
                    fprintf('\n  **********************************************************');
                    fprintf('\n  **********************************************************\n');
                                        
                    % Save settings in .mat file
                    inv_set.nmod(ix)=nmod(ix);
                    save(matfileinv,'-append','inv_set');
                    
                    export_gpdc=10; % Success !
                    delete(fullfile(dir_rep_ind,'*.report')); % Delete temp file
                    continue
                    
                catch
                    export_gpdc=export_gpdc+1; % Try again
                end
                delete(fullfile(dir_rep_ind,'*.txt')); % Delete temp file
            end
                    
            if export_gpdc==5 % Only 5 tries, after that go to next Xmid
                fprintf('\n  Unable to read reports after 5 attempts - Go to next Xmid\n\n');
                delete(fullfile(dir_rep_ind,'*.report')); % Delete temp file
                continue
            end
        end

        
        %% %% %%
        
        %%%% Read all models %%%%
        
        export_gpdc=0;
        fprintf('\n  Reading all calculated models\n');
        while export_gpdc<5
            try % Try reading reports (random crash probably due to memory lag)
                % Export calculated models from binary report files
                matgpdcreport(dir_rep_ind,1,inv_set.nrun(ix),nmodeinv(ix),nmaxmod(ix),wave);
                VSall=cell(nmaxmod(ix)*inv_set.nrun(ix),2); % Depth vs Vs cell array
                VPall=cell(nmaxmod(ix)*inv_set.nrun(ix),2); % Depth vs Vp cell array
                RHOall=cell(nmaxmod(ix)*inv_set.nrun(ix),2); % Depth vs Density cell array
                
                % Loop over all runs
                for n=1:inv_set.nrun(ix)
                    fprintf(['\n      Run ',num2str(n),'\n']);
                    vsfile=fullfile(dir_rep_ind,['vs',num2str(n),'.txt']);
                    vpfile=fullfile(dir_rep_ind,['vp',num2str(n),'.txt']);
                    rhofile=fullfile(dir_rep_ind,['rho',num2str(n),'.txt']);
                    [VSall(1+(n-1)*nmaxmod(ix):(n-1)*nmaxmod(ix)+nmaxmod(ix),:)]=...
                        readmodel(vsfile,nmaxmod(ix),n); % Read file
                    [VPall(1+(n-1)*nmaxmod(ix):(n-1)*nmaxmod(ix)+nmaxmod(ix),:)]=...
                        readmodel(vpfile,nmaxmod(ix),n); % Read file
                    [RHOall(1+(n-1)*nmaxmod(ix):(n-1)*nmaxmod(ix)+nmaxmod(ix),:)]=...
                        readmodel(rhofile,nmaxmod(ix),n); % Read file
                end
                
                % Read models
                dim=ndims(VSall{1});
                VS=cat(dim,VSall{:,1}); % All Vs
                VSbest=VS(:,(nmaxmod(ix)*(bestrun-1))+1); % Best Vs
                VS=VS(:,I);
                VP=cat(dim,VPall{:,1}); % All Vp
                VPbest=VP(:,(nmaxmod(ix)*(bestrun-1))+1); % Best Vp
                VP=VP(:,I);
                RHO=cat(dim,RHOall{:,1}); % All Rho
                RHObest=RHO(:,(nmaxmod(ix)*(bestrun-1))+1); % Best Rho
                RHO=RHO(:,I);
                ZALL=cat(dim,VSall{:,2}); % All depth
                ZALL(end,:)=maxdepth;
                Zbest=ZALL(:,(nmaxmod(ix)*(bestrun-1))+1); % Best depth
                THbest=diff(Zbest);
                THbest=round(THbest(1:2:end)*1e6)/1e6; % Best thickness
                ZALL=ZALL(:,I);
                THALL=diff(ZALL);
                THALL=round(THALL(1:2:end,:)*1e6)/1e6; % All thickness
                
                vsmax = ceil(max(VS(:)/10))*10;
                if isempty(vsMIN) || isempty(vsMAX)
                    vsMIN = 0;
                    vsMAX = vsmax;
                end
                
                nC(ix)=length(THALL(:,1)); % nb of layers
                
                % Get selected and rejected models
                VSin=VS(:,OKsort);
                VPin=VP(:,OKsort);
                VPout=VP(:,KOsort);
                RHOin=RHO(:,OKsort);
                RHOout=RHO(:,KOsort);
                Zin=ZALL(:,OKsort);
                THin=THALL(:,OKsort);
                VSout=VS(:,KOsort);
                Zout=ZALL(:,KOsort);
                
                % Save settings in .mat file
                inv_set.nC(ix)=nC(ix);
                inv_set.nmod(ix)=nmod(ix);
                inv_set.minmis(ix)=minmis(ix);
                save(matfileinv,'-append','inv_set');
                
                export_gpdc=10; % Success !
                delete(fullfile(dir_rep_ind,'*.report')); % Delete temp file
            catch
                export_gpdc=export_gpdc+1; % Try again
            end
            delete(fullfile(dir_rep_ind,'*.txt')); % Delete temp file
        end
        
        if export_gpdc==5 % Only 5 tries, after that go to next Xmid
            fprintf('\n  Unable to read reports after 5 attempts - Go to next Xmid\n\n');
            continue
        end
        
        %% %% %%
        
        %%%% Building average models %%%%
        
        if calcmod==1
            fprintf('\n  Building average models\n');
            
            %%% Mean layered model %%%
            VSmean=mean(VSin,dim);
            VPmean=mean(VPin,dim);
            RHOmean=mean(RHOin,dim);
            THmean=mean(THin,dim);
            Zmean=zeros(size(VSmean));
            Zmeantmp=cumsum(THmean);
            mm=0;
            for ll=2:2:length(VSmean)-1
                mm=mm+1;
                Zmean(ll:ll+1)=Zmeantmp(mm);
            end
            Zmean(end)=maxdepth;
            
            % Standard deviation layered model
            VSstd=std(VSin,0,dim);
            VPstd=std(VPin,0,dim);
            RHOstd=std(RHOin,0,dim);
            THstd=std(THin,0,dim);
            Zstd=zeros(size(VSstd));
            Zstdtmp=cumsum(THstd);
            mm=0;
            for ll=2:2:length(VSstd)-1
                mm=mm+1;
                Zstd(ll:ll+1)=Zstdtmp(mm);
            end
            Zstd(end)=maxdepth;
            
            % Max layered model
            VSmax=max(VSin,[],dim);
            VPmax=max(VPin,[],dim);
            RHOmax=max(RHOin,[],dim);
            THmax=max(THin,[],dim);
            Zmax=zeros(size(VSmax));
            Zmaxtmp=cumsum(THmax);
            mm=0;
            for ll=2:2:length(VSmax)-1
                mm=mm+1;
                Zmax(ll:ll+1)=Zmaxtmp(mm);
            end
            Zmax(end)=maxdepth;
            
            % Min layered model
            VSmin=min(VSin,[],dim);
            VPmin=min(VPin,[],dim);
            RHOmin=min(RHOin,[],dim);
            THmin=min(THin,[],dim);
            Zmin=zeros(size(VSmin));
            Zmintmp=cumsum(THmin);
            mm=0;
            for ll=2:2:length(VSmin)-1
                mm=mm+1;
                Zmin(ll:ll+1)=Zmintmp(mm);
            end
            Zmin(end)=maxdepth;
            
            %%% Mean smoothed model %%%
            
            % !!!!(to be optimized for constant density)!!!!
            Zi=Zin;
            Zi(2:2:end-2,:)=Zi(2:2:end-2,:)-dz;
            IVsR=zeros(nmod(ix),nZ);
            IVpR=IVsR; IrhoR=IVsR;
            Zi(end,:)=Zi(end,:)+dz*(Zi(end,:)==Zi(end-1,:));
            
            % Interpolate Vs, Vp and Rho
            IVs = arrayfun(@(ii)(interp1q(Zi(:,ii),VSin(:,ii),ZZ')), 1:size(VSin, 2), 'UniformOutput', false);
            IVs = cell2mat(IVs)';
            IVp = arrayfun(@(ii)(interp1q(Zi(:,ii),VPin(:,ii),ZZ')), 1:size(VPin, 2), 'UniformOutput', false);
            IVp = cell2mat(IVp)';
            Irho = arrayfun(@(ii)(interp1q(Zi(:,ii),RHOin(:,ii),ZZ')), 1:size(RHOin, 2), 'UniformOutput', false);
            Irho = cell2mat(Irho)';
            IVSmean=mean(IVs,1);
            IVPmean=mean(IVp,1);
            IRHOmean=mean(Irho,1);
            
            % Standard deviation smoothed model
            IVSstd=std(IVs,0,1);
            IVPstd=std(IVp,0,1);
            IRHOstd=std(Irho,0,1);
            
            % Max smoothed model
            IVSmax=max(IVs,[],1);
            IVPmax=max(IVp,[],1);
            IRHOmax=max(Irho,[],1);
            
            % Max smoothed model
            IVSmin=min(IVs,[],1);
            IVPmin=min(IVp,[],1);
            IRHOmin=min(Irho,[],1);
            
            %%% Ridge model %%%
            
            if ridgecalc==1
                % Velocity and density vectors for ridge search
                minVs=floor(min(min(IVs))/10)*10;
                maxVs=ceil(max(max(IVs))/10)*10;
                velocityS=minVs:20:maxVs;
                minVp=floor(min(min(IVp))/10)*10;
                maxVp=ceil(max(max(IVp))/10)*10;
                velocityP=minVp:50:maxVp;
                minrho=floor(min(min(Irho))/10)*10;
                maxrho=ceil(max(max(Irho))/10)*10;
                density=minrho:50:maxrho;
                for ii=1:length(velocityS)
                    flagg1=abs(velocityS(ii)-IVs);
                    if ii==1
                        flagg2=abs(velocityS(ii)-IVs);
                        IVsR(flagg1==flagg2)=velocityS(ii);
                    else
                        flagg2=abs(velocityS(ii-1)-IVs);
                        IVsR(flagg1<flagg2)=velocityS(ii);
                    end
                end
                for ii=1:length(velocityP)
                    flagg1=abs(velocityP(ii)-IVp);
                    if ii==1
                        flagg2=abs(velocityP(ii)-IVp);
                        IVpR(flagg1==flagg2)=velocityP(ii);
                    else
                        flagg2=abs(velocityP(ii-1)-IVp);
                        IVpR(flagg1<flagg2)=velocityP(ii);
                    end
                end
                for ii=1:length(density)
                    flagg1=abs(density(ii)-Irho);
                    if ii==1
                        flagg2=abs(density(ii)-Irho);
                        IrhoR(flagg1==flagg2)=density(ii);
                    else
                        flagg2=abs(density(ii-1)-Irho);
                        IrhoR(flagg1<flagg2)=density(ii);
                    end
                end
                NN=zeros(nZ,length(velocityS));
                NNP=zeros(nZ,length(velocityP));
                NNr=zeros(nZ,length(density));
                VSridge=zeros(1,nZ); VPridge=VSridge; RHOridge=VSridge;
                VSridge_min=VSridge; VSridge_max=VSridge;
                VPridge_min=VSridge; VPridge_max=VSridge;
                RHOridge_min=VSridge; RHOridge_max=VSridge;
                % Count number of model per cell
                a_std = 0.01; % Keep cells with at least(1-a_std)*std models
                for jj=1:nZ
                    NN(jj,:)=histc(IVsR(:,jj),velocityS);
                    VSridge(jj)=mean(velocityS(NN(jj,:)==max(NN(jj,:))));                   
                    VSridge_min(jj) = velocityS(find(NN(jj,:)>=a_std*std(NN(jj,:)),1,'first'));
                    VSridge_max(jj) = velocityS(find(NN(jj,:)>=a_std*std(NN(jj,:)),1,'last'));
                                        
                    NNP(jj,:)=histc(IVpR(:,jj),velocityP);
                    VPridge(jj)=mean(velocityP(NNP(jj,:)==max(NNP(jj,:))));
                    VPridge_min(jj) = velocityP(find(NNP(jj,:)>=a_std*std(NNP(jj,:)),1,'first'));
                    VPridge_max(jj) = velocityP(find(NNP(jj,:)>=a_std*std(NNP(jj,:)),1,'last'));
                    
                    NNr(jj,:)=histc(IrhoR(:,jj),density);
                    RHOridge(jj)=mean(density(NNr(jj,:)==max(NNr(jj,:))));
                    RHOridge_min(jj) = density(find(NNr(jj,:)>=a_std*std(NNr(jj,:)),1,'first'));
                    RHOridge_max(jj) = density(find(NNr(jj,:)>=a_std*std(NNr(jj,:)),1,'last'));
                end
                VSridge = median_filt(VSridge,9,1,length(VSridge));
                VSridge = mov_aver(VSridge',5,1,length(VSridge));
                VSridge_min = median_filt(VSridge_min,9,1,length(VSridge_min));
                VSridge_min = mov_aver(VSridge_min',5,1,length(VSridge_min));
                VSridge_max = median_filt(VSridge_max,9,1,length(VSridge_max));
                VSridge_max = mov_aver(VSridge_max',5,1,length(VSridge_max));
                NN(NN==0)=NaN;
                
                %% Plot ridgesearch results
                if nmod(ix)>1 %&& plotinvres==1
                    f2=plot_img_log(showplot,velocityS,ZZ,NN,flipud(autumn),...
                        1,1,1,fs,'Vs (m/s)',depthtitle,'Number of models',...
                        [0 vsmax],[dpMIN dpMAX],[1 max(NN(:))],vsticks,dticks,[],[],[],[],[0 0 24 18],[],0);
                    cm_saturation(0.5);
                    hold on
                    plot(VSridge,ZZ,'color','k','linewidth',1.5);
                    dashline(VSridge_min,ZZ,2,2,2,2,'color','k','linewidth',1.5);
                    dashline(VSridge_max,ZZ,2,2,2,2,'color','k','linewidth',1.5);
                    h=findall(gcf,'Type','Axes'); set(h,'FontSize',fs);

                    % Save figure
                    file_ridge=fullfile(dir_img_inv_mean,[num2str(XmidT(ix),xmidformat),...
                        '.mod1d.ridge.',imgform]);
                    save_fig(f2,file_ridge,imgform,imgres,1,1-testplot);
                    close(f2); figHandles = findall(0,'Type','figure');

                    if str2double(matrelease(1:4))>2014
                        figHandles=get(figHandles,'Number');
                    end
                    % Strange bug for first figure
                    if ismember(1,figHandles)==1
                        close(1);
                    end
                end
            end
            
            %%% Gaussian weighting models %%%
            
            if weightcalc==1 && nmod(ix)>2
                % Gaussian windowing for best model selection
                dgauss=0.1;
                meanmis=mean(misin);
                stdmis=std(misin);
                invers=1;
                step(1)=1000;
                flagg=2;
                while invers==1
                    t=min(misin):dgauss:max(misin)+dgauss;
                    gauss=(1/(stdmis*sqrt(2*pi)))*exp(-0.5.*((t-min(misin))/stdmis).^2);
                    coef=2*gauss*dgauss;
                    COEF=zeros(size(misin));
                    for gg=1:length(t)-1
                        flag1=find(misin>=t(gg) & misin<t(gg+1));
                        COEF(flag1)=coef(gg)/length(flag1);
                    end
                    step(flagg)=abs(sum(COEF)-1);
                    if step(flagg)>step(flagg-1) || flagg == 1000
                        invers=0;
                        dgauss=dgauss/0.75;
                        t=min(misin):dgauss:max(misin)+dgauss;
                        gauss=(1/(stdmis*sqrt(2*pi)))*exp(-0.5.*((t-min(misin))/stdmis).^2);
                        coef=2*gauss*dgauss;
                        COEF=zeros(size(misin));
                        for gg=1:length(t)-1
                            flag1=find(misin>=t(gg) & misin<t(gg+1));
                            COEF(flag1)=coef(gg)/length(flag1);
                        end
                    else
                        dgauss=dgauss*0.75;
                        flagg=flagg+1;
                    end
                end
                
                %%% Weighted layered model %%%
                
                THweight=sum(repmat(COEF,1,size(THin,1)).*THin',1);
                VSweight=sum(repmat(COEF,1,size(VSin,1)).*VSin',1)';
                VPweight=sum(repmat(COEF,1,size(VPin,1)).*VPin',1)';
                RHOweight=sum(repmat(COEF,1,size(RHOin,1)).*RHOin',1)';
                Zweight=zeros(size(VSweight));
                Zweighttmp=cumsum(THweight);
                mm=0;
                for ll=2:2:length(VSweight)-1
                    mm=mm+1;
                    Zweight(ll:ll+1)=Zweighttmp(mm);
                end
                Zweight(end)=maxdepth;
                
                %%% Weighted smoothed model %%%
                
                IVSweight=sum(repmat(COEF,1,size(IVs,2)).*IVs,1); % Weighted
                IVPweight=sum(repmat(COEF,1,size(IVp,2)).*IVp,1); % Weighted
                IRHOweight=sum(repmat(COEF,1,size(Irho,2)).*Irho,1); % Weighted
                
            end
            
            %%% Save velocity models %%%
            
            % Velocity model filenames
            % Best model
            nameDINbest=fullfile(dir_rep_ind,[num2str(XmidT(ix),xmidformat),...
                extens,'.Vms.best']); % Best model
            % Ridge
            nameDINridge=fullfile(dir_rep_ind,[num2str(XmidT(ix),xmidformat),...
                extens,'.Vms.ridge']); % Ridge model
            % Mean
            nameDINlay=fullfile(dir_rep_ind,[num2str(XmidT(ix),xmidformat),...
                extens,'.Vms.layered']); % Layered model
            nameDINsm=fullfile(dir_rep_ind,[num2str(XmidT(ix),xmidformat),...
                extens,'.Vms.smooth']); % Smooth model
            % Std
            nameDINlaystd=fullfile(dir_rep_ind,[num2str(XmidT(ix),xmidformat),...
                extens,'.VmsStd.layered']); % STD Layered model
            nameDINsmstd=fullfile(dir_rep_ind,[num2str(XmidT(ix),xmidformat),...
                extens,'.VmsStd.smooth']); % STD Smooth model
            % Min
            nameDINlaymin=fullfile(dir_rep_ind,[num2str(XmidT(ix),xmidformat),...
                extens,'.VmsMin.layered']); % Min Layered model
            nameDINsmmin=fullfile(dir_rep_ind,[num2str(XmidT(ix),xmidformat),...
                extens,'.VmsMin.smooth']); % Min Smooth model
            nameDINridgemin=fullfile(dir_rep_ind,[num2str(XmidT(ix),xmidformat),...
                extens,'.VmsMin.ridge']); % Min ridge model
            % Max
            nameDINlaymax=fullfile(dir_rep_ind,[num2str(XmidT(ix),xmidformat),...
                extens,'.VmsMax.layered']); % Max Layered model
            nameDINsmmax=fullfile(dir_rep_ind,[num2str(XmidT(ix),xmidformat),...
                extens,'.VmsMax.smooth']); % Max Smooth model
            nameDINridgemax=fullfile(dir_rep_ind,[num2str(XmidT(ix),xmidformat),...
                extens,'.VmsMax.ridge']); % Max ridge model
            % Weighted
            nameDINlayw=fullfile(dir_rep_ind,[num2str(XmidT(ix),xmidformat),...
                extens,'.Vws.layered']); % Weighted Layered model
            nameDINsmw=fullfile(dir_rep_ind,[num2str(XmidT(ix),xmidformat),...
                extens,'.Vws.smooth']); % Weighted Smooth model
            
            % Save best model
            dinsave(nameDINbest,THbest,VPbest(1:2:end),VSbest(1:2:end),RHObest(1:2:end));
            % Save ridge model
            if ridgecalc==1
                dinsave(nameDINridge,repmat(dz,1,length(ZZ)),VPridge,VSridge,RHOridge);
                dinsave(nameDINridgemin,repmat(dz,1,length(ZZ)),VPridge_min,VSridge_min,RHOridge_min);
                dinsave(nameDINridgemax,repmat(dz,1,length(ZZ)),VPridge_max,VSridge_max,RHOridge_max);
            end
            % Save mean models (layered, smoothed, layered smooth)
            dinsave(nameDINlay,THmean,VPmean(1:2:end),VSmean(1:2:end),RHOmean(1:2:end));
            dinsave(nameDINsm,repmat(dz,1,length(ZZ)),IVPmean,IVSmean,IRHOmean);
            dinsave(nameDINlaystd,THstd,VPstd(1:2:end),VSstd(1:2:end),RHOstd(1:2:end));
            dinsave(nameDINlaymin,THmin,VPmin(1:2:end),VSmin(1:2:end),RHOmin(1:2:end));
            dinsave(nameDINlaymax,THmax,VPmax(1:2:end),VSmax(1:2:end),RHOmax(1:2:end));
            dinsave(nameDINsmstd,repmat(dz,1,length(ZZ)),IVPstd,IVSstd,IRHOstd);
            dinsave(nameDINsmmin,repmat(dz,1,length(ZZ)),IVPmin,IVSmin,IRHOmin);
            dinsave(nameDINsmmax,repmat(dz,1,length(ZZ)),IVPmax,IVSmax,IRHOmax);
            if weightcalc==1 && nmod(ix)>2
                dinsave(nameDINlayw,THweight,VPweight(1:2:end),VSweight(1:2:end),RHOweight(1:2:end));
                dinsave(nameDINsmw,repmat(dz,1,length(ZZ)),IVPweight,IVSweight,IRHOweight);
            end
            
            % Plot and save average and weighted 1D models
            [f3,h0]=plot_curv(showplot,VSbest,Zbest,[],'-','k',[],1,1,0,fs,...
                'Vs (m/s)',depthtitle,[],[0 vsmax],[dpMIN dpMAX],[],[],[],...
                [],[],[],[0 0 24 18],[],0);
            hold on
            str0='Best model';
            h1=plot(VSmean,Zmean,'b-','linewidth',1.5);
            str1='Average layered model';
            h3=plot(IVSmean,ZZ,'r-','linewidth',1.5);
            str3='Average smooth model';
            if weightcalc==1 && nmod(ix)>2
                h4=plot(VSweight,Zweight,'c-','linewidth',1.5);
                str4='Weighted layered model';
                h6=plot(IVSweight,ZZ,'m-','linewidth',1.5);
                str6='Weighted smooth model';
            end
            if ridgecalc==1 && nmod(ix)>2
                h7=plot(VSridge,ZZ,'g-','linewidth',1.5);
                str7='Ridge model';
                dashline(VSridge_min,ZZ,2,2,2,2,'color','k','linewidth',1.5);
                dashline(VSridge_max,ZZ,2,2,2,2,'color','k','linewidth',1.5);
            end
            c=colorbar;
            if str2double(matrelease(1:4))>2014
                c.Label.String = 'Number of models';
                c.Label.Rotation = 270;
                c.Label.VerticalAlignment = 'Bottom';
            else
                c = cbhandle();
                cblabel('Number of models','Rotation', 270,'VerticalAlignment','Bottom');
            end
            if exist('NN','var')==1
                caxis([0 max(NN(:))]);
            end
            
            set(c,'visible','off');
            if (weightcalc==0 || (weightcalc==1 && nmod(ix)<=2)) && (ridgecalc==0 || (ridgecalc==1 && nmod(ix)<=2))
                h_legend=legend([h0,h1,h3],str0,str1,str3);
            elseif weightcalc==1 && ridgecalc==0
                h_legend=legend([h0,h1,h3,h4,h6],str0,str1,str3,str4,str6);
            elseif weightcalc==0 && ridgecalc==1
                h_legend=legend([h0,h1,h3,h7],str0,str1,str3,str7);
            else
                h_legend=legend([h0,h1,h3,h4,h6,h7],str0,str1,str3,str4,str6,str7);
            end
            set(h_legend,'FontSize',10,'linewidth',1,'location','northeast');
            hold off
            
            % Save figure
            file_mod=fullfile(dir_img_inv_mean,[num2str(XmidT(ix),xmidformat),'.mod1d.mean.',imgform]);
            save_fig(f3,file_mod,imgform,imgres,1,1-testplot);
            close(f3);
            % Concatenate ridge and average models figures
            if colnb>1
                colnb_tmp=2;
            else
                colnb_tmp=colnb;
            end
            filename_panel0=fullfile(dir_img_inv_mean,[num2str(XmidT(ix),xmidformat),...
                '.mod1d.mean.',imgform]);
            if exist('file_ridge','var')==1 && testplot==1 && nmod(ix)>1
                cat_img([file_mod,' ',file_ridge],imgform,colnb_tmp,'south',filename_panel0);
                if concat==1
                    delete(file_ridge);
                end
            end
        end
        
        %% %% %%
        
        %%%% Plot VS 2D pseudo-section %%%%
        
        if plot2dVS==1 && nmod(ix)>0 && calcmod==1 && Xlength>1
            % Get index corresponding to the investigation depth
            indf(ix)=round((maxdepth)/dz);
            % Look for topo index
            crit=abs(zround(ix)-depth);
            indi(ix)=find(crit==min(crit),1);
            if indf(ix)>length(depth)-indi(ix) || isnan(indf(ix))==1
                indf(ix)=length(depth)-indi(ix);
            end
            if indf(ix)>length(IVSmean)
                indf(ix)=length(IVSmean);
            end
            
            % Fill velocity matrix with average model
            layered(indi(ix):indi(ix)+indf(ix)-1,ix)=IVSmean(1:indf(ix));
            layered(layered==0)=NaN;
            % Plot pseudo-section
            f0=plot_img(10,XmidT,depth,layered,haxby(32),1,0,1,fs/2.5,...
                'X (m)','Elevation (m)','Vs (m/s)',[],...
                [floor(min(depth)/10)*10 ceil(max(depth)/10)*10],...
                [min(min(layered)) max(max(layered))],[],[],[],[],[],[],[],[],1);
            set(gca,'YDir','normal');
            drawnow;
        end
        
        %% %% %%
        
        %%%% Plot and save inversion results %%%%
        
        % Check number of colums for concatenation
        if colnb>nmodeinv(ix)+1
            colnb_tmp=nmodeinv(ix)+1;
        else
            colnb_tmp=colnb;
        end
        
        % Initialization of colormap
        if (plotinvres==1 || plotparam==1) && nmod(ix)>0
            % Create colormaps for selected and rejected models
            % Selected models (map2)
            if nmod(ix)==0
                col=[];
                youp=[];
            elseif nmod(ix)==1
                col=map2(1,:);
            else
                [col,youp]=createcolormap(misin,map2,Clogscale);
            end
            % Rejected models (map3)
            if nmod(ix)==nmaxmod(ix)*inv_set.nrun(ix)
                colout=[];
                youpout=[];
            elseif nmod(ix)==(nmaxmod(ix)*inv_set.nrun(ix))-1
                colout=map3(1,:);
            else
                [colout,youpout]=createcolormap(misout,map3,Clogscale);
            end
            nmod_total = inv_set.nrun(ix)*nmaxmod(ix);
            nmod_out = nmod_total - nmod(ix);
            if 2*std(misout)>min(misout)
                maxmisout=2*std(misout);
            else
                maxmisout=max(misout);
            end
            maxmisout=ceil(maxmisout*100)/100;
            colall=[colout;col];
        end
        
        %% %% %%
        
        if plotinvres==1 && nmod(ix)>0
            
            %% %% %%
            
            %%% Plot and save calculated dispersion curves %%%
            
            fprintf('\n  Saving calculated dispersion curves\n');
            for ip=1:nmodeinv(ix)
                m=modes(ip); % Mode number
                fprintf(['\n      Mode ',num2str(m),'\n']);
                tmp=D{ip}(I,:); % Read dispersion
                tmpin=tmp(OKsort,:); % Select selected models
                tmpout=tmp(KOsort,:); % Select rejected models
                tmpall=[tmpout;tmpin]; % Create matrix with rejected then selected models
                
                % Colorbar title
                str1=sprintf([' Accep. models ('...
                    num2str(nmod(ix)) ' / ' num2str(inv_set.nrun(ix)*nmaxmod(ix)) ')']);
                % Plot first dispersion curve
                [f1,~,~,~,c]=plot_curv(showplot,tmpall{1,1},1./tmpall{1,2},[],'-',colall(1,:),2,1,1,...
                    cbpos,fs,freqtitle_long,'Phase velocity (m/s)',str1,...
                    [fMIN fMAX],[VphMIN VphMAX],[],fticks,Vphticks,[],...
                    [],[],[0 0 24 18],[],0);
                hold on
                
                % Loop on all dispersion curves
                % Workaround for large number of plots
                if nmod_out > 5000
                    step_mod = 10;
                else
                    step_mod = 1;
                end
                % Plot rejected models
                for kk=2:step_mod:nmod_out
                    line(tmpall{kk,1}',1./tmpall{kk,2}','color',...
                        colall(kk,:),'linewidth',2); hold on;
                end
                % Workaround for large number of plots
                if nmod(ix) > 5000
                    step_mod = 3;
                else
                    step_mod = 1;
                end
                % Plot accepted models
                for ll=nmod_out+1:step_mod:nmod_total
                    line(tmpall{ll,1}',1./tmpall{ll,2}','color',...
                        colall(ll,:),'linewidth',2); hold on;
                end
                
                % Plot picked dispersion with errorbars
                plot(freqresamp{m+1},vresamp{m+1},'k+','linewidth',1.25,'markersize',4);
                if str2double(matrelease(1:4))>2014
                    h3=terrorbar(freqresamp{m+1},vresamp{m+1},deltaresamp{m+1},0.75,'units');
                    set(h3,'LineWidth',1.25,'Color',[0 0 0])
                else
                    h3=errorbar(freqresamp{m+1},vresamp{m+1},deltaresamp{m+1},'k+',...
                        'linewidth',1.25,'markersize',4);
                    xlimits=xlim;
                    tick_length=diff(xlimits)/100;
                    errorbar_tick(h3,tick_length,'units');
                end
                hold off
                colormap(map2);
                sizax=get(gca,'Position');
                set(c,'visible','off'); % Hide colorbar
                
                % Save figure
                filenamedsp=fullfile(dir_img_inv_1d,[num2str(XmidT(ix),xmidformat),...
                    '.calcdisp.M',num2str(m),'.',imgform]);
                save_fig(f1,filenamedsp,imgform,imgres,1,1-testplot);
                close(f1); clear f1;
            end
            % Get all filenames for all modes in one string
            filename_dspall=fullfile(dir_img_inv_1d,[num2str(XmidT(ix),xmidformat),...
                '.calcdisp.*.',imgform]);
            
            %% %% %%
            
            %%% Plot and save corresponding Vs models %%%
            
            % Get the 1D Vs average model to plot
            if plot1dVS==1
                if calcmod==1 % Get model calculated in current run
                    if strcmp(modeltype,'best')==1
                        VSplot=VSbest;
                        Zplot=Zbest;
                    elseif strcmp(modeltype,'layered')==1
                        if strcmp(avertype,'Vms')==1
                            VSplot=VSmean;
                            Zplot=Zmean;
                        else
                            VSplot=VSweight;
                            Zplot=Zweight;
                        end
                    elseif strcmp(modeltype,'smooth')==1
                        if strcmp(avertype,'Vms')==1
                            VSplot=IVSmean;
                            Zplot=ZZ;
                        else
                            if nmod(ix)>1
                                VSplot=IVSweight;
                                Zplot=ZZ;
                            else
                                VSplot=[];
                                Zplot=[];
                            end
                        end
                    elseif strcmp(modeltype,'ridge')==1
                        VSplot=VSridge;
                        Zplot=ZZ;
                    end
                else % Get model from file saved in previous run (obsolete)
                    % Velocity file name
                    filevel=fullfile(dir_rep_ind,[num2str(XmidT(ix),xmidformat),extens,'.',...
                        avertype,'.',modeltype]);
                    if exist(filevel,'file')==2
                        modvel=dlmread(filevel,'',1,0);
                        moddepth=[0;cumsum(modvel(:,1))];
                        if maxdepth>moddepth(end)
                            moddepth(end)=maxdepth;
                        else
                            modvel=modvel(moddepth<maxdepth,:);
                            moddepth=[0;cumsum(modvel(:,1))];
                        end
                        vssw=modvel(:,3);
                        
                        VSplot=[vssw;vssw(end)];
                        VSplotrep=repmat(VSplot,1,2)';
                        VSplot=VSplotrep(:)'; VSplot=VSplot(1:end-2)';
                        moddepthrep=repmat(moddepth,1,2)';
                        Zplot=moddepthrep(:)'; Zplot=Zplot(2:end-1)';
                    else
                        VSplot=[]; Zplot=[];
                    end
                end
            end
            
            % Plot all calculated models
            fprintf('\n  Saving calculated Vs models\n');
            tmpall=[VSout,VSin];
            zall=[Zout,Zin];
            
            % Colorbar titles for accepeted and rejected models
            str1=sprintf([' Accep. models ('...
                num2str(nmod(ix)) ' / ' num2str(inv_set.nrun(ix)*nmaxmod(ix)) ')']);
            str2=sprintf([' Rejec. models (' num2str((inv_set.nrun(ix)*nmaxmod(ix))-nmod(ix))...
                ' / ' num2str(inv_set.nrun(ix)*nmaxmod(ix)) ')']);
            % Plot first model
            [f4,~,~,~,c]=plot_curv(showplot,tmpall(:,1),zall(:,1),[],[],colall(1,:),2,1,1,cbpos,fs,...
                'Vs (m/s)',depthtitle,str1,[vsMIN vsMAX],[dpMIN dpMAX],[],vsticks,dticks,...
                [],[],[],[0 0 24 18],sizax,0);
            
            % Loop on all models
            % Workaround for large number of plots
            if nmod_out > 5000
                step_mod = 10;
            else
                step_mod = 1;
            end
            % Plot rejected models
            for kk=2:step_mod:nmod_out
                line(tmpall(:,kk)',zall(:,kk)','color',...
                    colall(kk,:),'linewidth',2); hold on;
            end
            % Workaround for large number of plots
            if nmod(ix) > 5000
                step_mod = 3;
            else
                step_mod = 1;
            end
            % Plot accepted models
            for ll=nmod_out+1:step_mod:nmod_total
                line(tmpall(:,ll)',zall(:,ll)','color',...
                    colall(ll,:),'linewidth',2); hold on;
            end
            
            colormap(map2);
            hold on
            % Plot 1d model
            if plot1dVS==1
                dashline(VSplot,Zplot,2,2,2,2,'color','k','linewidth',2);
            end
            hold off
            set(c,'visible','off'); % Hide colorbar
            
            % Save figure
            filename_modall=fullfile(dir_img_inv_1d,[num2str(XmidT(ix),xmidformat),...
                '.mod1d.all.',imgform]);
            save_fig(f4,filename_modall,imgform,imgres,1,1-testplot);
            close(f4); clear f4;
            
            filename_cbin=fullfile(dir_img_inv_1d,[num2str(XmidT(ix),xmidformat),...
                '.cbmisin.',imgform]);
            filename_cbout=fullfile(dir_img_inv_1d,[num2str(XmidT(ix),xmidformat),...
                '.cbmisout.',imgform]);
            
            % Save colorbars in separate file for nice display
            if nmod(ix)>1
                if cbpos==1
%                     cbticks=youp(1:3:end);
                    cbticks = logspace(log10(min(misin)),log10(max(misin)),6);
                else
%                     cbticks=youp(1:4:end);
                    cbticks = logspace(log10(min(misin)),log10(max(misin)),4);
                end
                f4 = plot_colorbar(showplot,[24, 1], cbpos, str1, map2, Clogscale, cbticks,...
                    [min(misin) max(misin)], fs, [0 0 24 18], sizax, 2);
                save_fig(f4,filename_cbin,imgform,imgres,1,0);
                close(f4); clear f4;
            end
            if nmod(ix)<inv_set.nrun(ix)*nmaxmod(ix)-1
                if cbpos==1
%                     cbticks=youpout(1:3:end);
                    cbticks = logspace(log10(min(misout)),log10(maxmisout),6);
                else
%                     cbticks=youpout(1:4:end);
                    cbticks = logspace(log10(min(misout)),log10(maxmisout),4);
                end
                f4 = plot_colorbar(showplot,[24, 1], cbpos, str2, map3, Clogscale, cbticks,...
                    [min(misout) maxmisout], fs, [0 0 24 18], sizax,1);
                save_fig(f4,filename_cbout,imgform,imgres,1,0);
                close(f4); clear f4;
            end
            
            %% %% %%
            
            %%% Concatenate dispersion curves and Vs models %%%
            
            if strcmp(imgform,'fig')~=1 && testplot==1
                fprintf('\n  Concatenate figures\n');

                if strcmp(imgform,'pdf')~=1 % Resize raster images
                    if isunix==1
                        [~,sizefig]=unix(['convert ',filename_dspall,' -format "%w,%h" info:']);
                        sizefig=str2num(sizefig); sizefig=sizefig(1,:);
                        [~,~]=unix(['convert ',filename_modall,' -gravity northeast -extent ',...
                            num2str(sizefig(1)),'x',num2str(sizefig(2)),' ',filename_modall]);
                    else
                        com1=['img_convert ',filename_dspall,' -format "%w,%h" info:'];
                        com1=strrep(com1,'\','/');
                        [~,sizefig]=unix(com1);
                        sizefig=str2num(sizefig); sizefig=sizefig(1,:);
                        com1=['img_convert ',filename_modall,' -gravity northeast -extent ',...
                            num2str(sizefig(1)),'x',num2str(sizefig(2)),' ',filename_modall];
                        com1=strrep(com1,'\','/');
                        [~,~]=unix(com1);
                    end
                end
                
                % Concatenate dispersion and 1D models
                filename_imgtmp=fullfile(dir_img_inv_1d,[num2str(XmidT(ix),xmidformat),...
                    '.imgtmp.',imgform]);
                cat_img([filename_dspall,' ',filename_modall],imgform,colnb_tmp,'east',filename_imgtmp,0);
                
                
                filename_cbtmp=fullfile(dir_img_inv_1d,[num2str(XmidT(ix),xmidformat),...
                    '.cb.',imgform]);
                filename_panel=fullfile(dir_img_inv_1d,[num2str(XmidT(ix),xmidformat),...
                    '.invresults.',imgform]);
                if  nmod(ix)>1 && nmod(ix)<inv_set.nrun(ix)*nmaxmod(ix)-1
                    % Concatenate colorbars
                    if (cbpos==1 && colnb_tmp==nmodeinv(ix)+1) || (cbpos==2 && colnb_tmp>1)
                        cat_img([filename_cbout,' ',filename_cbin],imgform,2,'center',filename_cbtmp,0);
                    elseif (cbpos==1 && colnb_tmp<nmodeinv(ix)+1) || (cbpos==2 && colnb_tmp==1)
                        cat_img([filename_cbout,' ',filename_cbin],imgform,1,'east',filename_cbtmp,0);
                    end
                    % Resize raster
                    if cbpos==2 && strcmp(imgform,'pdf')~=1
                        if isunix==1
                            [~,sizefig_cb]=unix(['convert ',filename_cbtmp,' -format "%h" info:']);
                            sizefig_cb=str2double(sizefig_cb);
                            [~,sizefig]=unix(['convert ',filename_imgtmp,' -format "%w" info:']);
                            sizefig=str2double(sizefig);
                            [~,~]=unix(['convert ',filename_cbtmp,' -gravity center -extent ',...
                                num2str(sizefig),'x',num2str(sizefig_cb),' ',filename_cbtmp]);
                        else
                            com1=['img_convert ',filename_cbtmp,' -format "%h" info:'];
                            com1=strrep(com1,'\','/');
                            [~,sizefig_cb]=unix(com1);
                            sizefig_cb=str2double(sizefig_cb);
                            com1=['img_convert ',filename_imgtmp,' -format "%w" info:'];
                            com1=strrep(com1,'\','/');
                            [~,sizefig]=unix(com1);
                            sizefig=str2double(sizefig);
                            com1=['img_convert ',filename_cbtmp,' -gravity center -extent ',...
                                num2str(sizefig),'x',num2str(sizefig_cb),' ',filename_cbtmp];
                            com1=strrep(com1,'\','/');
                            [~,~]=unix(com1);
                        end
                    end
                elseif nmod(ix)==1
                    movefile(filename_cbout,filename_cbtmp)
                elseif nmod(ix)>=inv_set.nrun(ix)*nmaxmod(ix)-1
                    movefile(filename_cbin,filename_cbtmp)
                end
                
                % Concatenate all
                if cbpos==1 && colnb_tmp==nmodeinv(ix)+1
                    cat_img([filename_imgtmp,' ',filename_cbtmp],imgform,2,'south',filename_panel);
                elseif cbpos==1 && colnb_tmp<nmodeinv(ix)+1
                    cat_img([filename_imgtmp,' ',filename_cbtmp],imgform,2,'center',filename_panel);
                elseif cbpos==2 && colnb_tmp==1
                    cat_img([filename_imgtmp,' ',filename_cbtmp],imgform,1,'east',filename_panel);
                else
                    cat_img([filename_imgtmp,' ',filename_cbtmp],imgform,1,'center',filename_panel);
                end
                if concat==1 % Delete single image files
                    if nmod(ix)>1 && nmod(ix)<inv_set.nrun(ix)*nmaxmod(ix)-1
                        delete(filename_cbin,filename_cbout,filename_cbtmp,filename_imgtmp,filename_modall,filename_dspall);
                    elseif nmod(ix)==1
                        delete(filename_cbout,filename_cbtmp,filename_imgtmp,filename_modall,filename_dspall);
                    elseif nmod(ix)>=inv_set.nrun(ix)*nmaxmod(ix)-1
                        delete(filename_cbin,filename_cbtmp,filename_imgtmp,filename_modall,filename_dspall);
                    end
                end
            end
        end
        
        %% %% %%
        
        %%%% Plot and save inversion parameters %%%%
        
        % Get selected parameters settings
        if plotparam==1 && nmod(ix)>0
            fprintf('\n  Saving inversion parameters\n');
            % First parameter
            if strcmp(param1,'Vs')==1
                axtit1='Vs (m/s)';
                p1=[VSout,VSin];
                p1=p1(1:2:end,:);
                xticks=vsticks;
                XMIN=[]; XMAX=[];
            elseif strcmp(param1,'Th')==1
                axtit1='Thickness (m)';
                p1=[Zout,Zin];
                p1=diff(p1);
                p1=p1(1:2:end,:);
                XMIN=[]; XMAX=[]; xticks=[];
            elseif strcmp(param1,'Dp')==1
                axtit1='Depth (m)';
                p1=[Zout,Zin];
                p1=p1(1:2:end,:);
                XMIN=[]; XMAX=[]; xticks=[];
            elseif strcmp(param1,'Vp')==1
                axtit1='Vp (m/s)';
                p1=[VPout,VPin];
                p1=p1(1:2:end,:);
                xticks=vpticks;
                XMIN=[]; XMAX=[];
            elseif strcmp(param1,'Dens')==1
                axtit1='Density (kg/m^3)';
                p1=[RHOout,RHOin];
                p1=p1(1:2:end,:);
                xticks=rhoticks;
                XMIN=[]; XMAX=[];
            end
            % Second parameter
            if strcmp(param2,'Vs')==1
                axtit2='Vs (m/s)';
                p2=[VSout,VSin];
                p2=p2(1:2:end,:);
                yticks=vsticks;
                YMIN=[]; YMAX=[];
            elseif strcmp(param2,'Th')==1
                axtit2='Thickness (m)';
                p2=[Zout,Zin];
                p2=diff(p2);
                p2=p2(1:2:end,:);
                YMIN=[]; YMAX=[]; yticks=[];
            elseif strcmp(param2,'Dp')==1
                axtit2='Depth (m)';
                p2=[Zout,Zin];
                p2=p2(1:2:end,:);
                YMIN=[]; YMAX=[]; yticks=[];
            elseif strcmp(param2,'Vp')==1
                axtit2='Vp (m/s)';
                p2=[VPout,VPin];
                p2=p2(1:2:end,:);
                yticks=vpticks;
                YMIN=[]; YMAX=[];
            elseif strcmp(param2,'Dens')==1
                axtit2='Density (kg/m^3)';
                p2=[RHOout,RHOin];
                p2=p2(1:2:end,:);
                yticks=rhoticks;
                YMIN=[];YMAX=[];
            end
            
            % Loop over all selected parameters
            for i=1:length(np1)
                for j=1:length(np2)
                    fprintf(['\n      ',param1,num2str(np1(i)),' vs ',param2,num2str(np2(j)),'\n']);
                    if exist('sizax','var')==0
                        sizax=[];
                    end
                    % Colorbar titles
                    str1=sprintf([' Accep. models ('...
                        num2str(nmod(ix)) ' / ' num2str(inv_set.nrun(ix)*nmaxmod(ix)) ')']);
                    str2=sprintf([' Rejec. models (' num2str((inv_set.nrun(ix)*nmaxmod(ix))-nmod(ix))...
                        ' / ' num2str(inv_set.nrun(ix)*nmaxmod(ix)) ')']);
                    
                    if length(unique(p1(np1(i),:)))==1
                        XMIN = unique(p1(np1(i),:)) - 0.1*unique(p1(np1(i),:));
                        XMAX = unique(p1(np1(i),:)) + 0.1*unique(p1(np1(i),:));
                    end
                    
                    if length(unique(p2(np2(j),:)))==1
                        YMIN = unique(p2(np2(j),:)) - 0.1*unique(p2(np2(j),:));
                        YMAX = unique(p2(np2(j),:)) + 0.1*unique(p2(np2(j),:));
                    end
                    
                    % Plot first model parameter
                    [f4,~,~,~,c]=plot_curv(showplot,p1(np1(i),:),p2(np2(j),:),[],'.','k',2,1,1,...
                        cbpos,fs,[axtit1,' [layer ',num2str(np1(i)),']'],[axtit2,' [layer ',num2str(np2(j)),']'],...
                        str1,[XMIN XMAX],[YMIN YMAX],[],vsticks,dticks,[],...
                        [],[],[0 0 24 18],sizax,0);
                    hold on
                    
                    % Loop on all models
                    % Workaround for large number of plots
                    if nmod_out > 5000
                        step_mod = 10;
                    else
                        step_mod = 1;
                    end
                    % Plot rejected models
                    for ii=2:step_mod:nmod_out
                        line(p1(np1(i),ii),p2(np2(j),ii),'color',colall(ii,:),'markerfacecolor',colall(ii,:),...
                            'markersize',5,'marker','o','linestyle','none');
                    end
                    % Workaround for large number of plots
                    if nmod(ix) > 5000
                        step_mod = 3;
                    else
                        step_mod = 1;
                    end
                    % Plot accepted models
                    for ii=nmod_out+1:step_mod:nmod_total
                        line(p1(np1(i),ii),p2(np2(j),ii),'color',colall(ii,:),'markerfacecolor',colall(ii,:),...
                            'markersize',5,'marker','o','linestyle','none');
                    end

                    hold off
                    colormap(map2);
                    sizax=get(gca,'Position');
                    set(c,'visible','off'); % Hide colorbar
                    
                    % Save figure
                    filename=fullfile(dir_img_inv_param,[num2str(XmidT(ix),xmidformat),...
                        '.param.',param1,num2str(np1(i)),param2,num2str(np2(j)),'.',imgform]);
                    save_fig(f4,filename,imgform,imgres,1,0);
                    close(f4); clear f4;
                    
                    filename_cbin=fullfile(dir_img_inv_param,[num2str(XmidT(ix),xmidformat),...
                        '.cbmisin.',imgform]);
                    filename_cbout=fullfile(dir_img_inv_param,[num2str(XmidT(ix),xmidformat),...
                        '.cbmisout.',imgform]);
                    
                    % Save colorbars
                    if i==1 && j==1
                        if nmod(ix)>1
                            if cbpos==1
                                %                     cbticks=youp(1:3:end);
                                cbticks = logspace(log10(min(misin)),log10(max(misin)),6);
                            else
                                %                     cbticks=youp(1:4:end);
                                cbticks = logspace(log10(min(misin)),log10(max(misin)),4);
                            end
                            f4 = plot_colorbar(showplot,[24, 1], cbpos, str1, map2, Clogscale, cbticks,...
                                [min(misin) max(misin)], fs, [0 0 24 18], sizax, 2);
                            save_fig(f4,filename_cbin,imgform,imgres,1,0);
                            close(f4); clear f4;
                        end
                        if nmod(ix)<inv_set.nrun(ix)*nmaxmod(ix)-1
                            if cbpos==1
                                %                     cbticks=youpout(1:3:end);
                                cbticks = logspace(log10(min(misout)),log10(maxmisout),6);
                            else
                                %                     cbticks=youpout(1:4:end);
                                cbticks = logspace(log10(min(misout)),log10(maxmisout),4);
                            end
                            f4 = plot_colorbar(showplot,[24, 1], cbpos, str2, map3, Clogscale, cbticks,...
                                [min(misout) maxmisout], fs, [0 0 24 18], sizax,1);
                            save_fig(f4,filename_cbout,imgform,imgres,1,0);
                            close(f4); clear f4;
                        end
                    end
                end
            end
            
            % Check number of column for concatenation
            if colnb>length(np1)*length(np2)
                colnb_tmp=length(np1)*length(np2);
            else
                colnb_tmp=colnb;
            end
            
            % Concatenation
            if testplot==1
                % Concatenate parameter plots
                filename_paramall=fullfile(dir_img_inv_param,[num2str(XmidT(ix),xmidformat),...
                    '.param.*.',imgform]);
                filename_imgtmp=fullfile(dir_img_inv_param,[num2str(XmidT(ix),xmidformat),...
                    '.imgtmp.',imgform]);
                cat_img(filename_paramall,imgform,colnb_tmp,'east',filename_imgtmp,0);
                
                % Concatenate colorbars
                filename_cbtmp=fullfile(dir_img_inv_param,[num2str(XmidT(ix),xmidformat),...
                    '.cb.',imgform]);
                if (cbpos==1 && colnb_tmp==length(np1)*length(np2)) || (cbpos==2 && colnb_tmp>1)
                    cat_img([filename_cbout,' ',filename_cbin],imgform,2,'center',filename_cbtmp,0);
                elseif (cbpos==1 && colnb_tmp<length(np1)*length(np2)) || (cbpos==2 && colnb_tmp==1)
                    cat_img([filename_cbout,' ',filename_cbin],imgform,1,'east',filename_cbtmp,0);
                end
                
                % Resize rasters
                if cbpos==2 && strcmp(imgform,'pdf')~=1
                    if isunix==1
                        [~,sizefig_cb]=unix(['convert ',filename_cbtmp,' -format "%h" info:']);
                        sizefig_cb=str2double(sizefig_cb);
                        [~,sizefig]=unix(['convert ',filename_imgtmp,' -format "%w" info:']);
                        sizefig=str2double(sizefig);
                        [~,~]=unix(['convert ',filename_cbtmp,' -gravity center -extent ',...
                            num2str(sizefig),'x',num2str(sizefig_cb),' ',filename_cbtmp]);
                    else
                        com1=['img_convert ',filename_cbtmp,' -format "%h" info:'];
                        com1=strrep(com1,'\','/');
                        [~,sizefig_cb]=unix(com1);
                        sizefig_cb=str2double(sizefig_cb);
                        com1=['img_convert ',filename_imgtmp,' -format "%w" info:'];
                        com1=strrep(com1,'\','/');
                        [~,sizefig]=unix(com1);
                        sizefig=str2double(sizefig);
                        com1=['img_convert ',filename_cbtmp,' -gravity center -extent ',...
                            num2str(sizefig),'x',num2str(sizefig_cb),' ',filename_cbtmp];
                        com1=strrep(com1,'\','/');
                        [~,~]=unix(com1);
                    end
                end
                
                % Concatenate all
                filename_param=fullfile(dir_img_inv_param,[num2str(XmidT(ix),xmidformat),...
                    '.param_',param1,num2str(np1(1)),'_',num2str(np1(end)),...
                    param2,num2str(np2(1)),'_',num2str(np2(end)),'.',imgform]);
                if cbpos==1 && colnb_tmp==length(np1)*length(np2)
                    cat_img([filename_imgtmp,' ',filename_cbtmp],imgform,2,'south',filename_param);
                elseif cbpos==1 && colnb_tmp<length(np1)*length(np2)
                    cat_img([filename_imgtmp,' ',filename_cbtmp],imgform,2,'center',filename_param);
                elseif cbpos==2 && colnb_tmp==1
                    cat_img([filename_imgtmp,' ',filename_cbtmp],imgform,1,'east',filename_param);
                else
                    cat_img([filename_imgtmp,' ',filename_cbtmp],imgform,1,'center',filename_param);
                end
                if concat==1 % Delete single image files
                    delete(filename_cbin,filename_cbout,filename_cbtmp,filename_imgtmp,filename_paramall);
                end
            end
        end
    end
    % End clock for single Xmid
    tend=toc(tstart);
    [time_string]=secs2hms(tend);
    if inversion==1 || calcmod==1 || plotinvres==1 || plotparam==1
        fprintf(['\n  Elapsed time: ',time_string,'\n']);
        fprintf('\n  **********************************************************');
        fprintf('\n  **********************************************************\n');
    end
end

%% %% %%

%%%%%% Plot and save post-inversion settings %%%%%%

if sum(isnan(nmodeinv))<Xlength && sum(isnan(nC))<Xlength
    fprintf('\n  Saving inversion QC\n');
    f5=figure(5);
    set(f5,'Units','centimeters');
    set(f5,'Position',[0,0,24,18])
    % Lambda max for each Xmid
    f11=subplot(5,1,1);
    plot(XmidT,lmaxpick,'k.','linewidth',1)
    xlabel('X (m)');ylabel('\lambda_{max}');
    xlim([xmin xmax]); ylim([0 10*ceil(max(lmaxpick)/10)]);
    set(gca,'TickDir','out','linewidth',1,'XMinorTick','on','YMinorTick','on');
    % Nb of modes for each Xmid
    f12=subplot(5,1,2);
    plot(XmidT,nmodeinv,'k.','linewidth',1)
    xlabel('X (m)');ylabel('Nb of modes');
    xlim([xmin xmax]); ylim([0 max(nmodeinv)+1]);
    set(gca,'TickDir','out','linewidth',1,'XMinorTick','on','YMinorTick','on');
    % Nb of layers for each Xmid
    f13=subplot(5,1,3);
    plot(XmidT,inv_set.nC,'k.','linewidth',1)
    set(gca,'YTick',0:2:50);
    xlabel('X (m)');ylabel('Nb of layers');
    xlim([xmin xmax]); ylim([0 max(inv_set.nC)+1]);
    set(gca,'TickDir','out','linewidth',1,'XMinorTick','on','YMinorTick','on');
    % Nb of models for each Xmid
    f14=subplot(5,1,4);
    semilogy(XmidT,inv_set.nmod,'k.','linewidth',1)
    xlabel('X (m)');ylabel('Accepted models');
    xlim([xmin xmax]); ylim([1 nmaxmod(ix)*max(inv_set.nrun(ix))]);
    set(gca,'Ytick',[10^0 10^1 10^2 10^3 10^4 10^5 ...
        10^6 10^7 10^8 10^9 10^10 10^11]);
    set(gca,'TickDir','out','linewidth',1,'XMinorTick','on','YMinorTick','on');
    % Lambda max for each Xmid
    f15=subplot(5,1,5);
    plot(XmidT,inv_set.minmis,'k.','linewidth',1)
    xlabel('X (m)');ylabel('Misfit min.');
    xlim([xmin xmax]); ylim([0 max(inv_set.minmis)+0.05]);
    set(gca,'TickDir','out','linewidth',1.5,'XMinorTick','on','YMinorTick','on');
    % Save figure
    filename=fullfile(dir_img_inv_mod,['C_SWIPinv_summary.',imgform]);
    save_fig(f5,filename,imgform,imgres,1);
    close(f5);
end

% Check elapsed time
Tend=toc(Tstart);
[time_string]=secs2hms(Tend);
fprintf(['\n  Total elapsed time: ',time_string,'\n']);
if length(Xmidselec)>1
    [time_string]=secs2hms(Tend/max(sum(nshot(Xmidselec,:)>0)));
    fprintf(['  Average time per Xmid: ',time_string,'\n']);
end