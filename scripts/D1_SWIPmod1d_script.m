%%% SURFACE-WAVE dispersion INVERSION & PROFILING (SWIP)
%%% MODULE D1 : SWIPmod1d.m
%%% S. Pasquet - V17.01.20
%%% SWIPmod1d.m plots observed and calculated dispersion for each Xmid
%%% It also plots 1D Vp, Vs, Vp/Vs and Poisson's ratio models

run('SWIP_defaultsettings')

if swip==0 && tomo==0 && user==0
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n   Select at least one model for forward calculation');
    fprintf('\n   Set either "swip", "tomo" or "user" to 1');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
    return
end

if plot1dcal==0 && plot1dmod==0
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n   Select at least one plot option');
    fprintf('\n   Set either "plot1dcal" or "plot1dmod" to 1');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
    return
end

% Initialization (same for D1 and D2)
if swip==1
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
else
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
dx=acquiparam.dx;
topo=acquiparam.topo;
nWmin=stackdisp.nWmin;
nWmax=stackdisp.nWmax;
xsca=pomega.xsca;
XmidT=xmidparam.XmidT; % Get Xmids
Xlength=length(XmidT); % Get Xmids
xmidformat=stackdisp.xmidformat;
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
zround=xmidparam.zround; % Get topography
if isempty(dpMAX)==1
    if swip==1
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
maxdepth=dpMAX;
depth=max(zround):-dz:min(zround)-maxdepth; % Depth vector with topo
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
    dir_img_inv_mod=fullfile(dir_img_inv,['models',extens]);
    dir_img_inv_1d=fullfile(dir_img_inv_mod,'1dmodels');
else
    dir_img_ind=fullfile(dir_img,'Usermodels');
    dir_img_inv_1d=fullfile(dir_img_ind,'1dmodels');
end
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

% Select velocity and STD models
if usevptomo==1 || tomo==1
    fprintf('\n  Select Vp model file\n');
    [filevel,pathvel]=uigetfile({'*.model;*.dat;*.xzv;*.txt'},'Select Vp model');
    Vpfile=fullfile(pathvel,filevel); % File with velocity (3 columns X,Z,Vp)
    if pathvel==0
        VpI=[]; VpItomo=[]; VsItomo=[];
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        fprintf('\n   No Vp model file selected - Ignore Vptomo');
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
    else
        try
            VpI=readtomo(Vpfile,0,XmidT,depth,xsca,vpaver,mean([nWmin,nWmax]),dx); % Read Vp tomo file
            VpItomo=readtomo(Vpfile,0,XmidT,depth,xsca); % Read Vp tomo file
        catch
            fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
            fprintf('\n   Invalid Vp model file - Ignore Vptomo');
            fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
            VpI=[]; VpItomo=[]; VsItomo=[];
        end
    end
    %     if plot1dstd==1
    %         [filevelstd,pathvelvelstd]=uigetfile({'*.model;*.dat;*.xzv;*.txt'},'Select Vp STD model');
    %         if pathvelvelstd==0
    %             VpIstd=[]; VpIstdtomo=[];
    %         else
    %             Vpstdfile=fullfile(pathvelstd,filevelstd); % File with velocity STD (3 columns X,Z,VpSTD)
    %             VpIstd=readtomo(Vpstdfile,0,XmidT,depth,xsca,vpaver,mean([nWmin,nWmax]),dx); % Read Vp STD tomo file
    %             VpIstdtomo=readtomo(Vpstdfile,0,XmidT,depth,xsca); % Read Vp STD tomo file
    %         end
    %     end
    if tomo==1 && isempty(VpItomo)==0
        fprintf('\n  Select Vs model file\n');
        [filevel,pathvel]=uigetfile({'*.model;*.dat;*.xzv;*.txt'},'Select Vs model');
        if pathvel==0
            VsItomo=[]; VpItomo=[];
            fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
            fprintf('\n   No Vs model file selected - Ignore Vstomo');
            fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
        else
            Vsfile=fullfile(pathvel,filevel); % File with velocity (3 columns X,Z,Vs)
            try
                VsItomo=readtomo(Vsfile,0,XmidT,depth,xsca); % Read Vs tomo file
                VpItomo(isnan(VsItomo)==1)=NaN;
                VsItomo(isnan(VpItomo)==1)=NaN;
                VsItomo=VsItomo/2;
            catch
                fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
                fprintf('\n   Invalid Vs model file - Ignore Vstomo');
                fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
                VsItomo=[]; VpItomo=[];
            end
            
        end
        %         if plot1dstd==1
        %             [filevelstd,pathvelvelstd]=uigetfile({'*.model;*.dat;*.xzv;*.txt'},'Select Vs STD model');
        %             if pathvelvelstd==0
        %                 VsIstdtomo=[];
        %             else
        %                 Vsstdfile=fullfile(pathvelstd,filevelstd); % File with velocity STD (3 columns X,Z,VsSTD)
        %                 VsIstdtomo=readtomo(Vsstdfile,0,XmidT,depth,xsca); % Read Vs STD tomo file
        %                 VpIstdtomo(isnan(VsIstdtomo)==1)=NaN;
        %                 VsIstdtomo(isnan(VpIstdtomo)==1)=NaN;
        %             end
        %         end
    end
else
    zinc=0:dz:maxdepth;
end

% Read user model
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

[testimgmgck,~]=unix('which montage');
[testpdfjam,~]=unix('which pdfjam');
testplot=((testpdfjam==0 && strcmp(imgform,'pdf')==1) || (testimgmgck==0 && strcmp(imgform,'pdf')==0 && strcmp(imgform,'fig')==0));
if concat == 0
    testplot = 0;
end

showplotinit=showplot;

fprintf('\n  **********************************************************');
fprintf('\n  **********************************************************\n');

for ix=Xmidselec
    if Xlength>1
        close all
        showplot=showplotinit;
    end
    if sum(nshot(ix,:))>=0
        fprintf(['\n  Xmid',num2str(ix),' = ',num2str(XmidT(ix),xmidformat),' m\n']);
    end
    if swip==1
        dir_rep_ind=[dir_rep_inv,'/',num2str(XmidT(ix),xmidformat),'_reports'];
    end
    dir_img_ind=[dir_img_inv_1d,'/mod1d_',num2str(XmidT(ix),xmidformat)];
    if exist(dir_img_ind,'dir')~=7
        mkdir(dir_img_ind);
    end
    
    if plot1dcal==1
        if plotflim==1
            flimsing=flim(ix);
        else
            flimsing=[];
        end
        dspfile_sum=fullfile(dir_dat,[num2str(XmidT(ix),xmidformat),'.sum.dsp']);
        if exist(dspfile_sum,'file')==2
            [dspmat,f,v]=dsp2dat(dspfile_sum,flip,0);
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
            f1=plot_curv(showplot,NaN,NaN,[],'.',[0 0 0],[],axetop,axerev,...
                0,fs,freqtitle_long,'Phase velocity (m/s)',[],...
                [fMIN fMAX],[VphMIN VphMAX],[],fticks,Vphticks,[],...
                [],[],[0 0 24 18],[]);
        end
        hold on
        hh=dashline(f,f,2,2,2,2,'color',[0 0 1],'linewidth',2);
        set(hh,'visible','off');
        if plotlamlim==1 && sampling==1
            dashline(f,f*max(resampvec),2,2,2,2,'color',[0 0 1],'linewidth',2);
        end
        if Flogscale==1
            set(gca,'xscale','log');
        end
        if plot1dobs==1
            if swip==0
                nametarg=fullfile(dir_targ,[num2str(XmidT(ix),xmidformat),'.target']);
            else
                nametarg=fullfile(dir_rep_ind,[num2str(XmidT(ix),xmidformat),'.target']);
                if exist(nametarg,'file')~=2 && sum(nshot(ix,:))>=0
                    if exist(nametarg(1:end-3),'file')==2
                        movefile(nametarg(1:end-3),nametarg);
                    else
                        fprintf(['\n  No dispersion picked for Xmid',num2str(ix),' = ',...
                            num2str(XmidT(ix),xmidformat),' m\n']);
                    end
                end
            end
            if exist(nametarg,'file')~=2
                if swip==1
                    npvc=0;
                else
                    pvcstruct=dir(fullfile(dir_pick,[num2str(XmidT(ix),xmidformat),'.*.pvc']));
                    npvc=length(pvcstruct);
                    if npvc==0 && sum(nshot(ix,:))>=0
                        fprintf(['\n  No dispersion picked for Xmid',num2str(ix),' = ',...
                            num2str(XmidT(ix),xmidformat),' m\n']);
                    end
                    modes=zeros(npvc,1)*NaN;
                    freqresamp=cell(npvc,1);
                    vresamp=freqresamp; deltaresamp=freqresamp;
                    for ip=1:npvc
                        pvcfile=pvcstruct(ip).name;
                        modes(ip)=str2double(pvcfile(end-4)); % Mode number
                        if modes(ip)>maxmodeinv(ix)
                            break
                        end
                        Vprev=load(fullfile(dir_pick,pvcfile));
                        % Resample in lambda or frequency
                        [freqresamp{modes(ip)+1},vresamp{modes(ip)+1},deltaresamp{modes(ip)+1}]=...
                            resampvel(Vprev(:,1),Vprev(:,2),Vprev(:,3),resampvec,sampling,1);
                    end
                end
            else
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
            
            for ip=1:npvc
                if mod(modes(ip),2)==0
                    col=pickcol1;
                else
                    col=pickcol2;
                end
                
                plot(freqresamp{modes(ip)+1},vresamp{modes(ip)+1},'.','linewidth',1.5,'markersize',9);
                if eb==1
                    if str2double(matrelease(1:4))>2014
                        han=terrorbar(freqresamp{modes(ip)+1},vresamp{modes(ip)+1},deltaresamp{modes(ip)+1},1,'units');
                        set(han,'LineWidth',1.5,'Color',col)
                    else
                        han=errorbar(freqresamp{modes(ip)+1},vresamp{modes(ip)+1},deltaresamp{modes(ip)+1},...
                            '.','Color',col,'linewidth',1.5,'markersize',9);
                        xlimits=xlim;
                        tick_length=diff(xlimits)/100;
                        errorbar_tick(han,tick_length,'units');
                    end
                end
            end
        end
    end
    if swip==1
        % Velocity file name
        filevel=fullfile(dir_rep_ind,[num2str(XmidT(ix),xmidformat),extens,'.',...
            avertype,'.',modeltype]);
        % Standard deviation file name
        filestd=fullfile(dir_rep_ind,[num2str(XmidT(ix),xmidformat),extens,'.',...
            'VmsStd.',modeltype]);
        if strcmp(modeltype,'best')==1
            filestd=fullfile(dir_rep_ind,[num2str(XmidT(ix),xmidformat),extens,'.',...
                'VmsStd.layered']);
        end
        if strcmp(modeltype,'ridge')==1
            filestd=fullfile(dir_rep_ind,[num2str(XmidT(ix),xmidformat),extens,'.',...
                'VmsStd.smlay']);
        end
        
        D=[];
        if exist(filevel,'file')==0 || exist(filestd,'file')==0
            if npvc>0 && sum(nshot(ix,:))>=0
                fprintf(['\n  No ',[num2str(XmidT(ix),xmidformat),extens,'.',...
                    avertype,'.',modeltype],' SWIP file for this Xmid\n']);
            end
            vpsw=[]; vssw=[]; rhosw=[];
            vpstd=[]; vsstd=[]; rhostd=[];
        else
            modvel=dlmread(filevel,'',1,0);
            moddepth=[0;cumsum(modvel(:,1))];
            modstd=dlmread(filestd,'',1,0);
            depthstd=[0;modstd(:,1)];
            
            if maxdepth>moddepth(end)
                moddepth(end)=maxdepth;
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
            
            depthstd(end)=0;
            vpstd=modstd(:,2);
            vsstd=modstd(:,3);
            rhostd=modstd(:,4);
            
            if usevptomo==0
                filedisp=[filevel,'.disp'];
            else
                filevel=[filevel,'_vptomo'];
                filedisp=[filevel,'_vptomo.disp'];
                vptomo=VpI(VpI(:,ix)>0,ix);
                if isempty(vptomo)~=1
                    vptomo=[vptomo(1);vptomo;vptomo(end)];
                    ztmp=dz.*ones(1,size(VpI(VpI(:,ix)>0,ix),1));
                    zinc=[0,cumsum(ztmp)];
                    if max(zinc)>maxdepth
                        vptomo=vptomo(1:length(zinc)+1);
                    elseif max(zinc)<maxdepth
                        zinc=0:dz:maxdepth;
                        vptomo2=vptomo(end)*ones(length(zinc)+1,1);
                        vptomo2(1:length(vptomo))=vptomo;
                        vptomo=vptomo2;
                    end
                    nlay=size(moddepth,1)-1;
                    [vpsw,~,~,vssw]=velresamp(zinc,vptomo,moddepth,vssw,0.1,0,0);
                    dinsave(filevel,thick,vpsw,vssw,rhosw);
                else
                    if sum(nshot(ix,:))>=0
                        fprintf(['\n  No Vp from tomo for Xmid',num2str(ix),' = ',...
                            num2str(XmidT(ix),xmidformat),' m\n']);
                    end
                    vpsw=[]; vssw=[]; rhosw=[];
                    vpstd=[]; vsstd=[]; rhostd=[];
                end
            end
            if plot1dcal==1 && isempty(vssw)==0
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
                %                 if usevptomo==1
                %                     unix(['rm -rf ',filevel]);
                %                 end
                delete(filedisp);
                if isempty(D)==0
                    for m=1:nmodemax
                        freqcal=D{m,1}; % Frequency
                        vcal=1./D{m,2}; % Frequency
                        plot(freqcal,vcal,'-','Color',[1 0 0],...
                            'linewidth',1.5,'markersize',10);
                    end
                else
                    fprintf('\n  Failed to compute theoretical dispersion from SWIP model\n');
                end
            end
        end
    end
    if tomo==1 && isempty(VpItomo)==0 && isempty(VsItomo)==0
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
            vpstd=[]; vsstd=[];
        end
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
                    plot(freqcal,vcal,'-','Color',[0 0.75 0],...
                        'linewidth',1.5,'markersize',10);
                end
            else
                fprintf('\n  Failed to compute theoretical dispersion from tomo model\n');
            end
        end
    end
    if user>0
        D=[];
        if user==1
            filevel=fullfile(dir_dat,[num2str(XmidT(ix),xmidformat),'.user']);
            filedisp=fullfile(dir_dat,[num2str(XmidT(ix),xmidformat),'.user.disp']);
            dinsave(filevel,[thkuser,thkuser(end)],vpuser,vsuser,rhouser);
        else
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
                    plot(freqcal,vcal,'-','Color',[0 0 1],...
                        'linewidth',1.5,'markersize',10);
                end
            else
                fprintf('\n  Failed to compute theoretical dispersion from user model\n');
            end
        end
    end
    
    % Get DOI (Lmax*doifact)
    if plotDOI==1
        DOI(ix)=lmaxpick(ix)*doifact;
    end
    % Get DOI with standard deviation
    if plotDOI==2
        if exist('vsstd','var')==1 && isempty(vsstd)==0
            flipmoddepth=flipud(moddepth);
            flipvsstd=flipud([vsstd;vsstd(end)]);
            if sum(vsstd)~=0 % Case more than one model is in the error bars
                indhsd=find(flipvsstd<stdMAX,1,'first');
                if indhsd~=1 % Case there are some VsSTD < stdMAX
                    hsdtmp=flipmoddepth(indhsd-1);
                elseif indhsd==1 % Case all VsSTD < stdMAX
                    hsdtmp=0;
                else % Case all VsSTD > stdMAX
                    hsdtmp=[];
                end
            else % Only one model => VsSTD=0
                hsdtmp=lmaxpick(ix)*doifact; % Fix higher limit to lmaxpick(ix)*doifact
            end
            if isempty(hsdtmp)==1 % Case all VsSTD > stdMAX
                hsdtmp=moddepth(find(vsstd<=median(vsstd),1,'first'));
            elseif hsdtmp==0 % Case all VsSTD < stdMAX
                hsdtmp=flipmoddepth(find(flipvsstd~=flipvsstd(1),1,'first'));
                if isempty(hsdtmp)==1
                    hsdtmp=flipmoddepth(1);
                end
            end
            if hsdtmp>lmaxpick(ix)*doifact
                hsdtmp=lmaxpick(ix)*doifact; % Fix higher limit to lmaxpick(ix)*doifact
            end
            DOI(ix)=hsdtmp;
        else
            DOI(ix)=maxdepth;
        end
    end
    if plotDOI==3
        hsdtmp=moddepth(find(vsstd<=median(vsstd),1,'last'));
        DOI(ix)=hsdtmp;
    end
    
    %%% PLOT DISPERSION IMAGES WITH CALCULATED DISPERSION CURVES %%%
    
    if plot1dcal==1
        if (swip==1 && isempty(vssw)==0) || user>0 || (tomo==1 && isempty(VpItomo)==0 && isempty(VsItomo)==0)
            fprintf('\n  Saving calculated dispersion\n');
            file1=fullfile(dir_img_ind,[num2str(XmidT(ix),xmidformat),'.xdisp.fwd.',imgform]);
            save_fig(f1,file1,imgform,imgres,1);
            if showplot==0
                close(f1);
            else
                showplot=showplot+1;
            end
        end
    end
    
    %%
    %%% PLOTS 1D MODELS %%%
    
    if plot1dmod==1
        if (swip==1 && isempty(vssw)==0) || user>0 || (tomo==1 && isempty(VpItomo)==0 && isempty(VsItomo)==0)
            fprintf('\n  Saving 1D models\n');
        end
        %%% PLOTS 1D VS MODELS %%%
        
        if swip==1 && isempty(vssw)==0
            % Plot Vs
            [VSplot,Zplot]=stair2plot(vssw,moddepth);
            f1=plot_curv(showplot,VSplot,Zplot,[],'-',[1 0 0],[],1,1,...
                0,fs,'Vs (m/s)',depthtitle,[],[vsMIN vsMAX],[dpMIN dpMAX],[],...
                vsticks,dticks,[],[],[],[0 0 24 18],[],0);
            hold on
            if plot1dstd==1 % Plot VsSTD on same plot
                depthstdcalc=[0;diff(vssw)<0;0];
                depthstdcalc(depthstdcalc==0)=-1;
                if errstd>0
                    vsstd=vssw*errstd/100;
                    depthstdup=sort(moddepth+depthstdcalc.*moddepth*errstd/100);
                    depthstddown=sort(moddepth-depthstdcalc.*moddepth*errstd/100);
                else
                    depthstdup=sort(moddepth+depthstdcalc.*depthstd);
                    depthstddown=sort(moddepth-depthstdcalc.*depthstd);
                end
                depthstdup(end)=maxdepth;
                depthstddown(end)=maxdepth;
                [~,Zplotup]=stair2plot(vsstd,depthstdup);
                [VSstdplot,Zplotdown]=stair2plot(vsstd,depthstddown);
                
                dashline(VSplot+VSstdplot,Zplotup,2,2,2,2,'color',[1,0,0],'linewidth',1.5);
                dashline(VSplot-VSstdplot,Zplotdown,2,2,2,2,'color',[1,0,0],'linewidth',1.5);
            end
            if plot1dvp==1 % Plot Vp on same plot
                [VPplot,Zplot]=stair2plot(vpsw,moddepth);
                plot(VPplot,Zplot,'m-','linewidth',1.5);
                xlabel('V (m/s)');
                %                 xlim([vsMIN vpMAX]); set(gca,'XTick',vpticks);
                if plot1dstd==1 % Plot VpSTD on same plot
                    depthstdcalc=[0;diff(vpsw)<0;0];
                    depthstdcalc(depthstdcalc==0)=-1;
                    if errstd>0
                        vpstd=vpsw*errstd/100;
                        depthstdup=sort(moddepth+depthstdcalc.*moddepth*errstd/100);
                        depthstddown=sort(moddepth-depthstdcalc.*moddepth*errstd/100);
                    else
                        depthstdup=sort(moddepth+depthstdcalc.*depthstd);
                        depthstddown=sort(moddepth-depthstdcalc.*depthstd);
                    end
                    depthstdup(end)=maxdepth;
                    depthstddown(end)=maxdepth;
                    depthstdup(end)=maxdepth;
                    depthstddown(end)=maxdepth;
                    [~,Zplotup]=stair2plot(vpstd,depthstdup);
                    [VPstdplot,Zplotdown]=stair2plot(vpstd,depthstddown);
                    
                    dashline(VPplot+VPstdplot,Zplotup,2,2,2,2,'color','m','linewidth',1.5);
                    dashline(VPplot-VPstdplot,Zplotdown,2,2,2,2,'color','m','linewidth',1.5);
                end
            end
        end
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
                plot(VSplotuser,Zplotuser,'b-','linewidth',1.5);
            else
                f1=plot_curv(showplot,VSplotuser,Zplotuser,[],'-',[0 0 1],[],1,1,...
                    0,fs,'Vs (m/s)',depthtitle,[],[vsMIN vsMAX],[dpMIN dpMAX],[],...
                    vsticks,dticks,[],[],[],[0 0 24 18],[],0);
                hold on
            end
            if plot1dstd==1 && errstd>0 % Plot VsSTD on same plot
                depthstdcalc=[0;diff(vsuserok)<0;0];
                depthstdcalc(depthstdcalc==0)=-1;
                vsstduser=vsuserok*errstd/100;
                depthstdup=sort(moddepthuser+depthstdcalc.*moddepthuser*errstd/100);
                depthstddown=sort(moddepthuser-depthstdcalc.*moddepthuser*errstd/100);
                depthstdup(end)=maxdepth;
                depthstddown(end)=maxdepth;
                [~,Zplotup]=stair2plot(vsstduser,depthstdup);
                [VSstduserplot,Zplotdown]=stair2plot(vsstduser,depthstddown);
                
                dashline(VSplotuser+VSstduserplot,Zplotup,2,2,2,2,'color','b','linewidth',1.5);
                dashline(VSplotuser-VSstduserplot,Zplotdown,2,2,2,2,'color','b','linewidth',1.5);
            end
            if plot1dvp==1 % Plot Vp on same plot
                [VPplotuser,Zplotuser]=stair2plot(vpuserok,moddepthuser);
                plot(VPplotuser,Zplotuser,'c-','linewidth',1.5);
                xlabel('V (m/s)');
                %                 xlim([vsMIN vpMAX]); set(gca,'XTick',vpticks);
                if plot1dstd==1 && errstd>0 % Plot VpSTD on same plot
                    depthstdcalc=[0;diff(vpuserok)<0;0];
                    depthstdcalc(depthstdcalc==0)=-1;
                    vpstd=vpuserok*errstd/100;
                    depthstdup=sort(moddepthuser+depthstdcalc.*moddepthuser*errstd/100);
                    depthstddown=sort(moddepthuser-depthstdcalc.*moddepthuser*errstd/100);
                    depthstdup(end)=maxdepth;
                    depthstddown(end)=maxdepth;
                    [~,Zplotup]=stair2plot(vpstd,depthstdup);
                    [VPstduserplot,Zplotdown]=stair2plot(vpstd,depthstddown);
                    
                    dashline(VPplotuser+VPstduserplot,Zplotup,2,2,2,2,'color','c','linewidth',1.5);
                    dashline(VPplotuser-VPstduserplot,Zplotdown,2,2,2,2,'color','c','linewidth',1.5);
                end
            end
        end
        if tomo==1 && isempty(VpItomo)==0 && isempty(VsItomo)==0
            % Plot Vs
            if (swip==1 && isempty(vssw)==0) || user==1
                plot(vstomo,cumsum(ztomo),'-','color',[0 0.75 0],'linewidth',1.5);
            else
                f1=plot_curv(showplot,vstomo,cumsum(ztomo),[],'-',[0 0.75 0],[],1,1,...
                    0,fs,'Vs (m/s)',depthtitle,[],[vsMIN vsMAX],[dpMIN dpMAX],[],...
                    vsticks,dticks,[],[],[],[0 0 24 18],[],0);
                hold on
            end
            if plot1dstd==1 % Plot VsSTD on same plot
                if errstd>0
                    dashline(vstomo+vstomo*errstd/100,cumsum(ztomo),2,2,2,2,'color',[0 0.75 0],'linewidth',1.5);
                    dashline(vstomo-vstomo*errstd/100,cumsum(ztomo),2,2,2,2,'color',[0 0.75 0],'linewidth',1.5);
                else
                    
                end
            end
            if plot1dvp==1 % Plot Vp on same plot
                plot(vptomo,cumsum(ztomo),'g-','linewidth',1.5);
                xlabel('V (m/s)');
                %                 xlim([vsMIN vpMAX]); set(gca,'XTick',vpticks);
                if plot1dstd==1 % Plot VpSTD on same plot
                    if errstd>0
                        dashline(vptomo+vptomo*errstd/100,cumsum(ztomo),2,2,2,2,'color','g','linewidth',1.5);
                        dashline(vptomo-vptomo*errstd/100,cumsum(ztomo),2,2,2,2,'color','g','linewidth',1.5);
                    else
                        
                    end
                end
            end
        end
        if (swip==1 && isempty(vssw)==0) || user>0 || (tomo==1 && isempty(VpItomo)==0 && isempty(VsItomo)==0)
            %             fprintf('\n  Saving Vs 1D model\n');
            % Plot DOI
            if plotDOI>0
                xL=get(gca,'XLim');
                han3=dashline(xL,[DOI(ix) DOI(ix)],2,2,2,2,'color','k','linewidth',1.5);
            end
            sizeax=get(gca,'Position');
            file_vs=fullfile(dir_img_ind,[num2str(XmidT(ix),xmidformat),'.xVs1D.',imgform]);
            save_fig(f1,file_vs,imgform,imgres,1,0);
            if showplot==0
                close(f1);
            else
                showplot=showplot+1;
            end
        end
        
        %%
        %%% PLOTS 1D VP MODELS %%%
        
        if swip==1 && isempty(vssw)==0
            % Plot Vp
            [VPplot,Zplot]=stair2plot(vpsw,moddepth);
            f1=plot_curv(showplot,VPplot,Zplot,[],'-',[1 0 0],[],1,1,...
                0,fs,'Vp (m/s)',depthtitle,[],[vpMIN vpMAX],[dpMIN dpMAX],[],...
                vpticks,dticks,[],[],[],[0 0 24 18],sizeax,0);
            hold on
            if plot1dstd==1 % Plot VpSTD on same plot
                depthstdcalc=[0;diff(vpsw)<0;0];
                depthstdcalc(depthstdcalc==0)=-1;
                if errstd>0
                    vpstd=vpsw*errstd/100;
                    depthstdup=sort(moddepth+depthstdcalc.*moddepth*errstd/100);
                    depthstddown=sort(moddepth-depthstdcalc.*moddepth*errstd/100);
                else
                    depthstdup=sort(moddepth+depthstdcalc.*depthstd);
                    depthstddown=sort(moddepth-depthstdcalc.*depthstd);
                end
                depthstdup(end)=maxdepth;
                depthstddown(end)=maxdepth;
                [~,Zplotup]=stair2plot(vpstd,depthstdup);
                [VPstdplot,Zplotdown]=stair2plot(vpstd,depthstddown);
                
                dashline(VPplot+VPstdplot,Zplotup,2,2,2,2,'color',[1,0,0],'linewidth',1.5);
                dashline(VPplot-VPstdplot,Zplotdown,2,2,2,2,'color',[1,0,0],'linewidth',1.5);
            end
        end
        if user>0
            % Plot Vp
            [VPplotuser,Zplotuser]=stair2plot(vpuserok,moddepthuser);
            if swip==1 && isempty(vssw)==0
                plot(VPplotuser,Zplotuser,'b-','linewidth',1.5);
            else
                f1=plot_curv(showplot,VPplotuser,Zplotuser,[],'-',[0 0 1],[],1,1,...
                    0,fs,'Vp (m/s)',depthtitle,[],[vpMIN vpMAX],[dpMIN dpMAX],[],...
                    vpticks,dticks,[],[],[],[0 0 24 18],sizeax,0);
                hold on
            end
            if plot1dstd==1 && errstd>0 % Plot VsSTD on same plot
                depthstdcalc=[0;diff(vpuserok)<0;0];
                depthstdcalc(depthstdcalc==0)=-1;
                vpstd=vpuserok*errstd/100;
                depthstdup=sort(moddepthuser+depthstdcalc.*moddepthuser*errstd/100);
                depthstddown=sort(moddepthuser-depthstdcalc.*moddepthuser*errstd/100);
                depthstdup(end)=maxdepth;
                depthstddown(end)=maxdepth;
                [~,Zplotup]=stair2plot(vpstd,depthstdup);
                [VPstduserplot,Zplotdown]=stair2plot(vpstd,depthstddown);
                
                dashline(VPplotuser+VPstduserplot,Zplotup,2,2,2,2,'color','b','linewidth',1.5);
                dashline(VPplotuser-VPstduserplot,Zplotdown,2,2,2,2,'color','b','linewidth',1.5);
            end
        end
        if tomo==1 && isempty(VpItomo)==0 && isempty(VsItomo)==0
            % Plot Vp
            if (swip==1 && isempty(vssw)==0) || user>0
                plot(vptomo,cumsum(ztomo),'-','color',[0 0.75 0],'linewidth',1.5);
            else
                f1=plot_curv(showplot,vptomo,cumsum(ztomo),[],'-',[0 0.75 0],[],1,1,...
                    0,fs,'Vp (m/s)',depthtitle,[],[vpMIN vpMAX],[dpMIN dpMAX],[],...
                    vpticks,dticks,[],[],[],[0 0 24 18],sizeax,0);
                hold on
            end
            if plot1dstd==1 % Plot VpSTD on same plot
                if errstd>0
                    plot(vptomo+vptomo*errstd/100,cumsum(ztomo),'--','color',...
                        [0 0.75 0],'linewidth',1.5);
                    plot(vptomo-vptomo*errstd/100,cumsum(ztomo),'--','color',...
                        [0 0.75 0],'linewidth',1.5);
                else
                    
                end
            end
        end
        if (swip==1 && isempty(vssw)==0) || user>0 || (tomo==1 && isempty(VpItomo)==0 && isempty(VsItomo)==0)
            %             fprintf('\n  Saving Vp 1D model\n');
            % Plot DOI
            if plotDOI>0
                xL=get(gca,'XLim');
                han3=dashline(xL,[DOI(ix) DOI(ix)],2,2,2,2,'color','k','linewidth',1.5);
            end
            file_vp=fullfile(dir_img_ind,[num2str(XmidT(ix),xmidformat),'.xVp1D.',imgform]);
            save_fig(f1,file_vp,imgform,imgres,1,0);
            if showplot==0
                close(f1);
            else
                showplot=showplot+1;
            end
        end
        
        %%
        %%% PLOTS 1D VP/VS MODELS %%%
        
        if swip==1 && isempty(vssw)==0
            % Plot Vp/Vs
            f1=plot_curv(showplot,VPplot./VSplot,Zplot,[],'-',[1 0 0],[],1,1,...
                0,fs,'Vp/Vs',depthtitle,[],[vpvsMIN vpvsMAX],[dpMIN dpMAX],[],...
                vpvsticks,dticks,[],[],[],[0 0 24 18],sizeax,0);
            hold on
            %         if plot1dstd==1 % Plot VsSTD on same plot
            %             depthstdcalc=[0;diff(vssw)<0;0];
            %             depthstdcalc(depthstdcalc==0)=-1;
            %             if errstd>0
            %                 vsstd=vssw*errstd/100;
            %                 depthstdup=sort(moddepth+depthstdcalc.*moddepth*errstd/100);
            %                 depthstddown=sort(moddepth-depthstdcalc.*moddepth*errstd/100);
            %             else
            %                 depthstdup=sort(moddepth+depthstdcalc.*depthstd);
            %                 depthstddown=sort(moddepth-depthstdcalc.*depthstd);
            %             end
            %             depthstdup(end)=maxdepth;
            %             depthstddown(end)=maxdepth;
            %             stairs(depthstdup,[vssw;vssw(end)]+[vsstd;vsstd(end)],...
            %                 'r--','linewidth',1.5);
            %             stairs(depthstddown,[vssw;vssw(end)]-[vsstd;vsstd(end)],...
            %                 'r--','linewidth',1.5);
            %         end
        end
        if user>0
            % Plot Vp/Vs
            if swip==1 && isempty(vssw)==0
                plot(VPplotuser./VSplotuser,Zplotuser,'b-','linewidth',1.5);
            else
                f1=plot_curv(showplot,VPplotuser./VSplotuser,Zplotuser,[],'-',[0 0 1],[],1,1,...
                    0,fs,'Vp/Vs',depthtitle,[],[vpvsMIN vpvsMAX],[dpMIN dpMAX],[],...
                    vpvsticks,dticks,[],[],[],[0 0 24 18],sizeax,0);
                hold on
            end
            %             if plot1dstd==1 && errstd>0 % Plot VsSTD on same plot
            %                 depthstdcalc=[0;diff(vpuserok)<0;0];
            %                 depthstdcalc(depthstdcalc==0)=-1;
            %                 vpstd=vpuserok*errstd/100;
            %                 depthstdup=sort(moddepthuser+depthstdcalc.*moddepthuser*errstd/100);
            %                 depthstddown=sort(moddepthuser-depthstdcalc.*moddepthuser*errstd/100);
            %                 depthstdup(end)=maxdepth;
            %                 depthstddown(end)=maxdepth;
            %                 stairs(depthstdup,[vpuserok;vpuserok(end)]+[vpstd;vpstd(end)],...
            %                     'b--','linewidth',1.5);
            %                 stairs(depthstddown,[vpuserok;vpuserok(end)]-[vpstd;vpstd(end)],...
            %                     'b--','linewidth',1.5);
            %             end
        end
        if tomo==1 && isempty(VpItomo)==0 && isempty(VsItomo)==0
            % Plot Vp/Vs
            if (swip==1 && isempty(vssw)==0) || user>0
                plot(vptomo./vstomo,cumsum(ztomo),'-','color',[0 0.75 0],'linewidth',1.5);
            else
                f1=plot_curv(showplot,vptomo./vstomo,cumsum(ztomo),[],'-',[0 0.75 0],[],1,1,...
                    0,fs,'Vp/Vs',depthtitle,[],[vpvsMIN vpvsMAX],[dpMIN dpMAX],[],...
                    vpvsticks,dticks,[],[],[],[0 0 24 18],sizeax,0);
                hold on
            end
            %             if plot1dstd==1 % Plot VpSTD on same plot
            %                 if errstd>0
            %                     plot(cumsum(ztomo),vptomo+vptomo*errstd/100,'b--','linewidth',1.5);
            %                     plot(cumsum(ztomo),vptomo-vptomo*errstd/100,'b--','linewidth',1.5);
            %                 else
            %
            %                 end
            %             end
        end
        if (swip==1 && isempty(vssw)==0) || user>0 || (tomo==1 && isempty(VpItomo)==0 && isempty(VsItomo)==0)
            %             fprintf('\n  Saving Vp/Vs 1D model\n');
            % Plot DOI
            if plotDOI>0
                xL=get(gca,'XLim');
                han3=dashline(xL,[DOI(ix) DOI(ix)],2,2,2,2,'color','k','linewidth',1.5);
            end
            file_vpvs=fullfile(dir_img_ind,[num2str(XmidT(ix),xmidformat),'.xVpVs1D.',imgform]);
            save_fig(f1,file_vpvs,imgform,imgres,1,0);
            if showplot==0
                close(f1);
            else
                showplot=showplot+1;
            end
        end
        
        %%
        %%% PLOTS 1D POISSON'S RATIO MODELS %%%
        
        if swip==1 && isempty(vssw)==0
            % Plot Poisson's ratio
            f1=plot_curv(showplot,poisson(VPplot,VSplot),Zplot,[],'-',[1 0 0],[],1,1,...
                0,fs,'Poisson',depthtitle,[],[poisMIN poisMAX],[dpMIN dpMAX],[],...
                poisticks,dticks,[],[],[],[0 0 24 18],sizeax,0);
            hold on
            %         if plot1dstd==1 % Plot VsSTD on same plot
            %             depthstdcalc=[0;diff(vssw)<0;0];
            %             depthstdcalc(depthstdcalc==0)=-1;
            %             if errstd>0
            %                 vsstd=vssw*errstd/100;
            %                 depthstdup=sort(moddepth+depthstdcalc.*moddepth*errstd/100);
            %                 depthstddown=sort(moddepth-depthstdcalc.*moddepth*errstd/100);
            %             else
            %                 depthstdup=sort(moddepth+depthstdcalc.*depthstd);
            %                 depthstddown=sort(moddepth-depthstdcalc.*depthstd);
            %             end
            %             depthstdup(end)=maxdepth;
            %             depthstddown(end)=maxdepth;
            %             stairs(depthstdup,[vssw;vssw(end)]+[vsstd;vsstd(end)],...
            %                 'r--','linewidth',1.5);
            %             stairs(depthstddown,[vssw;vssw(end)]-[vsstd;vsstd(end)],...
            %                 'r--','linewidth',1.5);
            %         end
        end
        if user>0
            % Plot Poisson's ratio
            if swip==1 && isempty(vssw)==0
                plot(poisson(VPplotuser,VSplotuser),Zplotuser,'b-','linewidth',1.5);
            else
                f1=plot_curv(showplot,poisson(VPplotuser,VSplotuser),Zplotuser,[],'-',[0 0 1],[],1,1,...
                    0,fs,'Poisson',depthtitle,[],[poisMIN poisMAX],[dpMIN dpMAX],[],...
                    poisticks,dticks,[],[],[],[0 0 24 18],sizeax,0);
                hold on
            end
            %             if plot1dstd==1 && errstd>0 % Plot VsSTD on same plot
            %                 depthstdcalc=[0;diff(vpuserok)<0;0];
            %                 depthstdcalc(depthstdcalc==0)=-1;
            %                 vpstd=vpuserok*errstd/100;
            %                 depthstdup=sort(moddepthuser+depthstdcalc.*moddepthuser*errstd/100);
            %                 depthstddown=sort(moddepthuser-depthstdcalc.*moddepthuser*errstd/100);
            %                 depthstdup(end)=maxdepth;
            %                 depthstddown(end)=maxdepth;
            %                 stairs(depthstdup,[vpuserok;vpuserok(end)]+[vpstd;vpstd(end)],...
            %                     'b--','linewidth',1.5);
            %                 stairs(depthstddown,[vpuserok;vpuserok(end)]-[vpstd;vpstd(end)],...
            %                     'b--','linewidth',1.5);
            %             end
        end
        if tomo==1 && isempty(VpItomo)==0 && isempty(VsItomo)==0
            % Plot Poisson's ratio
            if (swip==1 && isempty(vssw)==0) || user>0
                plot(poisson(vptomo,vstomo),cumsum(ztomo),'-','color',[0 0.75 0],'linewidth',1.5);
            else
                f1=plot_curv(showplot,poisson(vptomo,vstomo),cumsum(ztomo),[],'-',[0 0.75 0],[],1,1,...
                    0,fs,'Poisson',depthtitle,[],[poisMIN poisMAX],[dpMIN dpMAX],[],...
                    poisticks,dticks,[],[],[],[0 0 24 18],sizeax,0);
                hold on
            end
            %             if plot1dstd==1 % Plot VpSTD on same plot
            %                 if errstd>0
            %                     plot(cumsum(ztomo),vptomo+vptomo*errstd/100,'b--','linewidth',1.5);
            %                     plot(cumsum(ztomo),vptomo-vptomo*errstd/100,'b--','linewidth',1.5);
            %                 else
            %
            %                 end
            %             end
        end
        if (swip==1 && isempty(vssw)==0) || user>0 || (tomo==1 && isempty(VpItomo)==0 && isempty(VsItomo)==0)
            %             fprintf('\n  Saving Poisson''s ratio 1D model\n');
            % Plot DOI
            if plotDOI>0
                xL=get(gca,'XLim');
                han3=dashline(xL,[DOI(ix) DOI(ix)],2,2,2,2,'color','k','linewidth',1.5);
            end
            file_pois=fullfile(dir_img_ind,[num2str(XmidT(ix),xmidformat),'.xPoisson1D.',imgform]);
            save_fig(f1,file_pois,imgform,imgres,1,0);
            if showplot==0
                close(f1);
            else
                showplot=showplot+1;
            end
            if testplot==1
                filename_imgtmp=fullfile(dir_img_ind,[num2str(XmidT(ix),xmidformat),...
                    '.mod1d.tmp.',imgform]);
                filename_imgtmp2=fullfile(dir_img_ind,[num2str(XmidT(ix),xmidformat),...
                    '.mod1d.tmp2.',imgform]);
                filename_panel=fullfile(dir_img_ind,[num2str(XmidT(ix),xmidformat),...
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
