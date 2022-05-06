function quick_invert(inversion,plotinvres,plotparam,targetfile,paramfile,nrun,itmax,ns0,ns,nr,...
    nbest,outpoints,modeltype,imgform,imgres,concat,colnb,verbose)
%%% S. Pasquet - V17.03.29
%%% Quick inversion of dispersion data
%%% quick_invert(inversion,plotinvres,plotparam,targetfile,paramfile,nrun,itmax,ns0,ns,nr,...
%%%    calcmod,nbest,outpoints,modeltype,imgform,imgres,concat,colnb,verbose)

run('SWIP_defaultsettings')

% File name extensions
if nbest==0
    extens=['.bweb',num2str(outpoints)]; % Best within error bars
else
    extens=['.best',num2str(nbest)]; % Arbitrary nb
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

[testimgmgck,~]=unix('which montage');
[testpdfjam,~]=unix('which pdfjam');
testplot=((testpdfjam==0 && strcmp(imgform,'pdf')==1) || (testimgmgck==0 && strcmp(imgform,'pdf')==0 && strcmp(imgform,'fig')==0));
if concat == 0
    testplot = 0;
end

if isunix==1
    zipmethod=5;
else
    zipmethod=1;
end
tstart=tic;

%%
%%%%%% Inversion %%%%%%
if inversion==1
    if exist('targetfile','var')==0 || isempty(targetfile)==1
        [targetfile,targetpath]=uigetfile('*.target','Select target file');
        if targetfile==0
            fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
            fprintf('\n   Please select a target file');
            fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
            return
        end
    end
    nametarg=fullfile(targetpath,targetfile);
    
    if exist('paramfile','var')==0 || isempty(paramfile)==1
        [paramfile,parampath]=uigetfile('*.param','Select parameterization file');
        if paramfile==0
            fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
            fprintf('\n   Please select a parameterization file');
            fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
            return
        end
    end
    paramname=fullfile(parampath,paramfile);
    
    % Create report folder
    dir_rep_ind=fullfile(pwd,[targetfile(1:end-6),'report']);
    if exist(dir_rep_ind,'dir')~=7
        mkdir(dir_rep_ind);
    end
    
    % Run inversion
    status=matdinver(nametarg,paramname,nrun,itmax,ns0,ns,nr,dir_rep_ind,verbose);
    if status~=0
        return
    end
    copyfile(nametarg,dir_rep_ind);
    copyfile(paramname,dir_rep_ind);
    matzip(1,fullfile(dir_rep_ind,'*.report'),zipmethod,1);
else
    dir_rep_ind=uigetdir('./','Select inversion folder');
    if dir_rep_ind==0
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        fprintf('\n   Please select inversion folder');
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
        return
    else
        targstruct=dir(fullfile(dir_rep_ind,'*.target'));
        nametarg=fullfile(dir_rep_ind,targstruct.name);
        paramstruct=dir(fullfile(dir_rep_ind,'*.param'));
        paramname=fullfile(dir_rep_ind,paramstruct.name);
    end
end

if isempty(dpMAX)==1
    dpMIN = 0;
    [~,dpMAX]=param2mod(paramname);
end

%%

if calcmod==1
    
    %%
    %%%%%% Looking for the best models %%%%%%
    [~,targetfile,extension]=fileparts(nametarg);
    targetfile=[targetfile,extension];
    % Read target file to get picked dispersion curves
    [freqresamp,vresamp,deltaresamp,modes]=targ2pvc(nametarg);
    maxmodeinv=max(modes);
    nmaxmod=(itmax*ns)+ns0;
    nmodeinv=length(modes);
    for ii=1:nmodeinv
        lm(ii)=max(vresamp{ii}./freqresamp{ii});
    end
    lmaxpick=max(lm);
    
    if isempty(dpMAX)==1
        maxdepth=lmaxpick;
    else
        maxdepth=dpMAX;
    end
    
    fprintf('\n  Looking for models within the errorbars\n');
    % Export calculated dipsersion curves from binary report files
    if exist(fullfile(dir_rep_ind,'run_01.report'),'file')==2
        matzip(1,fullfile(dir_rep_ind,'*.report'),zipmethod,1);
    end
    matzip(0,fullfile(dir_rep_ind,'*.report.gz'),zipmethod,0);
    matgpdcreport(dir_rep_ind,0,nrun,maxmodeinv+1,nmaxmod,wave);
    if nbest>nmaxmod
        nbest=nmaxmod;
    end
    D=cell(nmodeinv,1);
    ok=cell(nmodeinv,1);
    bestmis=zeros(nrun,nmodeinv);
    misround=zeros(nrun*nmaxmod,nmodeinv);
    for ip=1:nmodeinv
        m=modes(ip);
        fprintf(['\n      Mode ',num2str(m),'\n']);
        minvelOK=vresamp{m+1}(isnan(vresamp{m+1})==0)-deltaresamp{m+1}(isnan(vresamp{m+1})==0);
        maxvelOK=vresamp{m+1}(isnan(vresamp{m+1})==0)+deltaresamp{m+1}(isnan(vresamp{m+1})==0);
        freqOK=freqresamp{m+1}(isnan(vresamp{m+1})==0);
        D{ip}=cell(nmaxmod*nrun,2); % Frequency vs slowness cell array
        dispselec=repmat(struct('modnum',[],'modok',[],...
            'misround',[],'nfreqsample',[]),1,nrun); % Models structure
        for n=1:nrun % Loop over all runs
            calcdisp=fullfile(dir_rep_ind,['best',num2str(n),'.M',num2str(m),'.txt']);
            [D{ip}(1+(n-1)*nmaxmod:(n-1)*nmaxmod+nmaxmod,:),dispselec(n)]=...
                readdisp(calcdisp,nmaxmod,n,nbest,outpoints,freqOK,minvelOK,maxvelOK); % Read file
            bestmis(n,ip)=dispselec(n).misround(1); % Get best model misfit
        end
        names=fieldnames(dispselec); % Get field names
        cellData=cellfun(@(f){vertcat(dispselec.(f))},names); % Collect field data into a cell array
        dispselec=cell2struct(cellData,names); % Convert the cell array into a structure
        ok{ip}=dispselec.modok; % Get selected models
        misround(:,ip)=dispselec.misround; % Get model misfits
    end
    delete(fullfile(dir_rep_ind,'best*.txt'));
    % Find best models indexes sorted by decreasing misfit
    bestrun=find(bestmis(:,1)==min(bestmis(:,1))); % Get best run nb to get best model
    mis=max(misround,[],2);
    [missort,I]=sortrows(mis,-1);
    if nbest==0
        test=sum(cell2mat(ok'),2);
        OKsort=test==nmodeinv;
        KOsort=test<nmodeinv;
        OKsort=OKsort(I);
        KOsort=KOsort(I);
    else
        OKsort=false(size(missort));
        OKsort(end-nbest+1:end)=true;
        KOsort=false(size(missort));
        KOsort(1:end-nbest)=true;
    end
    misin=missort(OKsort);
    misout=missort(KOsort);
    % Get nb of accepted models and misfit for the current Xmid
    nmod=size(misin,1);
    minmis=min(mis);
    fprintf(['\n  Minimum misfit = ',num2str(min(missort)),'\n']);
    if nbest==0
        fprintf(['\n  ',num2str(nmod),' models fit within the error bars\n']);
    else
        fprintf(['\n  The best ',num2str(nmod),' models have been selected\n']);
    end
    if nmod==0
        fprintf('\n  **********************************************************');
        fprintf('\n  **********************************************************\n');
    end
    
    %%
    %%%%%% Extraction of all models %%%%%%
    
    fprintf('\n  Reading all calculated models\n');
    % Export calculated models from binary report files
    matgpdcreport(dir_rep_ind,1,nrun,nmodeinv,nmaxmod,wave);
    VSall=cell(nmaxmod*nrun,2); % Depth vs Vs cell array
    VPall=cell(nmaxmod*nrun,2); % Depth vs Vs cell array
    RHOall=cell(nmaxmod*nrun,2); % Depth vs Vs cell array
    % Loop over all runs
    for n=1:nrun
        fprintf(['\n      Run ',num2str(n),'\n']);
        vsfile=fullfile(dir_rep_ind,['vs',num2str(n),'.txt']);
        vpfile=fullfile(dir_rep_ind,['vp',num2str(n),'.txt']);
        rhofile=fullfile(dir_rep_ind,['rho',num2str(n),'.txt']);
        [VSall(1+(n-1)*nmaxmod:(n-1)*nmaxmod+nmaxmod,:)]=...
            readmodel(vsfile,nmaxmod,n); % Read file
        [VPall(1+(n-1)*nmaxmod:(n-1)*nmaxmod+nmaxmod,:)]=...
            readmodel(vpfile,nmaxmod,n); % Read file
        [RHOall(1+(n-1)*nmaxmod:(n-1)*nmaxmod+nmaxmod,:)]=...
            readmodel(rhofile,nmaxmod,n); % Read file
    end
    delete(fullfile(dir_rep_ind,'*.txt'),fullfile(dir_rep_ind,'*.report'));
    % Read models
    dim=ndims(VSall{1});
    VS=cat(dim,VSall{:,1});
    VSbest=VS(:,(nmaxmod*(bestrun-1))+1);
    VS=VS(:,I);
    VP=cat(dim,VPall{:,1});
    VPbest=VP(:,(nmaxmod*(bestrun-1))+1);
    VP=VP(:,I);
    RHO=cat(dim,RHOall{:,1});
    RHObest=RHO(:,(nmaxmod*(bestrun-1))+1);
    RHO=RHO(:,I);
    ZALL=cat(dim,VSall{:,2});
    ZALL(end,:)=maxdepth;
    Zbest=ZALL(:,(nmaxmod*(bestrun-1))+1);
    THbest=diff(Zbest);
    THbest=round(THbest(1:2:end)*1e6)/1e6;
    ZALL=ZALL(:,I);
    THALL=diff(ZALL);
    THALL=round(THALL(1:2:end,:)*1e6)/1e6;
    
    nC=length(THALL(:,1)); % nb of layers
    ZZ=0:dz:maxdepth;
    nZ=length(ZZ);
    
    % Get accepted and rejected models
    VSin=VS(:,OKsort);
    VPin=VP(:,OKsort);
    VPout=VP(:,KOsort);
    RHOin=RHO(:,OKsort);
    RHOout=RHO(:,KOsort);
    Zin=ZALL(:,OKsort);
    THin=THALL(:,OKsort);
    VSout=VS(:,KOsort);
    Zout=ZALL(:,KOsort);
    
    %%
    %%%%%% Building average models %%%%%%
    
    if calcmod==1
        fprintf('\n  Building average models\n');
        
        %%%%%% Mean layered model %%%%%%
        
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
        
        %%%%%% Mean smoothed model %%%%%%
        
        % !!!!(to be optimized for constant density)!!!!
        Zi=Zin;
        Zi(2:2:end-2,:)=Zi(2:2:end-2,:)-dz;
        IVsR=zeros(nmod,nZ);
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
        
        %%%%%% Ridge model %%%%%%
        
        if ridgecalc==1
            % Velocity and density vectors for ridge search
            minVs=floor(min(min(IVs))/10)*10;
            maxVs=ceil(max(max(IVs))/10)*10;
            velocityS=minVs:25:maxVs;
            minVp=floor(min(min(IVp))/10)*10;
            maxVp=ceil(max(max(IVp))/10)*10;
            velocityP=minVp:25:maxVp;
            minrho=floor(min(min(Irho))/10)*10;
            maxrho=ceil(max(max(Irho))/10)*10;
            density=minrho:25:maxrho;
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
            % Count number of model per cell
            for jj=1:nZ
                NN(jj,:)=histc(IVsR(:,jj),velocityS);
                VSridge(jj)=velocityS(find(NN(jj,:)==max(NN(jj,:)),1,'first'));
                NNP(jj,:)=histc(IVpR(:,jj),velocityP);
                VPridge(jj)=velocityP(find(NNP(jj,:)==max(NNP(jj,:)),1,'first'));
                NNr(jj,:)=histc(IrhoR(:,jj),density);
                RHOridge(jj)=density(find(NNr(jj,:)==max(NNr(jj,:)),1,'first'));
            end
            NN(NN==0)=NaN;
            % Plot ridgesearch results
            if nmod>1
                [f2,~,~,~,c]=plot_img(showplot,velocityS,ZZ,NN,flipud(autumn),...
                    axetop,1,1,fs,'Vs (m/s)','Depth (m)','Number of models',...
                    [vsMIN vsMAX],[dpMIN dpMAX],[],vsticks,dticks,[],[],[],[],[0 0 24 18],[],0);
                hold on
                dashline(VSridge,ZZ,2,2,2,2,'color','k','linewidth',1.5);
                cpos=get(c,'position');
                set(c,'Position',[cpos(1)+0.04 cpos(2) cpos(3)*1.5 cpos(4)]);
                h=findall(gcf,'Type','Axes'); set(h,'FontSize',fs);
                % Save figure
                file_ridge=fullfile(dir_rep_ind,[targetfile(1:end-7),...
                    '.mod1d.ridge.',imgform]);
                save_fig(f2,file_ridge,imgform,imgres,1,1-testplot);
                close(f2); figHandles = findall(0,'Type','figure');
                
                matrelease=version('-release');
                if str2double(matrelease(1:4))>2014
                    figHandles=get(figHandles,'Number');
                end
                % Strange bug for first figure
                if ismember(1,figHandles)==1
                    close(1);
                end
            end
        end
        
        %%%%%% Mean layered smoothed model %%%%%%
        
        %         Zmeani=Zmean;
        %         Zmeani(2:2:end-2,:)=Zmeani(2:2:end-2,:)-dz;
        %         VSmeanI=interp1q(Zmeani,VSmean,ZZ');
        %         VPmeanI=interp1q(Zmeani,VPmean,ZZ');
        %         RHOmeanI=interp1q(Zmeani,RHOmean,ZZ');
        %         % Standard deviation layered smoothed model
        %         Zstdi=Zstd;
        %         Zstdi(2:2:end-2,:)=Zstdi(2:2:end-2,:)-dz;
        %         VSstdI=interp1q(Zmeani,VSstd,ZZ');
        %         VPstdI=interp1q(Zmeani,VPstd,ZZ');
        %         RHOstdI=interp1q(Zmeani,RHOstd,ZZ');
        
        %%%%%% Gaussian weighting models %%%%%%
        
        if weightcalc==1 && nmod>1
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
                if step(flagg)>step(flagg-1)
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
            
            %%%%%% Weighted layered model %%%%%%
            
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
            
            %%%%%% Weighted smoothed model %%%%%%
            
            IVSweight=sum(repmat(COEF,1,size(IVs,2)).*IVs,1); % Weighted
            IVPweight=sum(repmat(COEF,1,size(IVp,2)).*IVp,1); % Weighted
            IRHOweight=sum(repmat(COEF,1,size(Irho,2)).*Irho,1); % Weighted
            
            %%%%%% Weighted layered smoothed model %%%%%%
            
            %             Zweighti=Zweight;
            %             Zweighti(2:2:end-2,:)=Zweighti(2:2:end-2,:)-dz;
            %             VSweightI=interp1q(Zweighti,VSweight,ZZ');
            %             VPweightI=interp1q(Zweighti,VPweight,ZZ');
            %             RHOweightI=interp1q(Zweighti,RHOweight,ZZ');
        end
        
        %%%%%% Save velocity models %%%%%%
        
        % Velocity model filenames
        % Best model
        nameDINbest=fullfile(dir_rep_ind,[targetfile(1:end-7),...
            extens,'.Vms.best']); % Layered model
        % Ridge
        nameDINridge=fullfile(dir_rep_ind,[targetfile(1:end-7),...
            extens,'.Vms.ridge']); % Layered model
        % Mean
        nameDINlay=fullfile(dir_rep_ind,[targetfile(1:end-7),...
            extens,'.Vms.layered']); % Layered model
        nameDINsm=fullfile(dir_rep_ind,[targetfile(1:end-7),...
            extens,'.Vms.smooth']); % Smooth model
        %         nameDINsmlay=fullfile(dir_rep_ind,[targetfile(1:end-7),...
        %             extens,'.Vms.smlay']); % Smooth layered model
        % Std
        nameDINlaystd=fullfile(dir_rep_ind,[targetfile(1:end-7),...
            extens,'.VmsStd.layered']); % Layered model
        nameDINsmstd=fullfile(dir_rep_ind,[targetfile(1:end-7),...
            extens,'.VmsStd.smooth']); % Smooth model
        nameDINsmlaystd=fullfile(dir_rep_ind,[targetfile(1:end-7),...
            extens,'.VmsStd.smlay']); % Smooth layered model
        % Weighted
        nameDINlayw=fullfile(dir_rep_ind,[targetfile(1:end-7),...
            extens,'.Vws.layered']); % Layered model
        nameDINsmw=fullfile(dir_rep_ind,[targetfile(1:end-7),...
            extens,'.Vws.smooth']); % Smooth model
        %         nameDINsmlayw=fullfile(dir_rep_ind,[targetfile(1:end-7),...
        %             extens,'.Vws.smlay']); % Smooth layered model
        
        % Save best model
        dinsave(nameDINbest,THbest,VPbest(1:2:end),VSbest(1:2:end),RHObest(1:2:end));
        % Save ridge model
        if ridgecalc==1
            dinsave(nameDINridge,repmat(dz,1,length(ZZ)),VPridge,VSridge,RHOridge);
        end
        % Save mean models (layered, smoothed, layered smooth)
        dinsave(nameDINlay,THmean,VPmean(1:2:end),VSmean(1:2:end),RHOmean(1:2:end));
        dinsave(nameDINsm,repmat(dz,1,length(ZZ)),IVPmean,IVSmean,IRHOmean);
        %         dinsave(nameDINsmlay,repmat(dz,1,length(ZZ)),VPmeanI,VSmeanI,RHOmeanI);
        % Save std models (layered, smoothed, layered smooth)
        dinsave(nameDINlaystd,THstd,VPstd(1:2:end),VSstd(1:2:end),RHOstd(1:2:end));
        dinsave(nameDINsmstd,repmat(dz,1,length(ZZ)),IVPstd,IVSstd,IRHOstd);
        %         dinsave(nameDINsmlaystd,repmat(dz,1,length(ZZ)),VPstdI,VSstdI,RHOstdI);
        % Save weighted models (layered, smoothed, layered smooth)
        if weightcalc==1
            dinsave(nameDINlayw,THweight,VPweight(1:2:end),VSweight(1:2:end),RHOweight(1:2:end));
            dinsave(nameDINsmw,repmat(dz,1,length(ZZ)),IVPweight,IVSweight,IRHOweight);
            %             dinsave(nameDINsmlayw,repmat(dz,1,length(ZZ)),VPweightI,VSweightI,RHOweightI);
        end
        
        % Plot and save average and weighted models
        [f3,h0]=plot_curv(showplot,VSbest,Zbest,[],'-','k',[],1,1,0,fs,...
            'Vs (m/s)','Depth (m)',[],[vsMIN vsMAX],[dpMIN dpMAX],[],[],[],...
            [],[],[],[0 0 24 18],[],0);
        hold on
        str0='Best model';
        h1=plot(VSmean,Zmean,'b-','linewidth',1.5);
        str1='Average layered model';
        %             h2=plot(VSmeanI,ZZ,'g-','linewidth',1.5);
        %             str2='Average parameters and interpolation';
        h3=plot(IVSmean,ZZ,'r-','linewidth',1.5);
        str3='Average smooth model';
        if weightcalc==1 && nmod>1
            h4=plot(VSweight,Zweight,'c-','linewidth',1.5);
            str4='Weighted layered model';
            %                 h5=plot(VSweightI,ZZ*NaN,'g--','linewidth',1.5);
            %                 dashline(VSweightI,ZZ,2,2,2,2,'color','g','linewidth',1.5);
            %                 str5='Weighted parameters and interpolation';
            h6=plot(IVSweight,ZZ,'m-','linewidth',1.5);
            str6='Weighted smooth model';
        end
        if ridgecalc==1 && nmod>1
            h7=plot(VSridge,ZZ,'g-','linewidth',1.5);
            str7='Ridge model';
        end
        colorbar; cblabel('Number of models','Rotation', 270,'VerticalAlignment','Bottom');
        if exist('NN','var')==1
            caxis([0 max(NN(:))]);
        end
        set(cbhandle,'visible','off');
        if (weightcalc==0 || (weightcalc==1 && nmod==1)) && (ridgecalc==0 || (ridgecalc==1 && nmod==1))
            %                 h_legend=legend([h0,h1,h2,h3],str0,str1,str2,str3);
            h_legend=legend([h0,h1,h3],str0,str1,str3);
        elseif weightcalc==1 && ridgecalc==0
            %                 h_legend=legend([h0,h1,h2,h3,h4,h5,h6],str0,str1,str2,str3,str4,str5,str6);
            h_legend=legend([h0,h1,h3,h4,h6],str0,str1,str3,str4,str6);
        elseif weightcalc==0 && ridgecalc==1
            %                 h_legend=legend([h0,h1,h2,h3,h7],str0,str1,str2,str3,str7);
            h_legend=legend([h0,h1,h3,h7],str0,str1,str3,str7);
        else
            %                 h_legend=legend([h0,h1,h2,h3,h4,h5,h6,h7],str0,str1,str2,str3,str4,str5,str6,str7);
            h_legend=legend([h0,h1,h3,h4,h6,h7],str0,str1,str3,str4,str6,str7);
        end
        set(h_legend,'FontSize',10,'linewidth',1,'location','southwest');
        hold off
        % Save figure
        file_mod=fullfile(dir_rep_ind,[targetfile(1:end-7),'.mod1d.mean.',imgform]);
        save_fig(f3,file_mod,imgform,imgres,1,1-testplot);
        close(f3);
        
        if colnb>1
            colnb_tmp=2;
        else
            colnb_tmp=colnb;
        end
        filename_panel0=fullfile(dir_rep_ind,[targetfile(1:end-7),'.mod1d.mean.',imgform]);
        if exist('file_ridge','var')==1 && testplot==1 && nmod>1
            cat_img([file_mod,' ',file_ridge],imgform,colnb_tmp,'south',filename_panel0);
            if concat==1
                delete(file_ridge);
            end
        end
    end
    
    %%
    %%%%%% Plot and save inversion results %%%%%%
    
    if colnb>nmodeinv+1
        colnb_tmp=nmodeinv+1;
    else
        colnb_tmp=colnb;
    end
    
    if (plotinvres==1 || plotparam==1) && nmod>0
        % Create colormaps for accepted and rejected models
        % Selected models (color)
        if nmod==0
            col=[];
            youp=[];
        elseif nmod==1
            col=map2(1,:);
        else
            [col,youp]=createcolormap(misin,map2,Clogscale);
        end
        % Rejected models (grey)
        if nmod==nmaxmod*nrun
            colout=[];
            youpout=[];
        elseif nmod==(nmaxmod*nrun)-1
            colout=map3(1,:);
        else
            [colout,youpout]=createcolormap(misout,map3,Clogscale);
        end
        if 2*std(misout)>min(misout)
            maxmisout=2*std(misout);
        else
            maxmisout=max(misout);
        end
        maxmisout=ceil(maxmisout*100)/100;
        colall=[colout;col];
    end
    
    %%
    
    if plotinvres==1 && nmod>0
        
        %%
        %%%%%% Plot and save calculated dispersion curves %%%%%%
        
        fprintf('\n  Saving calculated dispersion curves\n');
        for ip=1:nmodeinv
            m=modes(ip); % Mode number
            fprintf(['\n      Mode ',num2str(m),'\n']);
            tmp=D{ip}(I,:); % Read dispersion
            tmpin=tmp(OKsort,:); % Select accepted models
            tmpout=tmp(KOsort,:); % Select rejected models
            tmpall=[tmpout;tmpin]; % Create matrix with rejected then accepted models
            
            str1=sprintf([' Accep. models (' num2str(nmod) ' / ' num2str(nrun*nmaxmod) ')']);
            
            f1=plot_curv(showplot,tmpall{1,1},1./tmpall{1,2},[],'-',colall(1,:),2,1,1,...
                cbpos,fs,'Frequency (Hz)','Phase velocity (m/s)',str1,...
                [fMIN fMAX],[VphMIN VphMAX],[],fticks,Vphticks,[],...
                [],[],[0 0 24 18],[],0);
            hold on
            for kk=2:nrun*nmaxmod
                line(tmpall{kk,1},1./tmpall{kk,2},'color',...
                    colall(kk,:),'linewidth',2);
            end
            
            plot(freqresamp{m+1},vresamp{m+1},'k+','linewidth',1.25,'markersize',4);
            if str2double(matrelease(1:4))>2014
                h3=terrorbar(freqresamp{m+1},vresamp{m+1},deltaresamp{m+1},0.75,'units');
                set(h3,'LineWidth',1.25,'Color',[0 0 0])
            else
                h3=errorbar(freqresamp{m+1},vresamp{m+1},deltaresamp{m+1},'k+',...
                    'linewidth',1.25,'markersize',4);
                errorbar_tick(h3,0.75,'units');
            end
            hold off
            colormap(map2);
            sizax=get(gca,'Position');
            set(cbhandle,'visible','off');
            
            % Save figure
            filenamedsp=fullfile(dir_rep_ind,[targetfile(1:end-7),...
                '.calcdisp.M',num2str(m),'.',imgform]);
            save_fig(f1,filenamedsp,imgform,imgres,1,1-testplot);
            close(f1); clear f1;
        end
        filename_dspall=fullfile(dir_rep_ind,[targetfile(1:end-7),...
            '.calcdisp.*.',imgform]);
        
        %%
        %%%%%% Plot and save calculated Vs models %%%%%%
        
        % Get the 1D Vs average model to plot
        if plot1dVS==1
            if calcmod==1
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
                        if nmod>1
                            VSplot=IVSweight;
                            Zplot=ZZ;
                        else
                            VSplot=[];
                            Zplot=[];
                        end
                    end
                    %                     elseif strcmp(modeltype,'smlay')==1
                    %                         if strcmp(avertype,'Vms')==1
                    %                             VSplot=VSmeanI;
                    %                             Zplot=ZZ;
                    %                         else
                    %                             VSplot=VSweightI;
                    %                             Zplot=ZZ;
                    %                         end
                elseif strcmp(modeltype,'ridge')==1
                    VSplot=VSridge;
                    Zplot=ZZ;
                end
            else
                % Velocity file name
                filevel=fullfile(dir_rep_ind,[targetfile(1:end-7),extens,'.',...
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
        
        str1=sprintf([' Accep. models (' num2str(nmod) ' / ' num2str(nrun*nmaxmod) ')']);
        str2=sprintf([' Rejec. models (' num2str((nrun*nmaxmod)-nmod) ' / ' num2str(nrun*nmaxmod) ')']);
        
        f4=plot_curv(showplot,tmpall,zall,[],[],colall,2,1,1,...
            cbpos,fs,'Vs (m/s)','Depth (m)',str1,...
            [vsMIN vsMAX],[dpMIN dpMAX],[],vsticks,dticks,[],...
            [],[],[0 0 24 18],sizax,0);
        colormap(map2);
        hold on
        if plot1dVS==1
            dashline(VSplot,Zplot,2,2,2,2,'color','k','linewidth',2);
        end
        hold off
        set(cbhandle,'visible','off');
        
        % Save figure
        filename_modall=fullfile(dir_rep_ind,[targetfile(1:end-7),...
            '.mod1d.all.',imgform]);
        save_fig(f4,filename_modall,imgform,imgres,1,1-testplot);
        close(f4); clear f4;
        
        filename_cbin=fullfile(dir_rep_ind,[targetfile(1:end-7),...
            '.cbmisin.',imgform]);
        filename_cbout=fullfile(dir_rep_ind,[targetfile(1:end-7),...
            '.cbmisout.',imgform]);
        
        % Save colorbars
        if nmod>1
            % Save colorbars
            if cbpos==1
                cbticks=youp(1:3:end);
            else
                cbticks=youp(1:4:end);
            end
            f4 = plot_colorbar(showplot,[24, 1], cbpos, str1, map2, Clogscale, cbticks,...
                [min(misin) max(misin)], fs, [0 0 24 18], sizax, 2);
            save_fig(f4,filename_cbin,imgform,imgres,1,0);
            close(f4); clear f4;
        end
        if nmod<nrun*nmaxmod-1
            if cbpos==1
                cbticks=youpout(1:2:end);
            else
                cbticks=youpout(1:3:end);
            end
            f4 = plot_colorbar(showplot,[24, 1], cbpos, str2, map3, Clogscale, cbticks,...
                [min(misout) maxmisout], fs, [0 0 24 18], sizax, 1);
            save_fig(f4,filename_cbout,imgform,imgres,1,0);
            close(f4); clear f4;
        end
        
        %%
        if strcmp(imgform,'fig')~=1 && testplot==1
            fprintf('\n  Saving panel figure\n');
            % Figure concatenation
            if strcmp(imgform,'pdf')~=1
                if isunix==1
                    [~,sizefig]=unix(['convert ',fullfile(dir_rep_ind,[targetfile(1:end-7),...
                        '.calcdisp.*.',imgform]),' -format "%w,%h" info:']);
                    sizefig=str2num(sizefig); sizefig=sizefig(1,:);
                    [~,~]=unix(['convert ',filename_modall,' -gravity northeast -extent ',...
                        num2str(sizefig(1)),'x',num2str(sizefig(2)),' ',filename_modall]);
                else
                    com1=['img_convert ',fullfile(dir_rep_ind,[targetfile(1:end-7),...
                        '.calcdisp.*.',imgform]),' -format "%w,%h" info:'];
                    com1=strrep(com1,'\','/');
                    [~,sizefig]=unix(com1);
                    sizefig=str2num(sizefig); sizefig=sizefig(1,:);
                    com1=['img_convert ',filename_modall,' -gravity northeast -extent ',...
                        num2str(sizefig(1)),'x',num2str(sizefig(2)),' ',filename_modall];
                    com1=strrep(com1,'\','/');
                    [~,~]=unix(com1);
                end
            end
            
            filename_imgtmp=fullfile(dir_rep_ind,[targetfile(1:end-7),...
                '.imgtmp.',imgform]);
            cat_img([filename_dspall,' ',filename_modall],imgform,colnb_tmp,'east',filename_imgtmp,0);
            
            filename_cbtmp=fullfile(dir_rep_ind,[targetfile(1:end-7),...
                '.cb.',imgform]);
            filename_panel=fullfile(dir_rep_ind,[targetfile(1:end-7),...
                '.invresults.',imgform]);
            
            if  nmod>1 && nmod<nrun*nmaxmod-1
                if (cbpos==1 && colnb_tmp==nmodeinv+1) || (cbpos==2 && colnb_tmp>1)
                    cat_img([filename_cbout,' ',filename_cbin],imgform,2,'center',filename_cbtmp,0);
                elseif (cbpos==1 && colnb_tmp<nmodeinv+1) || (cbpos==2 && colnb_tmp==1)
                    cat_img([filename_cbout,' ',filename_cbin],imgform,1,'east',filename_cbtmp,0);
                end
                
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
            elseif nmod==1
                movefile(filename_cbout,filename_cbtmp)
            elseif nmod>=nrun*nmaxmod-1
                movefile(filename_cbin,filename_cbtmp)
            end
            
            if cbpos==1 && colnb_tmp==nmodeinv+1
                cat_img([filename_imgtmp,' ',filename_cbtmp],imgform,2,'south',filename_panel);
            elseif cbpos==1 && colnb_tmp<nmodeinv+1
                cat_img([filename_imgtmp,' ',filename_cbtmp],imgform,2,'center',filename_panel);
            elseif cbpos==2 && colnb_tmp==1
                cat_img([filename_imgtmp,' ',filename_cbtmp],imgform,1,'east',filename_panel);
            else
                cat_img([filename_imgtmp,' ',filename_cbtmp],imgform,1,'center',filename_panel);
            end
            if concat==1
                if nmod>1 && nmod<nrun*nmaxmod-1
                    delete(filename_cbin,filename_cbout,filename_cbtmp,filename_imgtmp,filename_modall,filename_dspall);
                elseif nmod==1
                    delete(filename_cbout,filename_cbtmp,filename_imgtmp,filename_modall,filename_dspall);
                elseif nmod>=nrun*nmaxmod-1
                    delete(filename_cbin,filename_cbtmp,filename_imgtmp,filename_modall,filename_dspall);
                end
            end
        end
    end
    
    %%
    %%%%%% Plot and save inversion parameters %%%%%%
    
    % Get selected parameters settings
    if plotparam==1 && nmod>0
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
                str1=sprintf([' Accep. models ('...
                    num2str(nmod) ' / ' num2str(nrun*nmaxmod) ')']);
                str2=sprintf([' Rejec. models (' num2str((nrun*nmaxmod)-nmod)...
                    ' / ' num2str(nrun*nmaxmod) ')']);
                
                f4=plot_curv(showplot,p1(np1(i),:),p2(np2(j),:),[],'.','k',2,1,1,...
                    cbpos,fs,[axtit1,' [layer ',num2str(np1(i)),']'],[axtit2,' [layer ',num2str(np2(j)),']'],...
                    str1,[XMIN XMAX],[YMIN YMAX],[],vsticks,dticks,[],...
                    [],[],[0 0 24 18],sizax,0);
                hold on
                for ii=1:length(p1(np1(i),:))
                    line(p1(np1(i),ii),p2(np2(j),ii),'color',colall(ii,:),'markerfacecolor',colall(ii,:),...
                        'markersize',5,'marker','o','linestyle','none');
                end
                hold off
                colormap(map2);
                sizax=get(gca,'Position');
                set(cbhandle,'visible','off');
                drawnow
                
                % Save figure
                filename=fullfile(dir_rep_ind,[targetfile(1:end-7),...
                    '.param.',param1,num2str(np1(i)),param2,num2str(np2(j)),'.',imgform]);
                save_fig(f4,filename,imgform,imgres,1,0);
                close(f4); clear f4;
                
                filename_cbin=fullfile(dir_rep_ind,[targetfile(1:end-7),...
                    '.cbmisin.',imgform]);
                filename_cbout=fullfile(dir_rep_ind,[targetfile(1:end-7),...
                    '.cbmisout.',imgform]);
                
                % Save colorbars
                if i==1 && j==1
                    if nmod>1
                        % Save colorbars
                        if cbpos==1
                            cbticks=youp(1:3:end);
                        else
                            cbticks=youp(1:4:end);
                        end
                        f4 = plot_colorbar(showplot,[24, 1], cbpos, str1, map2, Clogscale,cbticks,...
                            [min(misin) max(misin)], fs, [0 0 24 18], sizax, 2);
                        save_fig(f4,filename_cbin,imgform,imgres,1,0);
                        close(f4); clear f4;
                    end
                    if nmod<nrun*nmaxmod-1
                        if cbpos==1
                            cbticks=youpout(1:2:end);
                        else
                            cbticks=youpout(1:3:end);
                        end
                        f4 = plot_colorbar(showplot,[24, 1], cbpos, str2, map3, Clogscale,cbticks,...
                            [min(misout) maxmisout], fs, [0 0 24 18], sizax, 1);
                        save_fig(f4,filename_cbout,imgform,imgres,1,0);
                        close(f4); clear f4;
                    end
                end
            end
        end
        
        if colnb>length(np1)*length(np2)
            colnb_tmp=length(np1)*length(np2);
        else
            colnb_tmp=colnb;
        end
        
        if testplot==1
            filename_paramall=fullfile(dir_rep_ind,[targetfile(1:end-7),...
                '.param.*.',imgform]);
            filename_imgtmp=fullfile(dir_rep_ind,[targetfile(1:end-7),...
                '.imgtmp.',imgform]);
            cat_img(filename_paramall,imgform,colnb_tmp,'east',filename_imgtmp,0);
            
            filename_cbtmp=fullfile(dir_rep_ind,[targetfile(1:end-7),...
                '.cb.',imgform]);
            if (cbpos==1 && colnb_tmp==length(np1)*length(np2)) || (cbpos==2 && colnb_tmp>1)
                cat_img([filename_cbout,' ',filename_cbin],imgform,2,'center',filename_cbtmp,0);
            elseif (cbpos==1 && colnb_tmp<length(np1)*length(np2)) || (cbpos==2 && colnb_tmp==1)
                cat_img([filename_cbout,' ',filename_cbin],imgform,1,'east',filename_cbtmp,0);
            end
            
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
            
            filename_param=fullfile(dir_rep_ind,[targetfile(1:end-7),...
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
            if concat==1
                delete(filename_cbin,filename_cbout,filename_cbtmp,filename_imgtmp,filename_paramall);
            end
        end
    end
end
tend=toc(tstart);
[time_string]=secs2hms(tend);
if inversion==1 || plotinvres==1 || plotparam==1
    fprintf(['\n  Elapsed time: ',time_string,'\n']);
    fprintf('\n  **********************************************************');
    fprintf('\n  **********************************************************\n');
end
end