%%% SURFACE-WAVE dispersion INVERSION & PROFILING (SWIP)
%%% MODULE B : SWIPparam.m
%%% S. Pasquet - V18.11.26
%%% SWIPparam.m creates the parameterization file required to invert with
%%% module C the dispersion curves picked in module A
%%% It can create either a single parameterization file for all the profile
%%% or create automatic parameterization for each Xmid based on a Vp model

%%%-------------------------%%%
%%% START OF INITIALIZATION %%%

run('SWIP_defaultsettings') % Read default settings

%%% END OF INITIALIZATION %%%
%%%-----------------------%%%

fprintf(['\n  Creating "type ',num2str(paramtype),'" parameterization\n']);
if paramtype==0 % User-defined parameterization
    if length(thmin)<nlay
        thmin = [thmin repmat(thmin(end),1,nlay-length(thmin))];
    end
    if length(thmax)<nlay
        thmax = [thmax repmat(thmax(end),1,nlay-length(thmax))];
    end
    
    paramname=mod2param(nlay,nsublay,thmin,thmax,shape,lvz,...
        Vpmin,Vpmax,Numin,Numax,Vsmin,Vsmax,Rhomin,Rhomax,...
        Vplink,Nulink,0,Rholink,paramname); % Create parameterization
    if isempty(paramname)==1
        return
    end
    parampath=fullfile(pwd,'file.param');
    if exist(parampath,'dir')==7
        movefile(paramname,parampath);
    end
    fprintf(['\n  Parameterization file saved in file.param as ',paramname,'\n\n']);
    
else % Semi-automatic parameterization
    dir_all=dir_create(0);
    if dir_all.dir_main==0
        return
    end
    % Get SWIP subfolder and parameters
    dir_dat=dir_all.dir_dat;
    dir_targ=dir_all.dir_targ;
    matstruct=dir(fullfile(dir_dat,'*.param.mat')); % .mat file containing subproject parameters
    matfile=fullfile(dir_dat,matstruct.name);
    load(matfile); % Read .mat file
    dx=acquiparam.dx;
    nWmin=stackdisp.nWmin;
    nWmax=stackdisp.nWmax;
    xsca=pomega.xsca;
    XmidT=xmidparam.XmidT; % Xmid positions
    Xlength=xmidparam.Xlength; % Number of Xmids
    xmidformat=stackdisp.xmidformat; % Format Xmid with correct number of decimals
    zround=xmidparam.zround; % Get topography
    
    if length(thmin) == 1
        thmin = repmat(thmin,1,nlay);
    elseif length(thmin)<nlay
        thmin = [thmin repmat(thmin(end),1,nlay-length(thmin))];
    end
    if length(thmax) == 1
        thmax = repmat(thmax,1,nlay);
    elseif length(thmax)<nlay
        thmax = [thmax repmat(thmax(end),1,nlay-length(thmax))];
    end
    thmean_vec = (thmin(:)+2*thmax(:))/3;
    thmean = mean(thmean_vec);
    lvz_start = lvz;
    
    maxdepth=sum(thmax(1:end));
    depth=max(zround):-dz:min(zround)-maxdepth; % Depth vector with topo
    % Select Xmids
    if exist('Xmidselec','var')~=1 || isempty(Xmidselec)==1
        Xmidselec=1:Xlength;
    end
    if max(Xmidselec)>Xlength
        Xmidselec=Xmidselec(Xmidselec<=Xlength);
    end
    nshot=xmidparam.nshot;
    
    wave=targopt.wave(1);
    resampvec=targopt.resampvec;
    sampling=targopt.sampling;
    lmaxpick=targopt.lmaxpick; 

    ZZ_emp = bsxfun(@minus,zround,resampvec'/depth_fac);
    vs_emp=zeros(size(ZZ_emp))*NaN;
    delta_emp = vs_emp;
    
    vph2dobs=zeros(length(resampvec),Xlength)*NaN;
    delta2dobs=zeros(length(resampvec),Xlength)*NaN;
    
    % Read refraction velocity model
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
    
    if filevel==0
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        fprintf('\n   Please select a Vp model file');
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
        return
    end
    Vpfile=fullfile(pathvel,filevel); % File with velocity (3 columns X,Z,Vp)
    try
        [VpI,XI,ZI]=readtomo(Vpfile,0,XmidT,depth,xsca,vpaver,[nWmin,nWmax],dx); % Read Vp tomo file
    catch
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        fprintf('\n   Invalid file - Please select a valid Vp model file');
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
        return
    end
    
%     depth2=max(zround)-cumsum(thmean_vec); % Depth vector with topo
    depth2 = [zround;repmat(zround,length(cumsum(thmean_vec)),1)-repmat(cumsum(thmean_vec),1,length(zround))];
    vpmat=zeros(size(depth2)).*NaN;
    vsmat=zeros(size(depth2)).*NaN;
    lvzmat=zeros(size(depth2)).*NaN;
    indf=zeros(Xlength,1); indi=indf;
    test_lvz = zeros(Xlength,1);
    
    if paramtype == 5
        linklvz = 1;
    end
    
    Vpmin_init = Vpmin;
    Vpmax_init = Vpmax;
    Vsmin_init = Vsmin;
    Vsmax_init = Vsmax;

    fprintf('\n  **********************************************************');
    fprintf('\n  **********************************************************\n');
    
    %%%%%% Loop over all Xmids %%%%%%
    for ix=Xmidselec
        if sum(nshot(ix,:))>=0
            vptomo=VpI(VpI(:,ix)>0,ix); % Get Vp non-NaN and non-zero values
            % Get index corresponding to the maximum HSD
%             indf(ix)=round((nlay-1)*thmean/thmean);
%             indf(ix) = nlay-1;
%             % Look for topo index
%             crit=abs(zround(ix)-depth2);
%             indi(ix)=find(crit==min(crit),1);
%             if indi(ix)+indf(ix)>length(depth2)
%                 indi(ix)=indi(ix)-1;
%             end
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
%                 thick=repmat(thmean(1),1,nlay);
                thick = thmean_vec';
                thinc = [0,cumsum(thick)];
                % Resample vptomo along parameterization thicknesses
                [vpmean,Vpmin,Vpmax]=velresamp(zinc,vptomo,thinc);
                diffmean=diff(vpmean);
                for i=1:length(vpmean)
                    if linklvz~=1 && i>1
                        if vpmean(i)< vpmean(i-1)
                            vpmean(i)=vpmean(i-1);
                            Vpmin(i)=Vpmin(i-1);
                            Vpmax(i)=Vpmax(i-1);
                        end
                    end
                    if vpmean(i)-Vpmin(i)<vfac*vpmean(i)
                        Vpmin(i)=vpmean(i)-vfac*vpmean(i);
                    end
                    if Vpmax(i)-vpmean(i)<vfac*vpmean(i)
                        Vpmax(i)=vpmean(i)+vfac*vpmean(i);
                    end
                end
                
                if linklvz==1
                    % Authorize Vs low velocity zone if one is present in Vp
                    if isempty(vpmean(diffmean<0))==0
                        lvz_tmp=zeros(1,nlay);
                        if length(Vplink)==1
                            lvz_tmp(find(diff(Vpmax)<0)+1)=1;
                        else
                            for i=1:nlay-1
                                if Vplink(i)==1 && diffmean(i)<0
                                    lvz_tmp(i)=1;
                                end
                            end
                        end
                    end
                end
                if paramtype==2 % Fixed Vp
                    deltav=lorentzerr(vpmean',cumsum(thmean_vec'),10,1);
                    Vpmax=vpmean+deltav;
                    Vpmin=vpmean-deltav;
                    if length(Vpmin_init) > 1
                        Vpmin(end) = Vpmin_init(end);
                        Vpmax(end) = Vpmax_init(end);
                    end
                elseif paramtype==3 % Fixed thickness
                    thmin=thmax;
                elseif paramtype==4 % Fixed Vp and thickness
                    Vpmax=vpmean;
                    Vpmin=vpmean;
                    thmin=thmax;
                elseif paramtype==5
                    if length(Vpmin_init) > 1
                        Vpmin(end) = Vpmin_init(end);
                        Vpmax(end) = Vpmax_init(end);
                    end
                end
%                 vpmat(indi(ix):indi(ix)+indf(ix),ix) = vpmean;
                    vpmat(:,ix) = [vpmean;vpmean(end)];
            else
                Vpmin = Vpmin_init;
                Vpmax = Vpmax_init;
%                 vpmat(indi(ix):indi(ix)+indf(ix),ix) = mean([Vpmin Vpmax]);
                vpmat(:,ix) = mean([Vpmin Vpmax]);
            end
            
            if paramtype == 5;
                nametarg=fullfile(dir_targ,[num2str(XmidT(ix),xmidformat),'.target']);
                if exist(nametarg,'file')==2
                    % Read target file to get picked dispersion curves
                    [freqresamp,vresamp,deltaresamp,modes]=targ2pvc(nametarg);
                    npvc=length(modes);
                    lmaxpicktmp=zeros(npvc,1);
                    for ip=1%:npvc
                        % Resample in lambda or frequency
                        [freqresamp{modes(ip)+1},vresamp{modes(ip)+1},deltaresamp{modes(ip)+1}]=...
                            resampvel(freqresamp{modes(ip)+1},vresamp{modes(ip)+1},...
                            deltaresamp{modes(ip)+1},resampvec,sampling,1);
                        lmaxpicktmp(ip)=max(vresamp{modes(ip)+1}./freqresamp{modes(ip)+1});
                        if ip == 1
                            meanveresamp = vresamp{ip};
                            meandeltaresamp = deltaresamp{ip};
                        else
                            meanveresamp = nanmin([meanveresamp;vresamp{ip}]);
                            meandeltaresamp = nanmean([meandeltaresamp;deltaresamp{ip}]);
                        end
                    end
                    lmaxpick(ix)=max(lmaxpicktmp);
                    vph2dobs(:,ix)=meanveresamp';
                    delta2dobs(:,ix)=meandeltaresamp';
                else
                    fprintf(['\n  No dispersion picked for Xmid',num2str(ix),' = ',...
                        num2str(XmidT(ix),xmidformat),' m\n']);
                    npvc=0;
                    continue
                end

                % Fill velocity matrix with average model
                vs_emp(:,ix) = vph2dobs(:,ix)/vel_fac;
                delta_emp(:,ix) = delta2dobs(:,ix)/vel_fac;
                indnan = find(~isnan(vs_emp(:,ix)),1,'last');
                if vs_emp(indnan,ix) < mean([Vsmin_init,Vsmax_init])
                    vs_emp(indnan+1:end,ix) = mean([Vsmin_init,Vsmax_init]);
                    delta_emp(indnan+1:end,ix) = (mean([Vsmin_init,Vsmax_init])*delta_emp(indnan,ix)/vs_emp(indnan,ix));
                else
                    vs_emp(indnan+1:end,ix) = vs_emp(indnan,ix);
                    delta_emp(indnan+1:end,ix) = delta_emp(indnan,ix);
                end
                
                vs_emp_tmp = vs_emp(vs_emp(:,ix)>0,ix); % Get Vs non-NaN and non-zero values
                delta_emp_tmp = delta_emp(vs_emp(:,ix)>0,ix); % Get delta non-NaN and non-zero values
                z_tmp = ZZ_emp(vs_emp(:,ix)>0,ix);
                
                if isempty(vs_emp_tmp)~=1
                    vs_emp_tmp = [vs_emp_tmp(1);vs_emp_tmp;vs_emp_tmp(end)];
                    delta_emp_tmp = [delta_emp_tmp(1);delta_emp_tmp;delta_emp_tmp(end)];
                    
                    depth_tmp = abs(z_tmp - zround(ix));
                    zinc = [0;depth_tmp];
                    if max(zinc) > maxdepth
                        vs_emp_tmp = vs_emp_tmp(1:length(zinc)+1);
                        delta_emp_tmp = delta_emp_tmp(1:length(zinc)+1);
                    elseif max(zinc) < maxdepth
                        zinc = zinc(zinc < maxdepth);
                        vs_emp_tmp2 = mean([Vsmin_init(end),Vsmax_init(end)])*ones(length(zinc)+1,1);
                        vs_emp_tmp2(1:length(vs_emp_tmp)-1) = vs_emp_tmp(1:length(vs_emp_tmp)-1);
                        vs_emp_tmp = vs_emp_tmp2;
                        delta_emp_tmp2 = (mean([Vsmin_init,Vsmax_init])*delta_emp_tmp(end)/vs_emp_tmp(end))*ones(length(zinc)+1,1);
                        delta_emp_tmp2(1:length(delta_emp_tmp)-1) = delta_emp_tmp(1:length(vs_emp_tmp)-1);
                        delta_emp_tmp = delta_emp_tmp2;
                    end
                    %                     thick = repmat(thmean(1),1,nlay);
                    thick = thmean_vec';
                    thinc = [0,cumsum(thick)];
                    
                    % Resample vs_emp along parameterization thicknesses
                    [vsmean,~,~] = velresamp(zinc,vs_emp_tmp,thinc);
                    [deltamean,deltamin,deltamax] = velresamp(zinc,delta_emp_tmp,thinc);
                    diffmean=diff(vsmean);
                    %                     Vsmin = vsmean-(deltamean + 8*vfac*deltamean);
                    %                     Vsmax = vsmean+(deltamean + 8*vfac*deltamean);
                    %                     Vsmin(Vsmin<Vsmin_init) = Vsmin_init;
                    %                     Vsmax(Vsmax>Vsmax_init) = Vsmax_init;
                    
                    % Authorize Vs low velocity zone if one is present in Vs_emp
                    
                    if lvz_start == 1
                        lvz=zeros(1,nlay);
                        if ~isempty(vsmean(diffmean<0))
                            ind_lvz = diff(vsmean')<0;
                            if any(ind_lvz(1:end-1))
                                test_lvz(ix) = 1;
                                lvz(find(ind_lvz)+1)=1;
                            end
                        end
                        if isempty(vs_emp_tmp(diff(vs_emp_tmp)<0))==0
%                             depth_lvz = zinc(diff(vs_emp_tmp)<-0.025*vs_emp_tmp(1:end-1));
                            depth_lvz = zinc(diff(vs_emp_tmp)<=-10 & diff(vs_emp_tmp)<-0.025*vs_emp_tmp(1:end-1));
                            for ii = 1:nlay-1
                                if any(depth_lvz>=thinc(ii) & depth_lvz<thinc(ii+1))
                                    lvz(ii+1) = 1;
                                    test_lvz(ix) = test_lvz(ix) + 1;
                                end
                            end
                        end
                    else
                        lvz=zeros(1,nlay);
                    end
                    lvzmat(:,ix) = [0;lvz'];
                    
                    % Get index corresponding to the maximum HSD
%                     indf(ix)=round((nlay-1)*thmean/thmean);
%                     indf(ix)=nlay-1;
%                     % Look for topo index
%                     crit=abs(zround(ix)-depth2);
%                     indi(ix)=find(crit==min(crit),1);
%                     if indi(ix)+indf(ix)>length(depth2)
%                         indi(ix)=indi(ix)-1;
%                     end
%                     vsmat(indi(ix):indi(ix)+indf(ix),ix) = vsmean;
                    
                    Vsmin = Vsmin_init;
                    Vsmax = Vsmax_init;
%                     vsmat(indi(ix):indi(ix)+indf(ix),ix) = mean([Vsmin Vsmax]);
                    vsmat(:,ix) = mean([Vsmin Vsmax]);
                end
            else
                vs_emp_tmp = [];
                Vsmin = Vsmin_init;
                Vsmax = Vsmax_init;
%                 vsmat(indi(ix):indi(ix)+indf(ix),ix) = mean([Vsmin Vsmax]);
                vsmat(:,ix) = mean([Vsmin Vsmax]);
            end
            
            if isempty(vptomo)==1 && ~isempty(vs_emp_tmp)
                fprintf(['\n  No Vp at Xmid',num2str(ix),' = ',num2str(XmidT(ix),xmidformat),' m - Use default parameterization\n']);
            end
            
            fprintf(['\n  Saving parameterization for Xmid',num2str(ix),' = ',num2str(XmidT(ix),xmidformat),' m \n']);
            
            % Save parameterization
            paramname=fullfile(dir_targ,[num2str(XmidT(ix),xmidformat),'.type',...
                num2str(paramtype),'.param']);
            paramname=mod2param(nlay,nsublay,thmin,thmax,1,lvz,...
                Vpmin,Vpmax,Numin,Numax,Vsmin,Vsmax,Rhomin,Rhomax,...
                Vplink,Nulink,0,Rholink,paramname); % Create parameterization
            if isempty(paramname)==1
                return
            end
        end
    end
    
    %%
    XX = repmat(XmidT,size(depth2,1),1);
    XX_temp = repmat(XmidT,size(ZZ_emp,1),1);
    if plot2dVP==1 && length(XmidT)>1
        f1=plot_img(1,XI,ZI,VpI,map5,1,0,1,14,'X (m)','Elevation (m)','Vp (m/s)',...
            [],[floor(min(depth)/10)*10 ceil(max(depth)/10)*10],[],...
            [],[],[],[],[],[],[0 20 24 12],[],1,3); drawnow;
        vplim = get(gca,'clim');
        
        f2=plot_img(2,XX,depth2,vpmat,map5,1,0,1,14,'X (m)','Elevation (m)','Vp (m/s)',...
            [],[floor(min(depth)/10)*10 ceil(max(depth)/10)*10],vplim,...
            [],[],[],[],[],[],[26 20 24 12],[],1,3);
    end
    
    if plot2dVS==1 &&  paramtype == 5 && length(XmidT)>1
        f1=plot_img(3,XX_temp,ZZ_emp,vs_emp,map5,1,0,1,14,'X (m)','Elevation (m)','Vs (m/s)',...
            [],[floor(min(depth)/10)*10 ceil(max(depth)/10)*10],[],...
            [],[],[],[],[],[],[0 0 24 12],[],1,3); drawnow;
        vslim = get(gca,'clim');
        
        f2=plot_img(4,XX,depth2,lvzmat,gray(2),1,0,1,14,'X (m)','Elevation (m)','LVZ',...
            [],[floor(min(depth)/10)*10 ceil(max(depth)/10)*10],[0 1],...
            [],[],[0 1],[],[],[],[26 0 24 12],[],1,3);
        
%         f3=plot_img(5,XX,depth2,lvzmat,map5,1,0,1,14,'X (m)','Elevation (m)','Vs (m/s)',...
%             [],[floor(min(depth)/10)*10 ceil(max(depth)/10)*10],[],...
%             [],[],[],[],[],[],[26 0 24 12],[],1,0);
    end
end