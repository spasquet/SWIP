%%% SURFACE-WAVE dispersion INVERSION & PROFILING (SWIP)
%%% MODULE B : SWIPparam.m
%%% S. Pasquet - V17.05.25
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
    maxdepth=(nlay-1)*thmax;
    depth=max(zround):-dz:min(zround)-maxdepth; % Depth vector with topo
    % Select Xmids
    if exist('Xmidselec','var')~=1 || isempty(Xmidselec)==1
        Xmidselec=1:Xlength;
    end
    if max(Xmidselec)>Xlength
        Xmidselec=Xmidselec(Xmidselec<=Xlength);
    end
    nshot=xmidparam.nshot;
    
    % Read refraction velocity model
    fprintf('\n  Select Vp model file\n');
    [filevel,pathvel]=uigetfile({'*.model;*.dat;*.xzv;*.txt'},'Select Vp model');
    if filevel==0
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        fprintf('\n   Please select a Vp model file');
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
        return
    end
    Vpfile=fullfile(pathvel,filevel); % File with velocity (3 columns X,Z,Vp)
    try
        [VpI,XI,ZI]=readtomo(Vpfile,0,XmidT,depth,xsca,vpaver,mean([nWmin,nWmax]),dx); % Read Vp tomo file
    catch
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        fprintf('\n   Invalid file - Please select a valid Vp model file');
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
        return
    end
    
    % Plot refraction model resampled along Xmids
    if plot2dVP==1
        f1=plot_img(1,XI,ZI,VpI,haxby(32),1,0,1,14,'X (m)',...
            'Altitude (m)','Vp (m/s)',[],[floor(min(depth)/10)*10 ...
            ceil(max(depth)/10)*10],[],[],[],[],[],[],[],[0 0 24 12],[],1);
    end
    
    depth2=max(zround):-thmax:min(zround)-maxdepth; % Depth vector with topo
    vpmat=zeros(length(depth2),Xlength).*NaN;
    indf=zeros(Xlength,1); indi=indf;
    
    fprintf('\n  **********************************************************');
    fprintf('\n  **********************************************************\n');
    
    %%%%%% Loop over all Xmids %%%%%%
    
    for ix=Xmidselec
        if sum(nshot(ix,:))>=0
            vptomo=VpI(VpI(:,ix)>0,ix); % Get Vp non-NaN and non-zero values
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
                thick=repmat(thmax(1),1,nlay);
                thinc=[0,cumsum(thick)];
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
                        lvz=zeros(1,nlay);
                        if length(Vplink)==1
                            lvz(find(diff(Vpmax)<0)+1)=1;
                        else
                            for i=1:nlay-1
                                if Vplink(i)==1 && diffmean(i)<0
                                    lvz(i)=1;
                                end
                            end
                        end
                    end
                end
                if paramtype==2 % Fixed Vp
                    Vpmax=vpmean;
                    Vpmin=vpmean;
                elseif paramtype==3 % Fixed thickness
                    thmin=thmax;
                elseif paramtype==4 % Fixed Vp and thickness
                    Vpmax=vpmean;
                    Vpmin=vpmean;
                    thmin=thmax;
                end
                % Get index corresponding to the maximum HSD
                indf(ix)=round((nlay-1)*thmax/thmax);
                % Look for topo index
                crit=abs(zround(ix)-depth2);
                indi(ix)=find(crit==min(crit),1);
                if indi(ix)+indf(ix)>length(depth2)
                    indi(ix)=indi(ix)-1;
                end
                vpmat(indi(ix):indi(ix)+indf(ix),ix)=vpmean;
                fprintf(['\n  Saving parameterization for Xmid',num2str(ix),' = ',num2str(XmidT(ix),xmidformat),' m \n']);
            else
                fprintf(['\n  No Vp from tomo at Xmid',num2str(ix),' = ',num2str(XmidT(ix),xmidformat),' m - Cannot create parameterization\n']);
            end
            % Save parameterization
            paramname=fullfile(dir_targ,[num2str(XmidT(ix),xmidformat),'.type',...
                num2str(paramtype),'.param']);
            paramname=mod2param(nlay,nsublay(1),thmin(1),thmax(1),1,lvz,...
                Vpmin,Vpmax,Numin(1),Numax(1),Vsmin(1),Vsmax(1),Rhomin(1),Rhomax(1),...
                Vplink(1),Nulink(1),0,Rholink(1),paramname); % Create parameterization
            if isempty(paramname)==1
                return
            end
        end
    end
    % Plot refraction model resampled along Xmids and parameterization
    if plot2dVP==1
        f2=plot_img(2,XmidT,depth2,vpmat,haxby(32),1,0,1,14,'X (m)',...
            'Altitude (m)','Vp (m/s)',[],[floor(min(depth)/10)*10 ...
            ceil(max(depth)/10)*10],[],[],[],[],[],[],[],[25 0 24 12],[],1);
    end
end