%%% TRAVELTIME TOMOGRAPHY PROCESSING
%%% MODULE A : TOMOpickfb.m
%%% S. Pasquet - V18.10.17
%%% TOMOpickfb.m allows to pick first breaks and plot seismograms

run('tomo_defaultsettings');

% Folder initialization
dir_start=pwd;
dir_img=fullfile(dir_start,'file.img');
if exist(dir_img,'dir')~=7
    mkdir(dir_img);
end
dir_img_seismo=fullfile(dir_img,'seismo');
if exist(dir_img_seismo,'dir')~=7
    mkdir(dir_img_seismo);
end
dir_xzv=fullfile(dir_start,'file.xzv');
if exist(dir_xzv,'dir')~=7
    mkdir(dir_xzv);
end
dir_mat=fullfile(dir_start,'file.mat');
if exist(dir_mat,'dir')~=7
    mkdir(dir_mat);
end
dir_pfb=fullfile(dir_start,'file.pfb');

% Read SU file
sustruct=dir(fullfile(dir_start,'*.su'));
if length(sustruct)>1
    sufile=uigetfile('*.su','Select SU file to use');
    if sufile==0
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!');
        fprintf('\n   Please select a SU file');
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
    end
else
    sufile=sustruct.name;
end
lineName=sufile(1:end-3);
acquiparam=get_acquiparam(sufile,[],1);
xsca=acquiparam.xsca;
dt=acquiparam.dt*1000;
fldr=sort(acquiparam.fldr);
NSx=length(fldr);
nbpicks=zeros(size(fldr))*NaN;
Gx=acquiparam.Gx;
Sx=acquiparam.Sx;
Gz=acquiparam.Gz;
Sz=acquiparam.Sz;
Offset=Gx-Sx;
Gxpick=[]; Gzpick=[];
Sxpick=[]; Szpick = [];
Tpick=[]; Offpick=[];
shotnumformat=['%0',num2str(order(NSx,10)+1),'i']; % Precision

if exist('dt_min','var')~=1 || isempty(dt_min)==1
    dt_min = dt;
end

% Get selected shots
if exist('Sxselec','var')~=1 || isempty(Sxselec)==1
    Sxselec=1:NSx;
end
if max(Sxselec)>NSx
    Sxselec=Sxselec(Sxselec<=NSx);
end

% Bandpass filter
if filt==1
    sufilefilt=[sufile,'.filt'];
    com1=sprintf('sufilter < %s f=%d,%d,%d,%d amps=1,1,0,0 > %s',...
        sufile,fcutlow,fcuthigh-taper,fcuthigh,fcuthigh+taper,sufilefilt);
    unix(com1);
    sufileOK=sufilefilt;
else
    sufileOK=sufile;
end

% Load seismogram or not
if pick>=1 || save_seis==1 || save_seis_pck==1
    loadseis=1;
else
    loadseis=0;
end

% Get auto limits for plots
if exist('xMIN','var')==0 || isempty(xMIN)==1 || exist('xMAX','var')==0 || isempty(xMAX)==1
    xMIN=min([Sx;Gx]); xMAX=max([Sx;Gx]);
end
if exist('offMIN','var')==0 || isempty(offMIN)==1 || exist('offMAX','var')==0 || isempty(offMAX)==1
    offMIN=min(Offset); offMAX=max(Offset);
end
if exist('offMAX_pick','var')==0 || isempty(offMAX_pick)==1
    offMAX_pick=max(Offset);
end

% Plot offset vs shot position diagram
if save_src_off==1 && plot_pos==1 && (pick==1 || pick==2) && showplot==1
    fig3=plot_scat(1,(Sx+Gx)/2,Offset,'k','c',markersize/2,1,map1,...
        axetop,0,1,fs,'CMP (m)','Offset (m)','Observed time (ms)',...
        [xMIN xMAX],[offMIN offMAX],[tMIN tMAX],xticks,offticks,tticks,[],[],[1 1 24 12],[],horex);
end

pfbmatfile=fullfile(dir_mat,[lineName '.pfb.mat']); % Name of file to open %
pfb_all = cell(NSx,1);
if exist(pfbmatfile,'file')==2
    load(pfbmatfile);
    if length(pfb_all) ~= NSx
        pfb_all = [pfb_all;cell(NSx-length(pfb_all),1)];
    end
else
    for ix=1:NSx
        pfbfile=fullfile(dir_pfb,[num2str(ix,shotnumformat),'.pfb.dat']); % Name of file to open %
        if exist(pfbfile,'file')==2
            aa=load(pfbfile);
            pfb_all{ix}.sx = aa(1,2);
            pfb_all{ix}.sz = aa(1,3);
            pfb_all{ix}.ix = aa(1,4);
            pfb_all{ix}.xpick = aa(2:end,2);
            pfb_all{ix}.zpick = aa(2:end,3);
            pfb_all{ix}.tpick = aa(2:end,4)*1000;
            pfb_all{ix}.tshift = 0;
        end
    end
    if exist(dir_pfb,'dir')==7
        rmdir(dir_pfb,'s');
    end
    save(pfbmatfile,'pfb_all');
end

if (pick==1 || pick==2) && showplot==1
    f3=plot_curv(3,[],NaN,NaN,'.-','r',1,0,0,0,fs+4,'Distance (m)','Traveltime (ms)',[],...
        [xMIN xMAX],[tMIN tMAX],[],xticks,tticks,[],[],[],[1 24 65 12],[]);
end

sx_prev = [];
hh=zeros(length(Sxselec),1);
i=0; % Shot number flag
while i<length(Sxselec)
    i=i+1; ix=Sxselec(i);
    if pick>=1 || save_seis==1 || save_seis_pck==1
        
        [seismomat,xseis,zseis,tseis,sx,sz]=shotselec(sufileOK,fldr(ix),xsca,loadseis,0,normalize); % load shot
        tseis=1000*tseis;
        
        if pick==0 && save_seis==1 || save_seis_pck==1
            fprintf(['\n  Processing shot located at ',num2str(sx),' m\n']);
        elseif pick>=1
            fprintf(['\n  Picking shot located at ',num2str(sx),' m\n']);
        end
        
        %         cor_fac = repmat(abs(xseis-sx),1,length(tseis));
        %         cor_fac(cor_fac==0) = 0.25;
        %         seismomat = (seismomat .* cor_fac);
        %         seismomat = seismomat./repmat(max(seismomat,[],2),1,length(tseis));
        %         seismomat = seismomat - repmat(mean(seismomat,2),1,length(tseis));
        
        if ~exist('tMAX_pick','var') || isempty(tMAX_pick)
            tMAX_seis_tmp = max(tseis);
            tMAX_pick = max(tseis);
        else
            tMAX_seis_tmp = tMAX_pick;
        end
        
        if ~exist('tMIN_pick','var') || isempty(tMIN_pick)
            tMIN_seis_tmp = min(tseis);
            tMIN_pick = min(tseis);
        else
            tMIN_seis_tmp = tMIN_pick;
        end
        if isempty(pfb_all{ix})==0
            tshift =  pfb_all{ix}.tshift;
            
            %             pfb_all{ix}.sx = sx;
            %             pfb_all{ix}.sz = sz;
            %             if length(pfb_all{ix}.xpick) == 96
            %                 pfb_all{ix}.xpick = xseis;
            %                 pfb_all{ix}.zpick = zseis;
            %             else
            %                 return
            %             end
            
        else
            tshift = 0;
        end
        if pick>=1
            ind_tseis=find(abs(tseis-tMAX_seis_tmp)==min(abs(tseis-tMAX_seis_tmp))); % Get index of tMAX_seis
            seismomat2=seismomat(:,1:ind_tseis);
            tseis2=tseis(1:ind_tseis); % Keep seismo down to tMAX_seis
            while min(diff(tseis2))<dt_min % Downsample seismogram to enhance pick speed
                tseis2=tseis2(1:2:end);
                seismomat2=seismomat2(:,1:2:end);
            end
        end
        % Get shot acquisition setup param.
        dx1 = abs(diff(xseis));
        dx = median(dx1);
        xmin=sx-offMAX_pick-dx/2;
        xmax=sx+offMAX_pick+dx/2;
        if xmin<=min(xseis)
            xmin=min(xseis)-dx/2;
        end
        if xmax>=max(xseis)
            xmax=max(xseis)+dx/2;
        end
        
        if pick >= 2
            fprintf('\n  Automatic picking\n');
            [seismomat_test, tpick_auto, xpick_auto] = autopickfb(seismomat2,xseis,tseis2,sx,1,pick);
            if autogain == 1
                seismomat2 = seismomat_test;
            end
            tpick_auto(xpick_auto == sx) = 0;
            if pick == 3
                xpick = xpick_auto';
                tpick = tpick_auto;
            else
                xpick = xpick_auto;
                tpick = tpick_auto';
            end
        elseif pick == 1
            if autogain == 1
                seismomat2 = autopickfb(seismomat2,xseis,tseis2,sx,0,0);
            end
        end
    end
    
    if pick>=1
        percnext=perc;
        scalnext=scal;
        while isempty(percnext)==0 || isempty(scalnext)==0
            
            if pick < 3
                [fig,a1]=plot_wiggle(2,polarity*seismomat2',xseis,tseis2,scalnext,clip,percnext,fs,'Distance (m)','Time (ms)',...
                    [xmin xmax],[tMIN_pick tMAX_pick],[],[],[1 1 45 25],[]); % Plot seismo
                hold on;
                x_aerial_1 = [sx;xseis(xseis>=sx)];
                x_aerial_2 = [xseis(xseis<=sx);sx];
                ha1 = plot(x_aerial_1,tshift + (1000/343)*(x_aerial_1-sx),'c');
                ha2 = plot(x_aerial_2,tshift + (-1000/343)*(x_aerial_2-sx),'c');
                ha3 = plot(sx,tshift,'r.','linewidth',2,'markersize',20);
                a2 = magnify2(fig,5,0.25); hold on; % Magnifier
                
                set(fig,'name',['Shot at ',num2str(sx),' m'],'numbertitle','off');
                set(gcf,'units','normalized','outerposition',[0.3 0.5 1 0.66]); % Half screen right
                %                         set(gcf,'units','normalized','outerposition',[0.3 0 1 1]); % Full screen right
                %             set(gcf,'units','normalized','outerposition',[-0.02 -0.02 0.5 0.5]); % % 1/4 screen left top
                %             set(gcf,'units','normalized','outerposition',[-0.02 -0.02 0.6 1]); % % 1/2 screen left
                if isempty(pfb_all{ix})==0
                    xpick_prev=pfb_all{ix}.xpick;
                    tpick_prev=pfb_all{ix}.tpick + tshift;
                    for iii=1:length(xpick_prev)
                        diff_xpick = abs(xseis-xpick_prev(iii));
                        ind_diff = find(diff_xpick==min(diff_xpick),1,'first');
                        xpick_prev(iii) = xseis(ind_diff);
                    end
                    if length(xpick_prev)~=length(unique(xpick_prev))
                        [a,b]=hist(xpick_prev,unique(xpick_prev));
                        ind_rep = b(a>1);
                        if length(find(ismember(xpick_prev,ind_rep)))>1
                            for iii=1:length(find(ismember(xpick_prev,ind_rep)))/2
                                ind_del = find(ismember(xpick_prev,ind_rep(iii)),1,'last');
                                xpick_prev(ind_del) = [];
                                tpick_prev(ind_del) = [];
                            end
                        end
                    end
                else
                    if pick == 2
                        xpick_prev=xpick';
                        tpick_prev=tpick';
                    else
                        xpick_prev=[];
                        tpick_prev=[];
                    end
                end
                
                if ~isempty(tpick_prev)
                    if showplot==1
                        if hh(i)~=0
                            try
                                delete(hh(i));
                            catch
                            end
                        end
                        % Plot picked traveltimes
                        figure(3);
                        if sx_prev~=sx
%                             close(3);
                            f3=plot_curv(3,[],NaN,NaN,'.-','r',1,0,0,0,fs+4,'Distance (m)','Traveltime (ms)',[],...
                                [xMIN xMAX],[tMIN tMAX],[],xticks,tticks,[],[],[],[1 24 65 12],[]);
                        end
                        hold on;
                        hh(i) = plot(xpick_prev,tpick_prev - tshift,'r.-','linewidth',2,'markersize',15);
                    end
                    
                    % Plot previous picks
                    figure(2); hold on;
                    line(xpick_prev,tpick_prev,'color','b','linestyle','none','marker','x','markersize',8,'linewidth',2,'parent',a1);
                    line(xpick_prev,tpick_prev,'color','b','linestyle','none','marker','x','markersize',12,'linewidth',2,'parent',a2);
                    
                    % Pick window
                    figure(2); hold on;
                    if pick == 1
                        [xpick,tpick,closefig,shotprev,percnext,tshift,clean_air,scalnext]=matpickfb(xseis,tseis,seismomat,xpick_prev',tpick_prev',tshift,sx);
                    elseif pick == 2
                        [xpick,tpick,closefig,shotprev,percnext,tshift,clean_air,scalnext]=matpickfb(xseis,tseis,seismomat,xpick',tpick',tshift,sx);
                    end
                else
                    % Pick window
                    [xpick,tpick,closefig,shotprev,percnext,tshift,clean_air,scalnext]=matpickfb(xseis,tseis,seismomat,[],[],tshift,sx);
                end
                
                if length(tpick) == 1 && tpick == 0 && isempty(xpick) == 0
                    xpick = xpick_prev';
                    tpick = tpick_prev' - tshift;
                end
                
                if exist('f3','var')==1 && showplot==1
                    figure(3); hold on;
                    if hh(i)~=0
                        delete(hh(i));
                        hh(i) = 0;
                    end
                end
                
                if isempty(xpick)==0
                    % Plot picked traveltimes
                    [xpick,I] = sort(xpick);
                    tpick = tpick(I);
                    hh(i) = plot(xpick,tpick,'k.-');
                end
                
                if closefig==0
                    close(fig)
                end
                
                if shotprev==1 && i==1
                    fprintf('\n  No previous shot - Stay on first shot\n');
                end
            else
                percnext = []; scalnext = [];
            end
            
            if ~isempty(tpick)
                if sum(tpick)~=0
                    % Save matfile with picks
                    pfb_all{ix}.sx = sx;
                    pfb_all{ix}.sz = sz;
                    pfb_all{ix}.ix = ix;
                    pfb_all{ix}.xpick = xpick';
                    
                    [~,ind_test] = ismember(xpick,xseis);
                    pfb_all{ix}.zpick = zseis(ind_test);
                    pfb_all{ix}.tpick = tpick';
                    pfb_all{ix}.tshift = tshift;
                    pfb_all{ix}.trace = find(ismember(xseis,xpick));
                    save(pfbmatfile,'pfb_all');
                end
            else
                pfb_all{ix}=[];
            end
        end
    end
    
    if save_seis==1 || save_seis_pck==1
        if isempty(tMAX_seis)==1
            tMAX_seis = max(tseis);
            tMIN_seis = min(tseis);
        end
        if autogain == 1
            ind_tseis=find(abs(tseis-tMAX_seis)==min(abs(tseis-tMAX_seis))); % Get index of tMAX_seis
            seismomat2=seismomat(:,1:ind_tseis);
            tseis2=tseis(1:ind_tseis); % Keep seismo down to tMAX_seis
            while min(diff(tseis2))<dt_min % Downsample seismogram to enhance pick speed
                tseis2=tseis2(1:2:end);
                seismomat2=seismomat2(:,1:2:end);
            end
            seismomat2 = autopickfb(seismomat2,xseis,tseis2,sx,0,0);
        else
            seismomat2 = seismomat;
            tseis2 = tseis;
        end
    end
    
    % Plot raw seismo
    if save_seis==1
        fig1=plot_wiggle(0,polarity*seismomat2',xseis,tseis2,scal,clip,perc,fs,'Distance (m)','Time (ms)',...
            [xMIN xMAX],[tMIN_seis tMAX_seis],xticks,tseisticks,[0 0 seismo_size(1) seismo_size(2)],[]);
        file=fullfile(dir_img_seismo,[num2str(ix,shotnumformat),...
            '_seismo.',imgform]);
        %         hold on
        %         plot(Gxpick(Sxpick == 0),Tpick(Sxpick==0),'r.','linewidth',2,'markersize',15);
        
        save_fig(fig1,file,imgform,imgres,1);
        close(fig1)
    end
    
    % Plot picked seismo
    if save_seis_pck==1 && isempty(pfb_all{ix})==0
        fig2=plot_wiggle(0,polarity*seismomat2',xseis,tseis2-tshift,scal,clip,perc,fs,'Distance (m)','Time (ms)',...
            [xMIN xMAX],[tMIN_seis tMAX_seis],xticks,tseisticks,[0 0 seismo_size(1) seismo_size(2)],[]);
        xpick=pfb_all{ix}.xpick;
        tpick=pfb_all{ix}.tpick;
        
        if err_pc == 1
            err_ok = err_val*(tpick)/1000;
            if err_ok > err_val_max/1000
                err_ok = err_val_max/1000;
            elseif err_ok < err_val_min/1000
                err_ok = err_val_min/1000;
            end
        else
            err_ok = err_val/1000;
        end
        
        hold on;
        line(xpick,tpick,'color','r','linestyle','none','marker','+','markersize',2.5,'linewidth',0.8);
        line(xpick,tpick+err_ok*1000,'color','r','linestyle','--','marker','none','markersize',8,'linewidth',0.8);
        line(xpick,tpick-err_ok*1000,'color','r','linestyle','--','marker','none','markersize',8,'linewidth',0.8);
        hold off
        file=fullfile(dir_img_seismo,[num2str(ix,shotnumformat),...
            '_seismo_pick.',imgform]);
        save_fig(fig2,file,imgform,imgres,1);
        close(fig2)
    end
    
    % Plot offset vs Sx diagram
    if (save_src_off==1 || save_picks==1) && isempty(pfb_all{ix})==0
        newGx=pfb_all{ix}.xpick(:);
        newGz=pfb_all{ix}.zpick(:);
        newSx=bsxfun(@plus,newGx*0,pfb_all{ix}.sx);
        newSz=bsxfun(@plus,newGx*0,pfb_all{ix}.sz);
        Gxpick=[Gxpick;newGx];
        Gzpick=[Gzpick;newGz];
        Sxpick=[Sxpick;newSx];
        Szpick=[Szpick;newSz];
        Offpick=[Offpick;newGx(:)-newSx];
        Tpick=[Tpick;pfb_all{ix}.tpick(:)];
        
        if (pick==1 || pick==2) && save_src_off==1
            if plot_pos==1
                figure(1); hold on;
            end
            [fig3,ax3]=plot_scat(1,(Sxpick+Gxpick)/2,Offpick,Tpick,marker,markersize,1,map1,...
                axetop,0,1,fs,'CMP (m)','Offset (m)','Observed time (ms)',...
                [xMIN xMAX],[offMIN offMAX],[tMIN tMAX],xticks,offticks,tticks,[],[],[0 0 24 12],[],horex);
        end
    end
    if pick>=1
        fprintf('\n  **********************************************************');
        fprintf('\n  **********************************************************\n');
        sx_prev = sx;
    end
    
    % Next shot to be displayed
    if pick>=1 && exist('shotprev','var')==1 && shotprev==1 && i~=1
        ix=ix-2;
        i=i-2;
    elseif pick>=1 && exist('shotprev','var')==1 && shotprev==1 && i==1
        ix=ix-1;
        i=i-1;
    elseif pick>=1 && exist('shotprev','var')==1 && shotprev==-1
        break
    end
end

if exist('f3','var')==1 && showplot==1
%     close(f3);
end

if exist('fig3','var')==1
    close(fig3);
end

pickfile = fullfile(dir_start,[sufile(1:end-3),'.dat']);
if save_picks == 1
    % Save pick file for tomography code
    fidb=fopen(pickfile,'wt'); % File opening %
    ii=0;
    sources = unique([Sxpick,Szpick],'rows');
    nsx = size(sources,1);
    
    for isx = 1:nsx
        sx_tmp = sources(isx,1);
        sz_tmp = sources(isx,2);
        
        xpick_tmp = [];
        zpick_tmp = [];
        tpick_tmp = [];
        
        for ix=Sxselec
            if isempty(pfb_all{ix})==0
                if ismember([pfb_all{ix}.sx,pfb_all{ix}.sz],[sx_tmp, sz_tmp],'rows')
                    xpick_tmp = [xpick_tmp(:); pfb_all{ix}.xpick(:)];
                    zpick_tmp = [zpick_tmp(:); pfb_all{ix}.zpick(:)];
                    tpick_tmp = [tpick_tmp(:); pfb_all{ix}.tpick(:)];
                end
            end
        end
        
        xpick_single = [];
        zpick_single = [];
        tpick_single = [];
        single_traces = unique([xpick_tmp zpick_tmp],'rows');
        for igx = 1:length(single_traces)
            xpick_single(igx) = single_traces(igx,1);
            zpick_single(igx) = single_traces(igx,2);
            tpick_single(igx) = mean(tpick_tmp(xpick_tmp == xpick_single(igx) & zpick_tmp == zpick_single(igx)));
        end
        
        ii=ii+1;
        fprintf(fidb,'%d\t%f\t%f\t%f\t%d\n',0,sx_tmp,sz_tmp,ii,0);
        nbpicks=length(tpick_single);
        
        for j=1:nbpicks
            fprintf(fidb,'%f\t%f\t%f\t%f\t%f\n',ii,xpick_single(j),zpick_single(j),tpick_single(j)/1000,0.2*tpick_single(j)/1000);
        end
    end
    fclose(fidb);
    save(pfbmatfile,'pfb_all');
end

% if save_picks == 1
%     % Save pick file for tomography code
%     fidb=fopen(pickfile,'wt'); % File opening %
%     ii=0;
%     for ix=Sxselec
%         if isempty(pfb_all{ix})==0
%             ii=ii+1;
%             fprintf(fidb,'%d\t%f\t%f\t%f\t%d\n',0,pfb_all{ix}.sx,pfb_all{ix}.sz,ii,0);
%             nbpicks=length(pfb_all{ix}.tpick);
% %                         [xpick_tmp,I] = sort(pfb_all{ix}.xpick);
% %                         zpick_tmp = pfb_all{ix}.zpick(I);
% %                         tpick_tmp = pfb_all{ix}.tpick(I);
%             xpick_tmp = pfb_all{ix}.xpick;
%             zpick_tmp = pfb_all{ix}.zpick;
%             tpick_tmp = pfb_all{ix}.tpick;
%             for j=1:nbpicks
%                 fprintf(fidb,'%f\t%f\t%f\t%f\t%f\n',ii,xpick_tmp(j),zpick_tmp(j),tpick_tmp(j)/1000,0.2*tpick_tmp(j)/1000);
%             end
%         end
%     end
%     fclose(fidb);
%     save(pfbmatfile,'pfb_all');
% end

if pick==0
    fprintf('\n  **********************************************************');
    fprintf('\n  **********************************************************\n');
end

if exist(pickfile,'file') == 2
    input = load(pickfile); % Traveltime file
    if ~isempty(input)
        [Tpick, Sxpick, Szpick, Gxpick, Gzpick, Offpick] = readpicks(input);
    end
end

%%
pickfile_GIMLI = fullfile(dir_start,[sufile(1:end-3),'.sgt']);
if save_picks == 1
    save_pygimli(pickfile_GIMLI,Gx,Gz,Sx,Sz,Gxpick,Gzpick,Sxpick,Szpick,Tpick,err_pc,err_val,err_val_min,err_val_max);
end

%%
test_surf = [];
if save_src_off==1 && ~isempty(input)
    array = check_depth_array(Gx,Gz,Sx,Sz);
    for i = 1:4
        if ~isempty(array{i})
            
            IGall = ismember([Gx Gz],array{i}.G_sing,'rows');
            ISall = ismember([Sx Sz],array{i}.S_sing,'rows');
            IGpick = ismember([Gxpick Gzpick],array{i}.G_sing,'rows');
            ISpick = ismember([Sxpick Szpick],array{i}.S_sing,'rows');
            if sum(IGpick & ISpick)>0 && (i==1 || i==2)
                test_surf = 1;
            end
            
            % Plot offset vs Sx diagram
            if plot_pos == 1
                fig3 = figure; set(fig3,'Units','centimeters','position',[1 1 24 12]);
                if showplot == 0
                    set(fig3,'visible','off');
                end
                plot_curv(fig3,(Gx(IGall & ISall) + Sx(IGall & ISall))/2,Offset(IGall & ISall),[],...
                    '.','k',1,axetop,0,1,fs,'CMP (m)','Offset (m)','Observed time (ms)',...
                    [xMIN xMAX],[offMIN offMAX],[tMIN tMAX],xticks,offticks,tticks,[],[],[1 1 24 12],[],horex);
                hold on;
                plot_scat(fig3,(Gxpick(IGpick & ISpick) + Sxpick(IGpick & ISpick))/2,Offpick(IGpick & ISpick),Tpick(IGpick & ISpick),...
                    marker,markersize,1,map1,axetop,0,1,fs,'CMP (m)','Offset (m)','Observed time (ms)',...
                    get(gca,'xlim'),get(gca,'ylim'),[tMIN tMAX],xticks,offticks,tticks,[],[],[1 1 24 12],[],horex);
                if showplot == 0
                    set(fig3,'visible','off');
                end
            else
                fig3 = plot_scat(showplot,(Gxpick(IGpick & ISpick) + Sxpick(IGpick & ISpick))/2,Offpick(IGpick & ISpick),Tpick(IGpick & ISpick),...
                    marker,markersize,1,map1,axetop,0,1,fs,'CMP (m)','Offset (m)','Observed time (ms)',...
                    [xMIN xMAX],[offMIN offMAX],[tMIN tMAX],xticks,offticks,tticks,[],[],[1 1 24 12],[],horex);
            end
            file = fullfile(dir_img,sprintf('%s_cmp_offset_obs%i_raw.%s',sufile(1:end-3),i,imgform));
            save_fig(fig3,file,imgform,imgres,1);
            if showplot >= 1
                showplot = showplot+1;
            else
                close(fig3)
            end
        end
    end
end

%%
% Plot pseudo-velocity model
if save_src_off==1 && ~isempty(Sxpick) && length(unique(Sxpick))>1 && ~isempty(test_surf)
    % Get shot acquisition setup param.
    dx1 = abs(diff(Gx));
    dx = median(dx1);
    if min([Gx;Sx])~=0
        xshift = -min([Gx;Sx]);
    else
        xshift = 0;
    end
    input(:,2) = input(:,2) + xshift; % Shift X values
    topo = unique(sortrows([array{1}.G_sing; array{1}.S_sing],1),'rows');
    topo_new = interp1(topo(:,1),topo(:,2),unique([Gx;Sx]),'linear','extrap');
    topo = [unique([Gx;Sx]) topo_new];
    topo(:,1) = topo(:,1) + xshift; % Shift X values
    maxdepth = max([max(abs(Offpick))/3 max(abs(Sz - Gz))+4*dx]);
    [xx,zz,xzv,nodes,MP,input] = genmesh(input,topo,dx/2,[dx/2 dx],maxdepth,dl,1,10,1,0);
    
    [F2, pseudo_X, pseudo_Z, pseudo_V] = tt2velinit_new(Tpick, Sxpick + xshift, Szpick, Gxpick + xshift, Gzpick, 0, 0.2, 2);
    
    %     [F2, pseudo_X, pseudo_Z, pseudo_V] = tt2velinit_slope(array, input, xshift, 0);
    [xzvfinal, xi, zi, vi_filt] = vel2xzvtomo(xx,xzv,F2,[],[],0);
    
    %         fig3=plot_scat([],pseudo_X,pseudo_Z,pseudo_V,marker,markersize,1,map1,...
    %         axetop,0,1,fs,'Distance (m)','Altitude (m)',[velocity,' (m/s)'],...
    %         [xMIN xMAX],[zMIN zMAX],[vMIN vMAX],xticks,zticks,vticks,[],[],[1 1 24 12],[],1);
    %     hold on
    %     plot(topo(:,1),topo(:,2),'k-','linewidth',2);hold off;
    %                         return
    f1=plot_img(showplot,xi-xshift,zi,vi_filt,map2,axetop,0,1,fs,'Distance (m)',...
        'Altitude (m)',[velocity,' (m/s)'],[xMIN xMAX],[zMIN zMAX],...
        [vMIN vMAX],xticks,zticks,vticks,[],[],vISO,[1 20 24 12],[],vertex,blocky);
    sizeax=get(gca,'Position');
    if plottopo==1
        hold on
        plot(topo(:,1)-xshift,topo(:,2),'k-','linewidth',2);
    end
    file1=fullfile(dir_img,[lineName,'_',velocity,'pseudo','.',imgform]);
    save_fig(f1,file1,imgform,imgres,1);
    if showplot >= 1
        showplot = showplot+1;
    end
    
    filexzv=fullfile(dir_xzv,[lineName,'_pseudomodel.xzv']);
    dlmwrite(filexzv,[xi(:)-xshift,zi(:),vi_filt(:)],',');
end
%%
if save_src_off==1 && ~isempty(input)
    for i = 1:4
        if ~isempty(array{i})
            f4=plot_curv(showplot,[],NaN,NaN,'.-','k',1,0,0,0,fs+4,'Distance (m)','Traveltime (ms)',[],...
                [xMIN xMAX],[tMIN tMAX],[],xticks,tticks,[],[],[],[1 24 30 10],[]);
            % Plot picked traveltimes
            hold on;
            for ix=Sxselec
                if ~isempty(pfb_all{ix})
                    IGpick = ismember([pfb_all{ix}.xpick(:) pfb_all{ix}.zpick(:)],array{i}.G_sing,'rows');
                    ISpick = ismember([pfb_all{ix}.sx pfb_all{ix}.sz],array{i}.S_sing,'rows');
                    if sum(IGpick)>0 && sum(ISpick)>0
                        xpick = pfb_all{ix}.xpick(:);
                        tpick = pfb_all{ix}.tpick(:);
                        plot(xpick,tpick,'k.-','linewidth',1.5,'markersize',10);
                    end
                end
            end
            hold off
            file = fullfile(dir_img,sprintf('%s_traveltimes_obs%i.%s',sufile(1:end-3),i,imgform));
            %             file=fullfile(dir_img,[sufile(1:end-3),'_traveltimes_obs.',imgform]);
            save_fig(f4,file,imgform,imgres,1);
            if showplot >= 1
                showplot = showplot+1;
            end
        end
    end
end
%%
% Remove filtered SU file
if filt==1
    unix(['rm -f ',sufilefilt]);
end