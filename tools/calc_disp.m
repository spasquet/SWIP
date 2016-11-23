function calc_disp(sufile,nWvec,dSmin,nray,fmin,fmax,vmin,vmax,flip,xsca,tsca,imgform,imgres)

%%% S. Pasquet - V16.11.22
% Quick extraction of dispersion image from seismogram
% calc_disp(sufile,nWvec,dSmin,nray,fmin,fmax,vmin,vmax,flip,xsca,tsca,imgform,imgres)

if exist('sufile','var')==0 || isempty(sufile)==1
    [sufile,supath]=uigetfile('*.su','Select seismogram file');
    if sufile==0
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        fprintf('\n   Please select a seismogram file');
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
        return
    end
else
    supath=fileparts(sufile);
end
fprintf('\n  Reading %s\n',sufile);

acquiparam=get_acquiparam(sufile,[]);

if exist('xsca','var')==0 || isempty(xsca)==1
    xsca=acquiparam.xsca;
end

run('SWIP_defaultsettings');

dSmax=dSmin;
dir_dat_xmid=supath;
dx=acquiparam.dx; % Mean inter-geophone spacing (m)
xmidformat=['%12.',num2str(log(xsca)/log(10)),'f']; % Precision

if dx==0
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n   Please provide inter-geophone spacing');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
    return
end
topo=acquiparam.topo; % Get topography (X,Z in m)
Gxsing=acquiparam.Gxsing; % Single geophones positions
Sxsing=acquiparam.Sxsing; % Single sources positions
if length(Sxsing)>1
    fprintf('\n  Please select single shot SU file\n');
    return
end

xmin=min(Gxsing); % Get starting X coordinate (m)
xmax=max(Gxsing); % Get ending X coordinate (m)
winsize=nWvec;
maxwinsize=(winsize-1)*dx;
nwin=length(winsize);

if mod(nWvec,2)==1
    XmidT=Gxsing(1+(nWvec-1)/2:dW:end-(nWvec-1)/2);
else
    XmidT=mean([Gxsing((nWvec)/2:dW:end-(nWvec)/2),...
        Gxsing(1+(nWvec)/2:dW:1+end-(nWvec)/2)],2);
end
XmidT=round(XmidT'*xsca)/xsca;
Xlength=length(XmidT);

for ix=1:Xlength
    j=0; % Stack flag
    for jw=1:nwin
        % Retrieve first and last geophone position for the current window
        if mod(nWvec,2)==1 % Non-even number of traces
            Gleft=Gxsing(find(Gxsing<XmidT(ix),(winsize(jw)-1)/2,'last'));
            Gright=Gxsing(find(Gxsing>XmidT(ix),(winsize(jw)-1)/2));
            ntr=length(Gleft)+length(Gright)+1;
            if ntr~=winsize(jw) && calc==1 % Check number of extracted traces
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
            if ntr~=winsize(jw) && calc==1 % Check number of extracted traces
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
        % Get nb of selected shot for the current window
        nshot(ix,jw)=length(Sselec);
        
        if nshot(ix,jw)>0
            fprintf(['\n  Xmid = ',num2str(XmidT(ix),xmidformat),' m \n']);
            fprintf(['\n  Shot at ',num2str(Sselec(1)),...
                ' m with nW = ',num2str(winsize(jw)),' and dS = ',num2str(dSmin),'\n']);
        end
        
        %%%%%% Loop over all selected shots %%%%%%
        
        for ks=1:nshot(ix,jw)
            % Windowing, muting and saving seismogram in .su file
            seismofile=fullfile(dir_dat_xmid,[sufile(1:end-3),'_',num2str(XmidT(ix),xmidformat),'.',...
                num2str(winsize(jw)),'.',num2str(Sselec(ks)),'.su']);
            [seismomat,xseis,tseis,ntr]=matwind(sufile,Sselec(ks),Gmin(ix,jw),Gmax(ix,jw),xsca,...
                winsize(jw),seismofile,0,mute,tmin1,tmin2,tmax1,tmax2);
            if ntr~=winsize(jw) % Check number of extracted traces
                fprintf('\n  Not enough traces - Go to next shot\n');
                nshot(ix,jw)=nshot(ix,jw)-1;
                continue
            end
            j=j+1; % Stack flag
            % P-Omega transform on seismogram and saving in .dsp file
            dspfile=fullfile(dir_dat_xmid,[sufile(1:end-3),'_',num2str(XmidT(ix),xmidformat),'.',...
                num2str(winsize(jw)),'.',num2str(Sselec(ks)),'.dsp']);
            [dspmat,f,v]=matpomegal(seismofile,1,nray,fmin,fmax,vmin,vmax,...
                flip,xsca,tsca,1,dspfile,0);
            unix('rm -f speclud');
            % Spectrogram calculation on seismogram and saving in .spec file
            specfile=fullfile(dir_dat_xmid,[sufile(1:end-3),'_',num2str(XmidT(ix),xmidformat),'.',...
                num2str(winsize(jw)),'.',num2str(Sselec(ks)),'.spec']);
            [specmat,fspec,xspec]=matspecfx(seismofile,xsca,specfile,0,0);
        end
    end
end

plot_disp(dspfile,0,Dlogscale,0,max(f(:)),imgform,imgres,flip,eb)
plot_spec(specfile,0,max(f(:)),imgform,imgres)
plot_seismo(seismofile,0,750,imgform,imgres)

end