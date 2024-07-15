
%%% S. Pasquet - V18.10.15

%% A_TOMOpickfb settings
% Main settings
% Pick first arrival (=1) or not (=)
if exist('pick','var')==0 || isempty(pick)==1
    pick=1;
end
% Save seismogram image (=1) or not (=)
if exist('save_seis','var')==0 || isempty(save_seis)==1
    save_seis=0;
end
% Save picked seismogram image (=1) or not (=)
if exist('save_seis_pck','var')==0 || isempty(save_seis_pck)==1
    save_seis_pck=0;
end
% Save source vs offset diagram (=1) or not (=)
if exist('save_src_off','var')==0 || isempty(save_src_off)==1
    save_src_off=1;
end
% Save picks in .dat file (=1) or not (=)
if exist('save_picks','var')==0 || isempty(save_picks)==1
    save_picks=1;
end

% Picking and filtering settings
% Time shift to correct from delayed triggering (ms)
if exist('tshift','var')==0 || isempty(tshift)==1
    tshift=0;
end
% Minimum time sample to display seismogram for picking (ms)
if exist('dt_min','var')==0 || isempty(dt_min)==1
    dt_min=[];
end
% Time correction from delayed triggering (=1) or not (=0)
if exist('time_cor','var')==0 || isempty(time_cor)==1
    time_cor=0;
end
% Use band pass filter (=1) or not (=0)
if exist('filt','var')==0 || isempty(filt)==1
    filt=0;
end
% Low freq. cut (Hz) (used if filt=1)
if exist('fcutlow','var')==0 || isempty(fcutlow)==1
    fcutlow=1;
end
% High freq. cut (Hz) (used if filt=1)
if exist('fcuthigh','var')==0 || isempty(fcuthigh)==1
    fcuthigh=100;
end
% Apodisation taper (Hz) (used if filt=1)
if exist('taper','var')==0 || isempty(taper)==1
    taper=5;
end

% Error settings
if ~exist('err_pc','var') || isempty(err_pc)
    err_pc = 1;
end
if ~exist('err_val','var') || isempty(err_val)
    if err_pc == 1
        err_val = 0.1;
    else 
        err_val = 1;
    end
end
if ~exist('err_val_min','var') || isempty(err_val_min)
    err_val_min = 1;
end
if ~exist('err_val_max','var') || isempty(err_val_max)
    err_val_max = 10;
end

% General display settings
% Fig. file format ('pdf', 'png', 'jpeg', 'tiff' or 'fig')
if exist('imgform','var')==0 || isempty(imgform)==1
    imgform='png';
end
% Fig. resolution (dpi) when saving as raster
if exist('imgres','var')==0 || isempty(imgres)==1
    imgres=500;
end
% Fig. font size (reduce by 40% with fs=40 => fs=16 in AI) (max=40)
if exist('fs','var')==0 || isempty(fs)==1
    fs=20;
end
% Plot Xaxis on top (=1) or bottom (=0)
if exist('axetop','var')==0 || isempty(axetop)==1
    axetop=1;
end
% Plot colorbar on the right (=1) or bottom (=2)
if exist('cb_pos','var')==0 || isempty(cb_pos)==1
    cb_pos=1;
end
% Size of seismograms
if exist('seismo_size','var')==0 || isempty(seismo_size)==1
    seismo_size=[24 8];
end

% Seismogram plot settings
% Min. time (ms)
if (exist('tMIN_seis','var')==0 || isempty(tMIN_seis)==1) || (exist('tMAX_seis','var')==0 || isempty(tMAX_seis)==1)
    tMIN_seis=[];
end
% Max. time (ms)
if (exist('tMIN_seis','var')==0 || isempty(tMIN_seis)==1) || (exist('tMAX_seis','var')==0 || isempty(tMAX_seis)==1)
    tMAX_seis=[];
end
% Time ticks (ms)
if exist('tseisticks','var')==0 || isempty(tseisticks)==1
    tseisticks=[];
end
% Max. offset (m)
if (exist('offMAX_pick','var')==0 || isempty(offMAX_pick)==1)
    offMAX_pick=[];
end
% Max. time (ms)
if (exist('tMAX_pick','var')==0 || isempty(tMAX_pick)==1)
    tMAX_pick=[];
end
% Clipping percentile (%)
if (exist('perc','var')==0 || isempty(perc)==1)
    perc=80;
end
% Amplitude scaling factor
if (exist('scal','var')==0 || isempty(scal)==1)
    scal=5;
end
% Clip (=1) or not (=0)
if (exist('clip','var')==0 || isempty(clip)==1)
    clip=1;
end
% Fill positive (=1) or negative (-1) amplitudes
if (exist('polarity','var')==0 || isempty(polarity)==1)
    polarity=-1;
end
% Automatic gain (=1) or not (=0)
if (exist('autogain','var')==0 || isempty(autogain)==1)
    autogain=0;
end
% Normalize (=1) or not (=0)
if (exist('normalize','var')==0 || isempty(normalize)==1)
    normalize=0;
end

% Source vs offset diagram settings
% Plot Xaxis on top (=1) or bottom (=0)
if exist('axetop','var')==0 || isempty(axetop)==1
    axetop=1;
end
% Plot traces positions (=1) or not (=0)
if exist('plot_pos','var')==0 || isempty(plot_pos)==1
    plot_pos=1;
end
% Colormap for traveltimes
if exist('map1','var')==0 || isempty(map1)==1
    map1=haxby(32);
end
% Marker type ('s'=square, 'o'=circle)
if exist('marker','var')==0 || isempty(marker)==1
    marker='s';
end
% Marker size
if exist('markersize','var')==0 || isempty(markersize)==1
    markersize=30;
end

% Min. X (m)
if (exist('xMIN','var')==0 || isempty(xMIN)==1) || (exist('xMAX','var')==0 || isempty(xMAX)==1)
    xMIN=[];
end
% Max. X (m)
if (exist('xMIN','var')==0 || isempty(xMIN)==1) || (exist('xMAX','var')==0 || isempty(xMAX)==1)
    xMAX=[];
end
% X ticks (m)
if exist('xticks','var')==0 || isempty(xticks)==1
    xticks=[];
end
% Min. offset (m)
if (exist('offMIN','var')==0 || isempty(offMIN)==1) || (exist('offMAX','var')==0 || isempty(offMAX)==1)
    offMIN=[];
end
% Max. offset (m)
if (exist('offMIN','var')==0 || isempty(offMIN)==1) || (exist('offMAX','var')==0 || isempty(offMAX)==1)
    offMAX=[];
end
% Offset ticks (m)
if exist('offticks','var')==0 || isempty(offticks)==1
    offticks=[];
end
% Min. traveltime (ms)
if (exist('tMIN','var')==0 || isempty(tMIN)==1) || (exist('tMAX','var')==0 || isempty(tMAX)==1)
    tMIN=[];
end
% Max. traveltime (ms)
if (exist('tMIN','var')==0 || isempty(tMIN)==1) || (exist('tMAX','var')==0 || isempty(tMAX)==1)
    tMAX=[];
end
% Traveltime ticks (ms)
if exist('tticks','var')==0 || isempty(tticks)==1
    tticks=[];
end

%% B_TOMOinvmc settings
% Main settings
% Run inversion (=1) or not (=0)
if exist('inversion','var')==0 || isempty(inversion)==1
    inversion=1;
end
% Save seismogram image (=1) or not (=0)
if exist('plot_all','var')==0 || isempty(plot_all)==1
    plot_all=1;
end
% Save picked seismogram image (=1) or not (=0)
if exist('plot_final','var')==0 || isempty(plot_final)==1
    plot_final=1;
end
% Show plots during inversion (=1) or not (=0)
if exist('showplot','var')==0 || isempty(showplot)==1
    showplot=0;
end

% Inversion settings (used if inversion=1)
% Cell width (m)
if exist('dx','var')==0 || isempty(dx)==1
    dx=1;
end
% Cell height (m)
if exist('dz','var')==0 || isempty(dz)==1
    dz=[dx 2*dx];
end
% Max. depth of model (m)
if exist('maxdepth','var')==0 || isempty(maxdepth)==1
    maxdepth=50;
end
% Raytracing node spacing (m)
if exist('dl','var')==0 || isempty(dl)==1
    dl=dx;
end
% Refining factor for average model update (i.e. dx=dx/ref_fac)
if exist('ref_fac','var')==0 || isempty(ref_fac)==1
    ref_fac=2;
end

% Starting model
if exist('nstart','var')==0 || isempty(nstart)==1
    nstart=1;
end
% Nb of model for MC starting model generation
if exist('nModels','var')==0 || isempty(nModels)==1
    nModels=50;
end
% Number of iteration per inversion
if exist('niter','var')==0 || isempty(niter)==1
    niter=15;
end
% Stop after nitlim iterations if RMS > rms1
if exist('nitlim','var')==0 || isempty(nitlim)==1
    nitlim=7;
end
% Stop if RMS <= frms
if exist('frms','var')==0 || isempty(frms)==1
    frms=2.5;
end
% Change regularization when rms1 is reached
if exist('rms1','var')==0 || isempty(rms1)==1
    rms1=8.5;
end
% Change regularization when rms2 is reached
if exist('rms2','var')==0 || isempty(rms2)==1
    rms2=4.5;
end
% Keep models below rmslim for final average model
if exist('rmslim','var')==0 || isempty(rmslim)==1
    rmslim=4;
end

% Random generation of initial model (=1) or not (=0)
if exist('mc_init','var')==0 || isempty(mc_init)==1
    mc_init=0;
end
% Random distribution (=1) or gaussian (=0) of starting models
if exist('rand_distrib','var')==0 || isempty(rand_distrib)==1
    rand_distrib=0;
end
% Min. starting gradient for initial model (m/s/m)
if exist('grad_min','var')==0 || isempty(grad_min)==1
    grad_min=50;
end
% Max. starting gradient for initial model (m/s/m)
if exist('grad_max','var')==0 || isempty(grad_max)==1
    grad_max=200;
end
% Min. starting velocity for initial model (m/s)
if exist('start_vel_min','var')==0 || isempty(start_vel_min)==1
    start_vel_min=100;
end
% Max. starting velocity for initial model (m/s)
if exist('start_vel_max','var')==0 || isempty(start_vel_max)==1
    start_vel_max=1000;
end
% Min. velocity for velocity update (m/s)
if exist('min_vel','var')==0 || isempty(min_vel)==1
    min_vel=100;
end
% Max. velocity for velocity update (m/s)
if exist('max_vel','var')==0 || isempty(max_vel)==1
    max_vel=4000;
end

% Use layered starting model (=1) or not (=0)
if exist('lay_mod','var')==0 || isempty(lay_mod)==1
    lay_mod=0;
end

% Use reference starting model (=1) or not (=0)
if exist('ref_mod','var')==0 || isempty(ref_mod)==1
    ref_mod=0;
end

% Random generation of regularisation (=1) or not (=0)
if exist('mc_regul','var')==0 || isempty(mc_regul)==1
    mc_regul=0;
end
%%% Regularization settings
%%% h# = regularization depth dependency
%%% invpar# = Tikhonov parameters [0th 1st_horizontal 1st_vertical 2nd_horizontal 2nd_vertical 0]
% Initial regularization settings
if exist('h0','var')==0 || isempty(h0)==1
    h0=@(rr)(rr);
end
if exist('invpar0','var')==0 || isempty(invpar0)==1
    invpar0=[40 200 200 200 200  0];
end
% Regularization settings when RMS<rms1
if exist('h1','var')==0 || isempty(h1)==1
    h1=@(rr)(rr.^3);
end
if exist('invpar1','var')==0 || isempty(invpar1)==1
    invpar1=[1 120 120 120 120 0];
end
% Regularization settings when RMS<rms1
if exist('h2','var')==0 || isempty(h2)==1
    h2=@(rr)(rr.^0.3);
end
if exist('invpar2','var')==0 || isempty(invpar2)==1
    invpar2=[6 80 80 80 80 0];
end

% Velocity model plot settings (used if plot_all=1 or plot_final=1)
% % Velocity type ('Vp' or 'Vs')
if exist('velocity','var')==0 || isempty(velocity)==1
    velocity='Vp';
end
% Vertical exageration
if exist('vertex','var')==0 || isempty(vertex)==1
    vertex=1;
end
% Blocky model
if exist('blocky','var')==0 || isempty(blocky)==1
    blocky=2;
end
% Plot topo profile on 2D sections (=1) or not (=0)
if exist('plottopo','var')==0 || isempty(plottopo)==1
    plottopo=1;
end
% Plot borehole (=1) or not (=0)
if exist('plotBH','var')==0 || isempty(plotBH)==1
    plotBH=0;
end
% Colormap for velocity
if exist('map2','var')==0 || isempty(map2)==1
    map2=haxby(32);
end
% Colormap for residuals
if exist('map3','var')==0 || isempty(map3)==1
    map3=polarmap(32);
end
% Colormap for standard deviation
if exist('map4','var')==0 || isempty(map4)==1
    map4=flipud(inferno(32));
end
% Masking type (0 no mask, 1 raymask, 2 stdmask1, 3 stdmask2)
if exist('mask','var')==0 || isempty(mask)==1
    mask=2;
end
% Keep masked values as NaN in .xzv files
if exist('keep_nans','var')==0 || isempty(keep_nans)==1
    keep_nans=1;
end
% Auto mask below mask_depth
if exist('mask_depth','var')==0 || isempty(mask_depth)==1
    mask_depth=[];
end

% Min altitude (m)
if (exist('zMAX','var')==0 || isempty(zMAX)==1) || (exist('zMIN','var')==0 || isempty(zMIN)==1)
    zMIN=[];
end
% Max altitude (m)
if (exist('zMAX','var')==0 || isempty(zMAX)==1) || (exist('zMIN','var')==0 || isempty(zMIN)==1)
    zMAX=[];
end
% Altitude ticks (m)
if exist('zticks','var')==0 || isempty(zticks)==1
    zticks=[];
end
% Min. velocity (m/s)
if (exist('vMIN','var')==0 || isempty(vMIN)==1) || (exist('vMAX','var')==0 || isempty(vMAX)==1)
    vMIN=[];
end
% Max velocity (m/s)
if (exist('vMIN','var')==0 || isempty(vMIN)==1) || (exist('vMAX','var')==0 || isempty(vMAX)==1)
    vMAX=[];
end
% Velocity ticks (m/s)
if exist('vticks','var')==0 || isempty(vticks)==1
    vticks=[];
end
% Min. velocity standard deviation (m/s)
if (exist('stdMIN','var')==0 || isempty(stdMIN)==1) || (exist('stdMAX','var')==0 || isempty(stdMAX)==1)
    stdMIN=[];
end
% Max. velocity standard deviation (m/s)
if (exist('stdMIN','var')==0 || isempty(stdMIN)==1) || (exist('stdMAX','var')==0 || isempty(stdMAX)==1)
    stdMAX=[];
end
% Velocity STD ticks (m/s)
if exist('stdticks','var')==0 || isempty(stdticks)==1
    stdticks=[];
end

% Min. residual (ms)
if (exist('residMIN','var')==0 || isempty(residMIN)==1) || (exist('residMAX','var')==0 || isempty(residMAX)==1)
    residMIN=[];
end
% Max. residual (ms)
if (exist('residMIN','var')==0 || isempty(residMIN)==1) || (exist('residMAX','var')==0 || isempty(residMAX)==1)
    residMAX=[];
end
% Residual ticks
if exist('residticks','var')==0 || isempty(residticks)==1
    residticks=[];
end

% Velocity isocontours (m/s)
if exist('vISO','var')==0 || isempty(vISO)==1
    vISO=[]; 
end
% STdVs isocontours (m/s)
if exist('stdISO','var')==0 || isempty(stdISO)==1
    stdISO=[]; 
end
% Horizontal exageration
if exist('horex','var')==0 || isempty(horex)==1
    horex=1;
end
% Mask over std_mask (%)
if exist('std_mask','var')==0 || isempty(std_mask)==1
    std_mask=10;
end
% Reverse X axis (=1) or not (=0)
if exist('reverse_x','var')==0 || isempty(reverse_x)==1
    reverse_x=0;
end