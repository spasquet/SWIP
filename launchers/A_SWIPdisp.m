clear all; clc; close all;

%%% SURFACE-WAVE dispersion INVERSION & PROFILING (SWIP)
%%% MODULE A : SWIPdisp.m
%%% S. Pasquet - V17.01.18
%%% SWIPdisp.m performs windowing and stacking of surface-wave dispersion
%%% It allows to pick dispersion curves and save dispersion, spectrogram and seismograms
%%% Comment line to use default settings (cf SWIP_defaultsettings.m)

%%% Option calc=1 requires one SU file containing all shot gathers
%%% Mandatory SU headers: fldr, tracf, gx, sx, ns, dt
%%% Recommended SU headers: gelev, selev, scalco, scalel, gdel

%%%-------------------------%%%
%%% START OF INITIALIZATION %%%

%%% Main settings
Xmidselec   = [];         % Select Xmids (comment or [] to select all)
calc        = 1;          % SWIP subproject: new (=1), existing (=0) or import disp. curves (=2)

%%% Windowing and stacking settings (used if calc=1)
nWvec       = 31;         % Vector of window sizes (nb of traces)
dW          = 1;          % Shift between two successive windows (nb of traces)
dSmin       = 0;          % Min. offset between window and shot (nb of traces)
dSmax       = 20;         % Max. offset between window and shot (nb of traces)
side        = 'B';        % Sources side (B=both, L=left, R=right)

%%% P-omega transform settings (used if calc=1)
fmin        = 0;          % Min. frequency for p-omega transform (Hz)
fmax        = 100;        % Max. frequency for p-omega transform (Hz)
nray        = 500;        % Number of velocity samples for p-omega transform
vmin        = 0;          % Min. velocity for p-omega transform (m/s)
vmax        = 1500;       % Max. velocity for p-omega transform (m/s)

%%% Filter and mute settings (used if calc=1)
filt        = 0;          % Use band pass filter (=1) or not (=0)
fcutlow     = 1;          % Low freq. cut (Hz) (used if filt=1)
fcuthigh    = 100;        % High freq. cut (Hz) (used if filt=1)
taper       = 5;          % Taper window (Hz) (used if filt=1)
mute        = 0;          % Use mute (=1) or not (=0)
tmin1       = 0;          % Mute before tmin1 at shortest offset (s) (used if mute=1)
tmin2       = 0;          % Mute before tmin2 at largest offset (s) (used if mute=1)
tmax1       = 1;          % Mute after tmax1 at shortest offset (s) (used if mute=1)
tmax2       = 1;          % Mute after tmax2 at largest offset (s) (used if mute=1)

%%% Dispersion picking settings
pick        = 0;          % Pick dispersion: manual (=1), auto (=2) or not (=0)
mappick     = parula(39); % Colormap for picked dispersion image
mappicklog  = 0;          % Pseudo-log colorscale for dispersion (=1) or linear (=0)
dvmin       = 2;          % Min. phase velocity sampling (m/s)
modeinit    = 0;          % First picked propagation mode (0 = fundamental)
pickstyle   = 1;          % Semi-automatic picking (=1) or manual (=0)
smoothpick  = 1;          % Smooth picking (=1) or not (=0)

%%% Dispersion curves sampling settings
target      = 0;          % Convert picks to dinver target (=1) or not (=0)
wave        = 'R';        % Surface-wave type ('R' => Rayleigh or 'L' => Love)
maxmodeinv  = [];         % Max. mode number in target file for inversion (empty to include all picked modes)
sampling    = 1;          % Sampling in wavelength (1) or frequency (0)
resampvec   = (1:1:50);   % Resampling vector (wavelength [m] or frequency [Hz])
freqlim     = 1;          % Min. frequency defined with amplitude threshold (=1) or not (=0)
specampmin  = 0.01;       % Amplitude threshold (between 0 and 1) (used if freqlim=1)

%%% Error settings (used if target=1)
err         = 1;          % No error (=0), Lorentz error (=1) or percentage error (=2)
nWfac       = 1;          % Adapt Lorentz error with nW=nWfac*mean([nWmin,nWmax]) (used if err=1)
minerrvel   = 15;         % Min. vel. error for Lorentz (m/s) (used if err=1)
maxerrrat   = 0.5;        % Max. vel. error ratio for Lorentz (used if err=1)
sigma       = 15;         % Standard deviation error (%) (used if err=2)

%%% Plot settings
plotdisp    = 1;          % Save stacked dispersion images (=1) or not (=0)
plotpckdisp = 1;          % Save stacked dispersion images with picked curves (=1) or not (=0)
% plotspec    = 1;          % Save spectrograms (=1) or not (=0)
% plotseismo  = 1;          % Save seismograms (=1) or not (=0)
% plotsingle  = 1;          % Save single images (=1) or not (=0)
% plotstkdisp = 1;          % Save intermediate stacked dispersion images (=1) or not (=0)
% plot1dobs   = 1;          % Plot picked dispersion on 1D single graph (=1) or not (=0)
plot2dobs   = 1;          % Plot picked dispersion on 2D pseudo-section (=1) or not (=0)
showplot    = 0;          % Show plots before saving (=1) or not (=0)

%%% Figure display and output settings
imgform     = 'png';      % Fig. file format ('pdf', 'png', 'jpeg', 'tiff' or 'fig')
imgres      = 500;        % Fig. resolution (dpi) when saving as raster
fs          = 20;         % Fig. font size
cbpos       = 1;          % Colorbar on the right (=1) or at the bottom (=2)

%%% Dispersion, spectrograms and seismograms settings
Dlogscale   = 0;          % Pseudo-log colorscale for dispersion (=1) or linear (=0)
Flogscale   = 0;          % Logscale frequency axis (=1) or linear (=0)
axetop      = 0;          % Plot Xaxis on top (=1) or bottom (=0)
axerev      = 0;          % Yaxis pointing down (=1) or up (=0)
cb_disp     = 1;          % Plot dispersion colorbar (=1) or not (=0)
plotflim    = 1;          % Plot low cut frequency (=1) or not (=0)
plotlamlim  = 1;          % Plot max. lambda defined by resampvec (=1) or not (=0)
eb          = 1;          % Plot dispersion curves with errorbars (=1) or not (=0)
pickcol1    = 'w';        % Picks color for even (0, 2,...) modes number (cf ColorSpec)
pickcol2    = 'w';        % Picks color for odd (1, 3,...) modes number (cf ColorSpec)

map0        = bone(39);             % Colormap for dispersion images and spectrograms
fMIN        = 0;                    % Min. frequency to display (Hz)
fMAX        = fmax;                 % Max. frequency to display (Hz)
% fticks      = (fMIN:20:fMAX);       % Frequency ticks (Hz)
VphMIN      = 0;                    % Min. phase velocity to display (m/s)
VphMAX      = vmax;                 % Max. phase velocity to display (m/s)
% Vphticks    = (VphMIN:250:VphMAX);  % Phase velocity ticks (m/s) 
tMIN        = [];                   % Min. time to display (ms)
tMAX        = [];                   % Max. time to display (ms)
% tticks      = (tMIN:100:tMAX);      % Time ticks (ms)

%%% Phase velocity pseudo-section settings
map1        = haxby(32);            % Colormap for phase velocity
xMIN        = [];                   % Min. X (m)
xMAX        = [];                   % Max. X (m)
% xticks      = (xMIN:40:xMAX);       % X ticks (m)
lamMIN      = [];                   % Min. wavelength (m)
lamMAX      = [];                   % Max. wavelength (m)
% lticks      = (lamMIN:20:lamMAX);   % Wavelength ticks (m)
vphMIN      = [];                   % Min. phase velocity (m/s)
vphMAX      = [];                   % Max. phase velocity (m/s)
% vphticks    = (vphMIN:200:vphMAX);  % Phase velocity ticks (m/s)
vphISO      = [];                   % Phase velocity isocontours (m/s)

%%% END OF INITIALIZATION %%%
%%%-----------------------%%%

run('A_SWIPdisp_script');