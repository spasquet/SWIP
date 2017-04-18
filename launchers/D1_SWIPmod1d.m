clear all; clc; close all;

%%% SURFACE-WAVE dispersion INVERSION & PROFILING (SWIP)
%%% MODULE D1 : SWIPmod1d.m
%%% S. Pasquet - V17.04.14
%%% SWIPmod1d.m plots observed and calculated dispersion for each Xmid
%%% It also plots 1D Vp, Vs, Vp/Vs and Poisson's ratio models

%%% Option swip=1 plots velocity models and calculates theoretical
%%% dispersion from SWIP inversion results (module C)
%%% Requires to select a subproject folder Wmin_max.dWx.dSmin_max.side
%%% and an inversion folder located in file.inv
%%% If usevptomo=1, dispersion is calculated with VP from refraction tomography
%%% Requires to select a 3-column ASCII file (X, Z, Vp)
%%% Option tomo=1 plots velocity models and calculates theoretical
%%% dispersion from refraction tomography files
%%% Requires to select a subproject folder Wmin_max.dWx.dSmin_max.side
%%% and two 3-column ASCII files (X, Z, Vp) and (X, Z, Vs)
%%% Option user=1 plots velocity models and calculates theoretical
%%% dispersion from user-defined velocity structure
%%% Requires to select a subproject folder Wmin_max.dWx.dSmin_max.side

%%% Comment line to use default settings (cf SWIP_defaultsettings.m)

%%%-------------------------%%%
%%% START OF INITIALIZATION %%%

%%% Main settings
Xmidselec  = [];    % Select Xmids (comment or [] to select all)
swip       = 1;     % Plot vel. models from SWIP (=1) or not (=0)
tomo       = 0;     % Plot vel. models from P- and SH-wave tomography (=1) or not (=0)
user       = 0;     % Plot vel. models defined by user below (=1) or with input file (=2)

%%% SWIP model settings (used if swip=1)
modeltype  = 5;     % Plotted 1D Vs model
%%% (1=best, 2=averaged layered, 3=average smooth, 4=weighted layered, 5=weighted smooth, 6=ridge)
nbest      = 0;     % Best models within error bars (=0) or nbest models (>0)
outpoints  = 0;     % No. of points allowed out of the error bars
usevptomo  = 0;     % Forward calc. with Vp from tomo. (=1) or from SW inversion (=0)
dz         = 0.2;   % Sampling in depth (m)

%%% User defined 1D model parameters (used if user=1)
vpuser     = [500,1850,2500,4000];  % Vp (m/s)
vsuser     = [200,925,925,1250];    % Vs (m/s)
rhouser    = [1800,1800,1800,1800]; % Rho (kg/m3)
thkuser    = [0.85,10,25];          % Thickness (m)

%%% Toggle plots
plot1dcal  = 1;     % Plot dispersion images with calculated dispersion curves (=1) or not (=0)
plot1dmod  = 1;     % Plot Vs, Vp, Vp/Vs and Poisson's ratio 1D models (=1) or not (=0)
showplot   = 0;     % Show plots before saving (=1) or not (=0)

%%% Figure display and output settings
imgform    = 'png'; % Fig. file format ('pdf', 'png', 'jpeg', 'tiff' or 'fig')
imgres     = 500;   % Fig. resolution (dpi) when saving as raster
fs         = 20;    % Fig. font size
concat     = 1;     % Save indiv. and merged figures (=2), merged figure (=1) or indiv. figures (=0)

%%% Dispersion curves and images settings (used if plot1dcal = 1)
nmodemax   = 5;     % No. of mode for forward calculation
Dlogscale  = 0;     % Pseudo-log colorscale for dispersion (=1) or linear (=0)
Flogscale  = 0;     % Logscale frequency axis (=1) or linear (=0)
axetop     = 0;     % Plot Xaxis on top (=1) or bottom (=0)
axerev     = 0;     % Yaxis pointing down (=1) or up (=0)
cb_disp    = 1;     % Plot dispersion colorbar (=1) or not (=0)
plotflim   = 0;     % Plot low cut frequency (=1) or not (=0)
plotlamlim = 0;     % Plot max. lambda defined by resampvec (=1) or not (=0)
plot1dobs  = 1;     % Plot picked dispersion curves (=1) or not (=0)
eb         = 1;     % Plot dispersion curves with errorbars (=1) or not (=0)
pickcol1   = 'w';   % Picks color for even (0, 2,...) modes number (cf ColorSpec)
pickcol2   = 'w';   % Picks color for odd (1, 3,...) modes number (cf ColorSpec)

map0       = bone(39);              % Colormap for dispersion image
fMIN       = [];                    % Min. frequency (Hz)
fMAX       = [];                    % Max. frequency (Hz)
% fticks     = (fMIN:20:fMAX);        % Frequency ticks (Hz)
VphMIN     = [];                    % Min. phase velocity (m/s)
VphMAX     = [];                    % Max. phase velocity (m/s)
% Vphticks   = (VphMIN:250:VphMAX);   % Phase velocity ticks (m/s)

%%% Vs, Vp, Vp/Vs and Poisson's ratio 1D models settings (used if plot1dmod = 1)
plot1dstd  = 1;     % Plot error enveloppe (=1) or none (=0)
errstd     = 0;     % Percentage error on velocity models (in %) (0 for STD from SWIP)
plotDOI    = 2;     % Plot DOI estimated from wavelength (=1), from VsSTD (=2) or not (=0)
doifact    = 0.66;  % DOI factor (DOI = max_wavelength*doifact) (used if plotDOI=1)
stdMAX     = 150;   % Max. STdVs (m/s) to determine DOI (used if plotDOI=2)
plot1dvp   = 0;     % Plot 1D Vp models on same graph with 1D Vs models (=1) or not (=0)

dpMIN      = [];                    % Min. depth (m)
dpMAX      = [];                    % Max. depth (m)
% dticks     = (dpMIN:10:dpMAX);      % Depth ticks (m)
vsMIN      = [];                    % Min. Vs (m/s)
vsMAX      = [];                    % Max. Vs (m/s)
% vsticks    = (vsMIN:500:vsMAX);     % Vs ticks (m/s)
vpMIN      = [];                    % Min. Vp (m/s)
vpMAX      = [];                    % Max. Vp (m/s)
% vpticks    = (vpMIN:1000:vpMAX);    % Vp ticks (m/s)
vpvsMIN    = [];                    % Min. Vp/Vs
vpvsMAX    = [];                    % Max. Vp/Vs
% vpvsticks  = (vpvsMIN:vpvsMAX);     % Vp/Vs ticks
poisMIN    = [];                    % Min. Poisson's ratio
poisMAX    = [];                    % Max. Poisson's ratio
% poisticks  = (poisMIN:0.1:poisMAX); % Poisson's ratio ticks

%%% END OF INITIALIZATION %%%
%%%-----------------------%%%

run('D1_SWIPmod1d_script');