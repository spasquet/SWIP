clear all; clc; close all;

%%% SURFACE-WAVE dispersion INVERSION & PROFILING (SWIP)
%%% MODULE D1 : SWIPmod1d.m
%%% S. Pasquet - V16.8.3
%%% SWIPmod1d.m plots observed and calculated dispersion for each Xmid
%%% It also plots 1D Vp, Vs, Vp/Vs and Poisson's ratio models
%%% Comment line to use default settings (cf defaultsettings.m)

%%%-------------------------%%%
%%% START OF INITIALIZATION %%%

Xmidselec  = [];    % Select Xmids (comment or [] to select all)

%%% Main settings
swip       = 1;     % Use vel. models from SWIP (=1) or not (=0)
tomo       = 0;     % Use vel. models from P- and SH-wave tomography (=1) or not (=0)
user       = 0;     % Use vel. models defined by user below (=1) or with gpdc file (=2)

%%% SWIP model settings (used if swip=1)
modeltype  = 2;     % Plotted 1D Vs model (1=best, 2=layered, 3=smooth, 4=ridge)
avertype   = 0;     % Mean model (=0) or weighted model (=1)
nbest      = 0;     % Best models within error bars (=0) or nbest models (>0)
outpoints  = 0;     % No. of points allowed out of the error bars
usevptomo  = 0;     % Use Vp from tomography (=1) or from SW inversion (=0)
vpaver     = 0;     % Average Vp below window (=1) or extract at Xmid (=0)
dz         = 0.2;   % Sampling in depth (m)

%%% User defined 1D model parameters (used if user=1)
vpuser     = [500,1850,2500,4000];  % Vp (m/s)
vsuser     = [200,925,925,1250];    % Vs (m/s)
rhouser    = [1800,1800,1800,1800]; % Rho (kg/m3)
thkuser    = [0.85,10,25];          % Thickness (m)

%%% Plot settings
plot1dcal  = 1;     % Plot dispersion images with calculated dispersion curves (=1) or not (=0)
plot1dmod  = 1;     % Plot Vs, Vp, Vp/Vs and Poisson's ratio 1D models (=1) or not (=0)
showplot   = 0;     % Show plots before saving (=1) or not (=0)

%%% General display settings
imgform    = 'png'; % Fig. file format ('pdf', 'png', 'jpeg', 'tiff' or 'fig')
imgres     = 500;   % Fig. resolution (dpi) when saving as raster
concat     = 1;     % Fig. concatenation panel (=1) or not (=0)

fs         = 20;    % Fig. font size

%%% Dispersion curves and images settings
nmodemax   = 5;     % No. of mode for forward calculation
Dlogscale  = 1;     % Log colorscale for dispersion (=1) or linear (=0)
Flogscale  = 0;     % Logscale frequency axis (=1) or linear (=0)
axetop     = 0;     % Plot Xaxis on top (=1) or bottom (=0)
axerev     = 0;     % Yaxis pointing down (=1) or up (=0)
cb_disp    = 0;     % Plot colorbar (=1) or not (=0)
plot1dobs  = 1;     % Plot picked dispersion curves (=1) or not (=0)
eb         = 1;     % Plot dispersion curves with errorbars (=1) or not (=0)
pickcol1   = 'w';   % Picks color for even (0, 2,...) modes number (cf ColorSpec)
pickcol2   = 'w';   % Picks color for odd (1, 3,...) modes number (cf ColorSpec)
plotflim   = 0;     % Plot low cut frequency (=1) or not (=0)

map0       = bone(39);              % Colormap for dispersion image
fMIN       = [];                    % Min. frequency (Hz)
fMAX       = [];                    % Max. frequency (Hz)
% fticks     = (fMIN:20:fMAX);        % Frequency ticks (Hz)
VphMIN     = [];                    % Min. phase velocity (m/s)
VphMAX     = [];                    % Max. phase velocity (m/s)
% Vphticks   = (VphMIN:250:VphMAX);   % Phase velocity ticks (m/s)

%%% Vs, Vp, Vp/Vs and Poisson's ratio 1D models settings
plot1dstd  = 1;     % Plot STD enveloppe (=1), % error (=2) or none (=0)
errstd     = 0;     % Percentage error on velocity models (in %) (0 for STD from SWIP)
plotDOI    = 2;     % Plot DOI estimated from wavelength (=1), from VsSTD (=2) or not (=0)
doifact    = 0.66;  % DOI factor (DOI = Lmax*fact) (used if plotDOI=1)
stdMAX     = 200;   % Max. STdVs (m/s) (used if plotDOI=2)
plot1dvp   = 0;     % Plot 1D Vp models with 1D Vs models (=1) or not (=0)

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
