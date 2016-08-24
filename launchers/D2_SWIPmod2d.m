clear all; clc; close all;

%%% SURFACE-WAVE dispersion INVERSION & PROFILING (SWIP)
%%% MODULE D2 : SWIPmod2d.m
%%% S. Pasquet - V16.6.28
%%% SWIPmod2d.m plots observed, calculated and residual 2D pseudo-sections
%%% It also plots Vp, Vs, Vp/Vs and Poisson's ratio 2D pseudo-sections
%%% Comment line to use default settings (cf defaultsettings.m)

%%%-------------------------%%%
%%% START OF INITIALIZATION %%%

Xmidselec  = [];    % Select Xmids (comment or [] to select all)

%%% Main settings
swip       = 1;     % Use vel. models from SWIP (=1) or from P- and SH-wave tomography (=2)

%%% SWIP model settings (used if swip=1)
modeltype  = 2;     % Plotted 1D Vs model (1=best, 2=layered, 3=smooth, 4=ridge)
avertype   = 0;     % Mean model (=0) or weighted model (=1)
nbest      = 0;     % Best models within error bars (=0) or nbest models (>0)
outpoints  = 0;     % No. of points allowed out of the error bars
usevptomo  = 0;     % Use Vp from tomography (=1) or from SW inversion (=0)
vpaver     = 0;     % Average Vp below window (=1) or extract at Xmid (=0)
dz         = 0.2;   % Sampling in depth (m)

%%% Plot and save settings
plot2dcal  = 1;     % Plot calculated Vphase 2D sections (=1) or not (=0)
plot2dmod  = 1;     % Plot models 2D sections (=1) or not (=0)
plot2dert  = 0;     % Plot electrical resistivity 2D section (=1) or not (=0)
showplot   = 0;     % Show plot before saving (=1) or not (=0)
savexzv    = 1;     % Save 2D models in .xzv ASCII file (=1) or not (=0)

%%% General display settings
imgform    = 'png'; % Fig. file format ('pdf', 'png', 'jpeg', 'tiff' or 'fig')
imgres     = 500;   % Fig. resolution when saving as raster
concat     = 1;     % Fig. concatenation panel (=1) or not (=0)

fs         = 20;    % Fig. font size
axetop     = 1;     % Plot Xaxis on top (=1) or bottom (=0)
cbpos      = 1;     % Colorbar on the right (=1) or at the bottom (=2)

%%% Phase velocity and residuals pseudo-section settings
map1       = haxby(32);              % Colormap for phase velocity
map4       = haxby(32);              % Colormap for phase velocity residuals
lamMIN     = [];                     % Min. wavelength (m)
lamMAX     = [];                     % Max. wavelength (m)
% lticks     = (lamMIN:20:lamMAX);     % Wavelength ticks (m)
vphMIN     = [];                     % Min. phase velocity (m/s)
vphMAX     = [];                     % Max. phase velocity (m/s)
% vphticks   = (vphMIN:200:vphMAX);    % Phase velocity ticks (m/s)
residMIN   = [];                     % Min. residual (m/s)
residMAX   = [];                     % Max. residual (m/s)
% residticks = (residMIN:20:residMAX); % Residual ticks (m/s)

%%% Vs, Vp, Vp/Vs, Poisson's ratio and Res. pseudo-section settings
blocky     = 0;     % Blocky (=0), smooth interp (=1) or smooth contour (=2) images
vertex     = 1;     % Vertical exageration
plottopo   = 1;     % Plot topo profile on 2D sections (=1) or not (=0)
plotDOI    = 0;     % Plot DOI estimated from wavelength (=1), from VsSTD (=2) or not (=0)
maskDOI    = 2;     % Mask below DOI from wavelength (=1), from VsSTD (=2) or not (=0)
doifact    = 0.66;  % DOI factor (DOI = Lmax*fact)
dpMAX      = [];    % Max. depth (m)

map5       = haxby(32);              % Colormap for Vp and Vs
map6       = haxby(32);              % Colormap for Vp/Vs and Poisson's ratio
map7       = haxby(32);              % Colormap for electrical resistivity

xMIN       = [];                     % Min. X (m)
xMAX       = [];                     % Max. X (m)
% xticks     = (xMIN:50:xMAX);         % X ticks (m)
zMIN       = [];                     % Min. altitude (m)
zMAX       = [];                     % Max. altitude (m)
% zticks     = (zMIN:20:zMAX);         % Altitude ticks (m)
vsMIN      = [];                     % Min. Vs (m/s)
vsMAX      = [];                     % Max. Vs (m/s)
% vsticks    = (vsMIN:500:vsMAX);      % Vs ticks (m/s)
vsISO      = [];                     % Vs isocontours (m/s)
vpMIN      = [];                     % Min. Vp (m/s)
vpMAX      = [];                     % Max. Vp (m/s)
% vpticks    = (vpMIN:1000:vpMAX);     % Vp ticks (m/s)
vpISO      = [];                     % Vp isocontours (m/s)
stdMIN     = [];                     % Min. StdVs (m/s)
stdMAX     = [];                     % Max. STdVs (m/s)
% stdticks   = (stdMIN:50:stdMAX);     % Vs STD ticks (m/s)
stdISO      = [];                     % STdVs isocontours (m/s)
vpvsMIN    = [];                     % Min. Vp/Vs
vpvsMAX    = [];                     % Max. Vp/Vs
% vpvsticks  = (vpvsMIN:vpvsMAX);      % Vp/Vs ticks
vpvsISO      = [];                     % Vp/Vs isocontours
poisMIN    = [];                     % Min. Poisson's ratio
poisMAX    = [];                     % Max. Poisson's ratio
% poisticks  = (poisMIN:0.1:poisMAX);  % Poisson's ratio ticks
poisISO      = [];                     % Poisson's ratio isocontours
resMIN     = [];                     % Min. resistivity (Ohm.m)
resMAX     = [];                     % Max. resistivity (Ohm.m)
% resticks   = [10 35 110 400];        % Resistivity ticks (Ohm.m)
resISO      = [];                     % Resistivity isocontours (Ohm.m)

%%% END OF INITIALIZATION %%%
%%%-----------------------%%%

run('D2_SWIPmod2d_script');
