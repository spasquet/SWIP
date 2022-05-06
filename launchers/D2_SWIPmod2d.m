clear all; clc; close all;

%%% SURFACE-WAVE dispersion INVERSION & PROFILING (SWIP)
%%% MODULE D2 : SWIPmod2d.m
%%% S. Pasquet - V20.02.17
%%% SWIPmod2d.m plots observed, calculated and residual pseudo-sections
%%% It also plots Vp, Vs, Vp/Vs, Poisson's ratio and auxiliary data 2D sections

%%% Option input_vel=1 plots velocity models and calculates theoretical
%%% dispersion from SWIP inversion results (module C)
%%% Requires to select a subproject folder Wmin_max.dWx.dSmin_max.side
%%% and an inversion folder located in file.inv
%%% If usevptomo=1, dispersion is calculated with VP from refraction tomography
%%% Requires to select a 3-column ASCII file (X, Z, Vp)
%%% Option input_vel=2 plots velocity models and calculates theoretical
%%% dispersion from refraction tomography files
%%% Requires to select a subproject folder Wmin_max.dWx.dSmin_max.side
%%% and two 3-column ASCII files (X, Z, Vp) and (X, Z, Vs - optional)
%%% Option input_aux=1 plots auxiliary data
%%% Requires to select a 3-column ASCII files (X, Z, auxdata)

%%% Comment line to use default settings (cf SWIP_defaultsettings.m)

%%%-------------------------%%%
%%% START OF INITIALIZATION %%%

%%% Main settings
Xmidselec  = [];    % Select Xmids (comment or [] to select all)
input_vel  = 1;     % Input vel. models from SWIP (=1) or from P- and SH-wave tomography (=2)
input_aux  = 0;     % Input auxiliary 2D data (=1) or not (=0)

%%% SWIP model settings (used if input_vel=1)
modeltype  = 3;     % Plotted 1D Vs model
%%% (1=best, 2=averaged layered, 3=average smooth, 4=weighted layered, 5=weighted smooth, 6=ridge)
nbest      = 0;     % Best models within error bars (=0) or nbest models (>0)
outpoints  = 0;     % No. of points allowed out of the error bars
usevptomo  = 0;     % Forward calc. with Vp from tomo. (=1) or from SW inversion (=0)
dz         = 0.2;   % Sampling in depth (m)

%%% Toggle plots and data save
plot2dcal  = 1;     % Plot Vphase and residuals 2D sections (=1) or not (=0)
plot2dmod  = 1;     % Plot models 2D sections (=1) or not (=0)
showplot   = 0;     % Show plot before saving (=1) or not (=0)
savexzv    = 1;     % Save 2D models in .xzv ASCII file (=1) or not (=0)

%%% Figure display and output settings
imgform    = 'png'; % Fig. file format ('pdf', 'png', 'jpeg', 'tiff' or 'fig')
imgres     = 200;   % Fig. resolution when saving as raster
fs         = 16;    % Fig. font size
concat     = 1;     % Save indiv. and merged figures (=2), merged figure (=1) or indiv. figures (=0)
cbpos      = 1;     % Colorbar on the right (=1) or at the bottom (=2)
axetop     = 1;     % Plot Xaxis on top (=1) or bottom (=0)

%%% Phase velocity and residuals pseudo-section settings (used if plot2dcal=1)
map1       = haxby(32);              % Colormap for phase velocity
map4       = polarmap(32);           % Colormap for phase velocity residuals

lamMIN     = [];                     % Min. wavelength (m)
lamMAX     = [];                     % Max. wavelength (m)
% lticks     = (lamMIN:20:lamMAX);     % Wavelength ticks (m)
vphMIN     = [];                     % Min. phase velocity (m/s)
vphMAX     = [];                     % Max. phase velocity (m/s)
% vphticks   = (vphMIN:200:vphMAX);    % Phase velocity ticks (m/s)
vphISO     = [];                     % Phase velocity isocontours (m/s)
residMIN   = [];                     % Min. residual (%)
residMAX   = [];                     % Max. residual (%)
% residticks = (residMIN:15:residMAX); % Residual ticks (%)

%%% 2D models settings (used if plot2dmod=1)
blocky     = 2;     % Blocky (=0), smooth interp (=1) or smooth contour (=2) images
vertex     = 1;     % Vertical exageration
plottopo   = 1;     % Plot topo profile on 2D sections (=1) or not (=0)
plotDOI    = 0;     % Plot DOI estimated from wavelength (=1), from VsSTD (=2) or not (=0)
maskDOI    = 2;     % Mask below DOI from wavelength (=1), from VsSTD (=2) or not (=0)
transpa    = 1;     % Transparency mask under DOI (=1) or not (=0)
doifact    = 0.66;  % DOI factor (DOI = Lmax*fact)
std_mask   = 400;   % VsSTD limit (m/s) to apply mask (used if plotDOI = 2)
dpMAX      = [];    % Max. depth (m)
plotiso    = 0;     % Plot specific isocontours on all plots (>0) or not (=0)
%%% Vp (=1), Vs (=2), StdVs (=3), Vp/Vs (=4), Poisson's ratio (=5), auxiliary data (=6)


map5       = haxby(32);              % Colormap for Vp and Vs
map6       = haxby(32);              % Colormap for Vp/Vs and Poisson's ratio
map7       = flipud(inferno(32));    % Colormap for VsSTD
map8       = haxby(32);              % Colormap for auxiliary data

xMIN       = [];                     % Min. X (m)
xMAX       = [];                     % Max. X (m)
% xticks     = (xMIN:40:xMAX);         % X ticks (m)
zMIN       = [];                     % Min. altitude (m)
zMAX       = [];                     % Max. altitude (m)
% zticks     = (zMIN:20:zMAX);         % Altitude ticks (m)
vsMIN      = [];                     % Min. Vs (m/s)
vsMAX      = [];                     % Max. Vs (m/s)
% vsticks    = (vsMIN:200:vsMAX);      % Vs ticks (m/s)
vsISO      = [];                     % Vs isocontours (m/s)
vpMIN      = [];                     % Min. Vp (m/s)
vpMAX      = [];                     % Max. Vp (m/s)
% vpticks    = (vpMIN:500:vpMAX);     % Vp ticks (m/s)
vpISO      = [];                     % Vp isocontours (m/s)
vpmask     = 0;                      % Mask Vp with SWIP mask (=1) or not (=0)
stdMIN     = [];                     % Min. StdVs (%) 
stdMAX     = [];                     % Max. STdVs (%) 
% stdticks   = (stdMIN:5:stdMAX);     % Vs STD ticks (%) 
stdISO     = [];                     % STdVs isocontours (%)
vpvsMIN    = [];                     % Min. Vp/Vs
vpvsMAX    = [];                     % Max. Vp/Vs
% vpvsticks  = (vpvsMIN:1:vpvsMAX);    % Vp/Vs ticks
vpvsISO    = [];                     % Vp/Vs isocontours
poisMIN    = [];                     % Min. Poisson's ratio
poisMAX    = [];                     % Max. Poisson's ratio
% poisticks  = (poisMIN:0.1:poisMAX);  % Poisson's ratio ticks
poisISO    = [];                     % Poisson's ratio isocontours
auxMIN     = [];                     % Min. auxiliary data
auxMAX     = [];                     % Max. auxiliary data
% auxticks   = [];                     % Auxiliary data ticks
auxISO     = [];                     % Auxiliary data isocontours
auxlogscal = 0;                      % Auxiliary data log scale (=1) or not (=0)
auxmask    = 0;                      % Mask auxiliary data with SWIP mask (=1) or not (=0)
auxtitle   = ' ';                    % Auxiliary data title
specISO    = [];                     % Specific isocontours

%%% END OF INITIALIZATION %%%
%%%-----------------------%%%

run('D2_SWIPmod2d_script');
