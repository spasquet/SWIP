clear all; clc; close all;

%%% SURFACE-WAVE dispersion INVERSION & PROFILING (SWIP)
%%% MODULE C : SWIPinv.m
%%% S. Pasquet - V16.9.15
%%% SWIPinv.m performs surface-wave inversion along the seismic profile
%%% and select best models for each Xmid to build a pseudo-2D Vs section
%%% Comment line to use default settings (cf defaultsettings.m)

%%%-------------------------%%%
%%% START OF INITIALIZATION %%%

%%% Main settings
Xmidselec   = [];         % Select Xmids (comment or [] to select all)
inversion   = 1;          % Run inversion (=1) or not (=0)
calcmod     = 1;          % Extract and build average models (=1) or not (=0)

%%% Inversion settings (used if inversion=1)
paramtype   = 0;          % Type of parameterization (=0,1,2,3 or 4) (cf Module B)
nrun        = 3;          % No. of run
itmax       = 150;        % No. of iteration per run
ns0         = 100;        % No. of starting models
ns          = 75;         % No. of models created at each iterations
nr          = 50;         % No. of previous models to build new sub-parameter space
verbose     = 0;          % Display inversion process (=1) or not (=0)

%%% Average models calculation settings (used if calcmod=1)
nbest       = 0;          % Select best models within error bars (=0) or nbest models (>0)
outpoints   = 0;          % No. of points allowed out of the error bars
dz          = 0.2;        % Depth sampling (m)

%%% Plot settings
plotinvres  = 1;          % Plot inversion results (=1) or not (=0)
plotparam   = 0;          % Plot inversion parameters (=1) or not (=0)
plot2dVS    = 1;          % Plot pseudo-2D Vs section (=1) or not (=0)
showplot    = 0;          % Show plots before saving (=1) or not (=0)

%%% General display settings
imgform     = 'png';      % Fig. file format ('pdf', 'png', 'jpeg', 'tiff' or 'fig')
imgres      = 500;        % Fig. resolution (dpi) when saving as raster
fs          = 20;         % Fig. font size
concat      = 1;          % Fig. concatenation panel (=1) or not (=0)
colnb       = 3;          % Number of columns for figure panels
cbpos       = 1;          % Colorbar on the right (=1) or at the bottom (=2)
Clogscale   = 1;          % Log colorscale (=1) or linear (=0)
map2        = hsv(16);    % Colormap for accepted models misfit
map3        = graycm(16); % Colormap for rejected models misfit

%%% Calculated dispersion curves display settings (used if plotinvres=1)
Flogscale   = 0;                   % Logscale frequency axis (=1) or linear (=0)
fMIN        = [];                  % Min. frequency (Hz)
fMAX        = [];                  % Max. frequency (Hz)
% fticks      = (fMIN:15:fMAX);      % Frequency ticks (Hz)
VphMIN      = [];                  % Min. phase velocity (m/s)
VphMAX      = [];                  % Max. phase velocity (m/s)
% Vphticks    = (VphMIN:500:VphMAX); % Phase velocity ticks (m/s)

%%% Calculated models display settings (used if plotinvres=1)
plot1dVS    = 1;                   % Plot final 1D Vs model (=1) or not (=0)
modeltype   = 5;                   % Plotted 1D Vs model
%%% (1=best, 2=averaged layered, 3=average smooth, 4=weighted layered, 5=weighted smooth, 6=ridge)

dpMIN       = [];                  % Min. depth (m) 
dpMAX       = [];                  % Max. depth (m) 
% dticks      = (dpMIN:10:dpMAX);    % Depth ticks (m) 
vsMIN       = [];                  % Min. Vs (m/s) 
vsMAX       = [];                  % Max. Vs (m/s) 
% vsticks     = (vsMIN:500:vsMAX);   % Vs ticks (m/s) 

%%% Parameter plot settings (used if plotparam=1)
param1      = 'Vs';                % First parameter to plot ('Vs', 'Th', 'Vp', 'Dens') 
param2      = 'Th';                % Second parameter to plot ('Vs', 'Th', 'Vp', 'Dens') 
np1         = (1:2);               % Vector of layer number for first parameter
np2         = (1:2);               % Vector of layer number for second parameter

%%% END OF INITIALIZATION %%%
%%%-----------------------%%%

run('C_SWIPinv_script');