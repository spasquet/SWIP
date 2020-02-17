clear all; clc; close all;

%%% SURFACE-WAVE dispersion INVERSION & PROFILING (SWIP)
%%% MODULE C : SWIPinv.m
%%% S. Pasquet - V20.02.17
%%% SWIPinv.m performs inversion of dispersion curves picked in module A
%%% and select best models for each Xmid to build a pseudo-2D Vs section
%%% It allows to plot all generated models and inversion parameters for
%%% each Xmid with misfit based colorbars for selected and rejected models

%%% Option inversion=1 starts a new inversion
%%% Requires to select a subproject folder Wmin_max.dWx.dSmin_max.side
%%% containing target files with dispersion curves created in module A
%%% Requires to run Module B beforehand to create parameterization
%%% Option inversion=0 works with previously run inversion
%%% Requires to select a subproject folder Wmin_max.dWx.dSmin_max.side
%%% and an inversion folder located in file.inv

%%% Option paramtype=0 uses a single user-defined parameterization
%%% Requires to select a .param file containing parameterization
%%% Option paramtype>0 uses a semi-automatic parameterization
%%% Requires to have .param files for each Xmid along with .target files

%%% Comment line to use default settings (cf SWIP_defaultsettings.m)

%%%-------------------------%%%
%%% START OF INITIALIZATION %%%

%%% Main settings
Xmidselec  = [];         % Select Xmids (comment or [] to select all)
inversion  = 1;          % Run new inversion (=1) or use existing (=0)

%%% Inversion settings (used if inversion=1)
paramtype  = 0;          % Type of parameterization (=0,1,2,3 or 4) (cf Module B)
nrun       = 2;          % No. of run
itmax      = 150;        % No. of iteration per run
ns0        = 100;        % No. of starting models
ns         = 75;         % No. of models created at each iterations
nr         = 50;         % No. of previous models to build new sub-parameter space
verbose    = 0;          % Display info during inversion (=1) or not (=0) (!! set to 1 to debug !!)

%%% Average models calculation settings
nbest      = 0;          % Select best models within error bars (=0) or nbest models (>0)
outpoints  = 0;          % No. of points allowed out of the error bars
dz         = 0.2;        % Depth sampling (m)

%%% Toggle plots
plotinvres = 0;          % Save inversion results (=1) or not (=0) (!! time consuming !!)
plotparam  = 0;          % Save inversion parameters (=1) or not (=0) (!! time consuming !!)
plot2dVS   = 1;          % Plot raw pseudo-2D Vs section during inversion (=1) or not (=0)
showplot   = 0;          % Show plots before saving (=1) or not (=0)

%%% Figure display and output settings
imgform    = 'png';      % Fig. file format ('pdf', 'png', 'jpeg', 'tiff' or 'fig')
imgres     = 250;        % Fig. resolution (dpi) when saving as raster
fs         = 20;         % Fig. font size
concat     = 1;          % Save indiv. and merged figures (=2), merged figure (=1) or indiv. figures (=0)
colnb      = 3;          % No. of columns when merging figures
cbpos      = 2;          % Colorbar on the right (=1) or at the bottom (=2)
Clogscale  = 1;                   % Log colorscale (=1) or linear (=0) for misfit
map2       = hsv(16);             % Colormap for accepted models misfit
map3       = graycm(16);          % Colormap for rejected models misfit

%%% Calculated dispersion curves display settings (used if plotinvres=1)
Flogscale  = 0;                   % Logscale frequency axis (=1) or linear (=0)
fMIN       = [];                  % Min. frequency (Hz)
fMAX       = [];                  % Max. frequency (Hz)
% fticks     = (fMIN:15:fMAX);      % Frequency ticks (Hz)
VphMIN     = [];                  % Min. phase velocity (m/s)
VphMAX     = [];                  % Max. phase velocity (m/s)
% Vphticks   = (VphMIN:500:VphMAX); % Phase velocity ticks (m/s)

plot1dVS   = 1;                   % Plot final 1D Vs model (=1) or not (=0)
modeltype  = 3;                   % 1D final Vs model to plot
%%% (1=best, 2=averaged layered, 3=average smooth, 4=weighted layered, 5=weighted smooth, 6=ridge)

dpMIN      = [];                  % Min. depth (m)
dpMAX      = [];                  % Max. depth (m)
% dticks     = (dpMIN:10:dpMAX);    % Depth ticks (m)
vsMIN      = [];                  % Min. Vs (m/s)
vsMAX      = [];                  % Max. Vs (m/s)
% vsticks    = (vsMIN:500:vsMAX);   % Vs ticks (m/s)

%%% Parameter plot settings (used if plotparam=1)
param1     = 'Vs';       % First parameter to plot ('Vs', 'Th', 'Vp', 'Dens')
param2     = 'Th';       % Second parameter to plot ('Vs', 'Th', 'Vp', 'Dens')
np1        = (1:2);      % Vector of layer number for first parameter
np2        = (1:2);      % Vector of layer number for second parameter

%%% END OF INITIALIZATION %%%
%%%-----------------------%%%

run('C_SWIPinv_script');
