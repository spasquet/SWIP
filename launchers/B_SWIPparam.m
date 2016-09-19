clear all; clc; close all;

%%% SURFACE-WAVE dispersion INVERSION & PROFILING (SWIP)
%%% MODULE B : SWIPparam.m
%%% S. Pasquet - V16.9.15
%%% SWIPparam.m creates the parameterization file for surface-wave inversion
%%% It can create either a global file for all the profile or create
%%% automatic parameterization for specific Xmid based on a velocity model
%%% Comment line to use default settings (cf defaultsettings.m)

%%%-------------------------%%%
%%% START OF INITIALIZATION %%%

%%% Main settings
Xmidselec = [];   % Select Xmids (comment or [] to select all)
paramname = [];   % Name of parameterization file (empty for autoname)
paramtype = 0;    % Type of parameterization
%%% 0 => default inversion (same parameterization for all Xmids)
%%% 1 => Vp and thicknesses vary in reduced range defined from velocity model
%%% 2 => Fixed Vp defined from velocity model - thicknesses vary with Vs
%%% 3 => Fixed thicknesses - Vp vary in reduced range defined from velocity model
%%% 4 => Fixed Vp and thicknesses defined from velocity model

%%% Parameter space settings
nlay      = 11;   % No. of layers (including half-space)

%%% The following settings can be either scalars (same parameters for all layers)
%%% or vectors with nlay elements (specific parameters for each layer)

nsublay   = 10;   % No. of sublayers per layer
thmin     = 0.5;  % Min. thickness per layer (m)
thmax     = 2.5;    % Max. thickness per layer (m)

lvz       = 0;    % Allow low velocity layer (=1) or not (=0)
shape     = 1;    % Shape of the velocity variation with depth
%%% (1='Uniform', 2='Linear', 3='LinearIncrease', 4='LinearDecrease', 5='PowerLaw')

Vsmin     = 10;   % Min. Vs (m/s)
Vsmax     = 2500; % Max. Vs (m/s)
Vpmin     = 10;   % Min. Vp (m/s)
Vpmax     = 5000; % Max. Vp (m/s)
Rhomin    = 2000; % Min. Rho (kg/m3)
Rhomax    = 2000; % Max. Rho (kg/m3)
Numin     = 0.1;  % Min. Poisson's ratio
Numax     = 0.5;  % Max. Poisson's ratio

Vplink    = 1;    % Vp linked (=1) or not (=0) to Vs
Rholink   = 1;    % Rho linked (=1) or not (=0) to Vs
Nulink    = 1;    % Nu linked (=1) or not (=0) to Vs

%%% Semi-automatic parameterization settings (used if paramtype=1,2,3 or 4)
plot2dVP  = 1;    % Plot imported Vp models (=1) or not (=0)
dz        = 0.2;  % Sampling in depth (m)
vfac      = 0.25; % Increase Vp range (Vpmin-vfac*Vpmin<v<Vpmax+vfac*Vpmax)
linklvz   = 0;    % Link Vs LVZ with Vp LVZ (=1) or not (=0)

%%% END OF INITIALIZATION %%%
%%%-----------------------%%%

run('B_SWIPparam_script');