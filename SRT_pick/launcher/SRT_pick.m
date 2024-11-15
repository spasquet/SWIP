clear all; clc; close all;

%%% TRAVELTIME TOMOGRAPHY PROCESSING
%%% MODULE A : TOMOpickfb.m
%%% S. Pasquet - V18.04.26
%%% TOMOpickfb.m allows to pick first breaks and plot seismograms

%%%-------------------------%%%
%%% START OF INITIALIZATION %%%

Sxselec       = [1:5:91];    % Select shots (comment or [] to select all)

%%% Main settings
pick          = 1;     % Pick first arrival (=1) or not (=)
save_seis     = 1;     % Save seismogram image (=1) or not (=)
save_seis_pck = 1;     % Save picked seismogram image (=1) or not (=)
save_src_off  = 1;     % Save source vs offset diagram (=1) or not (=)
save_picks    = 1;
showplot      = 0;

%%% Picking settings
% dt_min        = 0.25;  % Minimum time sample to display seismogram for picking (ms)
offMAX_pick   = [];   % Max. offset (m)
tMAX_pick     = 50;   % Max. time (ms)
scal          = 2;    % Amplitude scaling factor
clip          = 1;    % Clip (=1) or not (=0)
polarity      = -1;   % Fill positive (=1) or negative (-1) amplitudes
autogain      = 1;    % Automatic gain (=1) or not (=0)

%%% Error settings
err_pc        = 1;     % Use percentage error (=1) or absolute (=0)
err_val       = 0.05;  % Error value (in %/100 or ms)
err_val_min   = 0.5;   % Min absolute error when using percentage error (ms)
err_val_max   = 3;     % Max absolute error when using percentage error (ms)

%%% General display settings
imgform       = 'png'; % Fig. file format ('pdf', 'png', 'jpeg', 'tiff' or 'fig')
imgres        = 150;   % Fig. resolution (dpi) when saving as raster
fs            = 16;    % Fig. font size (reduce by 40% with fs=40 => fs=16 in AI) (max=40)

%%% Seismogram plot settings
tMIN_seis     = 0;        % Min. time (ms)
tMAX_seis     = 50;        % Max. time (ms)
% tseisticks    = (tMIN_seis:50:tMAX_seis);

%%% Source vs offset diagram settings
axetop        = 1;         % Plot Xaxis on top (=1) or bottom (=0)
horex         = 0.25;         % Horizontal exageration
plot_pos      = 1;         % Plot traces positions (=1) or not (=0)
marker        = 's';       % Marker type ('s'=square, 'o'=circle)
markersize    = 25;        % Marker size
map1          = haxby(32); % Colormap for traveltimes

% xMIN          = 0;        % Min. X (m)
% xMAX          = 270;        % Max. X (m)
% xticks        = (xMIN:50:xMAX);
offMIN        = [];        % Min. offset (m)
offMAX        = [];        % Max. offset (m)
% offticks      = (offMIN:25:offMAX);
% tMIN          = 0;        % Min. traveltime (ms)
% tMAX          = 100;        % Max. traveltime (ms)
% tticks        = (tMIN:25:tMAX);
vMIN          = [];        % Min. pseudo-velocity (m/s)
vMAX          = [];        % Max. pseudo-velocity (m/s)
% vticks      = (vMIN:250:vMAX);

%%% END OF INITIALIZATION %%%
%%%-----------------------%%%

run('SRT_pick_script');
