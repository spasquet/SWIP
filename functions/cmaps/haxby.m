function map = haxby(m)
%HAXBY  Haxby color map
%   HAXBY(M) returns an M-by-3 matrix containing a colormap with Haxby's
%   colors, commonly used for displaying bathymetry data.
%   HAXBY, by itself, is the same length as the current colormap.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(haxby)
%
%   Use
%             colormap(flipud(haxby))
%
%   for bathymetry data (positive downward).
%
%   Colormap is based on the colors used by W. F. Haxby's Gravity
%   field of World's oceans, 1985, developed for geoid and gravity maps.
%   The version used here is formed from a linear interpolation of
%   the GMT color table used by MB-System by David W. Caress and Dale N. Chayes.
%   <http://www.ldeo.columbia.edu/res/pi/MB-System>
%
%   See also HSV, GRAY, PINK, COOL, BONE, COPPER, FLAG, HOT
%   COLORMAP, RGBPLOT.

% Kelsey Jordahl
% Marymount Manhattan College
% Time-stamp: <Fri Oct 30 12:45:12 EDT 2009>

if nargin < 1, m = size(get(gcf,'colormap'),1); end
% mbm_grdplot Haxby color pallette
ncolors=31;
c=[ 10    0   121;    40   0   150;    20   5   175;   0	10	200;
    0	25	212;   0	40	224;   26	102	240;   13	129	248;
    25	175	255;   50	190	255;   68	202	255;   97	225	240; 106 235 225;
    124	235	200;    138	236	174;   172	245	168;   205	255	162;    223	245	141;
    240	236	121;    247	215	104;   255	189	87; 255	160	69; 244	117	75; 238	80	78;
    255	90	90; 255	124	124; 255 158	158	;245 179 174; 255	196	196; 255 215	215;
    255	230	230];
pp=1:(m-1)/(ncolors-1):m;
r=interp1(pp,c(:,1),1:m,'linear');
g=interp1(pp,c(:,2),1:m,'linear');
b=interp1(pp,c(:,3),1:m,'linear');
map=[r' g' b']/255;


