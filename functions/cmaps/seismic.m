function rgb = seismic(n)

% seismic(n) creates a colormap, ranging from dark blue via white to dark
%red.
%
% Nico Sneeuw
% Munich, 31/08/94

if nargin == 0, n = size(get(gcf,'colormap'),1); end

m    = ceil(n/3);
top  = ones(m,1);
bot  = zeros(m,1);
up   = (0:m-1)'/m;
down = flipud(up);

r    = [bot; up; 1; top; down];
g    = [bot; up; 1; down; bot];
b    = [up; top; 1; down; bot];
rgb  = [r g b];

% rgb-map has size 4m+1 now. The central part will be extracted.

xlarge = 4*m+1-n;
xblue  = round(xlarge/2);
xred   = xlarge - xblue;
rgb([1:xblue 4*m-xred+2:4*m+1],:) = [];
