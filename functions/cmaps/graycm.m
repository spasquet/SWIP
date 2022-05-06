function map=graycm(n)

%%% S. Pasquet - V16.6.30
% Gray colormap

map=repmat(linspace(0.2,0.8,n),3,1)';
map = flipud(map);
end