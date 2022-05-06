function pois = poisson(vp,vs)

%%% S. Pasquet - V16.11.18
% Compute Poisson's ratio
% pois = poisson(vp,vs)

pois=(0.5*(vp./vs).^2-1)./((vp./vs).^2-1);
pois(vp<vs)=0;
pois(pois<0)=0;
