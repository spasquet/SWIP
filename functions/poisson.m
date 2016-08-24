function pois = poisson(vp,vs)

%%% S. Pasquet - V16.3.10
% Compute Poisson's ratio

if vp<vs
    pois=0;
else
    pois=(0.5*(vp./vs).^2-1)./((vp./vs).^2-1);
end
