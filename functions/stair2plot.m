function [Vplot,Zplot]=stair2plot(V,moddepth)

%%% S. Pasquet - V16.11.18
% [Vplot,Zplot]=stair2plot(V,moddepth)

Vplot=[V;V(end)];
Vplotrep=repmat(Vplot,1,2)';
Vplot=Vplotrep(:)'; Vplot=Vplot(1:end-2)';
moddepthrep=repmat(moddepth,1,2)';
Zplot=moddepthrep(:)'; Zplot=Zplot(2:end-1)';

end