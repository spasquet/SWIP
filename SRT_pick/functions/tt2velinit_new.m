function [F2, pseudo_X, pseudo_Z, pseudo_V] = tt2velinit_new(Tpick, Sxpick, Szpick, Gxpick, Gzpick, extrap, c1, c2)

%%
array  = check_depth_array(Gxpick,Gzpick,Sxpick,Szpick);
 
IGpick = []; ISpick = [];
I_depth_surf = zeros(size(Gxpick));
I_depth_depth = zeros(size(Gxpick));

for i = [2 4]
    if ~isempty(array{i})
        IGpick = ismember([Gxpick Gzpick],array{i}.G_sing,'rows');
        ISpick = ismember([Sxpick Szpick],array{i}.S_sing,'rows');
        I_depth_surf = I_depth_surf | (IGpick & ISpick);
    end
end
if ~isempty(array{3})
    IGpick = ismember([Gxpick Gzpick],array{3}.G_sing,'rows');
    ISpick = ismember([Sxpick Szpick],array{3}.S_sing,'rows');
    I_depth_depth = I_depth_depth | (IGpick & ISpick);
end

a = 0;
b = (1-a)/2;

dist_X = (Gxpick - Sxpick);
dist_Z = (Gzpick - Szpick);
dist_AB = sqrt(dist_X.^2 + dist_Z.^2);

t_direct = (Tpick/1000);
slope_direct = t_direct./dist_AB;
slope_direct = slope_direct/nanmax(slope_direct(~isinf(slope_direct)));
% c = 0.1+0.15*10.^-(t_direct.*dist_AB);
c = (c1+(0.75-0.75*exp(-c2*(slope_direct).^2)))/(1+c1);

% plot_scat([],dist_AB,c,Tpick);

% % keyboard

% plot(dist_AB,c,'x')
% keyboard
% plot_scat([],Tpick,dist_AB,c)
% keyboard
% c = 0.25;
% figure;
% plot(c,Tpick);
% figure;
% plot(dist_AB.*c,abs(dist_AB)./(Tpick/1000));
% return
ind1 = (Gxpick>Sxpick & Gzpick>Szpick & ~I_depth_depth & ~I_depth_surf);
ind2 = (Gxpick>Sxpick & Gzpick<Szpick & ~I_depth_depth & ~I_depth_surf);
ind3 = (Gxpick<Sxpick & Gzpick>Szpick & ~I_depth_depth & ~I_depth_surf);
ind4 = (Gxpick<Sxpick & Gzpick<Szpick & ~I_depth_depth & ~I_depth_surf);

ind5 = (Gxpick>Sxpick & Gzpick>Szpick & I_depth_depth);
ind6 = (Gxpick>Sxpick & Gzpick<Szpick & I_depth_depth);
ind7 = (Gxpick<Sxpick & Gzpick>Szpick & I_depth_depth);
ind8 = (Gxpick<Sxpick & Gzpick<Szpick & I_depth_depth);

ind9 = (Gxpick>Sxpick & Gzpick>Szpick & I_depth_surf);
ind10 = (Gxpick>Sxpick & Gzpick<Szpick & I_depth_surf);
ind11 = (Gxpick<Sxpick & Gzpick>Szpick & I_depth_surf);
ind12 = (Gxpick<Sxpick & Gzpick<Szpick & I_depth_surf);

xA = Sxpick;
zA = Szpick;
xB = Gxpick;
zB = Gzpick;

alpha = (atan(dist_Z./dist_X));
beta = (atan((dist_AB.*c) ./ (dist_AB*b)));
gamma = (pi/2 - alpha - beta);

dist_AC = sqrt((dist_AB*b).^2 + (dist_AB.*c).^2);

% keyboard
xC = xA + dist_AC.*sin(-gamma);
xC(ind1) = xB(ind1) + dist_AC(ind1).*sin(-gamma(ind1));
xC(ind2) = xB(ind2) + dist_AC(ind2).*sin(-gamma(ind2));
xC(ind5 | ind6) = xA(ind5 | ind6) + dist_AC(ind5 | ind6).*sin(gamma(ind5 | ind6));
xC(ind7 | ind8) = xB(ind7 | ind8) + dist_AC(ind7 | ind8).*sin(gamma(ind7 | ind8));
xC(ind9) = xA(ind9) + 0.5 * abs(dist_X(ind9));
xC(ind10) = xA(ind10) + 0.5 * abs(dist_X(ind10));
xC(ind11) = xB(ind11) + 0.5 * abs(dist_X(ind11));
xC(ind12) = xB(ind12) + 0.5 * abs(dist_X(ind12));

% xC(~ind1 & ~ind2 & ~ind3 & ~ind4 & ~ind5 & ~ind6 & ~ind7 & ~ind8 & ~ind9 & ~ind10 & ~ind11 & ~ind12)

zC = zA - dist_AC.*cos(gamma);
zC(ind1) = zB(ind1) - dist_AC(ind1).*cos(gamma(ind1));
zC(ind2) = zB(ind2) - dist_AC(ind2).*cos(gamma(ind2));
zC(ind5 | ind6) = zA(ind5 | ind6) + dist_AC(ind5 | ind6).*cos(gamma(ind5 | ind6));
zC(ind7 | ind8) = zB(ind7 | ind8) + dist_AC(ind7 | ind8).*cos(gamma(ind7 | ind8));

zC(ind9) = zA(ind9) + 0.5 * abs(dist_Z(ind9));
zC(ind10) = zA(ind10) - 0.5 * abs(dist_Z(ind10));
zC(ind11) = zA(ind11) + 0.5 * abs(dist_Z(ind11));
zC(ind12) = zA(ind12) - 0.5 * abs(dist_Z(ind12));

dist_AC = sqrt((xC-xA).^2 + (zC-zA).^2);
dist_BC = sqrt((xC-xB).^2 + (zC-zB).^2);
diff_AC_BC = round(dist_AC - dist_BC);

xC(diff_AC_BC~=0) = NaN;
zC(diff_AC_BC~=0) = NaN;

% xC(I_depth_surf) = NaN;
% zC(I_depth_surf) = NaN;

% xE = xC + (dist_AB*(a/2)).*cos(alpha);
% zE = zC + (dist_AB*(a/2)).*sin(alpha);

% zC(I_depth_depth) = zA(I_depth_depth) + dist_AC(I_depth_depth).*cos(gamma(I_depth_depth));
% zE(I_depth_depth) = zC(I_depth_depth) - (dist_AB(I_depth_depth)*(a/2)).*sin(alpha(I_depth_depth));

pseudo_X = xC;
pseudo_Z = zC;
dist_ray = 2*dist_AC + dist_AB*a;
% dist_ray = dist_AB;

% a_vel = a_scal*abs(Offpick).^(1/10);
pseudo_V = dist_ray./(Tpick*1e-3);
if any(I_depth_surf)
    pseudo_V(I_depth_surf) = sqrt((xA(I_depth_surf) - xB(I_depth_surf)).^2 + (zA(I_depth_surf) - zB(I_depth_surf)).^2)./(Tpick(I_depth_surf)*1e-3);
end
a_vel = 1;
pseudo_V = a_vel.*pseudo_V;
%%
% figure(4);
% plot(Gxpick,Gzpick,'k.','markersize',10); xlim([0 190]); ylim([630 710]); hold on;
% plot(xA(20),zA(20),'m.','markersize',30); 
% plot(xB(20),zB(20),'m.','markersize',30); 
% plot(pseudo_X(20),pseudo_Z(20),'m.','markersize',30); 
% 
% keyboard
% 
% plot(xA(50),zA(50),'c.','markersize',20); 
% plot(xB(50),zB(50),'c.','markersize',20); 
% plot(pseudo_X(50),pseudo_Z(50),'c.','markersize',20); 
% 
% plot(xA(ind7(20)),zA(ind7(20)),'r.','markersize',15); 
% plot(xB(ind7(20)),zB(ind7(20)),'r.','markersize',15); 
% plot(pseudo_X(ind7(20)),pseudo_Z(ind7(20)),'r.','markersize',15); 
% 
% plot(xA(ind8(20)),zA(ind8(20)),'b.','markersize',15); 
% plot(xB(ind8(20)),zB(ind8(20)),'b.','markersize',15); 
% plot(pseudo_X(ind8(20)),pseudo_Z(ind8(20)),'b.','markersize',15); 
% 
% plot(Gxpick(end-15),Gzpick(end-15),'y.','markersize',20); 
% plot(Sxpick(end-15),Szpick(end-15),'y.','markersize',20); 
% plot(pseudo_X(end-15),pseudo_Z(end-15),'y.','markersize',20); 
% 
% plot(xA(1000),zA(1000),'k.','markersize',20); 
% plot(xB(1000),zB(1000),'k.','markersize',20); 
% plot(pseudo_X(1000),pseudo_Z(1000),'k.','markersize',20); 
% 
% hold off;


% close all;
% plot_scat([],pseudo_X,pseudo_Z,pseudo_V,[],[],[],[],[],[],1); axis equal
% hold on; plot(Gxpick,Gzpick,'r.','markersize',10); %xlim([0 190]); ylim([630 700]);

% plot_scat([],pseudo_X(ind1),pseudo_Z(ind1),pseudo_V(ind1),[],[],[],[],[],[],1); axis equal
% hold on; plot(Gxpick,Gzpick,'k.','markersize',10); %xlim([0 190]); ylim([630 700]);
% 
% plot_scat([],pseudo_X(ind2),pseudo_Z(ind2),pseudo_V(ind2),[],[],[],[],[],[],1); axis equal
% hold on; plot(Gxpick,Gzpick,'k.','markersize',10); %xlim([0 190]); ylim([630 700]);
% 
% plot_scat([],pseudo_X(ind3),pseudo_Z(ind3),pseudo_V(ind3),[],[],[],[],[],[],1); axis equal
% hold on; plot(Gxpick,Gzpick,'k.','markersize',10); %xlim([0 190]); ylim([630 700]);
% 
% plot_scat([],pseudo_X(ind4),pseudo_Z(ind4),pseudo_V(ind4),[],[],[],[],[],[],1); axis equal
% hold on; plot(Gxpick,Gzpick,'k.','markersize',10); %xlim([0 190]); ylim([630 700]);
% 
% plot_scat([],pseudo_X(ind5),pseudo_Z(ind5),pseudo_V(ind5),[],[],[],[],[],[],1); axis equal
% hold on; plot(Gxpick,Gzpick,'k.','markersize',10); %xlim([0 190]); ylim([630 700]);
% 
% % plot_scat([],pseudo_X(ind6),pseudo_Z(ind6),pseudo_V(ind6),[],[],[],[],[],[],1); axis equal
% % hold on; plot(Gxpick,Gzpick,'k.','markersize',10); %xlim([0 190]); ylim([630 700]);
% % 
% % plot_scat([],pseudo_X(ind7),pseudo_Z(ind7),pseudo_V(ind7),[],[],[],[],[],[],1); axis equal
% % hold on; plot(Gxpick,Gzpick,'k.','markersize',10); %xlim([0 190]); ylim([630 700]);
% 
% plot_scat([],pseudo_X(ind8),pseudo_Z(ind8),pseudo_V(ind8),[],[],[],[],[],[],1); axis equal
% hold on; plot(Gxpick,Gzpick,'k.','markersize',10); %xlim([0 190]); ylim([630 700]);

% plot_scat([],pseudo_X(ind9),pseudo_Z(ind9),pseudo_V(ind9),[],[],[],[],[],[],1); axis equal
% hold on; plot(Gxpick,Gzpick,'k.','markersize',10); %xlim([0 190]); ylim([630 700]);
% 
% plot_scat([],pseudo_X(ind10),pseudo_Z(ind10),pseudo_V(ind10),[],[],[],[],[],[],1); axis equal
% hold on; plot(Gxpick,Gzpick,'k.','markersize',10); %xlim([0 190]); ylim([630 700]);
% 
% plot_scat([],pseudo_X(ind11),pseudo_Z(ind11),pseudo_V(ind11),[],[],[],[],[],[],1); axis equal
% hold on; plot(Gxpick,Gzpick,'k.','markersize',10); %xlim([0 190]); ylim([630 700]);
% 
% plot_scat([],pseudo_X(ind12),pseudo_Z(ind12),pseudo_V(ind12),[],[],[],[],[],[],1); axis equal
% hold on; plot(Gxpick,Gzpick,'k.','markersize',10); %xlim([0 190]); ylim([630 700]);

% plot_scat([],pseudo_X(I_depth_depth),pseudo_Z(I_depth_depth),pseudo_V(I_depth_depth),[],[],[],[],[],[],1); axis equal
% hold on; plot(Gxpick,Gzpick,'k.','markersize',10); %xlim([0 190]); ylim([630 700]);
% 
% plot_scat([],pseudo_X(I_depth_surf),pseudo_Z(I_depth_surf),pseudo_V(I_depth_surf),[],[],[],[],[],[],1); axis equal
% hold on; plot(Gxpick,Gzpick,'k.','markersize',10); %xlim([0 190]); ylim([630 700]);

%%
% keyboard
ind_nan = isnan(pseudo_X) | isnan(pseudo_Z) | isnan(pseudo_V);
pseudo_X(ind_nan) = [];
pseudo_Z(ind_nan) = [];
pseudo_V(ind_nan) = [];

ind_inf = isinf(pseudo_X) | isinf(pseudo_Z) | isinf(pseudo_V);
pseudo_X(ind_inf) = [];
pseudo_Z(ind_inf) = [];
pseudo_V(ind_inf) = [];

if extrap == 1
    F2 = scatteredInterpolant(pseudo_X,pseudo_Z,pseudo_V,'natural','nearest');
else
    F2 = scatteredInterpolant(pseudo_X,pseudo_Z,pseudo_V,'natural','none');
end