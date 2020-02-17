function [gauss_weight, gauss_center_new] = stack_weight_lam_new(f,v,maxwinsize)

%%% S. Pasquet - V17.08.28
% Calculate wavelength based weights to stack multi-window dispersion
%
% [gauss_weight, gauss_center] = stack_weight_lam_new(f,v,maxwinsize)

a = 1.5; % increase a for larger gaussian (default 1.5)

[vv,ff] = meshgrid(v,f);
ll = vv./ff;

max_resamp_win = maxwinsize.*10.^(1./sqrt(0.5*maxwinsize));

% gauss_center = linspace(0,max(max_resamp_win),length(maxwinsize));
% diff_gauss_center = mean(diff(gauss_center));

% diff_gauss_center = mean(diff(max_resamp_win));
% gauss_center = max_resamp_win - [2*diff(max_resamp_win)/3 0];

diff_gauss_center = (max(max_resamp_win) - min(max_resamp_win))/length(max_resamp_win);
gauss_center = max_resamp_win;
% threshold = 0.01;

%%
% close all
gauss_weight = zeros(size(ll,1),size(ll,2),length(gauss_center));
% cm = haxby(length(gauss_center));
for i = 1:length(gauss_center)
    ll_uni = unique(ll);
        
%     gauss_tst = exp(-0.5*((gauss_center(i)-ll_uni)./diff_gauss_center*a).^2);
    
    sigma = 1/sqrt(2*pi);
    gauss_tst = 1/(sigma*sqrt(2*pi))*exp(-0.5*((ll_uni-gauss_center(i)+(3*diff_gauss_center*sigma))/(sigma*diff_gauss_center*a)).^2);
    gauss_tmp = 1/(sigma*sqrt(2*pi))*exp(-0.5*((ll-gauss_center(i)+(3*diff_gauss_center*sigma))/(sigma*diff_gauss_center*a)).^2);
    
%     ind_min = find(abs(gauss_tst(ll_uni>gauss_center(i)) - threshold) == min(abs(gauss_tst(ll_uni>gauss_center(i)) - threshold)),1);
%     ll_uni_up = ll_uni(ll_uni>gauss_center(i));
%     gauss_center(i) = gauss_center(i) - (ll_uni_up(ind_min) - max_resamp_win(i));
%     gauss_tmp = exp(-0.5*((gauss_center(i)-ll)./diff_gauss_center*a).^2);
    
    gauss_center_new(i) = gauss_center(i)-(3*diff_gauss_center*sigma);
    if i == length(gauss_center)
        gauss_tmp(ll>gauss_center_new(i)) = 1;
        gauss_tst(ll_uni>gauss_center_new(i)) = 1;
    elseif i == 1
        gauss_tmp(ll<gauss_center_new(i)) = 1;
        gauss_tst(ll_uni<gauss_center_new(i)) = 1;
    end
    
    gauss_weight(:,:,i) = gauss_tmp;
    
%     plot_img([],vv,ll,gauss_weight(:,:,i),[],1,1,1,[],[],[],[],[],[0 max(max_resamp_win)]);
%     
%     figure(10); plot(ll_uni,gauss_tst,'-','linewidth',3);  hold on;
%     xlim([0 max(max_resamp_win)]);
%     line([0 140],[threshold threshold]);
%     line([max_resamp_win(i) max_resamp_win(i)],[0 1]);
end



% hold off;